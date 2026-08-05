"""Robust CSV reader with automatic format detection (SFPPy / tools).

Handles common CSV variations encountered in heterogeneous tabular data:
- Multiple encodings (UTF-8 with/without BOM, Latin-1, CP1252, Mac Roman)
- Multiple delimiters (tab, semicolon, comma)
- Multiple quoting modes (minimal, all, none)
- Stray quotes, BOM characters, line ending variations

Author: Olivier Vitrac, PhD, HDR | olivier.vitrac@gmail.com | INRAE + Generative Simulation
"""

import csv
import io
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Optional, Union
import pandas as pd
import logging

logger = logging.getLogger(__name__)


@dataclass
class CSVReadResult:
    """Result of robust CSV reading attempt."""

    df: Optional[pd.DataFrame] = None
    success: bool = False
    encoding: Optional[str] = None
    delimiter: Optional[str] = None
    quoting: Optional[int] = None
    num_columns: int = 0
    num_rows: int = 0
    error: Optional[str] = None
    strategy_index: int = -1
    # Field-shift guard: fraction of sampled data rows whose RAW field count matches the
    # header width under the chosen (delimiter, quoting). 1.0 = every row aligns;
    # < 1.0 means some rows carry extra/missing fields (field-shift risk — the
    # mechanism that contaminates a column with a neighbour's value). Surfaced so
    # misalignment is visible instead of silently coerced by on_bad_lines='warn'.
    row_arity_consistency: float = 1.0

    def __repr__(self):
        if self.success:
            return (f"CSVReadResult(success=True, {self.num_rows}×{self.num_columns}, "
                   f"encoding={self.encoding}, delimiter={repr(self.delimiter)}, "
                   f"quoting={self.quoting})")
        else:
            return f"CSVReadResult(success=False, error={self.error})"


def read_csv_robust(
    filepath: Union[str, Path],
    min_columns: int = 1,
    verbose: bool = False,
    **pandas_kwargs
) -> CSVReadResult:
    """Read CSV with automatic format detection and fallback strategies.

    Tries multiple encoding/delimiter/quoting combinations until successful parse.

    Args:
        filepath: Path to CSV file
        min_columns: Minimum expected columns (validation threshold, default=1)
        verbose: Log all attempts
        **pandas_kwargs: Additional arguments passed to pd.read_csv()

    Returns:
        CSVReadResult with parsed DataFrame and metadata

    Strategy:
        1. Detect delimiter from file extension and content sniffing
        2. Try encodings in order: utf-8-sig, utf-8, latin-1, cp1252, mac_roman
        3. Try quoting modes: QUOTE_MINIMAL, QUOTE_NONE, QUOTE_ALL
        4. Validate: check column count, row count, no null column names
        5. Return first successful parse

    Example:
        >>> result = read_csv_robust("data.csv")
        >>> if result.success:
        ...     df = result.df
        ...     print(f"Loaded {result.num_rows} rows with {result.encoding}")
    """
    filepath = Path(filepath)

    if not filepath.exists():
        return CSVReadResult(
            success=False,
            error=f"File not found: {filepath}"
        )

    # Handle empty files gracefully
    if filepath.stat().st_size == 0:
        return CSVReadResult(
            df=pd.DataFrame(),
            success=True,
            encoding=None,
            delimiter=None,
            quoting=None,
            num_columns=0,
            num_rows=0,
            strategy_index=-3  # Special: empty file
        )

    # Encoding priority (most specific to most permissive)
    encodings = [
        'utf-8-sig',      # UTF-8 with BOM (strips BOM automatically)
        'utf-8',          # Standard UTF-8
        'latin-1',        # ISO-8859-1 (very permissive, never fails)
        'cp1252',         # Windows-1252 (common for Excel exports)
        'mac_roman',      # Common for legacy Mac exports
    ]

    # Delimiter detection (priority order)
    delimiters = _detect_delimiters(filepath)

    # Quoting modes (from most to least common)
    quoting_modes = [
        csv.QUOTE_MINIMAL,  # Standard CSV quoting
        csv.QUOTE_NONE,     # No quoting (handles quoted TSV)
        csv.QUOTE_ALL,      # Everything quoted
    ]

    # Build parsing strategies (cartesian product)
    strategies = []
    strategy_idx = 0

    for encoding in encodings:
        for delimiter in delimiters:
            for quoting in quoting_modes:
                strategies.append({
                    'index': strategy_idx,
                    'encoding': encoding,
                    'sep': delimiter,
                    'quoting': quoting,
                    'on_bad_lines': 'warn',
                    'engine': 'python' if quoting == csv.QUOTE_NONE else 'c',
                })
                strategy_idx += 1

    # First: Try special handling for "row-level quoted TSV" format
    # Format: "col1\tcol2\tcol3",\r\n (entire row is quoted with trailing comma)
    row_quoted_result = _try_row_level_quoted_tsv(filepath, min_columns, verbose, pandas_kwargs)
    if row_quoted_result.success:
        return row_quoted_result

    # Second: Try multi-space (fixed-width) format
    # Format: "Employee_ID     Name            Department" (2+ spaces between columns)
    multi_space_result = _try_multi_space_delimited(filepath, min_columns, verbose, pandas_kwargs)
    if multi_space_result.success:
        return multi_space_result

    # Try each standard strategy. Arity-gated selection: a well-formed
    # file returns on the first strategy whose rows align with the header (the
    # fast path = legacy first-success behaviour). Only a RAGGED first parse
    # re-ranks the remaining strategies by (multi-column, row-consistency), so
    # field-shift (wrong quoting / embedded delimiters) is caught and surfaced —
    # without ever failing a file that previously parsed.
    last_error = None
    HIGH_CONSISTENCY = 0.99
    best = None  # (rank_key, result)

    for strategy in strategies:
        try:
            # Merge with user-provided kwargs (user takes precedence)
            params = {**strategy, **pandas_kwargs}
            index = params.pop('index')

            if verbose:
                logger.info(f"Strategy {index}: encoding={params['encoding']}, "
                          f"sep={repr(params['sep'])}, quoting={params['quoting']}")

            # Attempt read
            df = pd.read_csv(filepath, **params)

            # Validation checks
            if len(df.columns) < min_columns:
                if verbose:
                    logger.debug(f"  ✗ Only {len(df.columns)} columns (min={min_columns})")
                continue

            # Check for null/empty column names
            if df.columns.isnull().any() or (df.columns == '').any():
                if verbose:
                    logger.debug(f"  ✗ Null or empty column names")
                continue

            # Post-processing: Clean column names
            df.columns = df.columns.str.strip()

            # Remove BOM from first column if present (fallback)
            if df.columns[0].startswith('\ufeff'):
                first_col = df.columns[0].lstrip('\ufeff')
                df.rename(columns={df.columns[0]: first_col}, inplace=True)

            # Remove stray quotes from column names
            df.columns = df.columns.str.strip('"').str.strip("'")

            # Raw per-row arity under this (delimiter, quoting).
            consistency = _row_arity_consistency(
                filepath, params['encoding'], params['sep'], params['quoting'])

            result = CSVReadResult(
                df=df,
                success=True,
                encoding=params['encoding'],
                delimiter=params['sep'],
                quoting=params['quoting'],
                num_columns=len(df.columns),
                num_rows=len(df),
                strategy_index=index,
                row_arity_consistency=consistency,
            )

            # Fast path: an aligned parse wins immediately (legacy behaviour).
            if consistency >= HIGH_CONSISTENCY:
                logger.info(f"✓ Successfully parsed: {result.num_rows} rows × "
                          f"{result.num_columns} columns "
                          f"(encoding={result.encoding}, delimiter={repr(result.delimiter)}, "
                          f"quoting={result.quoting}, strategy={index}, "
                          f"arity={consistency:.3f})")
                return result

            # Ragged parse: remember it, prefer multi-column then higher
            # consistency, and keep scanning for a better-aligned strategy.
            rank = (len(df.columns) >= 2, round(consistency, 3), len(df.columns))
            if best is None or rank > best[0]:
                best = (rank, result)
            if verbose:
                logger.debug(f"  ~ ragged parse (arity={consistency:.3f}), "
                             f"keeping as candidate")
            continue

        except Exception as e:
            last_error = str(e)
            if verbose:
                logger.debug(f"  ✗ Error: {e}")
            continue

    # No strategy reached the consistency bar — return the best-ranked parse we
    # obtained (never fail a file that parsed), and surface the misalignment so
    # downstream consumers see it instead of silent field-shift.
    if best is not None:
        result = best[1]
        logger.warning(
            "CSV parsed with row misalignment: %.1f%% of data rows match the "
            "header width (encoding=%s, delimiter=%r, quoting=%s, strategy=%d). "
            "Possible field-shift / embedded-delimiter contamination — verify the "
            "source delimiter and quoting.",
            result.row_arity_consistency * 100.0, result.encoding,
            result.delimiter, result.quoting, result.strategy_index)
        return result

    # All strategies failed
    return CSVReadResult(
        success=False,
        error=f"All {len(strategies)} parsing strategies failed. Last error: {last_error}"
    )


def _detect_delimiters(filepath: Path) -> list:
    """Detect likely delimiters from filename and content."""
    delimiters = []

    # Heuristic 1: File extension
    name_lower = filepath.name.lower()
    if '.tsv' in name_lower or 'tab' in name_lower:
        delimiters.append('\t')

    # Heuristic 2: Sniff first line (limit to 8KB for safety)
    try:
        with open(filepath, 'rb') as f:
            # Read first line (up to 8KB)
            sample = f.read(8192).decode('utf-8', errors='ignore')
            first_line = sample.split('\n')[0] if '\n' in sample else sample

            # Count occurrences
            tab_count = first_line.count('\t')
            semicolon_count = first_line.count(';')
            comma_count = first_line.count(',')

            # Delimiter priority based on counts
            counts = [
                (tab_count, '\t'),
                (semicolon_count, ';'),
                (comma_count, ','),
            ]

            # Sort by count (descending)
            counts.sort(reverse=True, key=lambda x: x[0])

            for count, delim in counts:
                if count > 0 and delim not in delimiters:
                    delimiters.append(delim)

    except Exception as e:
        logger.debug(f"Delimiter detection failed: {e}")

    # Fallback: standard order
    if not delimiters:
        delimiters = [',', '\t', ';']

    return delimiters


def _row_arity_consistency(
    filepath: Path,
    encoding: str,
    delimiter: str,
    quoting: int,
    sample: int = 500,
) -> float:
    """Fraction of sampled data rows whose RAW field count equals the header's.

    Parses the raw file with ``csv.reader`` under the SAME (delimiter, quoting)
    pandas used, so embedded/quoted delimiters are counted correctly. This is the
    one signal that catches field-shift: a row with an unquoted embedded delimiter
    carries an extra field, so its width != the header width, and the value under
    each later column belongs to its neighbour (a neighbouring-column contamination). pandas' ``on_bad_lines='warn'`` hides this by padding/
    truncating to header width — so the check must run on the raw stream.

    Returns 1.0 when it cannot assess (header-only, unreadable, single field) so
    a read error never penalises an otherwise-valid parse.
    """
    try:
        widths = []
        with open(filepath, 'r', encoding=encoding, newline='') as f:
            reader = csv.reader(f, delimiter=delimiter, quoting=quoting)
            for i, row in enumerate(reader):
                # skip wholly-blank lines (trailing newline artifacts)
                if not row or (len(row) == 1 and row[0] == ''):
                    continue
                widths.append(len(row))
                if len(widths) > sample:
                    break
        if len(widths) < 2:
            return 1.0
        header_w = widths[0]
        data = widths[1:]
        return sum(1 for w in data if w == header_w) / len(data)
    except Exception:
        return 1.0


def _detect_multi_space_format(filepath: Path) -> bool:
    """Detect if file uses multi-space (fixed-width) format.

    Returns True if the header line contains patterns of 2+ consecutive spaces
    that appear to be column separators (not tabs, commas, or semicolons).
    """
    import re

    try:
        with open(filepath, 'rb') as f:
            sample = f.read(4096).decode('utf-8', errors='ignore')
            # Strip BOM if present
            if sample.startswith('\ufeff'):
                sample = sample[1:]

            first_line = sample.split('\n')[0] if '\n' in sample else sample
            first_line = first_line.rstrip('\r')

            # Check for 2+ consecutive spaces pattern
            multi_space_matches = re.findall(r' {2,}', first_line)

            # Must have multiple instances of 2+ spaces
            if len(multi_space_matches) < 3:
                return False

            # Should NOT have significant tabs, commas, or semicolons
            tab_count = first_line.count('\t')
            comma_count = first_line.count(',')
            semicolon_count = first_line.count(';')

            # If we have more multi-space patterns than standard delimiters, it's likely fixed-width
            max_standard = max(tab_count, comma_count, semicolon_count)

            return len(multi_space_matches) > max_standard

    except Exception as e:
        logger.debug(f"Multi-space format detection failed: {e}")
        return False


def _try_multi_space_delimited(
    filepath: Path,
    min_columns: int,
    verbose: bool,
    pandas_kwargs: dict
) -> CSVReadResult:
    """Try parsing fixed-width files with multi-space delimiters.

    This handles files where columns are separated by 2+ consecutive spaces:
        Employee_ID     Name            Department
        8000001         John Smith      Engineering
        8000002         Jane Doe        Marketing

    Common in legacy fixed-width exports.

    Args:
        filepath: Path to CSV file
        min_columns: Minimum expected columns
        verbose: Enable verbose logging
        pandas_kwargs: Additional pandas arguments

    Returns:
        CSVReadResult (success=False if format doesn't match)
    """
    # First check if this looks like a multi-space format
    if not _detect_multi_space_format(filepath):
        return CSVReadResult(success=False, error="Not multi-space format")

    if verbose:
        logger.info("Detected multi-space (fixed-width) format, trying regex delimiter...")

    # Try multiple encodings
    encodings = ['utf-8-sig', 'utf-8', 'latin-1', 'cp1252']

    for encoding in encodings:
        try:
            # Use regex for 2+ spaces as delimiter
            df = pd.read_csv(
                filepath,
                sep=r'\s{2,}',  # 2 or more whitespace characters
                engine='python',  # Required for regex sep
                encoding=encoding,
                on_bad_lines='warn',
                **{k: v for k, v in pandas_kwargs.items()
                   if k not in ['sep', 'engine', 'encoding', 'on_bad_lines']}
            )

            # Validation
            if len(df.columns) < min_columns:
                if verbose:
                    logger.debug(f"  ✗ Only {len(df.columns)} columns (min={min_columns})")
                continue

            # Check for null/empty column names
            if df.columns.isnull().any() or (df.columns == '').any():
                if verbose:
                    logger.debug(f"  ✗ Null or empty column names")
                continue

            # Clean up column names
            df.columns = df.columns.str.strip()

            # Remove BOM from first column if present
            if df.columns[0].startswith('\ufeff'):
                first_col = df.columns[0].lstrip('\ufeff')
                df.rename(columns={df.columns[0]: first_col}, inplace=True)

            result = CSVReadResult(
                df=df,
                success=True,
                encoding=encoding,
                delimiter=r'\s{2,}',  # Regex pattern
                quoting=None,  # N/A for regex delimiter
                num_columns=len(df.columns),
                num_rows=len(df),
                strategy_index=-2  # Special strategy for multi-space
            )

            logger.info(f"✓ Multi-space delimited: {result.num_rows} rows × {result.num_columns} columns "
                       f"(encoding={result.encoding})")

            return result

        except Exception as e:
            if verbose:
                logger.debug(f"  ✗ Multi-space parse failed with {encoding}: {e}")
            continue

    return CSVReadResult(success=False, error="Multi-space format detected but parsing failed")


def _try_row_level_quoted_tsv(
    filepath: Path,
    min_columns: int,
    verbose: bool,
    pandas_kwargs: dict
) -> CSVReadResult:
    """Try parsing "row-level quoted TSV" format.

    This handles files where each entire row is wrapped in quotes with tabs inside:
        "col1\tcol2\tcol3",
        "val1\tval2\tval3",

    This format is common in some legacy/Excel exports.

    Args:
        filepath: Path to CSV file
        min_columns: Minimum expected columns
        verbose: Enable verbose logging
        pandas_kwargs: Additional pandas arguments

    Returns:
        CSVReadResult (success=False if format doesn't match)
    """
    try:
        # Read raw content with multiple encodings
        content = None
        detected_encoding = None

        for encoding in ['utf-8-sig', 'utf-8', 'latin-1', 'cp1252']:
            try:
                with open(filepath, 'r', encoding=encoding) as f:
                    content = f.read()
                detected_encoding = encoding
                break
            except UnicodeDecodeError:
                continue

        if content is None:
            return CSVReadResult(success=False, error="Could not decode file")

        # Normalize line endings
        content = content.replace('\r\n', '\n').replace('\r', '\n')
        lines = content.strip().split('\n')

        if len(lines) < 2:
            return CSVReadResult(success=False, error="Too few lines")

        # Check if first line matches pattern: starts with " and contains tabs
        first_line = lines[0]
        if not (first_line.startswith('"') and '\t' in first_line):
            return CSVReadResult(success=False, error="Not row-level quoted TSV format")

        if verbose:
            logger.info("Detected row-level quoted TSV format, preprocessing...")

        # Preprocess: strip outer quotes and trailing commas
        cleaned_lines = []
        for i, line in enumerate(lines):
            if not line.strip():
                continue

            # Strip outer quotes and trailing comma
            line = line.strip()
            if line.startswith('"'):
                line = line[1:]  # Remove leading quote
            if line.endswith('",'):
                line = line[:-2]  # Remove trailing ",
            elif line.endswith('"'):
                line = line[:-1]  # Remove trailing "

            # Handle escaped quotes within fields (Excel-style "" -> ")
            # But preserve double-quotes that are part of field values
            # This is tricky - we need to handle: ""Global Modern Services"," PLC""
            # Which should become: "Global Modern Services", PLC"
            # The doubled quotes at field boundaries need to be reduced

            cleaned_lines.append(line)

        # Join and parse as TSV
        cleaned_content = '\n'.join(cleaned_lines)

        # Parse with pandas
        df = pd.read_csv(
            io.StringIO(cleaned_content),
            sep='\t',
            quoting=csv.QUOTE_NONE,  # We already stripped outer quotes
            on_bad_lines='warn',
            **{k: v for k, v in pandas_kwargs.items()
               if k not in ['sep', 'quoting', 'encoding', 'on_bad_lines']}
        )

        # Validation
        if len(df.columns) < min_columns:
            return CSVReadResult(
                success=False,
                error=f"Only {len(df.columns)} columns (min={min_columns})"
            )

        # Clean up column names
        df.columns = df.columns.str.strip().str.strip('"').str.strip("'")

        # Clean up all string values: strip quotes from start/end
        for col in df.columns:
            if df[col].dtype == 'object':
                df[col] = df[col].apply(lambda x: _clean_quoted_value(x) if isinstance(x, str) else x)

        result = CSVReadResult(
            df=df,
            success=True,
            encoding=detected_encoding,
            delimiter='\t',
            quoting=csv.QUOTE_NONE,
            num_columns=len(df.columns),
            num_rows=len(df),
            strategy_index=-1  # Special strategy
        )

        logger.info(f"✓ Row-level quoted TSV: {result.num_rows} rows × {result.num_columns} columns "
                   f"(encoding={result.encoding})")

        return result

    except Exception as e:
        if verbose:
            logger.debug(f"Row-level quoted TSV detection failed: {e}")
        return CSVReadResult(success=False, error=str(e))


def _clean_quoted_value(value: str) -> str:
    """Clean up a value that may have leading/trailing quotes.

    Handles Excel-style escaped quotes ("") and field-level quoting.

    Excel CSV escaping rules:
    - "" inside a quoted field means a literal " character
    - ""text"" means the value is "text" (quotes included)

    Examples:
        '""Global Modern Services"," PLC""' -> '"Global Modern Services", PLC'
        'Regular' -> 'Regular' (unchanged)
        '""' -> '"' (single quote)
        '"simple"' -> 'simple' (simple quoting)
    """
    if not value:
        return value

    value = value.strip()

    # Handle Excel-style double-quote escaping: "" -> "
    # This is the standard CSV escaping for literal quote characters
    if '""' in value:
        value = value.replace('""', '"')

    # If value is now fully quoted (single quotes at edges), consider unquoting
    # This handles values like '"text"' -> 'text'
    if len(value) >= 2 and value.startswith('"') and value.endswith('"'):
        inner = value[1:-1]

        # DON'T unquote if inner contains quote-comma patterns
        # This indicates the quotes are part of the value (e.g., company names)
        # Pattern: "Company Name", LLC  -> the quotes ARE the value
        if '",' in inner or ',"' in inner:
            # Keep the value as-is (quotes are part of the name)
            pass
        # Only unquote if inner has NO quotes (simple quoted value)
        elif '"' not in inner:
            value = inner

    return value


def read_csv_simple(filepath: Union[str, Path], **kwargs) -> pd.DataFrame:
    """Simple wrapper that returns DataFrame or raises on failure.

    Convenience function for cases where you don't need result metadata.

    Args:
        filepath: Path to CSV file
        **kwargs: Arguments passed to read_csv_robust()

    Returns:
        DataFrame

    Raises:
        ValueError: If CSV parsing fails

    Example:
        >>> df = read_csv_simple("data.csv")
    """
    result = read_csv_robust(filepath, **kwargs)

    if not result.success:
        raise ValueError(f"Failed to parse CSV: {result.error}")

    return result.df


# ---------------------------------------------------------------------------
# CSV -> XLSX conversion (native .xlsx export)
# ---------------------------------------------------------------------------
# A comma/semicolon-delimited CSV can be misread by a locale-specific
# spreadsheet (e.g. French Excel splits on ';'), so a field that legitimately
# contains ';' or ',' spills into neighbouring columns on open. Writing a native
# .xlsx removes that ambiguity entirely: cell boundaries are explicit in the
# file, independent of the consumer's locale and with no text-to-columns step.
XLSX_MAX_ROWS = 1_048_576   # Excel/OOXML worksheet row limit (incl. header)
XLSX_MAX_COLS = 16_384      # Excel/OOXML worksheet column limit


def _safe_sheet_name(name: str) -> str:
    """Excel worksheet names: max 31 chars, and none of []:*?/\\ ."""
    cleaned = "".join(c for c in str(name) if c not in set('[]:*?/\\')).strip()
    return (cleaned or "data")[:31]


def to_xlsx(
    source: Union[str, Path, pd.DataFrame, CSVReadResult],
    xlsx_path: Optional[Union[str, Path]] = None,
    sheet_name: str = "data",
    min_columns: int = 1,
    verbose: bool = False,
    **read_kwargs,
) -> Path:
    """Load a CSV robustly and write a native ``.xlsx`` with correct cells.

    ``source`` may be a path to a CSV, an already-parsed :class:`CSVReadResult`,
    or a :class:`pandas.DataFrame`. When a CSV is read, the detected
    (encoding, delimiter, quoting) and the row-arity consistency are logged, so
    any embedded-delimiter contamination in the *source* is surfaced rather than
    silently propagated — the ``.xlsx`` faithfully mirrors the parsed cells.

    Args:
        source: CSV path, ``CSVReadResult``, or ``DataFrame``.
        xlsx_path: output path (default: the CSV path with a ``.xlsx`` suffix).
            Required when ``source`` is a ``DataFrame`` (no path to derive from).
        sheet_name: worksheet name (sanitised to Excel's 31-char / charset rules).
        min_columns: minimum expected columns when reading a CSV.
        verbose: verbose read logging.
        **read_kwargs: forwarded to :func:`read_csv_robust`.

    Returns:
        The output ``Path``.

    Raises:
        ValueError: on parse failure, a missing output path for a DataFrame, or
            a table exceeding the Excel worksheet limits.

    Example:
        >>> to_xlsx("per_group.csv")                      # -> per_group.xlsx
        >>> to_xlsx("data.csv", "out/data.xlsx", sheet_name="results")
    """
    if isinstance(source, pd.DataFrame):
        if xlsx_path is None:
            raise ValueError("xlsx_path is required when source is a DataFrame")
        df = source
    else:
        res = source if isinstance(source, CSVReadResult) else read_csv_robust(
            source, min_columns=min_columns, verbose=verbose, **read_kwargs)
        if not res.success:
            raise ValueError(f"robust CSV read failed: {res.error}")
        df = res.df
        logger.info(
            "to_xlsx: loaded %d×%d (encoding=%s, delimiter=%r, quoting=%s, arity=%.3f)",
            res.num_rows, res.num_columns, res.encoding, res.delimiter,
            res.quoting, res.row_arity_consistency)
        if res.row_arity_consistency < 0.999:
            logger.warning(
                "to_xlsx: source row-arity %.3f < 1.0 — possible embedded-delimiter "
                "contamination; the .xlsx mirrors the parsed cells as-is.",
                res.row_arity_consistency)
        if xlsx_path is None:
            xlsx_path = Path(source).with_suffix(".xlsx")

    xlsx_path = Path(xlsx_path)
    if len(df) + 1 > XLSX_MAX_ROWS:
        raise ValueError(
            f"{len(df)} data rows exceed the Excel worksheet limit "
            f"({XLSX_MAX_ROWS} incl. header); split into sheets or keep the CSV.")
    if len(df.columns) > XLSX_MAX_COLS:
        raise ValueError(
            f"{len(df.columns)} columns exceed the Excel limit ({XLSX_MAX_COLS}).")

    xlsx_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_excel(xlsx_path, index=False,
                sheet_name=_safe_sheet_name(sheet_name), engine="openpyxl")
    logger.info("to_xlsx: wrote %s (%.1f MB)", xlsx_path, xlsx_path.stat().st_size / 1e6)
    return xlsx_path


def _main(argv=None) -> int:
    """CLI: robustly convert one or more CSVs to native ``.xlsx``.

        python -m tools.robust_csv per_group.csv per_scenario.csv
        python -m tools.robust_csv data/*.csv --out-dir xlsx/
    """
    import argparse
    ap = argparse.ArgumentParser(
        description="Robust CSV reader + CSV->XLSX converter (SFPPy).")
    ap.add_argument("csv", nargs="+", help="CSV file(s) to convert to .xlsx")
    ap.add_argument("-o", "--out-dir", default=None,
                    help="output directory (default: alongside each CSV)")
    ap.add_argument("--sheet", default="data", help="worksheet name")
    ap.add_argument("-v", "--verbose", action="store_true")
    a = ap.parse_args(argv)
    logging.basicConfig(
        level=logging.INFO if a.verbose else logging.WARNING,
        format="%(levelname)s %(message)s")
    rc = 0
    for c in a.csv:
        try:
            out = (Path(a.out_dir) / (Path(c).stem + ".xlsx")) if a.out_dir else None
            path = to_xlsx(c, out, sheet_name=a.sheet, verbose=a.verbose)
            print(f"[ok] {c} -> {path}")
        except Exception as e:
            print(f"[fail] {c}: {e}", file=sys.stderr)
            rc = 1
    return rc


if __name__ == "__main__":
    sys.exit(_main())
