#!/usr/bin/env python3
"""
Synchronize thumbs_tiny/ with thumbs/ by downloading missing tiny thumbnails from PubChem.

This script:
1. Scans thumbs/ for all cached CIDs
2. Checks which CIDs are missing from thumbs_tiny/
3. Downloads missing tiny thumbnails from PubChem with rate limiting
4. Saves progress to allow resumption if interrupted

PubChem Image API:
- Tiny (100x100): https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{CID}/PNG?image_size=100x100

Usage:
    python sync_thumbs_tiny.py [--dry-run] [--limit N] [--delay SECONDS]

Author: Olivier Vitrac
Version: 1.0 (2026-01-09)
"""

import os
import sys
import time
import argparse
import requests
from pathlib import Path
from datetime import datetime

# Add project root for imports
SCRIPT_DIR = Path(__file__).parent
PROJECT_ROOT = SCRIPT_DIR.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

# Paths - script now lives in cache.PubChem/
CACHE_DIR = SCRIPT_DIR
THUMBS_DIR = CACHE_DIR / "thumbs"
THUMBS_TINY_DIR = CACHE_DIR / "thumbs_tiny"
PROGRESS_FILE = THUMBS_TINY_DIR / ".sync_progress.txt"

# PubChem API - use imgsrv.fcgi endpoint (more reliable than PUG REST)
# t=s for small (100x100), t=l for large
PUBCHEM_IMAGE_URL = "https://pubchem.ncbi.nlm.nih.gov/image/imgsrv.fcgi"
DEFAULT_DELAY = 0.17  # ~6 requests/second (PubChem allows up to 5/sec, but imgsrv is fast)
DEFAULT_TIMEOUT = 15  # seconds
MAX_RETRIES = 3
RETRY_BACKOFF = 2.0  # seconds, doubles each retry


def get_cached_cids(thumbs_dir: Path) -> set:
    """Get set of CIDs from existing thumbnail files."""
    cids = set()
    for f in thumbs_dir.glob("*.png"):
        try:
            cid = int(f.stem)
            cids.add(cid)
        except ValueError:
            continue
    return cids


def load_progress(progress_file: Path) -> set:
    """Load set of already-processed CIDs from progress file."""
    if not progress_file.exists():
        return set()
    try:
        with open(progress_file, 'r') as f:
            return {int(line.strip()) for line in f if line.strip().isdigit()}
    except Exception:
        return set()


def save_progress(progress_file: Path, processed: set):
    """Save processed CIDs to progress file."""
    with open(progress_file, 'w') as f:
        for cid in sorted(processed):
            f.write(f"{cid}\n")


def download_tiny_thumbnail(cid: int, output_path: Path, timeout: int = DEFAULT_TIMEOUT,
                            max_retries: int = MAX_RETRIES, quiet: bool = False) -> bool:
    """
    Download tiny thumbnail from PubChem with retry logic.

    Uses the imgsrv.fcgi endpoint with t=s for small (100x100) images.

    Returns True if successful, False otherwise.
    """
    params = {"cid": cid, "t": "s"}  # t=s for small thumbnail

    for attempt in range(max_retries):
        try:
            response = requests.get(PUBCHEM_IMAGE_URL, params=params, timeout=timeout)

            if response.status_code == 200:
                # Verify it's actually a PNG
                if response.content[:8] == b'\x89PNG\r\n\x1a\n':
                    with open(output_path, 'wb') as f:
                        f.write(response.content)
                    return True
                else:
                    if not quiet:
                        print(f"[WARN] CID {cid}: Response is not a PNG image")
                    return False

            elif response.status_code == 404:
                if not quiet:
                    print(f"[SKIP] CID {cid}: Not found on PubChem")
                return False

            elif response.status_code in (503, 429, 500):
                # Retriable errors - back off and retry
                if attempt < max_retries - 1:
                    backoff = RETRY_BACKOFF * (2 ** attempt)
                    if not quiet:
                        print(f"[RETRY] CID {cid}: HTTP {response.status_code}, waiting {backoff:.1f}s...")
                    time.sleep(backoff)
                    continue
                else:
                    if not quiet:
                        print(f"[FAIL] CID {cid}: HTTP {response.status_code} after {max_retries} retries")
                    return False
            else:
                if not quiet:
                    print(f"[FAIL] CID {cid}: HTTP {response.status_code}")
                return False

        except requests.Timeout:
            if attempt < max_retries - 1:
                backoff = RETRY_BACKOFF * (2 ** attempt)
                if not quiet:
                    print(f"[RETRY] CID {cid}: Timeout, waiting {backoff:.1f}s...")
                time.sleep(backoff)
                continue
            else:
                if not quiet:
                    print(f"[TIMEOUT] CID {cid}: Request timed out after {max_retries} retries")
                return False

        except requests.RequestException as e:
            if not quiet:
                print(f"[ERROR] CID {cid}: {e}")
            return False

    return False


def main():
    parser = argparse.ArgumentParser(
        description="Sync thumbs_tiny/ with thumbs/ by downloading from PubChem"
    )
    parser.add_argument(
        "--dry-run", "-n",
        action="store_true",
        help="Show what would be downloaded without actually downloading"
    )
    parser.add_argument(
        "--limit", "-l",
        type=int,
        default=0,
        help="Limit number of downloads (0 = no limit)"
    )
    parser.add_argument(
        "--delay", "-d",
        type=float,
        default=DEFAULT_DELAY,
        help=f"Delay between requests in seconds (default: {DEFAULT_DELAY})"
    )
    parser.add_argument(
        "--reset",
        action="store_true",
        help="Reset progress and re-download all missing thumbnails"
    )
    parser.add_argument(
        "--quiet", "-q",
        action="store_true",
        help="Suppress per-file output"
    )
    args = parser.parse_args()

    # Ensure directories exist
    THUMBS_TINY_DIR.mkdir(parents=True, exist_ok=True)

    # Get CIDs from thumbs/
    print(f"Scanning {THUMBS_DIR}...")
    source_cids = get_cached_cids(THUMBS_DIR)
    print(f"  Found {len(source_cids)} thumbnails in thumbs/")

    # Get CIDs already in thumbs_tiny/
    existing_tiny = get_cached_cids(THUMBS_TINY_DIR)
    print(f"  Found {len(existing_tiny)} thumbnails in thumbs_tiny/")

    # Load progress (CIDs we've already attempted)
    if args.reset and PROGRESS_FILE.exists():
        PROGRESS_FILE.unlink()
        print("  Progress reset")

    processed = load_progress(PROGRESS_FILE)
    print(f"  Progress file: {len(processed)} CIDs already processed")

    # Determine what needs to be downloaded
    # Skip CIDs that are already in thumbs_tiny/ OR have been processed (failed attempts)
    already_done = existing_tiny | processed
    to_download = sorted(source_cids - already_done)

    print(f"\nTo download: {len(to_download)} thumbnails")

    if args.limit > 0:
        to_download = to_download[:args.limit]
        print(f"  Limited to: {len(to_download)} thumbnails")

    if not to_download:
        print("\nNothing to download. All thumbnails are synchronized!")
        return 0

    if args.dry_run:
        print("\n[DRY RUN] Would download:")
        for cid in to_download[:20]:
            print(f"  - {cid}.png")
        if len(to_download) > 20:
            print(f"  ... and {len(to_download) - 20} more")
        return 0

    # Download loop
    print(f"\nStarting download (delay: {args.delay}s between requests)...")
    print(f"Estimated time: {len(to_download) * args.delay / 60:.1f} minutes\n")

    success_count = 0
    fail_count = 0
    start_time = time.time()

    for i, cid in enumerate(to_download, 1):
        output_path = THUMBS_TINY_DIR / f"{cid}.png"

        if not args.quiet:
            print(f"[{i}/{len(to_download)}] Downloading CID {cid}...", end=" ", flush=True)

        success = download_tiny_thumbnail(cid, output_path, quiet=args.quiet)

        if success:
            success_count += 1
            if not args.quiet:
                print("OK")
        else:
            fail_count += 1
            if not args.quiet:
                print("")  # newline after status messages

        # Mark as processed (whether success or failure)
        processed.add(cid)

        # Save progress every 50 downloads
        if i % 50 == 0:
            save_progress(PROGRESS_FILE, processed)
            elapsed = time.time() - start_time
            rate = i / elapsed if elapsed > 0 else 0
            remaining = (len(to_download) - i) / rate if rate > 0 else 0
            print(f"  [Progress] {i}/{len(to_download)} ({success_count} OK, {fail_count} failed) "
                  f"- ETA: {remaining/60:.1f} min")

        # Rate limiting
        if i < len(to_download):
            time.sleep(args.delay)

    # Final save
    save_progress(PROGRESS_FILE, processed)

    # Summary
    elapsed = time.time() - start_time
    print(f"\n{'='*60}")
    print(f"Sync complete!")
    print(f"  Downloaded: {success_count}")
    print(f"  Failed: {fail_count}")
    print(f"  Total time: {elapsed/60:.1f} minutes")
    print(f"  Rate: {len(to_download)/elapsed:.2f} downloads/second")

    # Final counts
    final_tiny = get_cached_cids(THUMBS_TINY_DIR)
    print(f"\nFinal status:")
    print(f"  thumbs/: {len(source_cids)} files")
    print(f"  thumbs_tiny/: {len(final_tiny)} files")
    print(f"  Coverage: {len(final_tiny)/len(source_cids)*100:.1f}%")

    return 0 if fail_count == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
