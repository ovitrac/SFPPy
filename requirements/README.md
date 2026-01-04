# SFPPy Requirements 🍏⏩🍎

Modular dependency files for different installation profiles.

## Installation Profiles

| Profile | Command | Description |
|---------|---------|-------------|
| **Core** | `pip install -r requirements/core.txt` | Minimal (patankar solver) |
| **Studio** | `pip install -r requirements/core.txt -r requirements/studio.txt` | Web GUI |
| **Survey** | `pip install -r requirements/core.txt -r requirements/survey.txt` | Batch processing |
| **Full** | `pip install -r requirements/full.txt` | Everything |
| **RAG** | `pip install -r requirements/core.txt -r requirements/rag.txt` | AI/LLM queries |
| **Dev** | `pip install -r requirements/core.txt -r requirements/dev.txt` | Development tools |

## Using pyproject.toml (Recommended)

Modern installation with optional dependencies:

```bash
# Core only
pip install -e .

# With Studio
pip install -e ".[studio]"

# With Survey
pip install -e ".[survey]"

# Full (Studio + Survey + dev)
pip install -e ".[full]"

# Everything including RAG
pip install -e ".[complete]"
```

## Dependency Files

### core.txt
Minimal scientific stack:
- numpy, scipy, matplotlib
- pandas, openpyxl
- pillow

### studio.txt
Web GUI dependencies:
- fastapi, uvicorn
- pydantic, jinja2
- requests, pyyaml

### survey.txt
Batch processing dependencies:
- Same as studio
- odfpy (ODS spreadsheets)
- aiofiles (async I/O)

### rag.txt
AI/LLM dependencies:
- llama-index
- sentence-transformers
- chromadb
- Requires Ollama: `ollama pull mistral`

### dev.txt
Development tools:
- pytest, pytest-cov
- black, isort, mypy
- flake8, bandit

### full.txt
Combines core + studio + survey (excludes RAG due to size)

## Conda Environment

For conda users:

```bash
# Create environment
conda create -n sfppy python=3.10 numpy scipy matplotlib pandas -y
conda activate sfppy

# Install remaining deps via pip
pip install -r requirements/studio.txt
```

Or use the environment.yml:

```bash
conda env create -f environment.yml
conda activate sfppy
```

## Version Constraints

All dependencies use minimum version constraints (`>=`) to allow flexibility while ensuring compatibility.

## Updating Dependencies

To upgrade all packages:

```bash
pip install --upgrade -r requirements/full.txt
```
