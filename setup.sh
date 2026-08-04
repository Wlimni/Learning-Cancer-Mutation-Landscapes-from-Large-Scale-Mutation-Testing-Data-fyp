#!/usr/bin/env bash
# Script Name : setup.sh
# Purpose     : One-time environment setup for the FYP notebooks -- creates a
#               Python venv, installs packages from requirements.txt, and
#               registers a Jupyter kernel so the notebooks have something to run on.
# Inputs      : requirements.txt (in the same directory as this script)
# Outputs     : .venv/ directory; a Jupyter kernel named "fyp-genie"
# Usage       : cd fyp && ./setup.sh
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"

PYTHON_BIN="${PYTHON_BIN:-python3.12}"
command -v "$PYTHON_BIN" >/dev/null 2>&1 || PYTHON_BIN=python3

if [ ! -d .venv ]; then
  echo "Creating virtual environment with $PYTHON_BIN..."
  "$PYTHON_BIN" -m venv .venv
fi

source .venv/bin/activate
pip install --upgrade pip -q
pip install -r requirements.txt -q

python -m ipykernel install --user --name fyp-genie --display-name "Python (fyp-genie)"

echo ""
echo "Setup complete. Open a notebook and select the 'Python (fyp-genie)' kernel."
