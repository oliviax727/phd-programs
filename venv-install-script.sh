#!/bin/bash

set -euo pipefail

rm -rf .venv

python3 -m venv .venv

# shellcheck disable=SC1091
source .venv/bin/activate
.venv/bin/python -m pip install --upgrade pip
.venv/bin/python -m pip install -r requirements.txt
.venv/bin/python -m ipykernel install --user --name=.venv
