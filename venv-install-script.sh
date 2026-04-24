#!/bin/bash
#find -name "*.py" -not -path "./.*" -exec sed -i '1i #!.venv/bin/python' {} ';'
python3.12 -m venv .venv
source .venv/bin/activate
python3 -m pip install --upgrade pip
python3 -m pip install -r requirements.txt
find -name "*.py" -not -path "./.*" -exec chmod u+x {} ';'