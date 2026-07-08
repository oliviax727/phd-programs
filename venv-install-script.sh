#!/bin/bash
#find -name "*.py" -not -path "./.*" -exec sed -i '1i #!.venv/bin/python' {} ';'
python3.12 -m venv .venv
source .venv/bin/activate
.venv/bin/python -m pip install --upgrade pip
.venv/bin/python -m pip install -r requirements.txt
.venv/bin/python -m pip install "git+https://github.com/OxfordSKA/OSKAR.git@master#egg=oskarpy&subdirectory=python"
ipython kernel install --user --name=.venv
#find -name "*.py" -not -path "./.*" -exec chmod u+x {} ';'