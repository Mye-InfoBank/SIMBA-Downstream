#!/bin/bash

sudo apt-get update -y && sudo apt-get install -y python3-pip python3-dev

python3 -m pip install -r requirements.txt

Rscript requirements.R