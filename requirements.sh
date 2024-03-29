#!/bin/bash

apt-get update -y && apt-get install -y python3-pip python3-dev && \
    python3 -m pip install --upgrade pip && \
    python3 -m pip install --user -r requirements.txt && \
    Rscript requirements.R