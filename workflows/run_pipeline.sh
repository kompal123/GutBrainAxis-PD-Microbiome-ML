#!/usr/bin/env bash
set -euo pipefail

export PYTHONPATH=src

python src/00_prepare_data.py
python src/01_differential_taxa.py
python src/02_train_ml.py
python src/03_gut_brain_story.py
