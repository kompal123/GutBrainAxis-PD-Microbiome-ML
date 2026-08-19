from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = PROJECT_ROOT / "data"
RAW_DIR = DATA_DIR / "raw"
PROCESSED_DIR = DATA_DIR / "processed"
RESULTS_DIR = PROJECT_ROOT / "results"
FIGURE_DIR = RESULTS_DIR / "figures"
TABLE_DIR = RESULTS_DIR / "tables"
DOCS_DIR = PROJECT_ROOT / "docs"
DOC_FIGURE_DIR = DOCS_DIR / "figures"
DOC_TABLE_DIR = DOCS_DIR / "tables"
MODEL_DIR = PROJECT_ROOT / "models"

SOURCE_WORKBOOK = RAW_DIR / "Source_Data_24Oct2022.xlsx"
ZENODO_SOURCE_URL = "https://zenodo.org/api/records/7246185/files/Source_Data_24Oct2022.xlsx/content"
RANDOM_STATE = 42
