from os.path import join
import platform

# Check if this is running on O2
IS_O2 = (platform.system() == "Linux")

# Directory / file constants
SRC_DIR = "src"
DATA_DIR = ("data" if not IS_O2 else "/n/data1/hms/dbmi/gehlenborg/lab/hd-vis-star")
RAW_DIR = join(DATA_DIR, "raw")
INTERMEDIATE_DIR = join(DATA_DIR, "intermediate")
PROCESSED_DIR = join(DATA_DIR, "processed")
