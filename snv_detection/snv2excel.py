import numpy as np
import pandas as pd
from snv_detector import SnvDetector
from pathlib import Path
import os

continents = [
    "Africa",
    "Europe",
    "North_America",
    "Oceania",
    "Asia",
]

DATA_DIR = "omicron2201"
MAX_CLUSTER_NUM = 10

if not Path(f"../{DATA_DIR}/output/{MAX_CLUSTER_NUM}/SNV_EXCEL").exists():
    os.mkdir(Path(f"../{DATA_DIR}/output/{MAX_CLUSTER_NUM}/SNV_EXCEL"))

for i in continents:
    # df = pd.read_csv(f"../CLUSTER_AND_ID/{i}.tsv", sep="\t")
    # cluster_threshhold = np.mean(df.groupby("labels").size())

    # cluster_threshhold = 50
    cluster_threshhold = 0
    d = SnvDetector(
        dir_path=f"../{DATA_DIR}/SNV/{MAX_CLUSTER_NUM}/CLUSTER/{i}",
        output_file=f"../{DATA_DIR}/output/{MAX_CLUSTER_NUM}/SNV_EXCEL/SNP_BASE_{i}_60_.xlsx",
        cluster_threshhold=cluster_threshhold)

    d.extract_snv(
        # dir_path=f"CLUSTER/{i}",
        # output_file=f"TEST_MEAN/SNP_BASE_{i}_60.xlsx",
        # threshold=60,
        threshhold=60,
    )