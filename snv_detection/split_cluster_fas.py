import sys

import pandas as pd
from tqdm import tqdm
from pathlib import Path
import os
from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor
import time

sys.path.append("..")

from utils.utils import get_logger

DATA_DIR = sys.argv[1]
MAX_CLUSTER_NUM = sys.argv[2]
CONTINENT = sys.argv[3]

HEADER = f"../{DATA_DIR}/ReceivedData/dataToFurukawa/Omicron.meta.id.fas.N1.id.meta"

FASTA_FILE = f"../{DATA_DIR}/ReceivedData/dataToFurukawa/Omicron.meta.id.fas.N1"

LOGGER = get_logger(f"../{DATA_DIR}/split_fas")
LOGGER.info(f"{DATA_DIR} CLUSTER_MAC: {MAX_CLUSTER_NUM} CONTINENT: {CONTINENT}")

START_TIME = time.time()


filepath = f"../{DATA_DIR}/output/{MAX_CLUSTER_NUM}/CLUSTER_AND_ID/{CONTINENT}.tsv"
outputdir = f"../{DATA_DIR}/SNV/{MAX_CLUSTER_NUM}/CLUSTER/{CONTINENT}"

cluster_and_id = pd.read_csv(filepath, sep="\t")
header = pd.read_csv(HEADER, sep="\t")
# print(header)

cluster_and_id = cluster_and_id.merge(header, left_on="Accession ID", right_on="Accession ID")
# print(cluster_and_id)
cluster_id = []
for i in range(0, max(cluster_and_id["labels"].values) + 1):
    cluster_id.append(cluster_and_id[cluster_and_id["labels"] == i]["Virus name"].values)

# print(cluster_id)
if not Path(f"../{DATA_DIR}/SNV/{MAX_CLUSTER_NUM}").exists():
    os.mkdir(Path(f"../{DATA_DIR}/SNV/{MAX_CLUSTER_NUM}"))
    os.mkdir(Path(f"../{DATA_DIR}/SNV/{MAX_CLUSTER_NUM}/CLUSTER/"))

if not Path(outputdir).exists():
    # continentまでのpath
    os.mkdir(Path(outputdir))

# def create_sep_list(seq):
    
with open(FASTA_FILE, "r") as r:
    # for i in tqdm(r.readlines()):
    #     # print(i)
    #     pass
    seqs = [f">{i}" for i in r.read().split(">") if i != ""]
    # seqs = [f">{i}" for i in tqdm(r.read().split(">")) if i != ""]

# def write_seq(seq, ids, w):
#     for id in ids:
#         if id in seq[:100]:
#             w.write(seq)
#     #         return True
#     # return False
#     return

seqs_h = [i.split("\n")[0] for i in seqs]

for i, ids in enumerate(cluster_id):
    with open(f"{outputdir}/CLUSTER{i}_ID.fas", "w+") as w:
        for i in range(len(seqs_h)):
        # for i in tqdm(range(len(seqs_h))):
            for id in ids:
                if id in seqs_h[i]:
                    w.write(seqs[i])

END_TIME = time.time()
delta = END_TIME - START_TIME
LOGGER.info(f"{CONTINENT} Done! 処理時間:{round(delta,3)}s\n")
# with ProcessPoolExecutor(max_workers=6) as executor:
#     for i, ids in enumerate(cluster_id):
#         with open(f"{outputdir}/CLUSTER{i}_ID.fas", "w+") as w:
#             for seq in tqdm(seqs):
#                 executor.submit(write_seq, seq, ids, w)
