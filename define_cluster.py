import os
import sys
import warnings
from pathlib import Path

import pandas as pd
from matplotlib import pyplot as plt

from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor
import time
 
sys.path.append("..")


from cluster_definition.cluster_definer import ClasterDefiner
from datareading.omicron import CreateOmicronHeaderDF, ShowImage
from utils.utils import get_logger

warnings.simplefilter("ignore")

# plt.rcParams["font.family"] = "Times New Roman"

DATA_COLUMNS = ["year", "month", "day"]
CLADE_COLUMNS = ["clade", "head2"]

JAPAN_HEADER_COLUMNS = ["head", "ID", "date"]
DELTA_HEADER_COLUMNS = ["head", "collection date", "submission date"]

ALL_DATA_HEADER_COLUMNS = [
    "head",
    "id",
    "continent",
    "country",
    "city",
    "host",
    "clade_head",
    "collection date",
]

MAX_CLUSTER_NUM = 300
DATA_DIR = "omicron2201"

LOGGER = get_logger(f"{DATA_DIR}/{MAX_CLUSTER_NUM}")
LOGGER.info(f"{DATA_DIR} CLUSTER_MAC: {MAX_CLUSTER_NUM}")

START_TIME = time.time()

# 行数
# pd.set_option("display.max_rows")
# pd.set_option("display.max_columns", 50)

def func_a(continent_df, blsom_weight, con):
    LOGGER.info(con)
    definer = ClasterDefiner(LOGGER)
    
    # value_countsだと同じ件数のものがランダムに並ぶのでindexと件数で並び替え
    agg_coordinate = continent_df["x y"].value_counts().reset_index().sort_values(by=["x y", "index"], ascending=[False, True])
    agg_coordinate = agg_coordinate.set_index('index')

    n = MAX_CLUSTER_NUM if len(agg_coordinate.index) > MAX_CLUSTER_NUM else len(agg_coordinate.index)
    _, top_n_coordinate = definer.data_from_top_n(agg_coordinate, n=n)
    LOGGER.info(len(top_n_coordinate))

    cluster_df = definer.df2cluster(continent_df, blsom_weight, top_n_coordinate, sigma=0.5)
    # definer.cluster_plot(cluster_df, (X, Y), output_file=f"{DATA_DIR}/output/{MAX_CLUSTER_NUM}/images/cluster/{con}/all_cluster.png")

    for i in cluster_df["labels"].unique():
        by3d = cluster_df[(cluster_df["labels"] == i)]
        ShowImage(X, Y).blsom_plot_3d(by3d, f"{DATA_DIR}/output/{MAX_CLUSTER_NUM}/images/cluster/{con}/3D/{i}.png")

    cluster_df.to_csv(f"{DATA_DIR}/output/{MAX_CLUSTER_NUM}/ANALYTICS_CLUSTER/{con}.tsv", sep="\t")
    cluster_df[["labels", "Accession ID"]].to_csv(f"{DATA_DIR}/output/{MAX_CLUSTER_NUM}/CLUSTER_AND_ID/{con}.tsv", "\t", index=False)

    END_TIME = time.time()
    delta = END_TIME - START_TIME
    LOGGER.info(f"{con} Done! 処理時間:{round(delta,3)}s")

    return


if __name__ in "__main__":
    file = f"{DATA_DIR}/ReceivedData/dataToFurukawa/Omicron.meta.id.fas.N1.id.meta"
    filename = f"{DATA_DIR}/ReceivedData/SOM/out_all"

    createdf = CreateOmicronHeaderDF(
        out_filename=filename, meta_filename=file, weight_file=f"{DATA_DIR}/ReceivedData/SOM/weight.100",
    )
    createdf.create_df()

    delta_blsom = createdf.delta_blsom
    X = createdf.X
    Y = createdf.Y
    blsom_weight = createdf.blsom_weight

    if DATA_DIR == "omicron2201":
        # ShowImage(delta_blsom, X, Y).show_color_by_continent()
        new_lin_file = f"{DATA_DIR}/ReceivedData/dataToFurukawa/Omicron.meta.id.fas.N1.lin"

        new_delta_header = pd.read_csv(new_lin_file, sep="\t")
        new_delta_header = new_delta_header.rename(columns={"lineage": "new_lineage"})
        delta_blsom = delta_blsom.merge(
            new_delta_header[["taxon", "new_lineage"]], right_on="taxon", left_on="Virus name", how="left",
        )

    delta_blsom.to_csv(f"{DATA_DIR}/output/ALL.tsv", sep="\t")

    # coordinates = delta_blsom["x y"].value_counts()

    # continent = ["Asia"]
    # continent = ["North_America", "Oceania", "Africa", "Europe"]
    # continent = ["Asia", "Oceania", "Africa"]
    continent = ["Asia", "North_America", "Oceania", "Africa", "Europe"]

    END_TIME = time.time()
    delta = END_TIME - START_TIME
    LOGGER.info(f"Read OK! 処理時間:{round(delta,3)}s")

    if not Path(f"{DATA_DIR}/output/{MAX_CLUSTER_NUM}").exists():
        os.mkdir(Path(f"{DATA_DIR}/output/{MAX_CLUSTER_NUM}"))
        os.mkdir(Path(f"{DATA_DIR}/output/{MAX_CLUSTER_NUM}/images/"))
        os.mkdir(Path(f"{DATA_DIR}/output/{MAX_CLUSTER_NUM}/images/cluster"))
        os.mkdir(Path(f"{DATA_DIR}/output/{MAX_CLUSTER_NUM}/ANALYTICS_CLUSTER"))
        os.mkdir(Path(f"{DATA_DIR}/output/{MAX_CLUSTER_NUM}/CLUSTER_AND_ID"))

    # with ThreadPoolExecutor(max_workers=8) as executor:            
    # with ProcessPoolExecutor(max_workers=6) as executor:
    for con in continent:
        if not Path(f"{DATA_DIR}/output/{MAX_CLUSTER_NUM}/images/cluster/{con}").exists():
            os.mkdir(Path(f"{DATA_DIR}/output/{MAX_CLUSTER_NUM}/images/cluster/{con}"))
            os.mkdir(Path(f"{DATA_DIR}/output/{MAX_CLUSTER_NUM}/images/cluster/{con}/3D"))
        continent_df = delta_blsom[delta_blsom["continent"] == con]
        # executor.submit(func_a, continent_df, blsom_weight, con)
        func_a(continent_df, blsom_weight, con)

    END_TIME = time.time()
    delta = END_TIME - START_TIME
    LOGGER.info(f"All Done! 処理時間:{round(delta,3)}s")
