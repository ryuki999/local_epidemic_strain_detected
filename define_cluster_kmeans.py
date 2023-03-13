import os
import warnings
from pathlib import Path

import pandas as pd
import time
from pathlib import Path
import copy
import pandas as pd
import random
import warnings

warnings.simplefilter('ignore')

from utils.utils import get_logger

from cluster_definition.cluster_definer import KmeansClusterDefiner
from datareading.omicron import CreateOmicronHeaderDF, Show3DImage, Show2DImage


MAX_CLUSTER_NUM = 10
DATA_DIR = "omicron2201"

LOGGER = get_logger(f"{DATA_DIR}/{MAX_CLUSTER_NUM}")
LOGGER.info(f"{DATA_DIR} CLUSTER_MAC: {MAX_CLUSTER_NUM}")

START_TIME = time.time()


from color import Color
import matplotlib.colors as mcolors
from PIL import ImageColor

random.seed(0)
colors = copy.copy(mcolors.CSS4_COLORS)
del colors["black"], colors["white"]

colors = [Color(color_name=i, hex=v, rgb=ImageColor.getcolor(v, "RGB")) for i, v in colors.items()]
colors = random.sample(colors, k=len(colors))


def plot(shower, df, classes, target, x, y, output_file):
    u_class = sorted(df["labels"].unique())
    rgb_color = {}
    for i, c in enumerate(u_class):
        colors[i].target = c
        rgb_color[c] = colors[c]
    s = shower(x, y, rgb_color, df)
    s.plot(title=len(classes), target=target, output_file=output_file)

if __name__ in "__main__":
    COMPARE_PASSED = True
    DATA_DIR = "omicron2201"
    # file = f"{DATA_DIR}/ReceivedData/dataToFurukawa/Omicron.meta.id.fas.N1.id.meta"
    # filename = f"{DATA_DIR}/ReceivedData/SOM/out_each_region"
    IMAGE_OUTPUT = f"{DATA_DIR}/output/images/"

    file = f"{DATA_DIR}/ReceivedData/dataToFurukawa/Omicron.meta.id.fas.N1.id.meta"
    filename = f"{DATA_DIR}/ReceivedData/SOM/out_all"

    createdf = CreateOmicronHeaderDF(
        out_filename=filename,
        meta_filename=file,
        weight_file=f"{DATA_DIR}/ReceivedData/SOM/weight.100",
    )
    createdf.create_df()

    delta_blsom = createdf.delta_blsom
    X = createdf.X
    Y = createdf.Y

    if COMPARE_PASSED:
        new_lin_file = f"{DATA_DIR}/ReceivedData/dataToFurukawa/Omicron.meta.id.fas.N1.lin"

        new_delta_header = pd.read_csv(new_lin_file, sep="\t")
        new_delta_header = new_delta_header.rename(columns={"lineage": "new_lineage"})

        delta_blsom = delta_blsom.merge(
            new_delta_header[["taxon", "new_lineage"]],
            right_on="taxon",
            left_on="Virus name",
            how="left",
        )

    delta_blsom.to_csv(f"{DATA_DIR}/ALL.tsv", sep="\t")

    delta_blsom = delta_blsom.merge(createdf.blsom_weight, right_index=True, left_on=["x y"])

    END_TIME = time.time()
    delta = END_TIME - START_TIME
    LOGGER.info(f"Read OK! 処理時間:{round(delta,3)}s")

    if not Path(f"{DATA_DIR}/output/{MAX_CLUSTER_NUM}").exists():
        os.mkdir(Path(f"{DATA_DIR}/output/{MAX_CLUSTER_NUM}"))
        os.mkdir(Path(f"{DATA_DIR}/output/{MAX_CLUSTER_NUM}/images/"))
        os.mkdir(Path(f"{DATA_DIR}/output/{MAX_CLUSTER_NUM}/images/cluster"))
        os.mkdir(Path(f"{DATA_DIR}/output/{MAX_CLUSTER_NUM}/ANALYTICS_CLUSTER"))
        os.mkdir(Path(f"{DATA_DIR}/output/{MAX_CLUSTER_NUM}/CLUSTER_AND_ID"))


    # 全部の重み
    definer = KmeansClusterDefiner(n_clusters=MAX_CLUSTER_NUM, cluster_range=5)
    w_cluster = definer.df2cluster(delta_blsom, createdf.blsom_weight)

    classes = sorted(w_cluster["labels"].unique())
    plot(shower=Show2DImage, df=w_cluster, classes=classes, target="labels", x=X, y=Y, output_file=f"{DATA_DIR}/output/{MAX_CLUSTER_NUM}/images/cluster/Cluster2D.png")

    continents = ["Europe", "North_America", "Asia", "Oceania"]
    for con in continents:
        if not Path(f"{DATA_DIR}/output/{MAX_CLUSTER_NUM}/images/cluster/{con}").exists():
            os.mkdir(Path(f"{DATA_DIR}/output/{MAX_CLUSTER_NUM}/images/cluster/{con}"))
            os.mkdir(Path(f"{DATA_DIR}/output/{MAX_CLUSTER_NUM}/images/cluster/{con}/2D"))
            os.mkdir(Path(f"{DATA_DIR}/output/{MAX_CLUSTER_NUM}/images/cluster/{con}/3D"))

        df = delta_blsom[delta_blsom["continent"] == con]
        definer = KmeansClusterDefiner(n_clusters=MAX_CLUSTER_NUM, cluster_range=3)
        w_cluster = definer.df2cluster(df, createdf.blsom_weight)
        plot(shower=Show2DImage, df=w_cluster, classes=classes, target="labels", x=X, y=Y, output_file=f"{DATA_DIR}/output/{MAX_CLUSTER_NUM}/images/cluster/{con}/Cluster2D_origin.png")

        w_cluster_5 = definer.extract_certain_range()
        print(sum(w_cluster_5["labels"].value_counts()))
        plot(shower=Show2DImage, df=w_cluster_5, classes=classes, target="labels", x=X, y=Y, output_file=f"{DATA_DIR}/output/{MAX_CLUSTER_NUM}/images/cluster/{con}/Cluster2D.png")


        for i in w_cluster_5["labels"].unique():
            by3d = w_cluster_5[w_cluster_5["labels"] == i]
            plot(shower=Show3DImage, df=by3d, classes=classes, target="labels", x=X, y=Y, output_file=f"{DATA_DIR}/output/{MAX_CLUSTER_NUM}/images/cluster/{con}/3D/{i}.png")
            plot(shower=Show2DImage, df=by3d, classes=classes, target="labels", x=X, y=Y, output_file=f"{DATA_DIR}/output/{MAX_CLUSTER_NUM}/images/cluster/{con}/2D/{i}.png")
        
        w_cluster_5 = w_cluster_5.drop(columns=[i for i in w_cluster_5.columns if "w" == str(i)[0]])
        w_cluster_5.to_csv(f"{DATA_DIR}/output/{MAX_CLUSTER_NUM}/ANALYTICS_CLUSTER/{con}.tsv", sep="\t")
        w_cluster_5[["labels", "Accession ID"]].to_csv(f"{DATA_DIR}/output/{MAX_CLUSTER_NUM}/CLUSTER_AND_ID/{con}.tsv", "\t", index=False)


    END_TIME = time.time()
    delta = END_TIME - START_TIME
    LOGGER.info(f"All Done! 処理時間:{round(delta,3)}s")
