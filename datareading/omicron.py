import sys
from collections import defaultdict
import re
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from tqdm import tqdm

sys.path.append("..")

import warnings

from datareading.blsom import  blsom_all_plot, weight_file2df

warnings.simplefilter("ignore")

plt.rcParams["font.family"] = "Times New Roman"

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


def blsom_outfile_to_df(
    filename: str,
):
    with open(f"{filename}", encoding="utf-8") as f:
        data = f.read().split("\n")
        max_x, max_y = map(int, data[0].split())

    header = pd.read_csv(filename, sep="\t", skiprows=2, skipfooter=2, names=["Virus name", "x y", "distance"])
    header["Virus name"] = header["Virus name"].str.replace(" >", "")
    header["x y"] = header["x y"].apply(lambda x: re.sub('[ 　]+', ' ', x))
    header["x y"] = header["x y"].apply(lambda x: x.lstrip())
    header["x y"] = header["x y"].apply(lambda x: x.rstrip())
    header = header.join(header["x y"].str.split(" ", expand=True)).rename(columns={0: "x", 1: "y"})
    return header, max_x, max_y


class CreateOmicronHeaderDF(object):
    def __init__(self, out_filename, meta_filename, weight_file=None):
        self.out_filename = out_filename
        self.meta_filename = meta_filename
        self.weight_file = weight_file

    def create_df(self):
        delta_header = pd.read_csv(self.meta_filename, sep="\t")
        delta_header = delta_header.iloc[delta_header["Virus name"].drop_duplicates().index]

        delta_header = delta_header.join(
            delta_header["Collection date"].str.split("-", expand=True).set_axis(DATA_COLUMNS, axis=1)
        )
        delta_header = delta_header.join(
            delta_header["Location"].str.split("/", expand=True).apply(lambda x: x.str.strip())
        ).rename(columns={0: "continent", 1: "country", 2: "city"})
        delta_header = delta_header.rename(columns={"Clade": "clade", "Pango lineage": "lineage"})    
        delta_header["continent"] = delta_header["continent"].str.replace(" ", "_")

        # outファイルの読み込み
        blsom, X, Y = blsom_outfile_to_df(self.out_filename)

        blsom["Virus name"] = blsom["Virus name"].str.split("|", expand=True).iloc[:,0]
        
        delta_blsom = blsom.merge(delta_header, left_on="Virus name", right_on="Virus name", how="left")

        if self.weight_file is not None:
            # weightファイルの読み込み
            blsom_weight = weight_file2df(filepath=self.weight_file)

        self.delta_blsom = delta_blsom
        self.X = X
        self.Y = Y
        self.blsom_weight = blsom_weight


class ShowImage(object):
    def __init__(self, x, y, dataframe=None):
        self.dataframe = dataframe
        self.x = x
        self.y = y

    def show_color_by_continent(self):
        target = "continent"
        plt.figure(figsize=(12, 8))
        img = blsom_all_plot(self.dataframe, target=target)
        #     plt.subplot(3, 2, i+1)
        #     ax = plt.gca()
        #     ax.axes.xaxis.set_visible(False)
        #     ax.axes.yaxis.set_visible(False)
        # plt.title(f"{con}")
        #     plt.title(f"{con}+{target}")
        plt.xlim(0, self.x)
        plt.ylim(self.y, 0)
        plt.xlabel("x")
        plt.ylabel("y")
        plt.imshow(img)
        plt.savefig(f"test.png")
        # plt.clf()

    def create_img(self):
        target = "continent"
        plt.figure(figsize=(12, 8))
        img = blsom_all_plot(self.dataframe, target=target)
        # plt.subplot(6, 8, num) # subplot(m, n, p): mは行, nは列, pは位置 → mは2や3だと上下離れてしまうので5にしています
        # plt.title(f"{con[:4]}+{mon}")
        plt.xlim(0, self.x)
        plt.ylim(self.y, 0)
        if img is not None:
            plt.imshow(img)
        else:
            plt.imshow(np.full((self.x + 1, self.y + 1, 3), [255, 255, 255]))

    def blsom_plot_3d(self, data, output_file):
        continent_color = {
            "Europe": "#FF1234",
            "North_America": "#00FF00",
            "Oceania": "#0000FF",
            "Asia": "#00cccc",
            "Africa": "#ff00ff",
            "South_America": "#ffcc00",
        }
        #     plt.figure(figsize=(124, 96))
        cor_count = {}
        for i in data.index:
        # for i in tqdm(data.index):
            key = f"{data['x'][i]} {data['y'][i]}"
            if key not in cor_count:
                cor_count[key] = defaultdict(int)
            cor_count[key][data["continent"][i]] += 1

        fig = plt.figure(figsize=(12, 8))
        ax = fig.add_subplot(111, projection="3d")

        for i in cor_count:
        # for i in tqdm(cor_count):
            x, y = map(int, i.split(" "))
            arr = sorted(cor_count[i].items(), key=lambda x: x[1])
            height = 0
            for j in arr:
                ax.bar3d(
                    x,
                    y,
                    cor_count[i][j[0]],
                    dx=0.5,
                    dy=0.5,
                    dz=-(cor_count[i][j[0]] - height),
                    color=continent_color[j[0]],
                )
                height = cor_count[i][j[0]]

        ax.set_xlabel("x")
        ax.set_ylabel("y")
        ax.set_zlabel("#Seqs")
        # ax.zaxis.set_rotate_label(False)
        # ax.set_zlabel("  #Seqs", rotation=0, fontsize=12)
        ax.set_xlim(0, self.x)
        ax.set_ylim(self.y, 0)
        #     plt.title(output_file.stem)
        # ax.set_zlim(0, 400)
        plt.savefig(output_file)
