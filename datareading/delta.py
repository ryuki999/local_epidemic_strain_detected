import sys
from collections import defaultdict

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from tqdm import tqdm

sys.path.append("..")

import warnings

from blsom import blsom_all_plot, weight_file2df

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
    filename: str, header_columns: list = None, tar="all",
):
    with open(f"{filename}", encoding="utf-8") as f:
        data = f.read().split("\n")
        header = []
        max_x, max_y = map(int, data[0].split())
        for i in range(2, len(data) - 3):
            line = data[i].split("\t")
            x, y = line[1].split()
            head = line[0].strip().split("|")
            if "Human" in data[i] or "human" in data[i] or tar == "delta":
                header.append(head + [x] + [y] + [line[2]])
        # print(header)
        header2 = pd.DataFrame(header, columns=header_columns)
    return header2, max_x, max_y
    
class CreateDeltaHeaderDF(object):
    def __init__(self, out_filename, meta_filename):
        self.out_filename = out_filename
        self.meta_filename = meta_filename

    def create_df(self):
        # metaファイルの読み込み
        delta_header = pd.read_csv(self.meta_filename, sep="\t")
        delta_header = delta_header.iloc[delta_header["Virus name"].drop_duplicates().index]
        delta_header["Virus name"] = delta_header["Virus name"].apply(lambda x: x.replace(" ", "_"))
        delta_header = delta_header.join(
            delta_header["Collection date"].str.split("-", expand=True).set_axis(DATA_COLUMNS, axis=1)
        )
        delta_header = delta_header.join(
            delta_header["Location"].str.split("/", expand=True).apply(lambda x: x.str.strip())
        ).rename(columns={0: "continent", 1: "country", 2: "city"})
        delta_header = delta_header.rename(columns={"Clade": "clade", "Pango lineage": "lineage"})
        delta_header["continent"] = delta_header["continent"].apply(lambda x: x.strip().replace(" ", "_"))

        # outファイルの読み込み
        blsom, X, Y = blsom_outfile_to_df(
            self.out_filename, DELTA_HEADER_COLUMNS + ["x", "y", "distance"], tar="delta",
        )
        blsom["head"] = blsom["head"].apply(lambda x: x.strip(">"))

        # weightファイルの読み込み
        weight = weight_file2df(filepath="weight.100")

        delta_blsom = blsom.merge(delta_header, left_on="head", right_on="Virus name")
        delta_blsom = weight.merge(delta_blsom, right_on=["x", "y"], left_on=["x", "y"], how="left")

        return delta_blsom, X, Y


class ShowImage(object):
    def __init__(self, dataframe, x, y):
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
        for i in tqdm(data.index):
            key = f"{data['x'][i]} {data['y'][i]}"
            if key not in cor_count:
                cor_count[key] = defaultdict(int)
            cor_count[key][data["continent"][i]] += 1

        fig = plt.figure(figsize=(12, 8))
        ax = fig.add_subplot(111, projection="3d")

        for i in tqdm(cor_count):
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
