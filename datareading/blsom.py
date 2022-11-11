import pickle
import random
import time
from collections import defaultdict

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from sklearn.preprocessing import LabelEncoder
from tqdm import tqdm

continent_color = {
    "Europe": "#FF1234",
    "North_America": "#00FF00",
    "Oceania": "#0000FF",
    "Asia": "#00cccc",
    "Africa": "#9966FF",
    "South_America": "#ffcc00",
}

targets = {
    "continent": {
        "Europe": "torch red",
        "North_America": "lime",
        "Oceania": "blue",
        #         "Asia": "sorbus",
        "Asia": "robin egg blue",
        "Africa": "fuchsia",
        #         "Africa": "violets are blue",
        "South_America": "yellow",
    },
    "continent4": {
        "Europe": "torch red",
        "North_America": "lime",
        "Oceania": "blue",
        #         "Asia": "sorbus",
        "Asia": "robin egg blue",
        "Africa": "fuchsia",
        #         "Africa": "violets are blue",
    },
    "Pango lineage": {
        "B.1.617.2": "fuchsia",
        "AY.4": "aqua",
        "AY.12": "blueviolet",
        "AY.3": "dark orange",
        "AY.6": "teal",
        "AY.9": "pink",
        "AY.5": "blue",
        "AY.7": "brown",
        "AY.11": "dark green",
    },  # AY3,1m AY.2,AY.10,AY.1,AY.8除外
    "month": {
        "11": "fuchsia",
        "12": "aqua",
        "01": "blueviolet",
        "02": "dark orange",
        "03": "teal",
        "04": "pink",
        "05": "blue",
        # "06": "",
        # "07": "brown",
        # "08": "dark green",
        # "09": "portage",
        # "10": "purple",
        # "11": "yellow",
    },
    "clade": {
        "G": "fuchsia",
        "GH": "aqua",
        "GR": "blueviolet",
        "GV": "dark orange",
        "GRY": "teal",
        "L": "pink",
        "O": "blue",
        "S": "brown",
        "V": "dark green",
    },
}


def create_color_code(n):
    targets = {}
    random.seed(0)
    while n > 0:
        r = random.randint(0, 255)
        g = random.randint(0, 255)
        b = random.randint(0, 255)
        if f"0x{hex(r)[2:]}{hex(g)[2:]}{hex(b)[2:]}" not in targets:
            targets[f"0x{hex(r)[2:]}{hex(g)[2:]}{hex(b)[2:]}"] = [r, g, b]
            n -= 1
    targets["white"] = [255, 255, 255]
    targets["black"] = [0, 0, 0]
    return targets


rgb_color = {
    "torch red": [255, 18, 52],
    "lime": [0, 255, 0],
    "blue": [0, 0, 255],
    "sorbus": [237, 107, 53],
    "fuchsia": [255, 0, 255],
    "yellow": [255, 255, 0],
    "brass": [181, 166, 66],
    "green": [0, 102, 0],
    "kelly green": [51, 204, 0],
    "aqua": [0, 255, 255],
    "portage": [153, 153, 255],
    "dark orange": [255, 127, 0],
    "orange red": [255, 69, 0],
    "sea green": [46, 139, 87],
    "bakers chocolate": [92, 51, 23],
    "barberry": [217, 217, 25],
    "dark green": [0, 100, 0],
    "indigo": [75, 0, 130],
    "deeppink": [255, 20, 147],
    "gray": [128, 128, 128],
    "purple": [128, 0, 128],
    "teal": [0, 128, 128],
    "sandy brown": [244, 164, 96],
    "brown": [143, 101, 82],
    "darkviolet": [148, 0, 211],
    "blueviolet": [138, 43, 226],
    "pink": [255, 192, 203],
    "violets are blue": [131, 102, 244],  # 9966FF
    "robin egg blue": [0, 204, 204],  # 00cccc
    "white": [255, 255, 255],
    "black": [0, 0, 0],
}



def weight_file2df(filepath):
    weight_dict = {}
    with open(filepath, "r") as r:
        datas = r.readlines()
        for i in range(1, len(datas) - 1, 2):
            data = datas[i+1].strip().split()
            weight_dict[datas[i].strip()] = np.array(list(map(float, data)))

    return weight_dict


def blsom_all_plot(data, rgb_color=None, target="clade"):
    if not rgb_color:
        rgb_color = create_color_code(1000)
    cor_count = {}
    le = LabelEncoder()

    labels = data[target]
    #     labels_id = le.fit_transform(labels)
    #     classes = le.classes_
    classes = np.array(list(set(labels)))
    # print(len(classes), classes)
    for i in data.index:
        key = f"{data['x'][i]} {data['y'][i]}"
        if key not in cor_count:
            cor_count[key] = defaultdict(int)

        cor_count[key][data[target][i]] += 1
    color_dict = {}
    # 画像データとして保存
    map_array = []
    target_c = []
    for cor in cor_count:
        count = 0
        x, y = map(int, cor.split(" "))
        map_array.append([x, y])

        for val in cor_count[cor]:
            if cor_count[cor][val] != 0:
                count += 1
                temp = val
        if count == 1:
            if target in targets.keys() and temp in targets[target].keys():
                target_c.append(rgb_color[targets[target][temp]])
            elif target not in targets.keys():
                #                 print(classes, temp, np.where(classes == temp))
                target_c.append(rgb_color[list(rgb_color.keys())[temp]])  # [np.where(classes == temp)[0][0]]]
                color_dict[temp] = list(rgb_color.keys())[np.where(classes == temp)[0][0]]
            else:
                target_c.append(rgb_color["white"])
        else:
            target_c.append(rgb_color["black"])
    # if color_dict:
    #     print(color_dict)
    if map_array == []:
        return
    sha = np.array(map_array).max(axis=0)
    img = np.full((sha[1] + 1, sha[0] + 1, 3), rgb_color["white"])

    for i, val in enumerate(map_array):
        img[val[1]][val[0]] = target_c[i]
    return img


def blsom_plot_3d(data, max_cors, output_file=None):
    continent_color = {
        "Europe": "#FF1234",
        "North_America": "#00FF00",
        "Oceania": "#0000FF",
        "Asia": "#00cccc",
        "Africa": "#ff00ff",
        "South_America": "#ffcc00",
    }
    X, Y = max_cors
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
    ax.set_xlim(0, X)
    ax.set_ylim(Y, 0)
    #     plt.title(output_file.stem)
    plt.title(data["labels"].unique()[0])
    # ax.set_zlim(0, 400)
    if output_file:
        plt.savefig(output_file)
    plt.show()
