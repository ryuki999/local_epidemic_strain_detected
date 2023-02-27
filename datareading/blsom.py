from collections import defaultdict

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from sklearn.preprocessing import LabelEncoder
from tqdm import tqdm
from color import create_color_code, targets


def weight_file2df(filepath):
    weight_dict = {}
    with open(filepath, "r") as r:
        datas = r.readlines()
        for i in range(1, len(datas) - 1, 2):
            data = datas[i+1].strip().split()
            weight_dict[datas[i].strip()] = np.array(list(map(float, data)))

    return weight_dict


def blsom_all_plot(data, rgb_color=None, target="clade"):

    cor_count = {}
    labels = data[target]

    classes = np.array(list(set(labels)))
    # print(len(classes), classes)
    for i in data.index:
        key = f"{data['x'][i]} {data['y'][i]}"
        if key not in cor_count:
            cor_count[key] = defaultdict(int)

        cor_count[key][data[target][i]] += 1

    if not rgb_color:
        rgb_color = create_color_code(1000)

        for i in range(len(classes)):
            v = rgb_color.pop(list(rgb_color.keys())[i])
            rgb_color[classes[i]] = v
    # print(rgb_color)

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
                target_c.append(rgb_color[temp].rgb)
                # target_c.append(rgb_color[targets[target][temp]].rgb)
            elif target not in targets.keys():
                target_c.append(rgb_color[temp].rgb)  # [np.where(classes == temp)[0][0]]]
                color_dict[temp] = list(rgb_color.keys())[np.where(classes == temp)[0][0]]
            else:
                target_c.append([255, 255, 255])
                # target_c.append(rgb_color["white"])
        else:
            target_c.append([0, 0, 0])
            # target_c.append(rgb_color["black"])
    # if color_dict:
    #     print(color_dict)
    if map_array == []:
        return
    sha = np.array(map_array).max(axis=0)
    img = np.full((sha[1] + 1, sha[0] + 1, 3), [255, 255, 255])

    for i, val in enumerate(map_array):
        img[val[1]][val[0]] = target_c[i]
    return img


def blsom_plot_3d(data, max_cors, output_file=None):
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
