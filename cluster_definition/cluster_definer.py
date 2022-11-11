import sys
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import time
sys.path.append("..")
from datareading.blsom import blsom_all_plot


class ClasterDefiner(object):
    def __init__(self, logger):
        self.logger = logger

    # def zscore(self, data, sigma=15, axis=None):
    #     x = data.values
    #     z = data.index
    #     xmean = x.mean(axis=axis, keepdims=True)
    #     xstd = np.std(x, axis=axis, keepdims=True)
    #     zscore = (x - xmean) / xstd
    #     print(zscore, xmean, xstd)
    #     score = zscore[x > xmean + sigma * xstd]
    #     return score, z[: len(score)]

    def data_from_top_n(self, data, n=5, axis=None):
        x = data.values
        z = data.index
        xmean = x.mean(axis=axis, keepdims=True)
        xstd = np.std(x, axis=axis, keepdims=True)
        zscore = (x - xmean) / xstd

        score = zscore[:n]
        return score, z[: n]

    def df2cluster(self, df, blsom_weight, center_cors, sigma=1.5, cluster_range=5):
        columns = list(df.columns) + ["dist", "labels"]
        eu2 = pd.DataFrame(columns=columns)
        # unique_labels = [eu[["x y"]== " ".join(i)]["labels"].unique()[0] for i in center_points]

        true_center_cors = []
        c = 0
        for _, center_point in enumerate(center_cors):
            self.logger.info(f"{c} {center_point} {len(eu2[eu2['x y'] == center_point].index)}")

            # 中心点が他のクラスタに含まれているときそのまま
            if len(eu2[eu2["x y"] == center_point].index) != 0:
                continue
            # 中心点が他のクラスタに含まれていないとき、新しいクラスタとして統合
            else:
                cen = blsom_weight[center_point]
                dist_list = [np.linalg.norm(cen - blsom_weight[com]) for com in df["x y"].values]

                df["dist"] = dist_list

                # omicron oceaniaの例のように 0 < -幾つかのようになり抽出データ0になる場合もある
                cluster = df[df["dist"] < sigma * np.std(df["dist"])]
                cluster["labels"] = c

                eu2 = pd.concat([eu2, cluster])
                true_center_cors.append((c, center_point))
                c += 1

        eu2 = self._remove_duplicates_data(eu2)
        # ±5に範囲を絞る
        eu3 = pd.DataFrame(columns=columns)
        for c, center_point in true_center_cors:
            base_x, base_y = map(int, center_point.split(" "))
            for x in range(base_x - cluster_range, base_x + cluster_range):
                for y in range(base_y - cluster_range, base_y + cluster_range):
                    fil_df = eu2[(eu2["labels"] == c) & (eu2["x y"] == f"{x} {y}")]
                    eu3 = pd.concat([eu3, fil_df])
        return eu3

    def _remove_duplicates_data(self, data):
        # 重複してるデータ
        duplicates_data = data[data.drop(columns=["dist", "labels"]).duplicated()]
        # w_columns = [i for i in duplicates_data.columns if str(i)[0] == "w"]
        unique_cors = duplicates_data["x y"].unique()

        filter_duplicates_data = pd.DataFrame()
        for cor in unique_cors:
            min_dist = min(duplicates_data[duplicates_data["x y"] == cor]["dist"])
            idx = duplicates_data[(duplicates_data["x y"] == cor) & (duplicates_data["dist"] == min_dist)].index

            virus_name = duplicates_data.loc[idx]
            filter_duplicates_data = pd.concat([filter_duplicates_data, virus_name])

        # 重複を削除したデータ
        drop_duplicates_data = data.loc[data.drop(columns=["dist", "labels"]).drop_duplicates(keep=False).index]
        drop_duplicates_data = pd.concat([drop_duplicates_data, filter_duplicates_data])
        return drop_duplicates_data

    def cluster_plot(self, data, max_cors, output_file):
        plt.figure(figsize=(32, 28))
        X, Y = max_cors
        for c in data["labels"].unique():
            img = None
            blsom_cluster = data[data["labels"] == c]
            blsom_cluster["labels"] = int(str(c).split("_")[0])
            target = "labels"
            plt.subplot(20, 5, c + 1)  # subplot(m, n, p): mは行, nは列, pは位置 → mは2や3だと上下離れてしまうので5にしています
            img = blsom_all_plot(blsom_cluster, target=target)
            plt.title(f"{c}")
            plt.xlim(0, X)
            plt.ylim(Y, 0)
            if img is not None:
                plt.imshow(img)
                con = str(data["continent"].iloc[0])
        plt.savefig(output_file)
