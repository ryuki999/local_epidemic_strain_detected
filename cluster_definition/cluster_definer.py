import sys
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from tqdm import tqdm
import time
from sklearn.cluster import KMeans

sys.path.append("..")
from datareading.blsom import blsom_all_plot


class KmeansTotalClusterDefiner(object):
    def __init__(self, n_clusters, kmeans_cluster_num, cluster_range, verbose=0, logger=None):
        # self.logger = logger
        self.n_clusters = n_clusters
        self.kmeans_cluster_num = kmeans_cluster_num
        # self.init = init
        self.verbose = verbose
        self.cluster_range = cluster_range
        self.class_labels = None
        self.df = None
        self.weight_labels = None

    def df2cluster(self, df, blsom_weight):
        init_points_coordinates = np.array(df[["x y"]].value_counts().reset_index().sort_values(by=[0, "x y"], ascending=[False, True])["x y"])
        # print(df[["x y"]].value_counts().reset_index().sort_values(by=[0, "x y"], ascending=[False, True]))
        init_points =  np.array([blsom_weight.loc[xy] for xy in tqdm(init_points_coordinates)])
        x = init_points
        
        kmeans = KMeans(n_clusters=self.kmeans_cluster_num, init=init_points[:self.kmeans_cluster_num], verbose=self.verbose)
        kmeans_model = kmeans.fit(x)
        self.kmeans_model = kmeans_model
        # print(kmeans_model.labels_)

        self.class_labels = np.unique(kmeans_model.labels_)
        kmeans_labels = pd.Series(kmeans_model.labels_, name="kmeans_labels")
        weight_labels = pd.Series(kmeans_labels, name="kmeans_labels")
        weight_labels.index = init_points_coordinates
        blsom = df.merge(weight_labels, left_on=["x y"], right_index=True)
        
        self.df = blsom
        self.weight_labels = weight_labels
        self.blsom_weight = blsom_weight
        return blsom

    def extract_certain_range(self):
        n = self.n_clusters if len(self.weight_labels.index) > self.n_clusters else len(self.weight_labels.index)
        center_points = self.weight_labels[:n]
        
        eu3 = []
        for i, (center_point, c) in tqdm(enumerate(zip(center_points.index, center_points.values))):
            # 中心点が他のクラスタに含まれているときそのまま
            if i > 0:
                temp = pd.concat(eu3)
            if i > 0 and len(temp[temp["x y"] == center_point].index) != 0:
                continue
            # 中心点が他のクラスタに含まれていないとき、新しいクラスタとして統合
            else:
                base_x, base_y = map(int, center_point.split(" "))
                fil_df = self.df[(self.df["kmeans_labels"] == c) &
                    (self.df["x"].astype(int) >= base_x - self.cluster_range) &
                    (self.df["x"].astype(int) < base_x + self.cluster_range) &
                    (self.df["y"].astype(int) >= base_y - self.cluster_range) &
                    (self.df["y"].astype(int) < base_y + self.cluster_range)]
                fil_df["labels"] = i
                eu3.append(fil_df)
        eu3 = pd.concat(eu3)
        # print(eu3[eu3.drop(columns=["labels"]).duplicated()])
        # return eu3

        drop_dup_subset_columns = [i for i in eu3.columns if i != "labels"]
        drop_duplicates_data = eu3.drop_duplicates(keep="first", subset=drop_dup_subset_columns)

        replace_dict = {}
        for i, c in enumerate(drop_duplicates_data["labels"].unique()):
            replace_dict[c] = i
 
        drop_duplicates_data["labels"] = drop_duplicates_data["labels"].map(replace_dict)
        print(drop_duplicates_data[drop_duplicates_data.drop(columns=["labels"]).duplicated(keep=False)])
        return drop_duplicates_data
# 845, 1410


class KmeansClusterDefiner(object):
    def __init__(self, n_clusters, cluster_range, verbose=0, logger=None):
        # self.logger = logger
        self.n_clusters = n_clusters
        # self.init = init
        self.verbose = verbose
        self.cluster_range = cluster_range
        self.class_labels = None
        self.df = None
        self.weight_labels = None

    def df2cluster(self, df, blsom_weight):
        """
        TODO: 初期点の決定の箇所がデータ件数の多い順になっていない
        """
        init_points_coordinates = np.array(df[["x y"]].value_counts().reset_index().sort_values(by=["x y"], ascending=[True])["x y"])
        init_points =  np.array([blsom_weight.loc[xy] for xy in tqdm(init_points_coordinates)])
        x = init_points
        
        kmeans = KMeans(n_clusters=self.n_clusters, init=init_points[:self.n_clusters], verbose=self.verbose)
        kmeans_model = kmeans.fit(x)
        self.kmeans_model = kmeans_model
        # print(kmeans_model.labels_)

        self.class_labels = np.unique(kmeans_model.labels_)
        labels = pd.Series(kmeans_model.labels_, name="labels")
        weight_labels = pd.Series(labels, name="labels")
        weight_labels.index = init_points_coordinates
        blsom = df.merge(weight_labels, left_on=["x y"], right_index=True)
        
        self.df = blsom
        self.weight_labels = weight_labels
        self.blsom_weight = blsom_weight
        return blsom
    
    def extract_certain_range(self):
        df = self.df
        center_points = []
        # for c in range(len(self.class_labels)):
        #     idx = self.weight_labels[self.weight_labels == c].index
        #     dist = np.linalg.norm(self.blsom_weight.loc[idx]-self.kmeans_model.cluster_centers_[c], axis=1)
        #     center_points.append(idx[np.argmin(dist)])

        for c in range(len(self.class_labels)):
            center_points.append(self.df[self.df["labels"]==c]["x y"].value_counts().sort_values(ascending=False).index[0])

        eu3 = []
        for c, center_point in tqdm(enumerate(center_points)):
            base_x, base_y = map(int, center_point.split(" "))
            fil_df = self.df[(self.df["labels"] == c) &
                (self.df["x"].astype(int) >= base_x - self.cluster_range) &
                (self.df["x"].astype(int) < base_x + self.cluster_range) &
                (self.df["y"].astype(int) >= base_y - self.cluster_range) &
                (self.df["y"].astype(int) < base_y + self.cluster_range)]
            eu3.append(fil_df)
        eu3 = pd.concat(eu3)

        return eu3

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
        eu3 = []
        for c, center_point in true_center_cors:
            base_x, base_y = map(int, center_point.split(" "))
            for x in range(base_x - cluster_range, base_x + cluster_range):
                for y in range(base_y - cluster_range, base_y + cluster_range):
                    fil_df = eu2[(eu2["labels"] == c) & (eu2["x y"] == f"{x} {y}")]
                    eu3.append(fil_df)
        eu3 = pd.concat(eu3)
        eu3.columns = columns
        return eu3

    def _remove_duplicates_data(self, data):
        # 重複してるデータ
        duplicates_data = data[data.drop(columns=["dist", "labels"]).duplicated()]
        # w_columns = [i for i in duplicates_data.columns if str(i)[0] == "w"]
        unique_cors = duplicates_data["x y"].unique()

        filter_duplicates_data = []
        for cor in unique_cors:
            min_dist = min(duplicates_data[duplicates_data["x y"] == cor]["dist"])
            idx = duplicates_data[(duplicates_data["x y"] == cor) & (duplicates_data["dist"] == min_dist)].index

            virus_name = duplicates_data.loc[idx]
            filter_duplicates_data.append(virus_name)

        filter_duplicates_data = pd.concat(filter_duplicates_data)
        filter_duplicates_data.columns = duplicates_data.columns

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
