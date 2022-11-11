from collections import defaultdict
import numpy as np
import pandas as pd

class Omicron2211PivotAdder(object):
    def __init__(self):
        self.threshholds = [i*0.01 for i in range(50, 101, 5)]
        # self.threshholds = [0.8]

    def _calculate_total_cluster_num(self, pivot):
        cluster_sum = defaultdict(int)
        for (index, values) in pivot.iterrows():
            cluster_sum[index[0]] += np.sum(values)

        return pd.Series(cluster_sum, name="Total")

    def _search_threshhold(self, pivot):
        match_cluster_nums, total_cluster_nums, extraction_ratios = [], [], []
        for threshhold in self.threshholds:
            count = sum([sum(values.drop(["Total", "cluster", "lineage2"]) / values["Total"] >= threshhold) for (_, values) in pivot.iterrows()])
            total = len(pivot.reset_index()["labels"].unique())

            match_cluster_nums += [count]
            total_cluster_nums += [total]
            extraction_ratios += [count / total]

        return pd.DataFrame(
            data = {"Threshhold": self.threshholds, "Match Cluster Num": match_cluster_nums, "Total Cluster Num": total_cluster_nums, "Extraction Ratio": extraction_ratios},
            columns=["Threshhold", "Match Cluster Num", "Total Cluster Num", "Extraction Ratio"])

    def _search_data_count_and_threshhold(self, pivot):
        max_num_in_cluster, total_data = [], []
        for c in pivot["cluster"].unique():
            max_num = np.max(pivot[pivot["cluster"] == c].drop(columns=["Total", "cluster", "lineage2"]).fillna(0).values)
            total_data += [pivot[pivot["cluster"] == c]["Total"].values[0]]
            max_num_in_cluster += [max_num]
        return pd.DataFrame(
            data = {"Cluster":pivot["cluster"].unique(), "MaxData":max_num_in_cluster, "TotalData": total_data,  "Percent": np.array(max_num_in_cluster)/np.array(total_data)},
            columns=["Cluster", "MaxData", "TotalData", "Percent"])
    
    @staticmethod
    def highlight_occupy_data(val):
        if val["Accession ID"] / val["Total"] >= 0.8:
            return ['background-color: #90ee90']*len(val)
        else:
            return ['background-color: white']*len(val) 

    @staticmethod
    def highlight_occupy_lineage_country(val):
        val_colors = []
        for i in val.index:
            # print(i)
            if i in ["cluster", "lineage2", "Total"]:
                val_colors.append('background-color: white')
                continue
            if val[i] / val["Total"] >= 0.8:
                val_colors[0] = 'background-color: #90ee90'
                val_colors[1] = 'background-color: #90ee90'
                # val_colors[2] = 'background-color: #90ee90'
                val_colors.append('background-color: #90ee90')
            else:
                val_colors.append('background-color: white')
        return val_colors

    @staticmethod
    def highlight_occupy_lineage(val):
        val_colors = []
        for i in val.index:
            # print(i)
            if i in ["cluster", "lineage2", "Total"]:
                val_colors.append('background-color: white')
                continue
            if val[i] / val["Total"] >= 0.8:
                val_colors[0] = 'background-color: #90ee90'
                val_colors[1] = 'background-color: #90ee90'
                val_colors.append('background-color: #90ee90')
            else:
                val_colors.append('background-color: white')
        return val_colors

    def add_pivot(self, writer, data):
        # 旧lineage別のデータ件数内訳
        pivot1 = pd.pivot_table(
            data, index=["labels", "lineage"], columns="country", values="Accession ID", aggfunc=len,
        )

        cluster_sum = self._calculate_total_cluster_num(pivot1)
        pivot1 = pivot1.merge(cluster_sum, left_on="labels", right_index=True)
        pivot1[["cluster", "lineage2"]] = pivot1.reset_index()[["labels", "lineage"]].values
        col = ["cluster", "lineage2"]
        col += [i for i in pivot1.columns if i not in col]
        pivot1 = pivot1.reindex(col, axis='columns')
        # pivot1.to_excel(writer, sheet_name="lineage", encoding="utf-8")
        pivot1.style.apply(self.highlight_occupy_lineage, axis=1).to_excel(writer, sheet_name="lineage", encoding="utf-8")

        # 新lineage別のデータ件数内訳
        # new_lineage = pd.pivot_table(
        #     data, index=["labels", "new_lineage"], values="Accession ID", aggfunc=len,
        # )
        # new_lineage = new_lineage.merge(cluster_sum, left_on="labels", right_index=True)
        # new_lineage['new_lineage2'] = new_lineage.reset_index()["new_lineage"].values
        # new_lineage = new_lineage.reindex(['new_lineage2', 'Accession ID', "Total"], axis='columns')
        # # new_lineage.to_excel(writer, sheet_name="new_lineage", encoding="utf-8")
        # new_lineage.style.apply(self.highlight_occupy_data, axis=1).to_excel(writer, sheet_name="new_lineage", encoding="utf-8")

       # country別のデータ件数内訳
        country = pd.pivot_table(
            data, index=["labels", "country"], values="Accession ID", aggfunc=len,
        )
        country = country.merge(cluster_sum, left_on="labels", right_index=True)
        country['country2'] = country.reset_index()["country"].values
        country = country.reindex(['country2', 'Accession ID', "Total"], axis='columns')
        # country.to_excel(writer, sheet_name="county", encoding="utf-8")
        country.style.apply(self.highlight_occupy_data, axis=1).to_excel(writer, sheet_name="county", encoding="utf-8")

        # 新lineage-country別のデータ件数内訳
        pivot2 = pd.pivot_table(
            data, index=["labels", "lineage"], columns="country", values="Accession ID", aggfunc=len,
        )
        pivot2 = pivot2.merge(cluster_sum, left_on="labels", right_index=True)
        pivot2[["cluster", "lineage2"]] = pivot2.reset_index()[["labels", "lineage"]].values
        col = ["cluster", "lineage2"]
        col += [i for i in pivot2.columns if i not in col]
        pivot2 = pivot2.reindex(col, axis='columns')
        # pivot2.to_excel(writer, sheet_name="new_lineage-country", encoding="utf-8")
        pivot2.style.apply(self.highlight_occupy_lineage_country, axis=1).to_excel(writer, sheet_name="new_lineage-country", encoding="utf-8")

        # pivot2の閾値探索
        result = self._search_threshhold(pivot2)
        result.to_excel(writer, sheet_name="search_threshhold", encoding="utf-8")

        # pivot2のデータ件数と閾値探索
        result = self._search_data_count_and_threshhold(pivot2)
        cor_max_data = np.corrcoef(result["MaxData"], result["Percent"])
        cor_total_data = np.corrcoef(result["TotalData"], result["Percent"])
        result.to_excel(writer, sheet_name="search_data_count_threshhold", encoding="utf-8")


class Omicron2201PivotAdder(object):
    def __init__(self):
        self.threshholds = [i*0.01 for i in range(50, 101, 5)]
        # self.threshholds = [0.8]

    def _calculate_total_cluster_num(self, pivot):
        cluster_sum = defaultdict(int)
        for (index, values) in pivot.iterrows():
            cluster_sum[index[0]] += np.sum(values)

        return pd.Series(cluster_sum, name="Total")

    def _search_threshhold(self, pivot):
        match_cluster_nums, total_cluster_nums, extraction_ratios = [], [], []
        for threshhold in self.threshholds:
            count = sum([sum(values.drop(["Total", "cluster", "lineage2", "new_lineage2"]) / values["Total"] >= threshhold) for (_, values) in pivot.iterrows()])
            total = len(pivot.reset_index()["labels"].unique())

            match_cluster_nums += [count]
            total_cluster_nums += [total]
            extraction_ratios += [count / total]

        return pd.DataFrame(
            data = {"Threshhold": self.threshholds, "Match Cluster Num": match_cluster_nums, "Total Cluster Num": total_cluster_nums, "Extraction Ratio": extraction_ratios},
            columns=["Threshhold", "Match Cluster Num", "Total Cluster Num", "Extraction Ratio"])

    def _search_data_count_and_threshhold(self, pivot):
        max_num_in_cluster, total_data = [], []
        for c in pivot["cluster"].unique():
            max_num = np.max(pivot[pivot["cluster"] == c].drop(columns=["Total", "cluster", "lineage2", "new_lineage2"]).fillna(0).values)
            total_data += [pivot[pivot["cluster"] == c]["Total"].values[0]]
            max_num_in_cluster += [max_num]
        return pd.DataFrame(
            data = {"Cluster":pivot["cluster"].unique(), "MaxData":max_num_in_cluster, "TotalData": total_data,  "Percent": np.array(max_num_in_cluster)/np.array(total_data)},
            columns=["Cluster", "MaxData", "TotalData", "Percent"])
    
    @staticmethod
    def highlight_occupy_data(val):
        if val["Accession ID"] / val["Total"] >= 0.8:
            return ['background-color: #90ee90']*len(val)
        else:
            return ['background-color: white']*len(val) 

    @staticmethod
    def highlight_occupy_lineage_country(val):
        val_colors = []
        for i in val.index:
            # print(i)
            if i in ["cluster", "lineage2", "new_lineage2", "Total"]:
                val_colors.append('background-color: white')
                continue
            if val[i] / val["Total"] >= 0.8:
                val_colors[0] = 'background-color: #90ee90'
                val_colors[1] = 'background-color: #90ee90'
                val_colors[2] = 'background-color: #90ee90'
                val_colors.append('background-color: #90ee90')
            else:
                val_colors.append('background-color: white')
        return val_colors

    @staticmethod
    def highlight_occupy_lineage(val):
        val_colors = []
        for i in val.index:
            # print(i)
            if i in ["cluster", "lineage2", "Total"]:
                val_colors.append('background-color: white')
                continue
            if val[i] / val["Total"] >= 0.8:
                val_colors[0] = 'background-color: #90ee90'
                val_colors[1] = 'background-color: #90ee90'
                val_colors.append('background-color: #90ee90')
            else:
                val_colors.append('background-color: white')
        return val_colors

    def add_pivot(self, writer, data):
        # 旧lineage別のデータ件数内訳
        pivot1 = pd.pivot_table(
            data, index=["labels", "lineage"], columns="country", values="Accession ID", aggfunc=len,
        )

        cluster_sum = self._calculate_total_cluster_num(pivot1)
        pivot1 = pivot1.merge(cluster_sum, left_on="labels", right_index=True)
        pivot1[["cluster", "lineage2"]] = pivot1.reset_index()[["labels", "lineage"]].values
        col = ["cluster", "lineage2"]
        col += [i for i in pivot1.columns if i not in col]
        pivot1 = pivot1.reindex(col, axis='columns')
        # pivot1.to_excel(writer, sheet_name="lineage", encoding="utf-8")
        pivot1.style.apply(self.highlight_occupy_lineage, axis=1).to_excel(writer, sheet_name="lineage", encoding="utf-8")

        # 新lineage別のデータ件数内訳
        new_lineage = pd.pivot_table(
            data, index=["labels", "new_lineage"], values="Accession ID", aggfunc=len,
        )
        new_lineage = new_lineage.merge(cluster_sum, left_on="labels", right_index=True)
        new_lineage['new_lineage2'] = new_lineage.reset_index()["new_lineage"].values
        new_lineage = new_lineage.reindex(['new_lineage2', 'Accession ID', "Total"], axis='columns')
        # new_lineage.to_excel(writer, sheet_name="new_lineage", encoding="utf-8")
        new_lineage.style.apply(self.highlight_occupy_data, axis=1).to_excel(writer, sheet_name="new_lineage", encoding="utf-8")

       # country別のデータ件数内訳
        country = pd.pivot_table(
            data, index=["labels", "country"], values="Accession ID", aggfunc=len,
        )
        country = country.merge(cluster_sum, left_on="labels", right_index=True)
        country['country2'] = country.reset_index()["country"].values
        country = country.reindex(['country2', 'Accession ID', "Total"], axis='columns')
        # country.to_excel(writer, sheet_name="county", encoding="utf-8")
        country.style.apply(self.highlight_occupy_data, axis=1).to_excel(writer, sheet_name="county", encoding="utf-8")

        # 新lineage-country別のデータ件数内訳
        pivot2 = pd.pivot_table(
            data, index=["labels", "lineage", "new_lineage"], columns="country", values="Accession ID", aggfunc=len,
        )
        pivot2 = pivot2.merge(cluster_sum, left_on="labels", right_index=True)
        pivot2[["cluster", "lineage2", "new_lineage2"]] = pivot2.reset_index()[["labels", "lineage", "new_lineage"]].values
        col = ["cluster", "lineage2", "new_lineage2"]
        col += [i for i in pivot2.columns if i not in col]
        pivot2 = pivot2.reindex(col, axis='columns')
        # pivot2.to_excel(writer, sheet_name="new_lineage-country", encoding="utf-8")
        pivot2.style.apply(self.highlight_occupy_lineage_country, axis=1).to_excel(writer, sheet_name="new_lineage-country", encoding="utf-8")

        # pivot2の閾値探索
        result = self._search_threshhold(pivot2)
        result.to_excel(writer, sheet_name="search_threshhold", encoding="utf-8")

        # pivot2のデータ件数と閾値探索
        result = self._search_data_count_and_threshhold(pivot2)
        cor_max_data = np.corrcoef(result["MaxData"], result["Percent"])
        cor_total_data = np.corrcoef(result["TotalData"], result["Percent"])
        result.to_excel(writer, sheet_name="search_data_count_threshhold", encoding="utf-8")

        # pivot4.style.apply(self._highlight_table, axis=1).to_excel(writer, sheet_name="county")

        # 新lineage別の平均値以上のクラスタのデータ件数内訳
        # cluster_num_avg = np.mean(cluster_sum)
        # pivot3 = pivot1[pivot1["Total"] >= cluster_num_avg]
        # pivot3.to_excel(writer, sheet_name="new_lineage平均")

        # cluster内の株数の平均値以上のクラスタのデータ件数内訳
        # pivot3 = pivot2[pivot2["Total"] >= cluster_num_avg]
        # pivot3.to_excel(writer, sheet_name="new_lineage-country平均")
        
        # pivot4の閾値探索
        # result = self._search_threshhold(pivot3)
        # result.to_excel(writer, sheet_name="search_threshhold平均")

        # cluster内の株数の50件以上のクラスタのデータ件数内訳
        # pivot4 = pivot2[pivot2["Total"] >= 50]
        # pivot4.to_excel(writer, sheet_name="new_lineage>50")

        # pivot4の閾値探索
        # result = self._search_threshhold(pivot4)
        # result.to_excel(writer, sheet_name="search_threshhold50")

        # pd.Series({
        #     "Cluster Num Avg": cluster_num_avg,
        # "Correlation To Max Data":cor_max_data,
        # "Correlation To Total Data":cor_total_data,
        # }).to_excel(writer, sheet_name="summary", header=None)


class DeltaPivotAdder(object):
    def __init__(self):
        pass

    def add_pivot(self, writer, data):
        pivot1 = pd.pivot_table(
            data, index=["labels", "lineage"], columns="country", values="Accession ID", aggfunc=len,
        )
        pivot1.to_excel(writer, sheet_name="pivot")