import re
from pathlib import Path

import numpy as np
import pandas as pd
from tqdm import tqdm
import time
import openpyxl

from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor
from codon import codon_table_v2


class SnvDetector(object):
    BASE = ["A", "T", "G", "C"]

    def __init__(self, cluster_threshhold, dir_path, output_file):
        self.cluster_threshhold = cluster_threshhold
        self.dir_path = dir_path
        self.output_file = output_file

    def extract_snv(self, threshhold):
        START_TIME = time.time()
        base = self._read_base()
        END_TIME = time.time()
        delta = END_TIME - START_TIME
        print(f"{delta}s")

        # 変異率が高いPOSのみを抽出する
        total_T = self._extract_high_mutation_rate_from_base(base, threshhold=threshhold)
        END_TIME = time.time()
        delta = END_TIME - START_TIME
        print(f"{delta}s")

        # 参照配列のコドンとSNP後のコドンを比較するためにSNP後のコドンを導出する
        total_snp = self._derive_codon_after_snp(base)
        END_TIME = time.time()
        delta = END_TIME - START_TIME
        print(f"{delta}s")

        columns = [
            "POS",
            "Reference Seq. (MN908947.3)",
            "遺伝子名：重複部位は後者の方を優先している(ORF7a/7b)",
            "遺伝子アミノ酸位置1",
            "遺伝子アミノ酸位置2",
            "[Comp]B.1.617.2",
            "Original",
            "Original_Amino",
            "Omicron(ALL_211209)"
        ]

        cluster_columns = [i for i in total_T.columns if i[0] == "C" and "_" not in i]
        result_columns = columns + [f"{c}{n}" for c in cluster_columns for n in ["", "_ALT", "_SNP_Codon", "_SNP_Amino", "_同義/非同義"]]

        result = total_T.merge(total_snp, left_on="POS", right_on="POS")
        result = result.reindex(result_columns, axis=1)
        result.to_excel(self.output_file, sheet_name="snv")

        # with pd.ExcelWriter(self.output_file, engine='openpyxl', mode="a") as writer:
        #     result.to_excel(writer, sheet_name=f"snv_unique")

        outs = []
        for c in range(self.cluster_num):
            out = result[columns+[i for i in result.columns if f"C{c}" == i.split("_")[0]]]
            out["クラスタ番号"] = f"C{c}"
            out = out[(out[f"C{c}"] >= out["Omicron(ALL_211209)"]) & (out[f"C{c}"] >= threshhold)]
            # out.to_excel(writer, sheet_name=f"C{c}")
            outs.append(pd.DataFrame(out.values, columns=columns+["Ratio", "SNV_BASE", "_SNV_Codon", "SNV_Amino", "同義/非同義", "クラスタ番号"]))
        outs = pd.concat(outs)
        outs["SNV"] = outs["Reference Seq. (MN908947.3)"] + outs["POS"].astype(str) + outs["SNV_BASE"]
        outs = outs.reindex(columns=["クラスタ番号"]+columns+["Ratio", "SNV_BASE", "_SNV_Codon", "SNV_Amino", "同義/非同義"])
        with pd.ExcelWriter(self.output_file, engine='openpyxl', mode="a") as writer:
            outs.to_excel(writer, sheet_name="snv_total")

        # total = self._arrange_excel(result, columns)
        # with pd.ExcelWriter(output_file, engine='openpyxl', mode="a") as writer:
        #     total.style.apply(self.highlight_occupy_snv, axis=1).to_excel(writer, sheet_name="snv_total")

    @staticmethod
    def highlight_occupy_snv(val):
        val_colors = []
        for i in val.index:
            # if i in ["cluster", "lineage2", "Total"]:
            #     val_colors.append('background-color: white')
            #     continue
            if i == "Ratio" and val["Ratio"] >= 80:
                # val_colors[0] = 'background-color: #90ee90'
                # val_colors[1] = 'background-color: #90ee90'
                val_colors.append('background-color: #90ee90')
            else:
                val_colors.append('background-color: white')
        return val_colors

    @staticmethod
    def _arrange_excel(result, columns):
        result_columns = columns + ["C", "Ratio", "SNV_Base", "SNV_Codon", "SMV_Amino", "同義/非同義", "SNV"]
        total = []
        for _, rows in result.iterrows():
            max_snv = rows[[i for i in rows.index if len(i) == 2]].astype(float)

            element = list(rows[columns].values)
            element += [max_snv.idxmax()]
            element += list(rows[[i for i in rows.index if max_snv.idxmax() == i[:2]]].values)
            element += [f"{rows['Reference Seq. (MN908947.3)']}{rows['POS']}{rows[f'{max_snv.idxmax()}_ALT']}"]
            total.append(pd.DataFrame([element], columns=result_columns))
        return pd.concat(total).reset_index(drop=True)

    def _read_base(self):
        """BASEファイルの読み込みとoutputファイルとの結合"""
        base = pd.read_excel("SNPcheck4PCRprimer211217_Omicrons.xlsx", skiprows=1)
        cluster_nums = []
        for filepath in Path(self.dir_path).glob("**/*.output"):
            cluster_num = re.search(r"(?<=CLUSTER)\d{1,2}", filepath.name).group()
            cluster = pd.read_csv(filepath, sep="\t", skiprows=3)

            if len(cluster.columns[10:]) < self.cluster_threshhold:
                continue
            
            # クラスタの全株数
            cluster["total"] = cluster.iloc[:, 10:].count(axis=1)
            # クラスタでSNVが起こっている株/クラスタの全株数
            cluster[f"C{cluster_num}"] = (cluster.iloc[:, 10:] == 1).sum(axis=1) * 100  / cluster["total"]
            cluster[f"C{cluster_num}_ALT"] = cluster["ALT"].apply(lambda x: x.split(",")[0])

            base = base.merge(
                cluster[["POS", f"C{cluster_num}", f"C{cluster_num}_ALT"]], left_on="POS", right_on="POS", how="outer",
            )
            base[f"C{cluster_num}"] = base[f"C{cluster_num}"].fillna(0)
        self.cluster_num = int(cluster_num)+1

        return base

    def _extract_high_mutation_rate_from_base(self, base, threshhold=None):
        """変異率が高いPOSのみを抽出する"""
        cluster_percent = base[[i for i in base.columns if i[0] == "C" and "_" not in i]]

        total = pd.DataFrame()
        for idx in cluster_percent.index:
            row = cluster_percent.iloc[idx]

            # 閾値処理
            if threshhold is None:
                total = pd.concat([total, base.iloc[idx]], axis=1)
            # 初期株群よりも高い変異率 & クラスタ唯一に起こっている SNVを抽出
            if threshhold is not None \
                and len(row[(row >= base.iloc[idx]["Omicron(ALL_211209)"])].index) > 0 \
                and len(row[row >= threshhold].index) == 1:
                total = pd.concat([total, base.iloc[idx]], axis=1)

        total_T = total.T.drop(columns=[i for i in total.T.columns if "T/F" in i])

        return total_T

        
    def _derive_codon_after_snp(self, base):
        """参照配列のコドンとSNP後のコドンを比較するためにSNP後のコドンを導出する"""
        total_snp = pd.DataFrame()
        # print(base.columns)
        c_nums = [
            re.search(r"(?<=C)\d{1,2}", i).group() for i in base.columns if len(i.split("_")) < 2 and i[0] == "C"
        ]
        base = base[(base["遺伝子アミノ酸位置1"].notna()) | (base["遺伝子アミノ酸位置2"].notna())]
        # base = base[(base["遺伝子アミノ酸位置1"].notna()) and (base["遺伝子アミノ酸位置2"].notna())]
        params = [(base[["POS", "Reference Seq. (MN908947.3)", f"C{c_num}_ALT", "遺伝子アミノ酸位置1", "遺伝子アミノ酸位置2", "遺伝子名：重複部位は後者の方を優先している(ORF7a/7b)"]].values, c_num) for _, c_num in enumerate(c_nums)]
        with ProcessPoolExecutor(max_workers=len(c_nums)) as executor:
            # print(f"CLUSTER{c_num}")
            results = executor.map(self._extract_codon, params)

        for i, adding_snp in enumerate(results):
            ### 大丈夫?
            if i != 0:
                adding_snp = adding_snp.drop(columns=["POS", "Original", "Original_Amino"])
            total_snp = pd.concat([total_snp, adding_snp], axis=1)
        return total_snp
        
    @staticmethod
    def _extract_codon(params):
        """コドンを導出する"""
        def _substitute_snv_base(Pos, codon, snv_codon, codon_pos):
            """SNV塩基を元の塩基と置換する"""
            list_codon = list(codon)
            list_snp_codon = list(snv_codon)
            list_codon[2-codon_pos] = list_snp_codon[2-codon_pos]
            amino = codon_table_v2[codon]
            snp_amino = ""
            if "".join(list_codon) in codon_table_v2:
                snp_amino = codon_table_v2["".join(list_codon)]
                
            if amino == snp_amino:
                synonymous_codon = "同義コドン"
            else:
                synonymous_codon = "非同義コドン"

            return pd.DataFrame(
                [[Pos - codon_pos, codon, "".join(list_codon), amino, snp_amino, synonymous_codon,]]
                # [[Pos - (codon_pos - 2), codon, "".join(list_codon), amino, snp_amino, synonymous_codon,]]
            )

        base_values, c_num = params
        codon = ""
        snv_codon = ""
        adding_snv = []
        BASE = ["A", "T", "G", "C", "*"]
        PreProtein = None
        for Pos, RefSeq, CnumALT, BasePos1, BasePos2, Protein in tqdm(base_values):
            if (type(BasePos1) != float and "@1" in BasePos1) or (type(BasePos2) != float and "@1" in BasePos2):
            # if PreProtein != Protein:
                codon = ""
                snv_codon = "" 

            codon += RefSeq
            if CnumALT is not np.nan and type(CnumALT) != float:
                snv_codon += CnumALT
            else:
                snv_codon += RefSeq

            if len(codon) >= 3:
                # if "112" in str(Pos):
                #     print(Pos, codon, snv_codon)
                # print(codon, snv_codon)
                new_snv0 = _substitute_snv_base(Pos, codon, snv_codon, 0)
                new_snv1 = _substitute_snv_base(Pos, codon, snv_codon, 1)
                new_snv2 = _substitute_snv_base(Pos, codon, snv_codon, 2)
                adding_snv.append(new_snv0)
                adding_snv.append(new_snv1)
                adding_snv.append(new_snv2)

                codon = ""
                snv_codon = ""
            PreProtein = Protein
        adding_snv = pd.concat(adding_snv)
        adding_snv = adding_snv.reset_index(drop=True)
        adding_snv.columns = [
            "POS",
            "Original",
            f"C{c_num}_SNP_Codon",
            "Original_Amino",
            f"C{c_num}_SNP_Amino",
            f"C{c_num}_同義/非同義",
        ]
        return adding_snv
