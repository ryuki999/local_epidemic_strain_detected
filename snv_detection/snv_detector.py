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

    def __init__(self, cluster_threshhold):
        self.cluster_threshhold = cluster_threshhold

    def extract_snv(self, dir_path, output_file, threshold):
        START_TIME = time.time()
        base = self._read_base(dir_path)
        END_TIME = time.time()
        delta = END_TIME - START_TIME
        print(f"{delta}s")

        # 変異率が高いPOSのみを抽出する
        total_T = self._extract_high_mutation_rate_from_base(base, threshold=threshold)
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
            "遺伝子アミノ酸位置",
            "[Comp]B.1.617.2",
            "Original",
            "Original_Amino",
        ]
        cluster_columns = [i for i in total_T.columns if i[0] == "C" and "_" not in i]
        columns += [f"{c}{n}" for c in cluster_columns for n in ["", "_ALT", "_SNP_Codon", "_SNP_Amino", "_同義/非同義"]]

        result = total_T.merge(total_snp, left_on="POS", right_on="POS")
        result = result.reindex(columns, axis=1)
        result.to_excel(output_file, sheet_name="snv")

        total = self._arrange_excel(result)
        with pd.ExcelWriter(output_file, engine='openpyxl', mode="a") as writer:
            total.style.apply(self.highlight_occupy_snv, axis=1).to_excel(writer, sheet_name="snv_total")

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
    def _arrange_excel(result):
        columns = [
            "POS",
            "Reference Seq. (MN908947.3)",
            "遺伝子名：重複部位は後者の方を優先している(ORF7a/7b)",
            "Original",
            "Original_Amino",
        ]
        result_columns = columns + ["C", "Ratio", "SNV_Base", "SNV_Codon", "SMV_Amino", "同義/非同義", "置換後塩基", "SNV"]
        total = []
        for i, rows in result.iterrows():
            max_snv = rows[[i for i in rows.index if len(i) == 2]].astype(float)
            BMatch = None
            for B in SnvDetector.BASE:
                if B in rows[f"{max_snv.idxmax()}_ALT"]:
                    BMatch = B
                    break

            element = list(rows[columns].values)
            element += [max_snv.idxmax()]
            element += list(rows[[i for i in rows.index if max_snv.idxmax() == i[:2]]].values)
            element += [BMatch,f"{rows['Reference Seq. (MN908947.3)']}{rows['POS']}{BMatch}"]
            total.append(pd.DataFrame([element], columns=result_columns))
        return pd.concat(total).reset_index(drop=True)

    def _read_base(self, dir_path: str):
        """BASEファイルの読み込みとoutputファイルとの結合"""
        base = pd.read_csv("SNP_BASE.csv")
        cluster_nums = []
        for filepath in Path(dir_path).glob("**/*.output"):
            cluster_num = re.search(r"(?<=CLUSTER)\d{1,2}", filepath.name).group()
            cluster = pd.read_csv(filepath, sep="\t", skiprows=3)

            if len(cluster.columns[10:]) < self.cluster_threshhold:
                continue

            cluster["count"] = (cluster.iloc[:, 10:] == 1).sum(axis=1)
            cluster["total"] = cluster.iloc[:, 10:].count(axis=1)
            cluster[f"C{cluster_num}"] = cluster["count"] * 100 / cluster["total"]
            cluster[f"C{cluster_num}_ALT"] = cluster["ALT"]
            base = base.merge(
                cluster[["POS", f"C{cluster_num}", f"C{cluster_num}_ALT"]], left_on="POS", right_on="POS", how="outer",
            )
            base[f"C{cluster_num}"] = base[f"C{cluster_num}"].fillna(0)

            cluster_nums.append(f"C{cluster_num}")

        return base

    def _extract_high_mutation_rate_from_base(self, base, threshold=None):
        """変異率が高いPOSのみを抽出する"""
        cluster_percent = base[[i for i in base.columns if i[0] == "C" and "_" not in i]]

        total = pd.DataFrame()
        for idx in cluster_percent.index:
            row = cluster_percent.iloc[idx]
            if threshold is None:
                total = pd.concat([total, base.iloc[idx]], axis=1)
            if threshold is not None and len(row[row >= threshold].index) == 1:
                total = pd.concat([total, base.iloc[idx]], axis=1)

        total_T = total.T.drop(columns=[i for i in total.T.columns if "T/F" in i])

        return total_T

        
    def _derive_codon_after_snp(self, base):
        """参照配列のコドンとSNP後のコドンを比較するためにSNP後のコドンを導出する"""
        total_snp = pd.DataFrame()
        c_nums = [
            re.search(r"(?<=C)\d{1,2}", i).group() for i in base.columns if len(i.split("_")) < 2 and i[0] == "C"
        ]
        base = base.dropna(subset=["遺伝子アミノ酸位置"])
        params = [(base[["POS", "Reference Seq. (MN908947.3)", f"C{c_num}_ALT"]].values, c_num) for _, c_num in enumerate(c_nums)]
        with ProcessPoolExecutor(max_workers=len(c_nums)) as executor:
            # print(f"CLUSTER{c_num}")
            results = executor.map(self._extract_codon, params)

        for i, adding_snp in enumerate(results):
            ### 大丈夫?
            if i != 0:
                adding_snp = adding_snp.drop(columns=["POS", "Original", "Original_Amino"])
            total_snp = pd.concat([total_snp, adding_snp], axis=1)
        return total_snp

        # for i, c_num in enumerate(c_nums):
        #     print(f"CLUSTER{c_num}")
        #     adding_snp = self._extract_codon(base, c_num)

        #     ### 大丈夫?
        #     if i != 0:
        #         adding_snp = adding_snp.drop(columns=["POS", "Original", "Original_Amino"])
        #     total_snp = pd.concat([total_snp, adding_snp], axis=1)
        # return total_snp

    @staticmethod
    def _extract_codon(params):
        """コドンを導出する"""
        def _substitute_snv_base(Pos, codon, snv_codon, codon_pos):
            """SNV塩基を元の塩基と置換する"""
            list_codon = list(codon)
            list_snp_codon = list(snv_codon)
            list_codon[codon_pos] = list_snp_codon[codon_pos]
            amino = codon_table_v2[codon]
            snp_amino = codon_table_v2["".join(list_codon)]
            if amino == snp_amino:
                synonymous_codon = "同義コドン"
            else:
                synonymous_codon = "非同義コドン"

            return pd.DataFrame(
                [[Pos - (codon_pos - 2), codon, "".join(list_codon), amino, snp_amino, synonymous_codon,]]
            )

        base_values, c_num = params
        codon = ""
        snv_codon = ""
        adding_snv = []
        BASE = ["A", "T", "G", "C"]
        
        for Pos, RefSeq, CnumALT in tqdm(base_values):
            codon += RefSeq
            flag = True
            if CnumALT is not np.nan and type(CnumALT) != float:
                    for B in BASE:
                        if B in CnumALT:
                            snv_codon += B
                            flag = False
                            break
            if flag:
                snv_codon += RefSeq

            if len(snv_codon) >= 3:

                new_snv0 = _substitute_snv_base(Pos, codon, snv_codon, 0)
                new_snv1 = _substitute_snv_base(Pos, codon, snv_codon, 1)
                new_snv2 = _substitute_snv_base(Pos, codon, snv_codon, 2)
                adding_snv.append(new_snv0)
                adding_snv.append(new_snv1)
                adding_snv.append(new_snv2)

                codon = ""
                snv_codon = ""
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
