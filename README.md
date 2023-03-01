# Sars-CoV-2の地域流行株の検出に関する研究

## Omicron株2022年01月取得データ解析
* 248,574件のOmicron株のデータ
  * 2022/1月までのデータ
  * 5連続塩基のオッズ比
  * Nが1%未満


## フォルダ構成
```
.
├── cluster_definition # クラスタ定義クラス
├── datareading # データ読み込みクラス
├── omicron2201
│   ├── ReceivedData # 解析データ
│   ├── SNV # クラスタ別SNV *.al / *.outputファイル
│   ├── output # クラスタ定義結果
├── snv_detection # SNV検出に関するファイル群
└── utils # ログ出力用間数群
└── # mainファイル群
```

## 解析対象データの詳細とフォルダの形式統一
* omicron2201が解析対象データにあたる
* fasta/blsomのoutファイル/データMetaファイルの名称統一
* outputフォルダに生成物を格納
* ReceivedDataにデータを格納

## BLSOMのプロットと画像化
`220216_cluster_defination.ipynb`を上からパラメータを変えて実行
* 生成物
  * `omicron2201/output_kmeans3_1bycluster/images`以下に全大陸データプロット結果と大陸別データプロット結果の図が出力される

## クラスタ定義
### 類似度マニュアル指定
* 定義方法
  * 学習されたBLSOMに大陸別にデータを分類し，格子点中のデータ件数の多い上位n件をクラスタの中心点の候補とする.
  * 中心点の候補と，各大陸のデータが分類されている格子点との距離を計算し，0.5σ以下の距離かつ中心点から±5の範囲の格子点をクラスタとして定義する.このとき，重なりのあるクラスタは統合させる.

`define_cluster.py`をパラメータを変えて実行する。

```
nohup python define_cluster.py &
```

パラメータは以下を変更する。
```
MAX_CLUSTER_NUM = 10 # 10/50/100/200/300
DATA_DIR = "omicron2201" # omicron2201/omicron2211
```

### K-meansクラスタに各1つ / ±3
* 定義方法
  * BLSOM上の全格子点中のデータ件数の多い上位n件をk-meansクラスタの中心点とする.
  * BLSOMの全ての重みをk-means法によりn個のクラスタに分割する.
  * クラスタの最もデータ件数が多い格子点から±3の範囲の格子点をクラスタとして定義する.


`define_cluster_kmeans.py`をパラメータを変えて実行する。

```
nohup python define_cluster_kmeans.py &
```

パラメータは以下を変更する。
```
MAX_CLUSTER_NUM = 10 # 10/50/100/200/300
DATA_DIR = "omicron2201" # omicron2201/omicron2211
```

### K-meansクラスタ横断 / ±3
* 定義方法
  * BLSOM上の全格子点中のデータ件数の多い上位n件をk-meansクラスタの中心点とする.
  * BLSOMの全ての重みをk-means法によりn個のk-meansクラスタに分割する.
  * K-meansクラスタに依らずにデータ件数が多い格子点m個をとり，これらを中心に±3の範囲の格子点をクラスタとして定義する.このとき，k-meansクラスタを跨ぐことを禁止し，重なりのあるクラスタは統合する.

`define_cluster_kmeans_total.py`をパラメータを変えて実行する。

```
nohup python define_cluster_kmeans_total.py &
```

パラメータは以下を変更する。
```
MAX_CLUSTER_NUM = 63000
KMEANS_CLUSTER_NUM = 10 # 10/50/100/200/300
DATA_DIR = "omicron2201" # omicron2201/omicron2211
```

### クラスタ定義結果の集計ファイルの出力
`arrange_data.ipynb`を上から実行する

* 生成物
  * omicron2201/output以下にクラスタ定義結果が出力される
  * クラスタ定義結果の集計ファイル(Excel)が出力される

## SNV解析
`snv_detection`フォルダに移動し，以下を実行する
```
sh split_cluster.sh
sh detection_snp.sh
python snv2excel.py
```

* 生成物
  * `omicron2201/output/10/SNV_EXCEL`以下にexcelファイルが出力される




<!-- ## omicron2211解析
### データ
* 2,414,620件のOmicron株のデータ
  * 2021/12月のデータ
  * 5連続塩基のオッズ比
  * Nが1%未満 -->