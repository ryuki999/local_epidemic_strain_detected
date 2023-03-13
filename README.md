# Sars-CoV-2の地域流行株の検出に関する研究
解析は以下の流れで行う.
* データの取得
* BLSOM解析と可視化結果の出力
* クラスタ定義と生成クラスタ集計ファイルの出力
* クラスタ別SNV集計ファイルの出力

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
│   ├── output # ユーザが確認する解析結果
├── snv_detection # SNV検出に関するファイル群
└── utils # ログ出力用間数群
└── # mainファイル群
```

## 解析対象データの詳細とフォルダの形式統一
* omicron2201が解析対象データにあたる(ReceivedData)
* fasta/blsomのoutファイル/データMetaファイルのフォーマット統一
* ユーザが確認する解析結果はすべてoutputフォルダに格納

## BLSOMのプロットと画像化
`220216_cluster_defination.ipynb`を上から，必要に応じてパラメータを変えて実行
* 生成物
  * `omicron2201/output/images`以下に全大陸データプロット結果と大陸別データプロット結果の図が出力される

## クラスタ定義
### BLSOMの重みを手動閾値でフィルタリング / ±3
* 定義方法
  * 学習されたBLSOMに大陸別にデータを分類し，格子点中のデータ件数の多い上位n件をクラスタの中心点の候補とする.
  * 中心点の候補と，各大陸のデータが分類されている格子点との距離を計算し，0.5σ以下の距離かつ中心点から±5の範囲の格子点をクラスタとして定義する.このとき，重なりのあるクラスタは統合させる.

`define_cluster.py`をパラメータを変えて実行する。

```
nohup python define_cluster.py &
```

パラメータは以下を変更する。
```
MAX_CLUSTER_NUM = 50 # 10/50/100/200/300
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
MAX_CLUSTER_NUM = 50 # 10/50/100/200/300
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
`cluster_analysis.ipynb`を上から実行する

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

このとき，`split_cluster.sh`ファイル内のpython scriptを実行する部分で，第一引数に解析対象データのあるフォルダ(omicron2201)を指定し，第二引数にクラスタ定義数(50)をする.また，第三引数では大陸名が指定される.
```
python split_cluster_fas.py omicron2201 50 ${i}
```

`detection_snp.sh`では，解析対象データのフォルダ(DATA_DIR)を指定し，クラスタ定義数(MAX_CLUSTER_NUM)を指定する.また，これは`snv2excel.py`でも同様の指定を行う.
```
DATA_DIR=omicron2201
MAX_CLUSTER_NUM=50
```

* 生成物
  * `omicron2201/output/50/SNV_EXCEL`以下にexcelファイルが出力される




<!-- ## omicron2211解析
### データ
* 2,414,620件のOmicron株のデータ
  * 2021/12月のデータ
  * 5連続塩基のオッズ比
  * Nが1%未満 -->