# Sars-CoV-2の地域流行株の検出に関する研究

## 解析の流れ
### 解析対象データの詳細とフォルダの形式統一
* omicron2201とomicron2211が解析対象データにあたる
* fasta/blsomのoutファイル/データMetaファイルの名称統一
* outputフォルダに生成物を格納
* ReceivedDataにデータを格納

### BLSOMのプロットと画像化
`220216_cluster_defination.ipynb`をパラメータを変えて実行

### クラスタ定義
`define_cluster.py`をパラメータを変えて実行する。

```
nohup python define_cluster.py &
```

パラメータは以下を変更する。
```python:define_cluster.py
MAX_CLUSTER_NUM = 10 # 10/50/100/200/300
DATA_DIR = "omicron2201" # omicron2201/omicron2211
```

### SNV解析
`snv_detection`フォルダ以下のプログラムを実行する。


## omicron2204解析
### データ
* 248,574件のOmicron株のデータ
  * 2021/12月のデータ
  * 5連続塩基のオッズ比
  * Nが1%未満


## omicron2211解析
### データ
* 2,414,620件のOmicron株のデータ
  * 2021/12月のデータ
  * 5連続塩基のオッズ比
  * Nが1%未満