
DATA_DIR=omicron2201
MAX_CLUSTER_NUM=10
LOG_FILE=../${DATA_DIR}/dec_snv_${MAX_CLUSTER_NUM}.log
# exec 1> >(awk '{print strftime("${DATA_DIR} CONTINENT: ${i} MAX_CLUSTER_NUM: ${MAX_CLUSTER_NUM}  } { fflush() } ' >>${LOG_FILE})
# exec 2> >(awk '{print strftime("[%Y-%m-%d %H:%M:%S]"),$0 } { fflush() } ' >>${LOG_FILE})

for i in Africa Asia Oceania North_America Europe
do
echo "${DATA_DIR} CONTINENT: ${i} MAX_CLUSTER_NUM: ${MAX_CLUSTER_NUM} Start!" | awk '{print strftime("[%Y-%m-%d %H:%M:%S]"), $0}' >>${LOG_FILE}
	for f in ../${DATA_DIR}/SNV/${MAX_CLUSTER_NUM}/CLUSTER/${i}/*.fas
	do
		# ファイルパスからディレクトリ名を抜いてファイル名を取得
		# f1=${f##*/}
		# ファイルパスから拡張子を抜いて取得
		f1=${f%.*}
		echo "${f1} MAFFT" | awk '{print strftime("[%Y-%m-%d %H:%M:%S]"), $0}' >> ${LOG_FILE}
		mafft --6merpair --thread -1 --keeplength  --quiet --addfragments ${f1}.fas MN908947_3.fasta > ${f1}.fas.al
		wait
		echo "${f1} SNP-SITES" | awk '{print strftime("[%Y-%m-%d %H:%M:%S]"), $0}' >> ${LOG_FILE}
		snp-sites -v -o ${f1}.output ${f1}.fas.al
	done
echo "${DATA_DIR} CONTINENT: ${i} MAX_CLUSTER_NUM: ${MAX_CLUSTER_NUM} Done!" | awk '{print strftime("[%Y-%m-%d %H:%M:%S]"), $0}' >>${LOG_FILE}
echo "" >> ${LOG_FILE}
done