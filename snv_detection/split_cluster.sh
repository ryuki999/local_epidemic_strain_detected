for i in Africa Asia Oceania North_America Europe
do
    python split_cluster_fas.py omicron2201 10 ${i}
    # python split_cluster_fas.py omicron2201 50 ${i}
    # python split_cluster_fas.py omicron2201 100 ${i}

    # python split_cluster_fas.py omicron2211 10 ${i}
    # python split_cluster_fas.py omicron2211 50 ${i}
    # python split_cluster_fas.py omicron2211 100 ${i}
done