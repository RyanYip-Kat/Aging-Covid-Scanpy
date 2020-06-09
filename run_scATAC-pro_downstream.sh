scATACP=/home/ye/Software/scATACP/scATAC-pro_1.1.1/scATAC-pro
config="configure_user.txt"

output=$1
mtx=$2  #"$path/matrix.mtx"

echo "### Clustering"
$scATACP -s clustering -i $mtx  -o $output -c $config

echo "### Motif analysis"
$scATACP -s motif_analysis  -i $mtx  -o $output -c $config

echo "### runCicero"
seurat_obj="$output/downstream_analysis/MACS2/FILTER/seurat_obj.rds"
$scATACP -s runCicero -i $seurat_obj -o $output -c $config

#echo "### footprint"
#$scATACP -s footprint  -i 0,1  -o $output -c $config

#echo "### runDA"
#$scATACP -s runDA  -i 0,1   -o $output -c $config

#$scATACP -s runGO  -i "$output/downstream_analysis/MACS2/FILTER/differential_accessible_features_0_vs_1.txt"   -o $output -c $config


#scATAC-pro -s visualize -i output/downstream_analysis/PEAK_CALLER/CELL_CALLER/VisCello_obj -c configure_user.txt

