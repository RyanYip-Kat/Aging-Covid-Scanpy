scATACP=/home/ye/Software/scATACP/scATAC-pro_1.1.1/scATAC-pro
output=$1
path=$2
config="configure_user.txt"
fragment="$path/outs/fragments.tsv.gz"
peak="$path/outs/peaks.bed"
bam="$path/outs/possorted_bam.bam"
# -------------------------
echo "### : convert10xbam"
$scATACP -s convert10xbam -i $bam -c $config -o $output

echo "### : call peak"
$scATACP -s call_peak  -i "$output/mapping_result/pbmc10k.positionsort.MAPQ30.bam" -c $config -o $output

echo "### : aggr_signal"
$scATACP -s aggr_signal -i "$output/mapping_result/pbmc10k.positionsort.MAPQ30.bam" -c $config -o $output

echo "### : get mtx"
$scATACP -s get_mtx  -i $fragment,$peak  -o $output -c $config

echo "### :qc_per_barcode"
$scATACP -s qc_per_barcode -i $fragment,$peak  -o $output -c $config

echo "### : call cell"
$scATACP -s call_cell -i "$output/raw_matrix/MACS2/matrix.mtx"  -o $output -c $config
# -------------------------

$scATACP -s get_bam4Cells  -i "$output/mapping_result/pbmc10k.positionsort.bam","$output/filtered_matrix/MACS2/FILTER/barcodes.txt" -c $config -o $output
