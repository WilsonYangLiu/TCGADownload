# Process `gencode.v22.annotation.csv.gz`, which was generated from `gencode.v22.annotation.gtf.gz` by using `Client.gtf.gtfClient` meathod 
# Extract gene
gzip -cd ../data/gencode.v22.annotation.csv.gz | awk '
BEGIN{FS=",";OFS=",";print "#gene_id,gene_name,Length"};
(NR>1&&$1=="gene"){
	print $4,$5,$3-$2+1
}' > ../data/gencode.v22.annotation.used4FPKM.csv

# Extract exon
gzip -cd ../data/gencode.v22.annotation.csv.gz | awk 'BEGIN{FS=",";OFS=","};(NR>1&&$1=="exon"){print $4,$5,$3-$2+1}' | sort | awk '
BEGIN{FS=",";OFS=",";gene="ENSG00000000003.13";name="TSPAN6";total=0;print "#gene_id,gene_name,Length"};
{
	if ($1==gene) {
		total=total+$3;}
	else {
		print gene,name,total;
		gene=$1;
		name=$2;
		total=$3;}
}
END{print gene,name,total}' >  ../data/gencode.v22.annotation.used4FPKM.exon.csv
