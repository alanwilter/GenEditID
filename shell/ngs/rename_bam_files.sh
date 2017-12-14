mkdir idbam
cd idbam
for bam in ../bam/*.bam
do
    file=$(basename $bam)
    barcode=$(echo $file | cut -d '.' -f 2).bam
    ln -s $bam $barcode
done
for bai in ../bam/*.bai
do
    file=$(basename $bai)
    barcode=$(echo $file | cut -d '.' -f 2).bai
    ln -s $bai $barcode
done
