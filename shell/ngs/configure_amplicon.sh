echo 'Creating pipeline directories'
mkdir amplicon_vardict
mkdir amplicon_gatk

echo "Setting project $1 into config files"
sed -i "s/GEPID/$1/g" config.vardict.xml
sed -i "s/GEPID/$1/g" config.gatk.xml

echo 'Making BAMs match the IDs'
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
