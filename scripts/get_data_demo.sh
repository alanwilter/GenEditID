#!/usr/bin/env bash
DIR=PROJECTS/demo
mkdir ${DIR}
cd ${DIR}
echo 'Data will be downloaded into' ${DIR}
set -x
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR909/005/SRR9091245/SRR9091245_1.fastq.gz -o SLX-15021.FLD0010.000000000-BJ8JD.s_1.r_1.fq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR909/005/SRR9091245/SRR9091245_2.fastq.gz -o SLX-15021.FLD0010.000000000-BJ8JD.s_1.r_2.fq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR909/002/SRR9091202/SRR9091202_1.fastq.gz -o SLX-15021.FLD0011.000000000-BJ8JD.s_1.r_1.fq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR909/002/SRR9091202/SRR9091202_2.fastq.gz -o SLX-15021.FLD0011.000000000-BJ8JD.s_1.r_2.fq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR909/009/SRR9091169/SRR9091169_1.fastq.gz -o SLX-15021.FLD0054.000000000-BJ8JD.s_1.r_1.fq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR909/009/SRR9091169/SRR9091169_2.fastq.gz -o SLX-15021.FLD0054.000000000-BJ8JD.s_1.r_2.fq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR909/000/SRR9091230/SRR9091230_1.fastq.gz -o SLX-15021.FLD0075.000000000-BJ8JD.s_1.r_1.fq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR909/000/SRR9091230/SRR9091230_2.fastq.gz -o SLX-15021.FLD0075.000000000-BJ8JD.s_1.r_2.fq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR909/001/SRR9091221/SRR9091221_1.fastq.gz -o SLX-15021.FLD0081.000000000-BJ8JD.s_1.r_1.fq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR909/001/SRR9091221/SRR9091221_2.fastq.gz -o SLX-15021.FLD0081.000000000-BJ8JD.s_1.r_2.fq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR909/005/SRR9091225/SRR9091225_1.fastq.gz -o SLX-15021.FLD0085.000000000-BJ8JD.s_1.r_1.fq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR909/005/SRR9091225/SRR9091225_2.fastq.gz -o SLX-15021.FLD0085.000000000-BJ8JD.s_1.r_2.fq.gz
set +x
echo 'Done'
