#!/usr/bin/env bash
DIR=PROJECTS/SRP199742_GEP000010
mkdir ${DIR}
cd ${DIR}
echo 'Data will be downloaded into' ${DIR}
set -x
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR914/000/SRR9141470/SRR9141470_1.fastq.gz -o SLX-15025.FLD0218.000000000-BNNFL.s_1.r_1.fq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR914/000/SRR9141470/SRR9141470_2.fastq.gz -o SLX-15025.FLD0218.000000000-BNNFL.s_1.r_2.fq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR914/007/SRR9141467/SRR9141467_1.fastq.gz -o SLX-15025.FLD0193.000000000-BNNFL.s_1.r_1.fq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR914/007/SRR9141467/SRR9141467_2.fastq.gz -o SLX-15025.FLD0193.000000000-BNNFL.s_1.r_2.fq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR914/009/SRR9141469/SRR9141469_1.fastq.gz -o SLX-15025.FLD0217.000000000-BNNFL.s_1.r_1.fq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR914/009/SRR9141469/SRR9141469_2.fastq.gz -o SLX-15025.FLD0217.000000000-BNNFL.s_1.r_2.fq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR914/008/SRR9141468/SRR9141468_1.fastq.gz -o SLX-15025.FLD0194.000000000-BNNFL.s_1.r_1.fq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR914/008/SRR9141468/SRR9141468_2.fastq.gz -o SLX-15025.FLD0194.000000000-BNNFL.s_1.r_2.fq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR914/006/SRR9141466/SRR9141466_1.fastq.gz -o SLX-15025.FLD0341.000000000-BNNFL.s_1.r_1.fq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR914/006/SRR9141466/SRR9141466_2.fastq.gz -o SLX-15025.FLD0341.000000000-BNNFL.s_1.r_2.fq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR914/004/SRR9141454/SRR9141454_1.fastq.gz -o SLX-15025.FLD0330.000000000-BNNFL.s_1.r_1.fq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR914/004/SRR9141454/SRR9141454_2.fastq.gz -o SLX-15025.FLD0330.000000000-BNNFL.s_1.r_2.fq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR914/005/SRR9141455/SRR9141455_1.fastq.gz -o SLX-15025.FLD0329.000000000-BNNFL.s_1.r_1.fq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR914/005/SRR9141455/SRR9141455_2.fastq.gz -o SLX-15025.FLD0329.000000000-BNNFL.s_1.r_2.fq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR914/001/SRR9141451/SRR9141451_1.fastq.gz -o SLX-15025.FLD0356.000000000-BNNFL.s_1.r_1.fq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR914/001/SRR9141451/SRR9141451_2.fastq.gz -o SLX-15025.FLD0356.000000000-BNNFL.s_1.r_2.fq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR914/000/SRR9141450/SRR9141450_1.fastq.gz -o SLX-15025.FLD0344.000000000-BNNFL.s_1.r_1.fq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR914/000/SRR9141450/SRR9141450_2.fastq.gz -o SLX-15025.FLD0344.000000000-BNNFL.s_1.r_2.fq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR914/009/SRR9141449/SRR9141449_1.fastq.gz -o SLX-15025.FLD0207.000000000-BNNFL.s_1.r_1.fq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR914/009/SRR9141449/SRR9141449_2.fastq.gz -o SLX-15025.FLD0207.000000000-BNNFL.s_1.r_2.fq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR914/005/SRR9141445/SRR9141445_1.fastq.gz -o SLX-15025.FLD0355.000000000-BNNFL.s_1.r_1.fq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR914/005/SRR9141445/SRR9141445_2.fastq.gz -o SLX-15025.FLD0355.000000000-BNNFL.s_1.r_2.fq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR914/008/SRR9141448/SRR9141448_1.fastq.gz -o SLX-15025.FLD0208.000000000-BNNFL.s_1.r_1.fq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR914/008/SRR9141448/SRR9141448_2.fastq.gz -o SLX-15025.FLD0208.000000000-BNNFL.s_1.r_2.fq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR914/006/SRR9141446/SRR9141446_1.fastq.gz -o SLX-15025.FLD0205.000000000-BNNFL.s_1.r_1.fq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR914/006/SRR9141446/SRR9141446_2.fastq.gz -o SLX-15025.FLD0205.000000000-BNNFL.s_1.r_2.fq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR914/002/SRR9141442/SRR9141442_1.fastq.gz -o SLX-15025.FLD0219.000000000-BNNFL.s_1.r_1.fq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR914/002/SRR9141442/SRR9141442_2.fastq.gz -o SLX-15025.FLD0219.000000000-BNNFL.s_1.r_2.fq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR914/004/SRR9141444/SRR9141444_1.fastq.gz -o SLX-15025.FLD0206.000000000-BNNFL.s_1.r_1.fq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR914/004/SRR9141444/SRR9141444_2.fastq.gz -o SLX-15025.FLD0206.000000000-BNNFL.s_1.r_2.fq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR914/003/SRR9141443/SRR9141443_1.fastq.gz -o SLX-15025.FLD0343.000000000-BNNFL.s_1.r_1.fq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR914/003/SRR9141443/SRR9141443_2.fastq.gz -o SLX-15025.FLD0343.000000000-BNNFL.s_1.r_2.fq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR914/009/SRR9141439/SRR9141439_1.fastq.gz -o SLX-15025.FLD0332.000000000-BNNFL.s_1.r_1.fq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR914/009/SRR9141439/SRR9141439_2.fastq.gz -o SLX-15025.FLD0332.000000000-BNNFL.s_1.r_2.fq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR914/008/SRR9141438/SRR9141438_1.fastq.gz -o SLX-15025.FLD0331.000000000-BNNFL.s_1.r_1.fq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR914/008/SRR9141438/SRR9141438_2.fastq.gz -o SLX-15025.FLD0331.000000000-BNNFL.s_1.r_2.fq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR914/007/SRR9141437/SRR9141437_1.fastq.gz -o SLX-15025.FLD0196.000000000-BNNFL.s_1.r_1.fq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR914/007/SRR9141437/SRR9141437_2.fastq.gz -o SLX-15025.FLD0196.000000000-BNNFL.s_1.r_2.fq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR914/000/SRR9141440/SRR9141440_1.fastq.gz -o SLX-15025.FLD0220.000000000-BNNFL.s_1.r_1.fq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR914/000/SRR9141440/SRR9141440_2.fastq.gz -o SLX-15025.FLD0220.000000000-BNNFL.s_1.r_2.fq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR914/006/SRR9141436/SRR9141436_1.fastq.gz -o SLX-15025.FLD0195.000000000-BNNFL.s_1.r_1.fq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR914/006/SRR9141436/SRR9141436_2.fastq.gz -o SLX-15025.FLD0195.000000000-BNNFL.s_1.r_2.fq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR914/002/SRR9141432/SRR9141432_1.fastq.gz -o SLX-15025.FLD0353.000000000-BNNFL.s_1.r_1.fq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR914/002/SRR9141432/SRR9141432_2.fastq.gz -o SLX-15025.FLD0353.000000000-BNNFL.s_1.r_2.fq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR914/003/SRR9141433/SRR9141433_1.fastq.gz -o SLX-15025.FLD0354.000000000-BNNFL.s_1.r_1.fq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR914/003/SRR9141433/SRR9141433_2.fastq.gz -o SLX-15025.FLD0354.000000000-BNNFL.s_1.r_2.fq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR914/007/SRR9141417/SRR9141417_1.fastq.gz -o SLX-15025.FLD0342.000000000-BNNFL.s_1.r_1.fq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR914/007/SRR9141417/SRR9141417_2.fastq.gz -o SLX-15025.FLD0342.000000000-BNNFL.s_1.r_2.fq.gz
set +x
echo 'Done'
