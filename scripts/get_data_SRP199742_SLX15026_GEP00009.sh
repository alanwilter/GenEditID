#!/usr/bin/env bash
DIR=PROJECTS/SRP199742_GEP00009
mkdir ${DIR}
cd ${DIR}
echo 'Data will be downloaded into' ${DIR}
set -x
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR914/001/SRR9141471/SRR9141471_1.fastq.gz -o SLX-15026.FLD0319.000000000-BJWVR.s_1.r_1.fq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR914/001/SRR9141471/SRR9141471_2.fastq.gz -o SLX-15026.FLD0319.000000000-BJWVR.s_1.r_2.fq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR914/003/SRR9141463/SRR9141463_1.fastq.gz -o SLX-15026.FLD0313.000000000-BJWVR.s_1.r_1.fq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR914/003/SRR9141463/SRR9141463_2.fastq.gz -o SLX-15026.FLD0313.000000000-BJWVR.s_1.r_2.fq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR914/005/SRR9141465/SRR9141465_1.fastq.gz -o SLX-15026.FLD0324.000000000-BJWVR.s_1.r_1.fq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR914/005/SRR9141465/SRR9141465_2.fastq.gz -o SLX-15026.FLD0324.000000000-BJWVR.s_1.r_2.fq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR914/004/SRR9141464/SRR9141464_1.fastq.gz -o SLX-15026.FLD0323.000000000-BJWVR.s_1.r_1.fq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR914/004/SRR9141464/SRR9141464_2.fastq.gz -o SLX-15026.FLD0323.000000000-BJWVR.s_1.r_2.fq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR914/002/SRR9141462/SRR9141462_1.fastq.gz -o SLX-15026.FLD0314.000000000-BJWVR.s_1.r_1.fq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR914/002/SRR9141462/SRR9141462_2.fastq.gz -o SLX-15026.FLD0314.000000000-BJWVR.s_1.r_2.fq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR914/000/SRR9141460/SRR9141460_1.fastq.gz -o SLX-15026.FLD0316.000000000-BJWVR.s_1.r_1.fq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR914/000/SRR9141460/SRR9141460_2.fastq.gz -o SLX-15026.FLD0316.000000000-BJWVR.s_1.r_2.fq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR914/001/SRR9141461/SRR9141461_1.fastq.gz -o SLX-15026.FLD0315.000000000-BJWVR.s_1.r_1.fq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR914/001/SRR9141461/SRR9141461_2.fastq.gz -o SLX-15026.FLD0315.000000000-BJWVR.s_1.r_2.fq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR914/009/SRR9141459/SRR9141459_1.fastq.gz -o SLX-15026.FLD0309.000000000-BJWVR.s_1.r_1.fq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR914/009/SRR9141459/SRR9141459_2.fastq.gz -o SLX-15026.FLD0309.000000000-BJWVR.s_1.r_2.fq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR914/008/SRR9141458/SRR9141458_1.fastq.gz -o SLX-15026.FLD0310.000000000-BJWVR.s_1.r_1.fq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR914/008/SRR9141458/SRR9141458_2.fastq.gz -o SLX-15026.FLD0310.000000000-BJWVR.s_1.r_2.fq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR914/007/SRR9141457/SRR9141457_1.fastq.gz -o SLX-15026.FLD0311.000000000-BJWVR.s_1.r_1.fq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR914/007/SRR9141457/SRR9141457_2.fastq.gz -o SLX-15026.FLD0311.000000000-BJWVR.s_1.r_2.fq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR914/006/SRR9141456/SRR9141456_1.fastq.gz -o SLX-15026.FLD0312.000000000-BJWVR.s_1.r_1.fq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR914/006/SRR9141456/SRR9141456_2.fastq.gz -o SLX-15026.FLD0312.000000000-BJWVR.s_1.r_2.fq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR914/002/SRR9141452/SRR9141452_1.fastq.gz -o SLX-15026.FLD0318.000000000-BJWVR.s_1.r_1.fq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR914/002/SRR9141452/SRR9141452_2.fastq.gz -o SLX-15026.FLD0318.000000000-BJWVR.s_1.r_2.fq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR914/003/SRR9141453/SRR9141453_1.fastq.gz -o SLX-15026.FLD0317.000000000-BJWVR.s_1.r_1.fq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR914/003/SRR9141453/SRR9141453_2.fastq.gz -o SLX-15026.FLD0317.000000000-BJWVR.s_1.r_2.fq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR914/007/SRR9141447/SRR9141447_1.fastq.gz -o SLX-15026.FLD0320.000000000-BJWVR.s_1.r_1.fq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR914/007/SRR9141447/SRR9141447_2.fastq.gz -o SLX-15026.FLD0320.000000000-BJWVR.s_1.r_2.fq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR914/001/SRR9141441/SRR9141441_1.fastq.gz -o SLX-15026.FLD0296.000000000-BJWVR.s_1.r_1.fq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR914/001/SRR9141441/SRR9141441_2.fastq.gz -o SLX-15026.FLD0296.000000000-BJWVR.s_1.r_2.fq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR914/004/SRR9141434/SRR9141434_1.fastq.gz -o SLX-15026.FLD0307.000000000-BJWVR.s_1.r_1.fq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR914/004/SRR9141434/SRR9141434_2.fastq.gz -o SLX-15026.FLD0307.000000000-BJWVR.s_1.r_2.fq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR914/005/SRR9141435/SRR9141435_1.fastq.gz -o SLX-15026.FLD0308.000000000-BJWVR.s_1.r_1.fq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR914/005/SRR9141435/SRR9141435_2.fastq.gz -o SLX-15026.FLD0308.000000000-BJWVR.s_1.r_2.fq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR914/001/SRR9141431/SRR9141431_1.fastq.gz -o SLX-15026.FLD0304.000000000-BJWVR.s_1.r_1.fq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR914/001/SRR9141431/SRR9141431_2.fastq.gz -o SLX-15026.FLD0304.000000000-BJWVR.s_1.r_2.fq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR914/009/SRR9141429/SRR9141429_1.fastq.gz -o SLX-15026.FLD0306.000000000-BJWVR.s_1.r_1.fq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR914/009/SRR9141429/SRR9141429_2.fastq.gz -o SLX-15026.FLD0306.000000000-BJWVR.s_1.r_2.fq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR914/008/SRR9141428/SRR9141428_1.fastq.gz -o SLX-15026.FLD0305.000000000-BJWVR.s_1.r_1.fq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR914/008/SRR9141428/SRR9141428_2.fastq.gz -o SLX-15026.FLD0305.000000000-BJWVR.s_1.r_2.fq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR914/000/SRR9141430/SRR9141430_1.fastq.gz -o SLX-15026.FLD0303.000000000-BJWVR.s_1.r_1.fq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR914/000/SRR9141430/SRR9141430_2.fastq.gz -o SLX-15026.FLD0303.000000000-BJWVR.s_1.r_2.fq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR914/007/SRR9141427/SRR9141427_1.fastq.gz -o SLX-15026.FLD0300.000000000-BJWVR.s_1.r_1.fq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR914/007/SRR9141427/SRR9141427_2.fastq.gz -o SLX-15026.FLD0300.000000000-BJWVR.s_1.r_2.fq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR914/005/SRR9141425/SRR9141425_1.fastq.gz -o SLX-15026.FLD0302.000000000-BJWVR.s_1.r_1.fq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR914/005/SRR9141425/SRR9141425_2.fastq.gz -o SLX-15026.FLD0302.000000000-BJWVR.s_1.r_2.fq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR914/006/SRR9141426/SRR9141426_1.fastq.gz -o SLX-15026.FLD0299.000000000-BJWVR.s_1.r_1.fq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR914/006/SRR9141426/SRR9141426_2.fastq.gz -o SLX-15026.FLD0299.000000000-BJWVR.s_1.r_2.fq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR914/003/SRR9141423/SRR9141423_1.fastq.gz -o SLX-15026.FLD0322.000000000-BJWVR.s_1.r_1.fq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR914/003/SRR9141423/SRR9141423_2.fastq.gz -o SLX-15026.FLD0322.000000000-BJWVR.s_1.r_2.fq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR914/004/SRR9141424/SRR9141424_1.fastq.gz -o SLX-15026.FLD0301.000000000-BJWVR.s_1.r_1.fq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR914/004/SRR9141424/SRR9141424_2.fastq.gz -o SLX-15026.FLD0301.000000000-BJWVR.s_1.r_2.fq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR914/000/SRR9141420/SRR9141420_1.fastq.gz -o SLX-15026.FLD0321.000000000-BJWVR.s_1.r_1.fq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR914/000/SRR9141420/SRR9141420_2.fastq.gz -o SLX-15026.FLD0321.000000000-BJWVR.s_1.r_2.fq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR914/009/SRR9141419/SRR9141419_1.fastq.gz -o SLX-15026.FLD0294.000000000-BJWVR.s_1.r_1.fq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR914/009/SRR9141419/SRR9141419_2.fastq.gz -o SLX-15026.FLD0294.000000000-BJWVR.s_1.r_2.fq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR914/008/SRR9141418/SRR9141418_1.fastq.gz -o SLX-15026.FLD0293.000000000-BJWVR.s_1.r_1.fq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR914/008/SRR9141418/SRR9141418_2.fastq.gz -o SLX-15026.FLD0293.000000000-BJWVR.s_1.r_2.fq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR914/002/SRR9141422/SRR9141422_1.fastq.gz -o SLX-15026.FLD0298.000000000-BJWVR.s_1.r_1.fq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR914/002/SRR9141422/SRR9141422_2.fastq.gz -o SLX-15026.FLD0298.000000000-BJWVR.s_1.r_2.fq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR914/001/SRR9141421/SRR9141421_1.fastq.gz -o SLX-15026.FLD0297.000000000-BJWVR.s_1.r_1.fq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR914/001/SRR9141421/SRR9141421_2.fastq.gz -o SLX-15026.FLD0297.000000000-BJWVR.s_1.r_2.fq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR914/006/SRR9141416/SRR9141416_1.fastq.gz -o SLX-15026.FLD0295.000000000-BJWVR.s_1.r_1.fq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR914/006/SRR9141416/SRR9141416_2.fastq.gz -o SLX-15026.FLD0295.000000000-BJWVR.s_1.r_2.fq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR914/005/SRR9141415/SRR9141415_1.fastq.gz -o SLX-15026.FLD0290.000000000-BJWVR.s_1.r_1.fq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR914/005/SRR9141415/SRR9141415_2.fastq.gz -o SLX-15026.FLD0290.000000000-BJWVR.s_1.r_2.fq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR914/004/SRR9141414/SRR9141414_1.fastq.gz -o SLX-15026.FLD0289.000000000-BJWVR.s_1.r_1.fq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR914/004/SRR9141414/SRR9141414_2.fastq.gz -o SLX-15026.FLD0289.000000000-BJWVR.s_1.r_2.fq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR914/003/SRR9141413/SRR9141413_1.fastq.gz -o SLX-15026.FLD0292.000000000-BJWVR.s_1.r_1.fq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR914/003/SRR9141413/SRR9141413_2.fastq.gz -o SLX-15026.FLD0292.000000000-BJWVR.s_1.r_2.fq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR914/002/SRR9141412/SRR9141412_1.fastq.gz -o SLX-15026.FLD0291.000000000-BJWVR.s_1.r_1.fq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR914/002/SRR9141412/SRR9141412_2.fastq.gz -o SLX-15026.FLD0291.000000000-BJWVR.s_1.r_2.fq.gz
set +x
echo 'Done'
