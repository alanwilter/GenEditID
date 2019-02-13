# NGS pipeline

## Sanger pipe line

* Paper : https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-016-2563-z

* Github link: https://github.com/richysix/Crispr

* Aligner used - BBmap

* Pipeline is only useful for calling Indels (no SNVs)

* Dindel is used for indel calling

* Dindel link : https://sites.google.com/site/keesalbers/soft/dindel

* Last update of Dindel is 2010

* We find it difficult to install Dindel

* Rich was able to install Dindel, but we were never able to run it.

Therefore we decided not to proceed with sanger pipeline.


## Alignment against amplicons or whole genome?

* Whole genome
    *   slower
    *   Can remove reads from unidentified off target regions

* Amplicon
    *   Faster
    *   Needs to create references every time
    *   Needs to map the amplicon and genomics coordinates
    *   Forces to reds to map, even if read is not from that region


We decided initially we use both options and compare.
