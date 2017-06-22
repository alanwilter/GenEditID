# AmpliconSeq Analysis Pipeline Package

author: Matt Eldridge, Cancer Research UK Cambridge Institute, University of Cambridge

## Introduction

AmpliconSeq is an analysis pipeline package that implements a workflow for
calling single nucleotide variants (SNVs) and indels in targeted amplicon
sequencing data. The pipeline offers a choice of variant callers, including
GATK HaplotypeCaller, GATK UnifiedGenotyper, FreeBayes and VarDict, and carries
out post-calling filtering to select high-confidence variants. In addition,
alignment and target coverage metrics are computed and compiled into a summary
report. Variants are annotated using Ensembl Variant Effect Predictor with
details about consequence on protein coding, accession numbers for known
variants and associated minor allele frequencies from the 1000 Genomes project.

The AmpliconSeq pipeline has the following features:

* Choice of variant callers: GATK HaplotypeCaller, GATK UnifiedGenotyper, FreeBayes and VarDict
* Alignment and coverage QC report using metrics calculated by Picard tools
* Annotation of variants using Ensembl Variant Effect Predictor
* Support for overlapping amplicon targets by partitioning reads prior to variant calling
* Support for calling and filtering low allele fraction SNVs, e.g. for circulating tumour
DNA in plasma samples with allele fractions down to 0.1%, by fitting probability distributions
to model background noise
* *Specific calling* of known mutations
* Support for Slurm and Platform LSF job schedulers for running on a high-performance computing cluster
* Accompanying visualization tool for viewing and assessing SNV calls

The AmpliconSeq pipeline was developed in collaboration with James Brenton's
group at Cancer Research UK Cambridge Institute, particularly with the support
of Anna Piskorz, Tedi Goranova and Geoff McIntyre. The figure below is a
schematic of the AmpliconSeq workflow; note that the pipeline does not peform
alignment of sequences to a reference genome but this is included in the diagram
to indicate that this is a necessary step in the bioinformatic analysis prior
to running the pipeline.

The workflow is executed using the workflow engine developed by Richard Bowers
at the Cancer Research UK Cambridge Institute (CRUK-CI) and distributes the
computational workload across multiple processors either on a powerful
multi-processor server or on a high-performance compute cluster.

The workflow system is a Java application that separates the mechanics of
running sequences of dependent jobs from the specifics of a given pipeline
specification. It offers a number of advantages for constructing pipelines
including the ability to restart pipelines that fail during execution without
re-running tasks that completed successfully and allowing pipelines to run on
either a single multi-processor machine or a compute cluster by simply changing
a run-time option. Job scheduling systems supported by the workflow engine
include Platform LSF, TORQUE and SLURM. The workflow system also handles memory
allocation for jobs, including setting the maximum heap size for Java programs,
e.g. Picard and GATK, and can attempt to re-run jobs that fail through having
insufficient memory by incrementally increasing the memory allocation up to a
specified maximum value.

The workflow system defines XML schema for task definitions, pipeline
construction (describing how the tasks are wired together) and also for the
configuration file used for a specific workflow instance or run. Task and
pipeline definitions are stored in the tasks and pipelines subdirectories
respectively within the AmpliconSeq pipeline installation directory and should
not normally be modified. Sample configuration files for the GATK, FreeBayes and
VarDict variant callers can be found within the config subdirectory.

-----------

## Installation

### Pre-requisites

The following software needs to be installed:

* Unix Bourne and Bash shells
* Java SE 8 or above
* Genome Analysis Toolkit (GATK) version 3.6 or above
* FreeBayes (optional)
* VarDict (Java version only, optional) version 1.4.10 or above
* R including the packages Nozzle.R1, base64, scales, forcats, tidyr, dplyr, ggplot2, fitdistrplus
* R packages shiny, DT and highcharter for R/Shiny visualization tool (optional)
* Perl
* Ensembl Variant Effect Predictor and dependent Perl packages (File::Copy::Recursive, Archive::Zip, DBI)

The workflow system on which the AmpliconSeq pipeline package is built creates
shell scripts for some of the tasks so both Bourne and Bash shells need to be
available.

GATK is required regardless of which variant caller has been selected as it is
used in filtering variants and calculating SNV metrics.

The workflow system and GATK are both Java applications requiring a suitable
Java Runtime Environment. You can test the availability and version Java of on
your system using the following command:

```bash
java -version
```

Recent versions of the HTSJDK library and Picard tools package, including the
version distributed and used in this pipeline require Java 8, i.e. version 1.8
or above when running the above command.

Support for VarDict is limited to the Java version, VarDictJava; the Perl version
has not been tested. The installation directory contains scripts for postprocessing
the output from VarDict taken from the VarDict GitHub repository (see the licensing
section below). These have been tested with VarDictJava version 1.4.10 and may not
work with earlier or later releases. These location of these scripts can be changed
in the configuration file to use updated versions compatible with later releases of
VarDict.

Report(s) are generated using R. The pipeline expects the Rscript executable
to be available on the path and that the dependent packages listed above are
pre-installed.

Ensembl Variant Effect Predictor (VEP) is used to annotate variants with the
consequence on protein coding, details of known variants and associated minor
allele frequencies from the 1000 Genomes project. This is run in offline mode
so additionally the cache files for the reference genome against which the
sequence data were aligned are required.

### Installing the AmpliconSeq pipeline package

The AmpliconSeq pipeline is available in pre-packaged form and can be unpacked
as follows (note that the version number should be substituted in the following
command).

```bash
tar zxf ampliconseq-pipeline-X.Y-distribution.tar.gz
```

### Installing the GATK package

The GATK package cannot be redistributed within this package under the terms
of Broad Institute Software License Agreement. Instead this needs to be
downloaded and unpacked separately.

By default the pipeline contained in this package is configured to use the
GATK jar file called gatk.jar in the installation directory. The GATK jar
file (usually named GenomeAnalysisTK.jar) can be moved or copied to the
installation directory and renamed to gatk.jar. Alternatively create a
symbolic link called gatk.jar to the location of the GATK jar file or
configure the gatkClasspath variable within the pipeline configuration file
to be the full path for the GATK jar file.

-----------

## Running the pipeline

### Configuring the pipeline

1. Copy a sample configuration file from the config subdirectory of the AmpliconSeq installation to a working directory.

2. Modify the configuration settings as appropriate, particularly variables corresponding to input file(s). For details of the configuration options available for a pipeline, refer to the comments on each setting within the corresponding sample configuration file and see the Configuration section below.

#### Running the pipeline

Run the pipeline using the run-pipeline command within the bin subdirectory of the AmpliconSeq installation, e.g.

```
/path_to/ampliconseq/bin/run-pipeline --mode local config.xml
```

where config.xml is the configuration file and the mode option is used to
instruct the workflow engine to run in local mode, i.e. not using a job
scheduler. Note that the --mode command line option is not required but if
set will override the value specified in the configuration file.

To run on a compute cluster using the LSF scheduler:

```bash
/path_to/ampliconseq/bin/run-pipeline --mode lsf config.xml
```

Similarly, to run on a compute cluster using the Slurm scheduler:

```bash
/path_to/ampliconseq/bin/run-pipeline --mode slurm config.xml
```

To view other options that can be specified at run-time:

```bash
/path_to/ampliconseq/bin/run-pipeline --help
```

-----------

## Configuration

The pipeline configuration file is an XML file in which a number of parameters
or variables can be set for a given run. Examples for each variant caller supported
are given in the config subdirectory. These can be copied to the run directory and
modified as appropriate.

There is a comment for each setting within the configuration file that
describes the purpose of that configuration variable and should help in
deciding whether this needs to be set specifically for a run. An example
is given below for setting the variant caller to be used.

```
        <!-- the variant caller (one of HaplotypeCaller, UnifiedGenotyper,
             FreeBayes or VarDict)
        -->
        <variantCaller>VarDict</variantCaller>
```

Configuration options relating to the run-time mode or scheduler can be found
in the executionConfiguration section of the configuration file. More details
are available in the comments in the sample configuration files and in the
workflow system user guide (available separately).

Options that are specific to the AmpliconSeq workflow being run are defined in
a variables block.

An identifier for a pipeline run must be specified in the runId variable within
the configuration file, e.g.

```
    <variables>
        <runId>MDE123</runId>
        ...
    </variables>
```

This is used in creating job names and in constructing names for temporary,
intermediate and final output files.

The ${...} notation is used to substitute variable values within pipeline
definitions and the configuration file. For example it is possible to define
variables using the values of other variables, e.g.

```
    <variables>
        <runId>MDE123</runId>
        <outputFilePrefix>${runId}.</outputFilePrefix>
        ...
    </variables>
```

The ${install} variable is a special variable that refers to the installation
directory. An example of its use is in instructing the workflow engine which
pipeline to run within the configuration file, e.g.

```
    <pipeline>${install}/pipelines/ampliconseq.xml</pipeline>
```

The working directory for a pipeline run is set by default to be the current
working directory in which the run-pipeline command was executed using another
special property, @{user.dir}

```
        <work>@{user.dir}</work>
```

-----------

## Pipeline inputs

The primary inputs to the pipeline are a sample metadata file that provides
sample details for each sequencing library or dataset and a corresponding set
of BAM files, containing aligned sequence reads, one for each dataset or
library.

#### Samples file

The samples file associates information about each of the DNA samples with the
corresponding sequencing dataset/library identifier or barcode. This is a
tab-delimited file with a column containing dataset identifiers with heading ID
(or one of the synonyms Barcode or Dataset ID) and any number of additional
columns.

Other sample details that, if present, will be extracted and included in the
variant output file, include the sample name or identifier (Sample, Sample ID
or Sample name column), sample type (Sample type), and patient identifier
(Patient identifier or Patient ID). Columns can be in any order and matching of
column headings is case-insensitive and allows for some flexibility in terms of
use of spaces or underscore characters, e.g. sample_name, SampleID, etc.

The dataset identifier or barcode usually corresponds to a single sequencing
library. There may be multiple libraries for each DNA sample, as in the
following example where there are duplicate libraries for each sample. Replicate
libraries for a single sample are recognized by sharing a common sample name.
Variant details for replicates are combined into a single row in the collated
variant output files.

An example samples file containing dulicate libraries for each sample is given
below.

```
ID        SAMPLE
FLD0013   JBLAB-2063
FLD0014   JBLAB-2150
FLD0015   JBLAB-2151
FLD0016   JBLAB-2137
FLD0017   JBLAB-199
FLD0037   JBLAB-2063
FLD0038   JBLAB-2150
FLD0039   JBLAB-2151
FLD0040   JBLAB-2137
FLD0041   JBLAB-199
```

The samples file is configured with the sampleSheet variable.

```
        <sampleSheet>${work}/samples.txt</sampleSheet>
```

#### BAM files

The BAM files must be contained within a single directory which is specified
using the bamDir variable in the configuration file, e.g. to use the
subdirectory of the current working directory named bam_files

```
        <bamDir>${work}/bam_files</bamDir>
```

See Configuration section below for more details on using variables and,
specifically, the current working directory within pipelines.

The pipeline expects the BAM files to be named using the dataset identifier (or
barcode) in the samples file as a prefix. For example, in the samples file given
above, BAM files named FLD0013.bam, FLD0014,bam, etc. are expected.

The suffix can be specified in the configuration file by altering the
inputBamFileSuffix variable, e.g. for BAM files such as FLD0013.bwa.bam

```
        <bamSuffix>bwa.bam</bamSuffix>
```

Where BAM files do not follow the required naming scheme or are distributed
across multiple directories, symbolic links can be created in another
directory to avoid moving or renaming files.

Dataset identifiers are used in the naming of intermediate and output files and
for jobs created and run by the workflow manager.

#### Target and amplicon intervals files

Two interval files are required, one specifying the locations of the amplicons,
the other containing the location of the target within each amplicon (i.e.
without the ends corresponding to primer sequences).

These must be in Picard-style interval format and should contain unique
identifiers for each target amplicon. The following is an excerpt from
an example amplicon file.

```
@SQ SN:chr17  LN:81195210
chr17 7572850 7573030 + TP53_E00001757276_1
chr17 7573859 7574054 + TP53_E00001728015_1
```

Each target amplicon should be represented once in the targets file and once in
the amplicons file using the same identifier, so the targets corresponding to
the amplicons in the previous example might be as follows.

```
@SQ SN:chr17  LN:81195210
chr17 7572870 7573010 + TP53_E00001757276_1
chr17 7573881 7574036 + TP53_E00001728015_1
```

The target and amplicon intervals files are configured with the
targetIntervals and ampliconIntervals variables.

```
        <ampliconIntervals>${work}/amplicons.txt</ampliconIntervals>

        <targetIntervals>${work}/targets.txt</targetIntervals>
```

#### Reference data

The reference genome sequence file needs to be specified using the
referenceSequence variable. This should be an indexed FASTA file.

```
        <referenceSequence>/reference_data/GRCh37.fa</referenceSequence>
```

Similarly there a number of variables for setting the location of the
Ensembl VEP cache that will need to be configured.

```
        <!-- the offline cache directory used by Ensembl Variant Effect
             Predictor
        -->
        <variantEffectPredictorCacheDirectory>~/.vep</variantEffectPredictorCacheDirectory>

        <!-- the species annotation data to be used by Ensembl Variant Effect
             Predictor
        -->
        <variantEffectPredictorSpecies>homo_sapiens</variantEffectPredictorSpecies>

        <!-- the genome assembly to be used by Ensembl Variant Effect Predictor
        -->
        <variantEffectPredictorAssembly>GRCh37</variantEffectPredictorAssembly>

        <!-- the reference sequence FASTA file used by Ensembl Variant Effect
             Predictor
        -->
        <variantEffectPredictorReferenceSequence>${referenceSequence}</variantEffectPredictorReferenceSequence>
```

In the above, Ensembl VEP is configured to use the same reference sequence
file as previously specified for use by other tasks in the pipeline. If using
a different version of the reference genome to that used by Ensembl, e.g. one
that uses UCSC-style naming of chromosomes (chrM, chr1, chr2, etc.) then the
variantEffectPredictorReferenceSequence should be set to the version downloaded
from Ensembl as part of the VEP cache to avoid some issues with HGVS annotations
affecting the cDNA and protein effect columns in the output variant spreadsheet.

### Optional inputs

There are a number of optional input files relating to specific calling of
known variants, blacklisting variants and specifying which SNVs to use in
creating an allele fraction table. More details are given in the Notes section
below.

-----------

## Pipeline outputs

The pipeline produces separate SNV and indel VCF files for each dataset along
with alignment summary and target coverage metrics files, including an HTML
report. The SNV and indel calls are collated for all datasets in which details
for replicates for each sample are combined into a single entry or row that also
includes annotations from Ensembl Variant Effect Predictor; these collated
results are produced as both tab-separated text files and an Excel spreadsheet.

Output files can be written to a specific directory by setting the outputDir
variable and can be given a prefix using the outputFilePrefix variable. Note
that a separator character ('.') is not added after the prefix when constructing
output file names so should be added if required; this allows for prefixes of
the form '${runId}_'.

For example to prefix all output files with the run identifier, set the
outputFilePrefix variable as follows.

```
        <outputFilePrefix>${runId}.</outputFilePrefix>
```

See the Configuration section below for more details on the run identifier.

The following are the main outputs of the AmpliconSeq pipeline (note that the
file names may be prefixed as described above).

File                           | Description
-------------------------------|-----------------------------------------------------------------------------
alignment_coverage_report.html | HTML report containing alignment summary metrics and target coverage metrics
amplicon_coverage.xlsx         | Excel spreadsheet containing coverage for each library and amplicon pairing
read_counts.txt                | Table of read counts for each base at every target position for all libraries
snv.txt                        | Collated SNV calls for all sample/patient replicate libraries
snv.filtered.txt               | Filtered SNV calls that exclude low-confidence calls not marked for specific calling
snv_allele_fraction.xlsx       | Excel spreadsheet containing allele fractions for selected SNVs
indel.txt                      | Collated indel calls for all sample/patient replicate libraries
indel.filtered.txt             | Filtered indel calls that exclude low-confidence calls not marked for specific calling
variants.xlsx                  | Excel spreadsheet containing collated, filtered SNV and indel calls


### Variants tables

The variants.xlsx Excel spreadsheet is the principal output from the pipeline
and includes two tables for the SNV and indel calls.

In creating this spreadsheet, variants from the VCF files (see below) are
collated for all libraries along with Ensembl VEP annotations into the snv.txt
and indel.txt tab-delimited files.

These are then filtered to exclude low-confidence and blacklisted variants that
have not been listed for specific calling and blacklisted variants to create the
snv.filtered.txt and indel.filtered.txt tab-delimited files, from which the
spreadsheet is created.

The SNV tabular files, snv.txt and snv.filtered.txt, can be read into the
visualization tool packaged with this pipeline (see the section on low allele
fraction SNVs below for more details).

### VCF files

The pipeline generates indexed SNV and indel VCF files for each library which
can amount to a large number of files. The output directory for these VCF files
defaults to a subdirectory named vcf within the output directory (specified using
the outputDir variable) but can be configured to be a different directory using
the outputVcfDir variable.

### Coverage table

In addition the alignment and coverage HTML report, an Excel spreadsheet is
generated that contains the mean depth covered by each amplicon within each
dataset in tabular form. This uses conditional formatting to highlight
library/amplicon pairs with low coverage; the coverageThreshold setting is
used to set a threshold below which depth of coverage will be highlighted.

### SNV allele fraction table

Allele fractions for a predetermined set of SNVs are calculated and can help
in identifying sample swaps. A tab-delimited file containing the SNVs to be
used can be specified by configuring the snvsForAlleleFractionTable variable
(see note below). If this is not specified the union of all medium and
high-confidence SNVs called in all samples is used.

Where replicate libraries for samples have been included in the run, those
SNVs where there is a substantial difference in the allele fractions between
replicates are highlighted and the number of SNVs flagged are tabulated in
the Flagged column. The maximum tolerated difference in allele fraction can
be set using the maxAlleleFractionDifference configuration variable.

### Read counts file

The read counts file contains the number of reads for every allele at each
target location. This is used in modelling the background noise and filtering
of SNV calls, particularly for low-frequency calling.

The read counts file is also used for depth and allele fraction estimates
where these are not given by the variant caller, e.g. for specific calling
of SNVs that were not called.

The visualization tool uses the read counts file to display scatter plots
of allele fractions for selected SNVs or any given substitution to help
with assessment of whether a variant is real (see the section on low allele
fraction SNVs below for more details).

The way in which overlapping read pairs contribute to allele counts can be
configured using the countOverlappingReadPairsMode setting, including how to
handle positions where the base called in the two reads disagree. Other
configuration settings for the calculation of allele read counts and SNV metrics
are available, including minimum mapping and base qualities; see the sample
configuration file(s) for details.

-----------

## Notes

### Replicate libraries for each sample

The annotated variant output file contains a single row for each variant within
each sample in which the details for each replicate are collated, i.e. each row
contains multiple columns for the quality, depth and allele fraction
for each replicate library. The number of columns for each of these values
depends on the maximum number of replicates for any sample. There is a potential
issue with this if there are a larger number of replicates for control samples,
e.g. water. To avoid an ungainly output file the control libraries can be given
different sample names in the samples file, e.g.

```
ID       SAMPLE
FLD0336  NTC1(water)
FLD0048  NTC2(water)
FLD0072  NTC3(water)
FLD0024  NTC4(water)
FLD0096  NTC5(water)
FLD0312  NTC6(water)
```

### Filters

SNV and indel filtering is carried out using the GATK VariantFiltration tool.
Filters are specified using the snvFilterExpression and indelFilterExpression
variables and take the form of logical expressions. Multiple filters can be applied
in a comma-separated list.

There are also snvFilterName and indelFilterName variables that give names to each
of the filter expressions; these names are included in the output VCF files and are
the names used in the filter columns of the output variant spreadsheet.

An example configuration of filters for the GATK variant callers is given
below. The attributes used in these expressions, e.g. QD, FD, MQ, etc., are
values computed by the variant caller and contained in the VCF file, or can
be additional metrics computed by the pipeline and added to the VCF file.

```
        <!-- the names of filters used in selecting high confidence SNV calls;
             can contain multiple filter names in a comma-separated list where
             each has a corresponding filter expression in the
             snvFilterExpression variable
        -->
        <snvFilterName>QD,FS,MQ,MQRankSum,PositionNoise,DatasetNoise</snvFilterName>

        <!-- the filter expression to use for selecting high confidence SNV
             calls; can contain multiple expressions in a comma-separated list
             where each has a corresponding name in the snvFilterName variable
        -->
        <snvFilterExpression><![CDATA[QD < 2.0,FS > 60.0,MQ < 40.0,MQRankSum < -12.5,VariantAlleleFraction < PositionSubstitutionAFThreshold,VariantAlleleFraction < DatasetSubstitutionAFThreshold]]></snvFilterExpression>

        <!-- the names of filters used in selecting high confidence indel calls;
             can contain multiple filter names in a comma-separated list where
             each has a corresponding filter expression in the
             indelFilterExpression variable
        -->
        <indelFilterName>QD,FS</indelFilterName>

        <!-- the filter expression to use for selecting high confidence indel
             calls; can contain multiple expressions in a comma-separated list
             where each has a corresponding name in the indelFilterName variable
        -->
        <indelFilterExpression><![CDATA[QD < 2.0,FS > 200.0]]></indelFilterExpression>
```

There should be an equal number of names and expressions in the snvFilterName
and snvFilterExpression variables, similarly an equal number in the
indelFilterName and indelFilterExpression variables.

Note that the <![CDATA[...]]> construct is used to allow characters such as '<'
and '>' within filter expressions that would otherwise be interpreted as opening
or closing elements as part of the XML syntax.

Filters are specified using Java Expression Language (JEXL) syntax in the
snvFilterExpression and indelFilterExpression variables. See the GATK
documentation for more details about using JEXL expressions for filtering
variants.

### Confidence level for variants called in replicate libraries

The variants in the collated variant output files are given a confidence level
that can be high, medium or low. Factors in determining the confidence include
the number of replicates in which the call is made without being filtered and
whether a minimum depth threshold was reached for the replicate library in order
to make the call. The depth or coverage threshold is specified in the
configuration file using the coverageThreshold variable; by default this is set
to 100.

Confidence | Criteria
-----------|----------------------------------------
High       | Call passes filters in all replicates and depth in each is not below the minimum coverage threshold.
Medium     | Call passes filters in at least one replicate with a depth that is not below the minimum coverage threshold.
Low        | Calls which don't pass filters in any replicate or for which there is insufficient coverage.

Low-confidence variants are excluded from the variant spreadsheet, variants.xlsx
file, unless listed for specific calling. The tab-delimited files snv.txt and
indel.txt include all variants called including those with low confidence.

### Depths and allele fractions

Depths and allele fractions for variants called in one or more but not all
replicates are included in the output variant tables but where a call was not
made these are computed from the read counts that are generated using GATK.
For replicates for which a call is made, depths and allele fractions will
preferentially be taken from the VCF file output by the variant caller. If the
variant caller doesn't call the variant for a replicate or doesn't output the
depth or allele fraction for a call, the values included in the final variants
table will instead be those computed from the read counts.

### Specific variant calling

There is support for "specific calling" for SNVs, i.e. the depth and allele
fraction for specified SNVs will be included in the output variant spreadsheet
regardless of whether those variants were called or filtered as low-confidence
calls. This is potentially useful for longitudinal tracking of patient samples
in which the allele fractions of known mutations can be monitored over time.

A tabular file containing variants for specific calling can be configured by
setting the specificCallingVariants variable. If specific calling is not
required this variable can be left unset.

```
        <specificCallingVariants>${work}/known_variants.txt</specificCallingVariants>
```

The file should be tab-delimited and contain columns for the sample, chromosome,
position, reference base and alternate allele as shown in the following example.

```
Sample   Chromosome  Position  Ref  Alt
mde1012  chr17       7577538   C    T
mde1025  chr17       7577594   A    T
mde1079  chr17       7577538   C    T
mde1200  chr17       7578406   C    T
```

The output variant spreadsheet, variants.xlsx, includes a column called
Specific that allows for filtering within Excel of the variants specified
for specific calling.

In the current release there is only limited support for specific calling of
indels. If indels are included in the specificCallingVariants file, these
will be included in the indels worksheet in the variants spreadsheet with a
"yes" in the Specific column and Ensembl VEP annotations but depth and allele
fractions will only be available if the indel was called.

### Blacklisting variants

For some target amplicon panels it may be desirable to blacklist certain
variants, e.g. where those variants are very common in the population at
large. This can be done by creating a 4-column tab-delimited file with
columns for the chromosome, position, reference base and alternate allele.

The blacklist file is configured by setting the blacklistedVariants variable.
If blacklisting is not needed this variable can be left unset.

```
        <blacklistedVariants>${work}/blacklist.txt</blacklistedVariants>
```

Both SNVs and indels can be included in the blacklist file, an example of
which is given below.

```
Chromosome  Position  Ref  Alt
chr17       7577644   C    G
chr17       7578210   T    C
chr17       7579472   G    C
chr17       7579920   A    AC
```

Variants in the blacklist are excluded from the final variants spreadsheet
unless these have also been specified for specific calling.

### SNV allele fraction table for identifying sample swaps

The SNV allele fraction spreadsheet can help in identifying possible sample
swaps. This contains allele fractions for every sample for a predetermined
set of SNVs, for example common SNPs from the 1000 Genomes project lying
within the target amplicons.

The SNV file is configured by setting the snvsForAlleleFractionTable
variable and is a tab-delimited file containing 4 columns for the
chromosome, position, reference base and alternate allele.

```
        <snvsForAlleleFractionTable>${work}/1000_genomes_snps.txt</snvsForAlleleFractionTable>
```

The following snippet shows the format of the SNV file.

```
Chromosome  Position  Ref  Alt
chr7        55248926  T    C
chr7        55248930  C    T
chr7        55248965  C    T
chr7        55248991  C    G
chr7        55249063  G    A
```

Allele fractions are computed using GATK and may differ from the allele
fractions computed by the variant caller; note that the variant caller will
not necessarily call the selected SNVs for which allele fractions are to be
calculated.

Potential sample swaps may be indicated if there are several sizeable differences
in the allele fractions between pairs of sample replicates. The variant allele
fraction is highlighted if it differs from that for a replicate by more than a
tolerated amount specified using the maxAlleleFractionDifference setting and
where the depth is above a minimum value specified with the coverageThreshold
variable.

### Differing depths reported by various tools within the pipeline

Three different measures of the depth of sequencing are reported in the various
output files for this pipeline.

The Picard CollectTargetedPcrMetrics tool gives the mean target coverage for
each amplicon target region defined in the targets file (see above). The depth
computed by Picard is given in the target coverage output file and used in the
HTML report.

The pipeline segregates reads by the target amplicon and, depending on the
settings for the maximum distance of the alignment start or end position to the
amplicon start or end position (maxDistanceFromAmpliconEnd setting), i.e. how
close the alignment has to be to the expected location of reads, and whether both
ends of a read pair need to be anchored within this distance
(requireBothEndsAnchored setting), the depth computed will be lower than that
given by Picard. In addition, where there are overlapping amplicons Picard will
count the reads in the overlapping region twice when calculating the average
depth for each of the overlapping amplicon. The depth given in the amplicon
coverage file and spreadsheet is that computed for reads assigned to each
amplicon separately.

Finally, the depths reported in the collated variant output file(s) are taken
from the VCF files produced by the variant caller. Depending on the variant
caller, these may be the number of reads used to call variants and may exclude
reads filtered because of low mapping quality or base quality at the variant
position.

Generally, the depth from the Picard metrics will be higher than the depth
computed by the pipeline for reads assigned to separate amplicons, and this
will likely be higher than the number of reads used to call variants.

-----------

## Low allele fraction SNVs

The AmpliconSeq pipeline supports calling of SNVs with low allele fractions through
modelling the background noise in the sequence data. This allows for calling SNVs from
circulating tumour DNA in plasma samples with allele fractions as low as 0.1%, although
this is only possible for certain substitution types in which the background noise levels
are low. For example, with Illumina sequencers C>G, G>C, A>C and T>G variants are more
amenable to calling at low allele fractions than A>G, T>C, C>T and G>A.

FreeBayes and VarDict can call variants with very low allele fractions by setting
the minimumAlleleFraction configuration variable, e.g. to 0.1% as shown below.
The GATK variant callers do not call variants with allele fractions below around 10 - 15%.

```
        <!-- the minimum allele fraction for calling variants
             (applies to FreeBayes and VarDict)
        -->
        <minimumAlleleFraction>0.001</minimumAlleleFraction>
```

### Background noise filters

A Beta probability distribution is fitted to the distribution of allele fractions for
all datasets in the run at each target position and for each of the three possible
substitutions at that position. This is used to derive an allele fraction threshold
below which there would be low confidence in a SNV call, using the quantile corresponding
to probability, p = 0.9999. Fitting the distribution for each amplicon target position
takes account of both the background noise associated with the substitution type and
the position within the amplicon, in effect identifying and accounting for noisy positions.

The allele fraction threshold is added to the VCF file and then can be used during
filtering by adding a filter expression as shown in the following snippet from the
configuration file. Note that there would usually be other filters; see the Filters
sections above for more details.

```
        <!-- the names of filters used in selecting high confidence SNV calls;
             can contain multiple filter names in a comma-separated list where
             each has a corresponding filter expression in the
             snvFilterExpression variable
        -->
        <snvFilterName>PositionNoise</snvFilterName>

        <!-- the filter expression to use for selecting high confidence SNV
             calls; can contain multiple expressions in a comma-separated list
             where each has a corresponding name in the snvFilterName variable
        -->
        <snvFilterExpression><![CDATA[VariantAlleleFraction < PositionSubstitutionAFThreshold,]]></snvFilterExpression>
```

Similarly, probability distributions are fitted separately for each substitution
type within each sample dataset using all positions for which that substitution is
possible. This models the background noise at the library/dataset level, so higher
allele fraction thresholds are obtained for noisy libraries.

The sample configuration files provided in the pipeline distribution in the config
subdirectory contain filter expressions that are suitable for each of supported
variant callers and include position and dataset filters.

### R/Shiny visualization tool

The AmpliconSeq package contains a visualization tool written in R and using
the Shiny framework. To start the application, navigate on the command line to
the R/shiny subdirectory within the installation and run the following.

```
Rscript start_shiny_server.R
```

This should open the application within a web browser. Brief instructions for
using the applications are given below.

*Warning:* this application should be considered to a beta release and
occasionally freezes, in which case either the web page should be refreshed
and the data reloaded or, failing that, the application may need to be
restarted.

1. From the *Read counts* tab, load the read counts file, read_counts.txt
(with specified prefix), produced as one of the output files by the AmpliconSeq
pipeline. The read counts file is large and it can take a few seconds before
the table is populated.

2. From the *Locations* tab, select a target position and alternate allele
to show the allele fractions for all libraries/datasets within the run as a
scatter plot. If there are multiple overlapping target amplicons at that position
these are selected using a drop down menu and viewed separately. Superimposed on
the scatter plot is a box and whiskers plot. The allele fraction threshold below
which SNV calls should be filtered is shown as a red line.

3. Within the *Locations* tab, select the *Density plot* tab to view the allele
fractions in the form of a kernel density plot. Some filtering of the data points
with the highest allele fractions and those with values of zero is performed to
aid with the fitting of a probability distribution. This filtered density is also
shown as well as the fitted probability distribution. The parameters for filtering
data points to exclude from fitting can be modified using the dialog on the left
hand side. Also, there is a choice of probability distribution that can be fitted
with the Normal, Log Normal and Beta distributions available.

4. The *Cullen and Frey* graph shows the degree of kurtosis and skewness for the
allele fractions for the selected position and substitution in the context of a
number of theoretical distributions. In most cases, the beta distribution seems
a good choice and is the distribution used in fitting by the AmpliconSeq pipeline.

5. The scatter and density plots are interactive, supporting zooming, selection
of data points and display of tooltips that provide a summary for each data point.
Selecting a data point provides details about the substitution in a table below
the scatter plot. When SNV calls are also loaded (see below) the points in the
scatter plot are categorized according to call status; clicking on the category
in the legend toggles the display of points within that category.

6. From the *SNVs* tab, load the SNVs called by the AmpliconSeq pipeline. Either
the full set of SNV calls (snv.txt) or the filtered calls (snv.filtered.txt) that
are included in the final variants spreadsheet can be loaded. The table is populated
with details of the SNV calls including the filters, if any, that were applied
and the confidence of the call. The table is interactive and supports sorting by
clicking on a column and filtering using the search box above the table or the
individual column search boxes below the table. Selecting a SNV call in this
table will update the selections in the Locations and Libraries tabs.

7. The *Libraries* tab is very similar to the *Locations* tab but shows the
allele fractions for all substitutions within a selected library/dataset.

-----------

## Useful links

[Genome Analysis Toolkit (GATK)](https://www.broadinstitute.org/gatk)

[FreeBayes](https://github.com/ekg/freebayes)

[VarDict](https://github.com/AstraZeneca-NGS/VarDictJava)

-----------

## License

The AmpliconSeq pipeline is made available under the MIT License (see separate
LICENSE file).

### VarDict copyright and license

The scripts teststrandbias.R and var2vcf_valid.pl from the VarDict package are
redistributed in the scripts subdirectory under the MIT license (see
LICENSE.VarDict in the licenses subdirectory).

### Highcharts license

The visualization tool uses the Highcharts JavaScript library via the highcharter
R package. This is free for non-profit organizations. For commercial use, licensing
terms and conditions can be found at [www.highcharts.com](https://www.highcharts.com).

-----------

[workflow]: ampliconseq_portrait.svg "Schematic overview of the AmpliconSeq workflow"
[scatter_plot]: scatter_plot.png "Allele fraction scatter plot"
[density_plot]: density_plot.png "Allele fraction density plot"
[cullen_frey_plot]: cullen_frey_plot.png "Cullen and Frey graph"
[snv_table]: snv_table.png "Table of SNV calls"
