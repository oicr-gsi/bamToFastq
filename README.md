# bamToFastq

Given aligned reads in bam format, this workflow will backextract to generate fastq files based on the readgroups.  By default, the files are named based on the readgroup IDs, but a custom naming scheme based on other readgroup fields can be provided.  This takes the form of mixing text with {FD} symbols, where FD is a readgroup FD (SM, PU, etc)

## Overview

## Dependencies

* [java 8](https://www.java.com/en/download/manual.jsp)
* [samtools 1.9](https://github.com/samtools/samtools/archive/0.1.19.tar.gz)
* [picard 2.21.2](https://broadinstitute.github.io/picard/command-line-overview.html)
* [python 3.6](https://www.python.org/downloads/)


## Usage

### Cromwell
```
java -jar cromwell.jar run bamToFastq.wdl --inputs inputs.json
```

### Inputs

#### Required workflow parameters:
Parameter|Value|Description
---|---|---
`bamFile`|File|A BAM file with one or more readgroups


#### Optional workflow parameters:
Parameter|Value|Default|Description
---|---|---|---
`fileNaming`|String|"{ID}"|The naming scheme for the extracted FASTQs


#### Optional task parameters:
Parameter|Value|Default|Description
---|---|---|---
`examineBam.prefix`|String|"output"|String prepended to flagstat file
`examineBam.memory`|Int|24|Memory allocated for this job
`examineBam.timeout`|Int|12|Time in hours before task timeout
`examineBam.modules`|String|"samtools/0.1.19"|Required environment modules
`nameCheck.modules`|String|""|Required environment modules
`nameCheck.memory`|Int|24|Memory allocated for this job
`nameCheck.timeout`|Int|12|Time in hours before task timeout
`backExtract.modules`|String|"picard/2.21.2"|Required environment modules
`backExtract.memory`|Int|24|Memory allocated for this job
`backExtract.timeout`|Int|12|Time in hours before task timeout
`renameFastqs.modules`|String|""|Required environment modules
`renameFastqs.memory`|Int|24|Memory allocated for this job
`renameFastqs.timeout`|Int|12|Time in hours before task timeout


### Outputs

Output | Type | Description
---|---|---
`flagStat`|File|A TXT file containing flag information about the BAM file
`fastq`|Array[File]?|one or more fastq files as determined by the readgroup information in the bam file


## Commands
 This section lists command(s) run by bamToFastq workflow
 
 * Running bamToFastq
 
 bamToFastq generates R1 and R2 fastqs for BAM files with one or more readgroups.
 
 ### Extract readgroup and flagstatinfo from bam file
 Readgroup information is needed to rename files.  Flagstat information is provided for comparison to the fastq ouput.
 
 ```
 samtools flagstat ~{bamFile} > ~{prefix}.flagstat.txt
 samtools view -H ~{bamFile} | grep -G "^@RG" > readgroups.tsv
 
 ```
 
 ### Verify that fields in the requested naming scheme are found in the bam file.
 By default, the fastqs are named using the read IDs. To modify names, readgroup fields can be indicated to use their values, mixing with text strings. The input parameter fileNaming should be formatted like this:
 
 ```
 "{tag1}_{tag2}_text"
 eg.
 "{PU}_{PL}". -> {PU}_{PL}_Rn.fastq.gz
 "{LB}_{PU}_test01" -> {LB}_{PU}_test01_Rn.fastq.gz
 
 where Rn is the read number (R1, R2)
 ```
 
 The file naming scheme will be declared invalid if:
 - It results in non-unique file names
 - It contains an invalid RG tag
 
 One should inspect the bam header readgroups prior to choosing a naming scheme.
 
 ```
 inline python code available in the WDL file
 
 ```
 ### Backextract fastq files from bam file per readgroup.
 ```
 java -Xmx20g -jar $PICARD_ROOT/picard.jar SamToFastq \
      INPUT=~{bamFile} \
      RG_TAG="ID" \
      OUTPUT_DIR=. \
      OUTPUT_PER_RG=true \
      COMPRESS_OUTPUTS_PER_RG=true \
      NON_PF=true \
      RE_REVERSE=true \
      VALIDATION_STRINGENCY=LENIENT
 
 ```
 
 ### Rename fastqs based on naming scheme.
 the default naming scheme is to use the RG IDs.  If an alternate naming pattern was provided then files will be renamed
 
 ```
 inline python code available in the WDL file
 
 ```
 ## Support

For support, please file an issue on the [Github project](https://github.com/oicr-gsi) or send an email to gsi@oicr.on.ca .

_Generated with generate-markdown-readme (https://github.com/oicr-gsi/gsi-wdl-tools/)_
