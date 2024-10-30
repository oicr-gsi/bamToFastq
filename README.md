# bamToFastq

Given aligned reads in bam format, this workflow will backextract to generate fastq files based on the readgroups.  By default, the files are named based on the readgroup IDs, but a custom naming scheme based on other readgroup fields can be provided.  This takes the form of mixing text with {FD} symbols, where FD is a readgroup FD (SM, PU, etc)

## Overview

## Dependencies

* [java 8](https://www.java.com/en/download/manual.jsp)
* [samtools 1.16.1](https://github.com/samtools/samtools/releases/tag/1.16)
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
`outputFilePrefix`|String|prefix to use to identify bam summary files


#### Optional workflow parameters:
Parameter|Value|Default|Description
---|---|---|---
`finalOutputDirectory`|String|""|Optional final output directory. Copy out fastq files using a task
`fileNaming`|String|"{ID}"|The naming scheme for the extracted FASTQs


#### Optional task parameters:
Parameter|Value|Default|Description
---|---|---|---
`examineBam.memory`|Int|24|Memory allocated for this job
`examineBam.timeout`|Int|12|Time in hours before task timeout
`examineBam.modules`|String|"samtools/1.16.1"|Required environment modules
`nameCheck.modules`|String|"samtools/1.16.1"|Required environment modules
`nameCheck.memory`|Int|24|Memory allocated for this job
`nameCheck.timeout`|Int|12|Time in hours before task timeout
`unsortBam.memory`|Int|24|Memory allocated for this job
`unsortBam.timeout`|Int|12|Time in hours before task timeout
`unsortBam.modules`|String|"samtools/1.16.1"|Required environment modules
`backExtractByRG.modules`|String|"picard/2.21.2"|Required environment modules
`backExtractByRG.memory`|Int|24|Memory allocated for this job
`backExtractByRG.timeout`|Int|96|Time in hours before task timeout
`haveCustomDirRename.modules`|String|""|Required environment modules
`haveCustomDirRename.memory`|Int|24|Memory allocated for this job
`haveCustomDirRename.timeout`|Int|12|Time in hours before task timeout
`haveCustomDirReview.modules`|String|"fastqc/0.11.9"|Required environment modules
`haveCustomDirReview.memory`|Int|24|Memory allocated for this job
`haveCustomDirReview.timeout`|Int|12|Time in hours before task timeout
`copyOutFastq.memory`|Int|24|Memory allocated for this job
`copyOutFastq.timeout`|Int|12|Time in hours before task timeout
`composeLog.memory`|Int|4|Memory allocated for this job
`composeLog.timeout`|Int|2|Time in hours before task timeout
`noCustomDirRename.modules`|String|""|Required environment modules
`noCustomDirRename.memory`|Int|24|Memory allocated for this job
`noCustomDirRename.timeout`|Int|12|Time in hours before task timeout
`noCustomDirReview.modules`|String|"fastqc/0.11.9"|Required environment modules
`noCustomDirReview.memory`|Int|24|Memory allocated for this job
`noCustomDirReview.timeout`|Int|12|Time in hours before task timeout
`summarize.memory`|Int|24|Memory allocated for this job
`summarize.timeout`|Int|12|Time in hours before task timeout


### Outputs

Output | Type | Description | Labels
---|---|---|---
`bamFileFlagstat`|File|A TXT file containing flag information about the BAM file|vidarr_label: bamFileFlagstat
`fastqc`|Array[File]?|An optional array of FastQC report files|vidarr_label: fastqc
`fastqs`|Array[File]?|An optional array of fastq.gz files|vidarr_label: fastqs
`summary`|File|Summary file generated by processing fastQC and samstat report data|vidarr_label: summary
`copyLog`|File?|log file from copy out task, provisioned if a custom output dir requested|vidarr_label: copyLog


## Commands

This section lists command(s) run by bamToFastq workflow
 
 * Running bamToFastq
 
Given aligned reads in bam format, this workflow will backextract to generate fastq files based on the readgroups. By default, the files are named based on the readgroup IDs, but a custom naming scheme based on other readgroup fields can be provided. This takes the form of mixing text with {FD} symbols, where FD is a readgroup FD (SM, PU, etc)
 
 
### Un-sort bam (sort by read name)
```
     samtools sort --threads 8 -n ~{bamFileToSort} -O bam -o unsorted.bam
```
 
### Examine Bam: collect stats with samtools
 
This runs samtools, collecting read group data and stats
 
```
     set -euo pipefail
     samtools view -H ~{bamFile} | grep -G "^@RG" > readgroups.tsv
     samtools stats --threads 8 ~{bamFile} > ~{prefix}.samstats.txt
```
 
### nameCheck: ensure we do not have collision among output files
 
Use read groups from the input bam file and naming schema so that the names for the
outputs are unique.
 
```
         
         samtools view -H ~{bamFile} | grep -G "^@RG" > readgroups.tsv
         python3<<CODE
 
         import difflib
         import re
         with open(r"readgroups.tsv",'r') as rg_data:
             rg = rg_data.read()
        
         fileNaming = "~{fileNaming}"
         
         def validate(valid="true"):
             isValid = open("isValid.txt", 'w')
             isValid.write(valid)
             isValid.close()
         
         validate()
 
             # Check No. 1:
             # Comparing RG tags from readgroups file and fileNaming
             # Exits if a tag is missing
 
         inputs = re.findall(r"[^{\}]+(?=})", fileNaming)
         tags = list(set(re.findall(r"[^\t:]+(?=:)", rg)))
 
         if all(i in tags for i in inputs) == False:  
             print("FAILED: Missing input tag")
             d = difflib.Differ()
             diff = d.compare(tags, inputs)
             errorMessage = '\n'.join(diff)
             validate("false")
             raise ValueError("Input tag does not exist:" + '\n' + errorMessage)
 
             # Check No. 2:
             # Predicting modified file name:
             # Exits if non-unique names are created
 
         rgArray = []
         dictArray = []
 
         for line in rg.splitlines():
             data = line[1:].split()
             rgArray.append(data)
 
         for row in rgArray:
             del row[0]
             for i in range(len(row)):
                 k = row[i].split(":")
                 row[i] = k
             dictArray.append({row[i][0]: row[i][1] for i in range(len(row))})
 
         fastqNames = []
         rgData = {}
 
         for j in range(len(dictArray)):
             idData = []
             for readNum in [1, 2]:
                 newName = fileNaming.format_map(dictArray[j])
                 predictedFileName = f"{newName}_{readNum}.fastq.gz"
                 if predictedFileName in fastqNames:
                     validate("false")
                     raise ValueError("File name results in non-unique names")
                 else:
                     fastqNames.append(predictedFileName)
                     idData.append(predictedFileName)
             rgData[dictArray[j]["ID"]]=idData
 
         print(rgData)
         
         CODE
 
```
 
### Backextract fastq files from bam file per readgroup.
 
We use picard tools for this (SamToFastq)
 
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
 
             ls *.fastq.gz > outfilenames            
 
```
 
### Rename Fastqs according to selected naming schema
 
the default naming scheme is to use the RG IDs. If an alternate naming pattern was provided then files will be renamed
 
```
             python3<<CODE
 
             import os
             import ast
             import re
 
             rgData = ast.literal_eval("~{rgData}")
             f = "~{sep=' ' rawFastqs}"
             fastqs = f.split()
 
             for fastq in fastqs:
                path = os.getcwd()
                ## get the filename, without the extension
                fastqID = re.sub(".fastq.gz","",os.path.basename(fastq))
                # determines read number
                readNum = int(fastqID[-1])
                ### get rid of the last two characters _N
                fastqID = fastqID[:-2]
                ### get the new name from the rgData dictionary
                newName = rgData[fastqID][(readNum-1)]
                ### format the file name
                formattedFileName = f"{path}/{newName}"
                ### rename the files
                os.rename(fastq, formattedFileName)
             
             CODE
 
             ls *.fastq.gz > outfilenames
 
```
 
### Generate fastqc reports for each fastq file
  
```
             ## run fastqc, this will generate an html report and zipped files, all in the working directory
             fastqc ~{fastq} --outdir .
             ### locate the data file in the zipped output
             datafile=`unzip -l ~{id}_fastqc.zip | grep fastqc_data | sed 's/.* //'`
             ### pull the text data out to a file
             unzip -p ~{id}_fastqc.zip $datafile > ~{id}_fastqc_data.txt
 
             ls *fastqc* > outfilenames
             ls *fastqc_data* > datafilenames
```
 
### Summarize samstats
 
```
           ### samstats
           bamtotal=`cat ~{samstats} | grep ^SN | cut -f 2- | grep "raw total sequences:" | cut -f2`
           bampaired=`cat ~{samstats} | grep ^SN | cut -f 2- | grep "reads paired:" | cut -f2`
           lnavg=`cat ~{samstats}| grep ^SN | cut -f 2- | grep "average length:" | cut -f2`
           lnmax=`cat ~{samstats}| grep ^SN | cut -f 2- | grep "maximum length:" | cut -f2`
           insertsize=`cat ~{samstats} | grep ^SN | cut -f 2- | grep "insert size average:" | cut -f2`
 
           echo "INPUT bam file stats" > "~{id}.summary.txt"
           echo -e "bam total\t$bamtotal" >> "~{id}.summary.txt"
           echo -e "bam paired\t$bampaired" >> "~{id}.summary.txt"
           echo -e "bam mean readlength\t$lnavg" >> "~{id}.summary.txt"
           echo -e "bam max readlength\t$lnmax" >> "~{id}.summary.txt"
           echo -e "bam insertsize\t$insertsize" >> "~{id}.summary.txt"
 
           echo ""  >> "~{id}.summary.txt"
           echo "fastq ReadCounts"  >> "~{id}.summary.txt"
           totalreads=0
           files="~{sep=' ' fastqc}"
           for f in $files
           do
             id=`basename $f "_fastqc_data.txt"`
             count=`cat $f | grep "Total Sequences" | cut -f2`
             echo -e "$id\t$count"  >> "~{id}.summary.txt"
             totalreads=$(($totalreads + $count)) 
           done
           echo -e "fastq total\t$totalreads"  >> "~{id}.summary.txt"
 
           echo ""  >> "~{id}.summary.txt"
           echo "fastq Sequence Length Distributions"  >> "~{id}.summary.txt"
           for f in $files
           do
             id=`basename $f "_fastqc_data.txt"`	
             echo $id  >> "~{id}.summary.txt"
             cat $f | grep "Sequence Length Distribution" -A 999999999 | grep "Sequence Duplication Levels" -B 999999999 | grep -v ">>" >> "~{id}.summary.txt"
           done
 
```
 
### Copy out files if output directory specified
 
```
           set -euxo pipefail
           if [[ -e ~{outputDir} && -d ~{outputDir} ]]; then
             cp ~{Fastq} ~{outputDir}
             echo "File ~{basename(Fastq)} copied to ~{outputDir}"
           else
             echo "Final Output Directory was not configured, ~{basename(Fastq)} was not provisioned properly"
           fi
```
 
### Compose log file from messages produced by copy out shard
 
```
           python3 <<CODE
           import json
           import re
           m_lines = re.split(",", "~{sep=',' messages}")
           with open("copy_out.log", "w") as log:
               for line in m_lines:
                   log.write(line + "\n")
           log.close()
           CODE
```
## Support

For support, please file an issue on the [Github project](https://github.com/oicr-gsi) or send an email to gsi@oicr.on.ca .

_Generated with generate-markdown-readme (https://github.com/oicr-gsi/gsi-wdl-tools/)_
