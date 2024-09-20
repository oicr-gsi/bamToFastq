# bamToFastq

Generating R1 and R2 from bam files

## Overview
![](./bamToFastq.png?raw=true "Workflow diagram")

## Dependencies

* [java 8](https://www.java.com/en/download/manual.jsp)
* [samtools 1.14](https://github.com/samtools/samtools/archive/0.1.19.tar.gz)
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
`fileNaming`|String|The naming scheme for the extracted FASTQs


#### Optional workflow parameters:
Parameter|Value|Default|Description
---|---|---|---
`finalOutputDirectory`|String|""|Optional final output directory. Copy out fastq files using a task


#### Optional task parameters:
Parameter|Value|Default|Description
---|---|---|---
`countFlags.prefix`|String|"output"|String prepended to flagstat file
`countFlags.memory`|Int|24|Memory allocated for this job
`countFlags.timeout`|Int|12|Time in hours before task timeout
`countFlags.modules`|String|"samtools/1.14"|Required environment modules
`nameCheck.memory`|Int|24|Memory allocated for this job
`nameCheck.timeout`|Int|12|Time in hours before task timeout
`backExtract.modules`|String|"picard/2.21.2"|Required environment modules
`backExtract.memory`|Int|24|Memory allocated for this job
`backExtract.timeout`|Int|12|Time in hours before task timeout
`haveCustomDirRename.modules`|String|""|Required environment modules
`haveCustomDirRename.memory`|Int|24|Memory allocated for this job
`haveCustomDirRename.timeout`|Int|12|Time in hours before task timeout
`copyOutFastq.memory`|Int|24|Memory allocated for this job
`copyOutFastq.timeout`|Int|12|Time in hours before task timeout
`composeLog.memory`|Int|4|Memory allocated for this job
`composeLog.timeout`|Int|2|Time in hours before task timeout
`noCustomDirRename.modules`|String|""|Required environment modules
`noCustomDirRename.memory`|Int|24|Memory allocated for this job
`noCustomDirRename.timeout`|Int|12|Time in hours before task timeout


### Outputs

Output | Type | Description | Labels
---|---|---|---
`flagStat`|File|A TXT file containing flag information about the BAM file|vidarr_label: flagStat
`modFastqs`|Array[File]?|An optional array of fastq.gz files|vidarr_label: modFastqs
`copyLog`|File?|log file from copy out task, provisioned if a custom output dir requested|vidarr_label: copyLog


## Commands
 
This section lists command(s) run by bamToFastq workflow
 
* Running bamToFastq
 
bamToFastq generates R1 and R2 fastqs for BAM files with one or more readgroups.
 
### Pulling out readgroups and flag data from bam file
 
```
     module load samtools
     samtools flagstat ~{bamFile} > ~{prefix}.flagstat.txt
     samtools view -H ~{bamFile} | grep -G "^@RG" > readgroups.tsv
```
 
### Checking name validity and creating new filenames.
 
By default, the fastqs are named using the read IDs. To use RG tags in the file names of the generated fastqs, input paramenter filNaming should be formatted like this:
 
```
 "{tag1}_{tag2}"
 "{PU}_{PL}"
 "{LB}_{PU}_test01"
```
(fileNaming follows Python 3.9 string-mapping formatting)
 
The file naming scheme will be declared invalid if:
 - It results in non-unique file names
 - It contains an invalid RG tag
 
```
         python3<<CODE
 
         import difflib
         import re
         with open(r"~{readGroups}",'r') as rg_data:
             rg = rg_data.read()
        
         fileName = "~{fileName}"
         
         def validate(valid="true"):
             isValid = open("isValid.txt", 'w')
             isValid.write(valid)
             isValid.close()
         
         validate()
 
             # Check No. 1:
             # Comparing RG tags from readgroups file and fileName
             # Exits if a tag is missing
 
         inputs = re.findall(r"[^{\}]+(?=})", fileName)
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
                 newName = fileName.format_map(dictArray[j])
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
 
```
     java -Xmx20g -jar $PICARD_ROOT/picard.jar SamToFastq \
     INPUT=~{bamFile} \
     RG_TAG="ID" \
     OUTPUT_DIR=. \
     OUTPUT_PER_RG=true \
     NON_PF=true \
     RE_REVERSE=true \
     VALIDATION_STRINGENCY=LENIENT
 
```
### Rename fastqs based on naming scheme.
 
```
     python3<<CODE
 
     import os
     import ast
 
     rgData = ast.literal_eval("~{rgData}")
 
     f = "~{sep=' ' rawFastqs}"
     fastqs = f.split()
 
     for fastq in fastqs:
 	path = os.getcwd()
 	fastqID = os.path.basename(fastq)[:-8] # determines ID
 	readNum = int(fastq[-7]) # determines read number
 	newName = rgData[fastqID][(readNum-1)]
 	formattedFileName = f"{path}/{newName}"
 	os.rename(fastq, formattedFileName)
     
     CODE
 
     ls *.fastq.gz > outfilenames
 
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
