version 1.0

workflow bamToFastq {
    input {
        File bamFile
        String fileNaming = "{ID}"
        String outputFilePrefix
    }

    parameter_meta {
        bamFile: "A BAM file with one or more readgroups"
        fileNaming: "The naming scheme for the extracted FASTQs"

    }

    meta {
        author: "Murto Hilali, Lawrence Heisler"
        email: "mhilali@oicr.on.ca, lheisler@oicr.on.ca"
        description: "Given aligned reads in bam format, this workflow will backextract to generate fastq files based on the readgroups.  By default, the files are named based on the readgroup IDs, but a custom naming scheme based on other readgroup fields can be provided.  This takes the form of mixing text with {FD} symbols, where FD is a readgroup FD (SM, PU, etc)"
        dependencies: [
            {
                name: "java/8",
                url: "https://www.java.com/en/download/manual.jsp"
            },
            {
                name: "samtools/1.16.1",
                url: "https://github.com/samtools/samtools/releases/tag/1.16"
            },
            {
                name: "picard/2.21.2",
                url: "https://broadinstitute.github.io/picard/command-line-overview.html"
            },
            {
                name: "python/3.6",
                url: "https://www.python.org/downloads/"
            }
        ]
        output_meta: {
            bamFileFlagstat: "A TXT file containing flag information about the bam file",
            fastq: "one or more fastq files as determined by the readgroup information in the bam file"
        }
    }
    
    # call function to assess content with flagstat and pull out the readgroups
    call examineBam {
        input:
            bamFile = bamFile,
            prefix = outputFilePrefix
    }

    ### check on the rg names before doing anything
    ### by default this will check that the ID tag is present, which it should be
    ### if fileRenaming has other values, it will look to ensure those tags are present
    call nameCheck {
           input:
              readGroups = examineBam.readGroups,
              fileName = fileNaming
    }

    
    ### proceed with backextraction only if the default ID is indicating for filename
    ###. or if the namecheck passed
    if (nameCheck.valid == true){
        call unsortBam{
            input:
               bamFileToSort = bamFile
        }
        call backExtractByRG { 
            input:
                bamFile = unsortBam.bamFile
        }
        call backExtract { 
            input:
                bamFile = unsortBam.bamFile
        }
    }

    if (fileNaming != "{ID}"){
      call renameFastqs { 
        input:
          rawFastqs = backExtractByRG.rawFastqs,
          rgData = nameCheck.rgData
      }
    }

    ### fastq reports
    scatter( fq in select_first([renameFastqs.modFastqs,backExtractByRG.rawFastqs]) ) {
       String fqid = basename(fq, ".fastq.gz")
       call reviewFastq {
          input:
            fastq = fq,
            id = fqid	   
       }
    }

    output {
      File flagStats = examineBam.samflagstats
      File samStats = examineBam.samstats 
      Array[File] fastq = select_first([renameFastqs.modFastqs,backExtractByRG.rawFastqs])
      Array[Array[File]] fastqc = reviewFastq.fastqcFiles
    } 
}


task unsortBam {
        input {
            File bamFileToSort

            Int memory = 24
            Int timeout = 12
            String modules = "samtools/1.16.1"
        }

        parameter_meta {
            bamFile: "A BAM file with one or more readgroups"
            modules: "Required environment modules"
            memory: "Memory allocated for this job"
            timeout: "Time in hours before task timeout"
        }

        command <<<
            samtools sort -n ~{bamFileToSort} -O bam -o unsorted.bam
        >>>

        runtime {
            modules: "~{modules}"
            memory: "~{memory}G"
            timeout: "~{timeout}"
        }

        output {
            File bamFile = "unsorted.bam"
        }

        meta {
            output_meta: {
                bamFile: "unsorted bam file (sorted by readname)"
            }
        }

}


task examineBam {
        input {
            File bamFile
            String prefix

            Int memory = 24
            Int timeout = 12
            String modules = "samtools/1.16.1"

        }

        parameter_meta {
            bamFile: "A BAM file with one or more readgroups"
            prefix: "String prepended to flagstat file"
            modules: "Required environment modules"
            memory: "Memory allocated for this job"
            timeout: "Time in hours before task timeout"
        }

        command <<<
            samtools flagstat ~{bamFile} > ~{prefix}.flagstat.txt
            samtools stats ~{bamFile} > ~{prefix}.samstats.txt
            samtools view -H ~{bamFile} | grep -G "^@RG" > ~{prefix}.readgroups.tsv
        >>>

        runtime {
            modules: "~{modules}"
            memory: "~{memory}G"
            timeout: "~{timeout}"
        }

        output {

            File samflagstats = "~{prefix}.flagstat.txt"
            File samstats = "~{prefix}.samstats.txt"
            File readGroups = "~{prefix}.readgroups.tsv"

        }

        meta {
            output_meta: {
                samflagstats: "A TXT file containing flag information about the BAM file",
                samstats : "A TXT file generated by samtools stats",
                readGroups: "A TSV file containing information about the merged BAM file"
            }
        }
}

task nameCheck {
        input {
            File readGroups
            String fileName

            String modules = ""
            Int memory = 24
            Int timeout = 12
        }

        parameter_meta {
            readGroups: "A TXT file containing information about the merged BAM file"
            fileName: "The naming scheme for the FASTQ files"
            modules: "Required environment modules"
            memory: "Memory allocated for this job"
            timeout: "Time in hours before task timeout"

        }

        command <<<

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

        >>>

        runtime {
            modules: "~{modules}"
            memory: "~{memory}G"
            timeout: "~{timeout}"
        }

        output {
            String rgData = read_string(stdout())
            Boolean valid = read_boolean("isValid.txt")

        }

        meta {
            output_meta: {
                rgData: "A Python dictionary containing RG data",
                valid: "A boolean value to determine workflow continuation"
            }
        }
}

task backExtractByRG {
        input {
            File bamFile

            String modules = "picard/2.21.2"
            Int memory = 24
            Int timeout = 12
        }

        parameter_meta {
            bamFile: "A BAM file with one or more readgroups"
            modules: "Required environment modules"
            memory: "Memory allocated for this job"
            timeout: "Time in hours before task timeout"
        }

        command <<<

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

        >>>

        runtime {
            modules: "~{modules}"
            memory: "~{memory}G"
            timeout: "~{timeout}"
        }

        output {

            #Array[File]? rawFastqs = glob("*.fastq.gz")
            Array[File]? rawFastqs = read_lines("outfilenames")
        }

        meta {
            output_meta: {
                rawFastqs: "An array of FASTQ files"
            }
        }
}



task backExtract {
        input {
            File bamFile

            String modules = "picard/2.21.2"
            Int memory = 24
            Int timeout = 12
        }

        parameter_meta {
            bamFile: "A BAM file with one or more readgroups"
            modules: "Required environment modules"
            memory: "Memory allocated for this job"
            timeout: "Time in hours before task timeout"
        }

        command <<<

            java -Xmx20g -jar $PICARD_ROOT/picard.jar SamToFastq \
            INPUT=~{bamFile} \
            FASTQ=R1.fastq \
            SECOND_END_FASTQ=R2.fastq \
            UNPAIRED_FASTQ=unpaired.fastq \
            NON_PF=true \
            RE_REVERSE=true \
            VALIDATION_STRINGENCY=LENIENT

        >>>

        runtime {
            modules: "~{modules}"
            memory: "~{memory}G"
            timeout: "~{timeout}"
        }

        output {
            File Read1 = "R1.fastq"
            File Read2 = "R2.fastq"
            File Unpaired = "unpaired.fastq"
        }

        meta {
            output_meta: {
                rawFastqs: "An array of FASTQ files"
            }
        }
}




task renameFastqs {
        input {
            Array[File]? rawFastqs
            String rgData
            String modules = ""
            Int memory = 24
            Int timeout = 12
        }

        parameter_meta {
            rawFastqs: "An array of FASTQ files"
            rgData: "A Python array of dictionaries containing RG data"
            modules: "Required environment modules"
            memory: "Memory allocated for this job"
            timeout: "Time in hours before task timeout"
        }

        command <<<
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

        >>>

        runtime {
            modules: "~{modules}"
            memory: "~{memory}G"
            timeout: "~{timeout}"
        }

        output {

            Array[File]? modFastqs = read_lines("outfilenames")

        }

        meta {
            output_meta: {
                modFastqs: "FASTQs renamed in accordance to input"

            }
        }
}




task reviewFastq {
        input {
            File fastq
            String id
            String modules = "fastqc/0.11.9"
            Int memory = 24
            Int timeout = 12
        }

        parameter_meta {
            fastq: "the fastq file to review"
            id: "the expected id for the report, which should be the basename of the fastq file" 
            modules: "Required environment modules"
            memory: "Memory allocated for this job"
            timeout: "Time in hours before task timeout"
        }

        command <<<
            ## run fastqc, this will generate an html report and zipped files, all in the working directory
            fastqc ~{fastq} --outdir .
            ### locate the data file in the zipped output
            datafile=`unzip -l ~{id}_fastqc.zip | grep fastqc_data | sed 's/.* //'`
            ### pull the text data out to a file
            unzip -p ~{id}_fastqc.zip $datafile > ~{id}_fastqc_data.txt

            ls *fastqc* > outfilenames
        >>>

        runtime {
            modules: "~{modules}"
            memory: "~{memory}G"
            timeout: "~{timeout}"
        }

        output {
            #File fastqcReport = "~{id}_fastqc.html"
            #File fastqcZip = "~{id}_fastqc.zip"
            #File fastqcData = "~{id}_fastqc_data.txt"
            Array[File] fastqcFiles = read_lines("outfilenames")
        }

        meta {
            output_meta: {
             fastqcFiles: "FASTQC output, html report, txt data file and zipped content"
            }
        }
}
