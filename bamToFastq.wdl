version 1.0

workflow bamToFastq {
    input {
        File bamFile
        String fileNaming = "{ID}"
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
                name: "samtools/1.9",
                url: "https://github.com/samtools/samtools/archive/0.1.19.tar.gz"
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
            bamFile = bamFile
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
        call backExtract { 
            input:
                bamFile = bamFile
        }
    }

    if (fileNaming != "{ID}"){
      call renameFastqs { 
        input:
          rawFastqs = backExtract.rawFastqs,
          rgData = nameCheck.rgData
      }
    }


    output {
      File flagStat = examineBam.bamFileFlagstat
      Array[File]? fastq = select_first([renameFastqs.modFastqs,backExtract.rawFastqs])
    } 
}

task examineBam {
        input {
            File bamFile
            String prefix = "output"

            Int memory = 24
            Int timeout = 12
            String modules = "samtools/0.1.19"

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
            samtools view -H ~{bamFile} | grep -G "^@RG" > readgroups.tsv
        >>>

        runtime {
            modules: "~{modules}"
            memory: "~{memory}G"
            timeout: "~{timeout}"
        }

        output {

            File bamFileFlagstat = "~{prefix}.flagstat.txt"
            File readGroups = "readgroups.tsv"

        }

        meta {
            output_meta: {
                bamFileFlagstat: "A TXT file containing flag information about the BAM file",
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
            module unload cromwell ## <- only needed for local testing
            module unload java ## <- only needed for local testing
            module load picard

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
