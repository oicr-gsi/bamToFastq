version 1.0

workflow bamToFastq {
    input {
        File bamFile
        String fileNaming
    }

    parameter_meta {
        bamFile: "A BAM file with one or more readgroups"
        fileNaming: "The naming scheme for the extracted FASTQs"

    }

    meta {
        author: "Murto Hilali"
        email: "mhilali@oicr.on.ca"
        description: "QC workflow to assess UMI components"
        dependencies: [
        ]
        output_meta: {
            bamFileFlagstat: "A TXT file containing flag information about the BAM file",
            readGroups: "A TXT file containing information about the merged BAM file"
        }
    }

    call countFlags {
        input:
            bamFile = bamFile
    }

    call nameCheck {
        input:
            readGroups = countFlags.readGroups,
            fileName = fileNaming
    }

    if (nameCheck.valid == true){
        call backExtract { 
            input:
                bamFile = bamFile
        }
    }

    call renameFastqs { 
        input:
            rawFastqs = backExtract.rawFastqs,
            rgData = nameCheck.rgData
    }

    output {
        File flagStat = countFlags.bamFileFlagstat
        File readGroups = countFlags.readGroups
        Array[File]? rawFastqs = backExtract.rawFastqs
        Array[File]? modFastqs = renameFastqs.modFastqs
        String rgData = nameCheck.rgData

    } 
}

task countFlags {
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

            module load samtools
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
                readGroups: "A TXT file containing information about the merged BAM file"
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
                rgData: "A Python dictionary containing RG data"
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
            modules: "Required environment modules"
            memory: "Memory allocated for this job"
            timeout: "Time in hours before task timeout"
        }

        command <<<
            module unload cromwell ## tmp for local testing
            module unload java ## tmp for local testing
            module load picard

            java -Xmx20g -jar $PICARD_ROOT/picard.jar SamToFastq \
            INPUT=~{bamFile} \
            RG_TAG="ID" \
            OUTPUT_DIR=. \
            OUTPUT_PER_RG=true \
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

            Array[File]? rawFastqs = glob("*.fastq")

        }

        meta {
            output_meta: {
                rawFastqs: "A file array of FASTQ files"
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
            rawFastqs: "A file array of FASTQ files"
            rgData: "A Python array of dictionaries containing RG data"
            modules: "Required environment modules"
            memory: "Memory allocated for this job"
            timeout: "Time in hours before task timeout"
        }

        command <<<
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

        >>>

        runtime {
            modules: "~{modules}"
            memory: "~{memory}G"
            timeout: "~{timeout}"
        }

        output {

            Array[File]? modFastqs = glob("*.fastq.gz")

        }

        meta {
            output_meta: {
                modFastqs: "FASTQs renamed in accordance to input"

            }
        }
}
