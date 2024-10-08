version 1.0

workflow bamToFastq {
    input {
        File bamFile
        String fileNaming
        String finalOutputDirectory = ""
    }

    parameter_meta {
        bamFile: "A BAM file with one or more readgroups"
        fileNaming: "The naming scheme for the extracted FASTQs"
        finalOutputDirectory: "Optional final output directory. Copy out fastq files using a task"
    }

    meta {
        author: "Murto Hilali"
        email: "mhilali@oicr.on.ca"
        description: "Generating R1 and R2 from bam files"
        dependencies: [
            {
                name: "java/8",
                url: "https://www.java.com/en/download/manual.jsp"
            },
            {
                name: "samtools/1.14",
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
    flagStat: {
        description: "A TXT file containing flag information about the BAM file",
        vidarr_label: "flagStat"
    },
    modFastqs: {
        description: "An optional array of fastq.gz files",
        vidarr_label: "modFastqs"
    },
    copyLog: {
        description: "log file from copy out task, provisioned if a custom output dir requested",
        vidarr_label: "copyLog"
    }
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

    if (finalOutputDirectory != "") {
      call renameFastqs as haveCustomDirRename {
        input:
            rawFastqs = backExtract.rawFastqs,
            rgData = nameCheck.rgData
      }

      scatter(fq in select_first([haveCustomDirRename.modFastqs])) {
        call copyOutFastq {
           input:
               Fastq = fq,
               outputDir = finalOutputDirectory 
        }
      }

      call composeLog {
        input:
            messages = copyOutFastq.message
      }

    } 

    if (finalOutputDirectory == "") {
      call renameFastqs as noCustomDirRename {
        input:
            rawFastqs = backExtract.rawFastqs,
            rgData = nameCheck.rgData
      }
    }

    output {
        File flagStat = countFlags.bamFileFlagstat
        Array[File]? modFastqs = noCustomDirRename.modFastqs
        File? copyLog = composeLog.log
    }
}

task countFlags {
        input {
            File bamFile
            String prefix = "output"
            Int memory = 24
            Int timeout = 12
            String modules = "samtools/1.14"

        }

        parameter_meta {
            bamFile: "A BAM file with one or more readgroups"
            prefix: "String prepended to flagstat file"
            modules: "Required environment modules"
            memory: "Memory allocated for this job"
            timeout: "Time in hours before task timeout"
        }

        command <<<
            set -euo pipefail
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
            Int memory = 24
            Int timeout = 12
        }

        parameter_meta {
            readGroups: "A TXT file containing information about the merged BAM file"
            fileName: "The naming scheme for the FASTQ files"
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

task copyOutFastq {
        input {
            File Fastq
            String outputDir
            Int memory = 24
            Int timeout = 12
        }

        parameter_meta {
            Fastq: "FASTQ file to copy"
            outputDir: "Oputput directory"
            memory: "Memory allocated for this job"
            timeout: "Time in hours before task timeout"
        }

        command <<<
          set -euxo pipefail
          if [[ -e ~{outputDir} && -d ~{outputDir} ]]; then
            cp ~{Fastq} ~{outputDir}
            echo "File ~{basename(Fastq)} copied to ~{outputDir}"
          else
            echo "Final Output Directory was not configured, ~{basename(Fastq)} was not provisioned properly"
          fi
        >>>

        runtime {
            memory: "~{memory}G"
            timeout: "~{timeout}"
        }

        output {
            String message = read_string(stdout())
        }

        meta {
            output_meta: {
                message: "Message from the copy task"
            }
        }
}

task composeLog {
        input {
            Array[String] messages
            Int memory = 4
            Int timeout = 2
        }

        parameter_meta {
            messages: "log messages from copyOut task"
            memory: "Memory allocated for this job"
            timeout: "Time in hours before task timeout"
        }

        command <<<
          python3 <<CODE
          import json
          import re
          m_lines = re.split(",", "~{sep=',' messages}")
          with open("copy_out.log", "w") as log:
              for line in m_lines:
                  log.write(line + "\n")
          log.close()
          CODE
        >>>

        runtime {
            memory: "~{memory}G"
            timeout: "~{timeout}"
        }

        output {
            File log = "copy_out.log"
        }

        meta {
            output_meta: {
                log: "Log file with all messages from copyOut task"
            }
        }


}
