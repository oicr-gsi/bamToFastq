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

    call backExtract { 
        input:
            bamFile = bamFile
    }

    call renameFastqs { 
        input:
            fileNaming = fileNaming,
            readGroups = countFlags.readGroups,
            rawFastqs = backExtract.rawFastqs
    }

    #call task { 
    #    input:
    #}

    #call repeatedTask as taskName {
    #    input:
    #}

    output {
        File flagStat = countFlags.bamFileFlagstat
        File readGroups = countFlags.readGroups
        Array[File] rawFastqs = backExtract.rawFastqs
        Array[File] modFastqs = renameFastqs.modFastqs

    } 
}

task countFlags {
        input {
            File bamFile
            Int memory = 24
            Int timeout = 12
            String modules = "samtools/0.1.19"
            String prefix = "output"
        }

        parameter_meta {
            bamFile: "A BAM file with one or more readgroups"
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

            #rawFastqs=$(ls -d $PWD/*.fastq)
            #printf '%s\n' ${rawFastqs[*]}

        >>>

        runtime {
            modules: "~{modules}"
            memory: "~{memory}G"
            timeout: "~{timeout}"
        }

        output {

            #Array[File] rawFastqs = read_lines(stdout())
            Array[File] rawFastqs = glob("*.fastq")

        }

        meta {
            output_meta: {
                rawFastqs: "A file array of FASTQ files"
            }
        }
}

task renameFastqs {
        input {
            Array[File] rawFastqs
            String fileNaming
            File readGroups
            String modules = ""
            Int memory = 24
            Int timeout = 12
        }

        parameter_meta {
            rawFastqs: "A file array of FASTQ files"
            fileNaming: "The naming scheme for the FASTQ files"
            modules: "Required environment modules"
            memory: "Memory allocated for this job"
            timeout: "Time in hours before task timeout"
        }

        command <<<
            python3 <<CODE
            import csv
            import os
            import re
            import difflib

            ## open readroups.tsv file
            rg_tsv = open(r"~{readGroups}")
            rg = csv.reader(rg_tsv, delimiter="\t")
            
            f = "~{sep=' ' rawFastqs}"
            fastqs = f.split()
            fileName = "~{fileNaming}"

            ## uses EG tag dictionaries to move 
            ## and rename files into current directory

            def rename(fastq, fileName, rg_dict):
            # check if input fileName is valid
                tags = list(rg_dict.keys())
                inputs = re.findall(r"[^{\}]+(?=})", fileName)

                if all(i in tags for i in inputs) == False:  
                    print("FAILED: Missing input tag")
                    d = difflib.Differ()
                    diff = d.compare(tags, inputs)
                    print('\n'.join(diff))
                    exit()

                # rename files
                path = os.getcwd()
                readNum = fastq[-7] # determines r1 or r2
                newName = fileName.format_map(rg_dict)
                formattedFileName = f"{path}/{newName}_{readNum}.fastq.gz"
                os.rename(fastq, formattedFileName)

            ## creating an array of dictionaries
            ## each element contains K-V pair for RG tags
            ## [probably could be cleaner]
            dict_array = []
            for row in rg:
                del row[0]
                for i in range(len(row)):
                    k = row[i].split(":")
                    row[i] = k
                dict_array.append({row[i][0]: row[i][1] for i in range(len(row))})

            ## searches through array of fastqs for ID matches
            ## and renames them accordingly
            for i in range(len(dict_array)):
                pattern = str(dict_array[i]["ID"])
                for fastq in fastqs:
                    if pattern in fastq:
                        rename(fastq=fastq, fileName=fileName, rg_dict=dict_array[i])
                        #one flaw: will delete files if two modFastqs have the same name
            CODE

        >>>

        runtime {
            modules: "~{modules}"
            memory: "~{memory}G"
            timeout: "~{timeout}"
        }

        output {

            Array[File] modFastqs = glob("*.fastq.gz")

        }

        meta {
            output_meta: {
                modFastqs: "FASTQs renamed in accordance to input"

            }
        }
}
