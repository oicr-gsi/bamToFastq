version 1.0

workflow bamToFastq {
    input {
        String bamFile
    }

    parameter_meta {
        bamFile: "A BAM file with one or more readgroups"

    }

    meta {
        author: "Murto Hilali"
        email: "mhilali@oicr.on.ca"
        description: "QC workflow to assess UMI components"
        dependencies: [
        ]
        output_meta: {
        }
    }

    call task { 
        input:
    }

    call repeatedTask as taskName {
        input:
    }

    output {
    } 
}

task countFlags {
        input {
            File bamFile
            String modules = ""
            Int memory = 24
            Int timeout = 12
            String modules = "samtools/1.9"
        }

        parameter_meta {
            bamFile: "A BAM file with one or more readgroups"
            modules: "Required environment modules"
            memory: "Memory allocated for this job"
            timeout: "Time in hours before task timeout"
        }

        command <<<

            module load samtools
            samtools flagstat ~{bamFile} > ~{bamFile}.flagstat.txt
            samtools view -H ~{bamFile} | grep \"\^\@RG\/" > readgroups.txt



        >>>

        runtime {
            modules: "~{modules}"
            memory: "~{memory}G"
            timeout: "~{timeout}"
        }

        output {

            File bamFileFlagstat = "~{bamFile}.flagstat.txt"
            File readGroups = "readgroups.txt"


        }

        meta {
            output_meta: {
                bamFileFlagstat: "~{bamFile}.flagstat.txt"
            }
        }
}

task backExtract {
        input {
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

            java -Xmx20g -jar $PICARD_ROOT/picard.jar SamToFastq
            
            INPUT=$bamfile
            OUTPUT_PER_RG=true    ### this will generate a set of files for each readgroup
            OUTPUT_DIR=.
            NON_PF=true
            RE_REVERSE=true
            VALIDATION_STRINGENCY=LENIENT

        >>>

        runtime {
            modules: "~{modules}"
            memory: "~{memory}G"
            timeout: "~{timeout}"
        }

        output {

        }

        meta {
            output_meta: {
            }
        }
}

task task {
        input {
            String modules = ""
            Int memory = 24
            Int timeout = 12
        }

        parameter_meta {
            modules: "Required environment modules"
            memory: "Memory allocated for this job"
            timeout: "Time in hours before task timeout"
        }

        command <<<

        >>>

        runtime {
            modules: "~{modules}"
            memory: "~{memory}G"
            timeout: "~{timeout}"
        }

        output {

        }

        meta {
            output_meta: {
            }
        }
}