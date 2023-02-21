## 2023-02-21
- Added outputFilePrefix as a required parameter, this will be used to better name bam summary files, expected value is the name of the bam file or another appropriate identifier
- Added samtools stats to the examineBam task; both flagstat and stats files will be provisioned out
- Added unsortBam task which will resort the bam by the read name. This will avoid the fastq files being sorted by their previous alignment positions
- there are now two backExtranction tasks
--- backExtract : this will extract all into a single file or pair of files, and unpaired reads into a separate file
--- backExtractByRG : the previous backExtract function, which uses readgroup information to generate multiple fastq pairs
----- both may be used, as backExtractByRG does not generate unpaired reads. Currently both run, but only the byRg results are carried forward
- Added reviewFastq task which runs Fastqc on fastq files, pulling out the html report, the zipped archive and also the data table, which may be used in a not yet available accounting task. This is called with a scatter/gather across all the fastq files that are generated


## 2023-02-19
- Incomplete work from coop term, completed now for current use. Modifications were required to have this functioning properly
- added COMPRESS_OUTPUTS_PER_RG to the backExtract SamToFastq command
- set the default fileNaming scheme to "{ID}"
- file renaming task is skipped if fileNaming is the default "{ID}" as this is what the command will do
- added a selectFirst in later tasks, to use either the renamed fastq or the original fastq in renaming is not needed
