# Examples using Picard Tools

To see all of the tools that Picard tools go to its **[website](https://broadinstitute.github.io/picard/)**

Here are some common examples of use.

Picard Tools runs on the java language, so the beginning of every command will call the java (jar) file, followed by the specific Picard Tools command:

```
java -Xmx2g -jar /path/to/jarfile/picard.jar picardCommand
```

the -Xmx2g part above tells java how much computing resources to provide to the program. For bigger operations requiring more RAM, for example, you can bump that up to -Xmx4g or -Xmx8g. Do not increase too much without checking with someone.

## Convert SAM/BAM file to fastq

This command will extract all of the read information (sequence and quality) from your SAM or BAM file and output as fastq format. Though the command name is SamToFastq, either SAM or BAM file can be used.

```
java -Xmx2g -jar /path/to/jarfile/picard.jar SamToFastq \
INPUT=inputBamFile \
FASTQ=outputFastqFile.fq \
```

