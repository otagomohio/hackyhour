# Examples using samtools

Here are some examples of what you can do with samtools

## Getting help
samtools does not use a '-h' or '--help' argument to get help. Instead, just type samtools at the prompt:

```
samtools
```

This will give you a list of commands, organised by category. To get more details about any of these individual commands, just type samtools followed by the name of the command, e.g.

```
samtools view
```

## Sorting and Indexing

Many programs require that a BAM file is sorted or indexed before using. This is to speed up operations, but in some cases will make the files smaller.

To sort the file with defaults

```
samtools sort -@ 4 -O bam inputBamFile.bam > sortedBamFile.bam
```

the -@ argument instructs samtools how many processors to use. This is optional. The -O is the output file format (SAM, BAM, or CRAM)

The default sorting is by leftmost coordinates. To sort by read name use the -n argument:

```
samtools sort -@ 4 -O bam -n inputBamFile.bam > sortedBamFile.bam
```

Notice that there is no -i or -f argument to indicate the input file. It is just added after the command and options. 

To index a BAM file, the input file must be sorted first, as above. 

```
samtools index sortedBamFile.bam
```

Note that the output file does not need to be specified. Samtools will create an output index file with a .bai suffix on the end:

sortedBamFile.bam.bai

If you add the -c argument then a .csi index (used for longer reference sequences) will be added instead. This option is good to use if you have very large chromosomes (~512 Mbp).

<br/><br/>

## samtools view

The view command name is misleading; there are actually many things you can do with this command. The primary purpose is to extract reads that map to specific chromosomes. The default is to output to the terminal, but adding the linux redirect arrow ('>') puts it in a file. You have to add the -b argument to keep it in bam format (default output is sam). 

```
samtools view -b inputBamFile.bam Chromosome1 > outputBamChr1only.bam
```

Note that the Chromosome1 has to be the same name as in the header part of the file (e.g. SN:Chromosome1 LN:1000), which is the same as the sequence name in the genome or reference fasta file (e.g. >Chromosome1).

Any sequence in the reference file can be used, as long as the names match. You can use the same command to reads mapping to multiple regions, just put a space between each name:

```
samtools view -b inputBamFile.bam contig1 contig2 contig3 > outputBamW3contigs.bam
```

I have used this command to extract reads mapping to several hundred contigs at once, and it is very fast. 

Because you can specify the output format, this command is handy to convert your bam file to sam (so you can look at it).

```
samtools view -h inputBamFile.bam > outputSamFile.sam
```

The -h option will include the header component in the output sam file. If you do not include this, then just the read alignment information will be output. 

You can also use the samtools view command to filter your bam file. The following command will filter out any unmapped reads:

```
samtools view -b -F4 inputBamFile.bam > outputFilteredBamFile.bam
```





<br/><br/>


