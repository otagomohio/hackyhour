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









