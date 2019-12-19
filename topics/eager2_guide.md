# Using the Eager2 ancient DNA pipeline

*Note: this is for running the Eager2 pipeline on the departmental server (boros)*

The Eager2 pipeline is now set up on the boros server. The default settings are working and many alternate steps should work, however, I have not tested everything. At the moment, the PMD tools is throwing an error, but MapDamage (default) is working. I will keep working on it, and appreciate any and all feedback so that I can keep improving the pipeline. Below are instructions for getting started, a couple of example commands, and some tips to deal with some of the pipeline's idiosyncracies.

## Loading the Eager software

```
module load nextflow
module load conda
source activate nf-core-eager-2.0.6
```

Now you should be ready to go. If you want to run any of the pipeline components separately from the main pipeline, then just do the last two commands, and any of the programs in Eager2 should be working and of the same version. Note, when updates come along, I will add these and keep the old ones in case you need to go back. It is possible with Eager2 to specify a version (as long as it is in the system).

## Example Commands

The Eager2 pipeline uses the [**Nextflow**](https://www.nextflow.io/) workflow framework to run the pipeline (e.g. ''nextflow run''), so that is the start of all commands. After this comes the all the Eager-specific parameters.

All parameter options are described on the Eager2 Usage website:

[https://github.com/nf-core/eager/blob/master/docs/usage.md](https://github.com/nf-core/eager/blob/master/docs/usage.md)

I will continue to add additional examples for other use-cases. Feel free to contribute your own (especially if they work), to the **Issues** tab of this GitHub page. Also, if there is something specific you need, please add this to **Issues**.

## Examples

The following are examples for specific cases. Refer to the above link for more details on adjusting parameters.

### Mapping to human mitochondrion with no reference index provided

In the following command, a path to the reference MT is provided, but the parameter --saveReference is added, so all index files from this run will be saved and can be used for subsequent runs

```
 nextflow run nf-core/eager --pairedEnd \
  --reads '/path/to/data/PREFIX*_R{1,2}_001.fastq.gz' \
  --fasta '/path/to/mito/reference/rCRS_NC_012920.fasta' \
  --outdir '/path/to/output/human/mt_tests/eager_humanMT_test1_noRef' \
  --max_memory '128.GB' --max_cpus 8 -name 'eager_humanMT_test1_noRef' \
  --max_time '8.h' --saveReference
```

### Mapping to human mitochondrion with existing reference index

The saved references from the previous run can be utilised in subsequent runs.

```
nextflow run nf-core/eager --pairedEnd \
  --reads '/path/to/data/PREFIX*_R{1,2}_001.fastq.gz' \
  --fasta '/path/to/mito/reference/rCRS_NC_012920.fasta' \
  --bwa_index '/path/to/bwa_index/folder/' \
  --seq_dict '/path/to/reference/files/rCRS_NC_012920.dict' \
  --fasta_index '/path/to/reference/files/rCRS_NC_012920.fasta.fai' \
  --outdir '/path/to/output/human/mt_tests/eager_humanMT_test2_wRef' \
  --max_memory '128.GB' --max_cpus 8 -name 'eager_humanMT_test2_wRef' \
  --max_time '8.h'
```

### Example for running pipeline on modern data

The main differences here are to turn off collapsing reads and damage calculation, as well as using the bwa mem algorithm, which is more standard and probably better for modern data.

```
nextflow run nf-core/eager --pairedEnd \
  --reads '/path/to/data/PREFIX*_R{1,2}_001.fastq.gz' \
  --fasta '/path/to/mito/reference/REF_MT_SEQ.fasta' \
  --skip_collapse  \
  --saveReference  \
  --bwamem  \
  --skip_damage_calculation \
  --outdir '/path/to/output/human/mt_tests/eager_mod_humanMT_test1_noRef' \
  --max_memory '128.GB' --max_cpus 8 -name 'eager_mod_humanMT_test1_noRef' \
  --max_time '8.h'
```

### Tips for a better run

If you look at the Usage page, you will see one of the first options is ''-profile''. Do not add this for now. I have set it up to run from a conda environment. I am working on setting up a Singularity profile (similar to Docker) but for now it is running fine as per above guidelines.

Other things to note:

* Notice that most parameter arguments are enclosed in quotes. It looks like all but straight numbers need quotes. Pay attention to this if you are exploring other parameters

* The reads entry can be a little hard to understand. Note that both pairs (if you are using paired end reads) are listed in the R{1,2} part of the file. Pay attention to the suffix part of the file. I will add more examples soon so hopefully this will be more clear.

* If you use an existing bwa index, note that you only give the path to the folder as there are several files included in a bwa index. Put these files in their own folder, otherwise the program gets confused. 

* For the genome reference (''-fasta'' parameter above), the bwa index was not recognized if the suffix '.fa' or '.fna' was used. Just changing the name to '.fasta' seemed to fix this. 



I will be adding more content as testing continues. Let me know if you have any questions.





