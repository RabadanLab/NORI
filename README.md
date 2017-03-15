# NORI

## Running the NORI pipeline

#### _Rabadan lab_

#### _2017-02-19_

This vignette will demonstrate the usage of the NORI package.

Before starting, make sure `bedtools` and `cpat.py` are located on your system path. You will be alerted if they are not found on the system path when they are attempted to be used.

This package requires `data.table`, `doParallel`, and `Rsubread`.

# Standard operation

To run the NORI pipeline, simply supply the paths for the input GTF file, the input BAM files, and the reference genomes against which to compare:

```
    # load package
    library(NORI)

    # define paths
    gtf_path <- "./merged.gtf"                          # input GTF
    ens_path <- "./Homo_sapiens.GRCh37.75.gtf"          # Ensembl reference
    refseq_path <- "./RefSeq-GRCh37-hg19.bed"           # RefSeq reference
    pseudo_path <- "./Human74.bed"                      # Pseudogene build
    ref_fa <- "./hg19.fa"                               # Reference assembly for CPAT
    bam_paths <- list.files(path = ".",
        pattern = "[.]bam$", full.names = TRUE)     # character vector of paths to BAM files

    # run NORI
    NORI(gtf_path, ens_path, refseq_path, pseudo_path, ref_fa_path, bam_paths)

```

NORI saves the filtered BED file in the current working directory by default. For convenience, NORI also saves a separate BED file at each step as well as a log file.

Note that at this time the Ensembl input must be in GTF format. The RefSeq input may be in BED (preferred) or GTF format.

# Optional arguments

There are a number of optional arguments that may be passed. The most important arguments will be covered here. See `?NORI` for more information.

## Output file paths

File paths for each step are returned invisibly and can be captured by storing the output as a variable

```
    lnc_list <- NORI(gtf_path, ens_path, refseq_path, pseudo_path, ref_fa, bam_paths)
```

## Output and logs

Outputs default to the current working directory. To change the location of the output files and the logs, simply specify them as follows:

```
    NORI(gtf_path, ens_path, refseq_path, pseudo_path, ref_fa, bam_paths,
        output_dir = "./output", logs_path = "./output/logs.txt")
```

## Skipping steps

Any step can easily be skipped as follows:

```
    # skip Pseudogene step
    NORI(gtf_path, ens_path, refseq_path, pseudo_path, ref_fa, bam_paths, filter_pseudogenes = FALSE)
```

… and all previously completed steps can be skipped by specifying `overwrite = FALSE`

```
    # skip existing files
    NORI(gtf_path, ens_path, refseq_path, pseudo_path, ref_fa, bam_paths, overwrite = FALSE)
```
# Step-specific options

## CPAT

By default, NORI uses the recommended CPAT cutoff value of 0.364\. You may specify your own cutoff as follows:

```
    NORI(gtf_path, ens_path, refseq_path, pseudo_path, ref_fa, bam_paths,
        cpat_cutoff = 0.370)
```

For convenience, the package includes the prebuilt human hexamer frequency table and training model for CPAT. These files can be located using the following commands:

```
    logit_model <- system.file("extdata", "Human_train.RData", package = "NORI")
    hexamer_dat <- system.file("extdata", "Human_Hexamer.tab", package = "NORI")
```

`NORI` automatically uses these files by default. You may pass your own frequency table and training model by simply specifying the paths to your files:

```
    NORI(gtf_path, ens_path, refseq_path, pseudo_path, ref_fa, bam_paths,
        logit_model = "./logit_model.RData", hexamer_dat = "./hexamer_dat.tab")
```

See the [CPAT documentation](http://rna-cpat.sourceforge.net/) for more information.

## RPKM

The threshold for RPKM defaults to <span class="math inline">\(0.05\cdot n\)</span>, where <span class="math inline">\(n\)</span> is the number of samples.

This threshold can easily be changed by passing the desired value:

```
    num_samples <- length(bam_paths)
    new_threshold <- 0.01 * num_samples

    NORI(gtf_path, ens_path, refseq_path, pseudo_path, ref_fa, bam_paths,
        rpkm_threshold = new_threshold)
```

<div id="computing-rpkm-separately-using-a-sun-grid-engine-cluster" class="section level3">

### Computing RPKM separately using a Sun Grid Engine cluster

The RPKM step calculates read counts in series using the featureCounts function from the Rsubread package. Although the number of threads can be specified as follows:

```
    NORI(gtf_path, ens_path, refseq_path, pseudo_path, ref_fa, bam_paths, nthreads = 4)
```

… computing the read counts in this manner can be slow with a large number of BAM files.

If you have access to a Sun Grid Engine cluster and would instead prefer to parallelize the RPKM step by submitting individual read counts jobs to the cluster (and, optionally, aggregating the output in parallel afterwards using the doParallel package), you may do so as follows (note that this requires you to have `Subread` on your system path):

```
    # run the NORI pipeline up to the RPKM step
    # we do not need to pass bam_paths since RPKM is not calculated here
    lnc_list <- NORI(gtf_path, ens_path, refseq_path, pseudo_path, ref_fa, filter_rpkm = FALSE)
```

Next, submit jobs to the Sun Grid Engine cluster to calculate read counts in parallel.

```
    # choose output directory for read counts
    read_counts_output_dir <- "./read_counts"

    # select last BED path to be used for read counts
    penultimate <- lnc_list[length(lnc_list)]

    # calculate read counts
    read_counts_sge <- function(penultimate, bam_paths, read_counts_output_dir) {
```

(see `?read_counts_sge` for optional arguments)

After the SGE jobs complete and the read counts have been calculated, calculate RPKM and filter according to chosen threshold.

```
    # after the jobs complete, get other components of RPKM and filter
    #     according to passed or default threshold
    num_cores <- 4 # optional
    rpkm_threshold <- 0.05 * length(bam_paths) # optional

    # select only read counts outputs (not summaries)
    # adjust the contents of list.files() according to your files
    read_counts_paths <- list.files(read_counts_output_dir, pattern = "TCGA", full.names = TRUE)
    read_counts_paths <- read_counts_paths[!grepl(".summary", read_counts_paths)]

    # finally, aggregate outputs to determine RPKM and filter according to chosen threshold
    rpkm <- calculate_rpkm_sge(read_counts_paths, num_cores = num_cores, write_rpkm_to_file = TRUE)
    final <- filter_rpkm(rpkm, penultimate, rpkm_threshold)
```
