# About this dataset

This dataset is a holdover from a time when this repo was intended for
visualizing seekdeep datasets. You can find the seekdeep repo here:
https://github.com/bailey-lab/seekdeep_illumina_snakemake

Because I didn't have a well-established publicly available amplicon panel at
the time, I subsampled from a publicly available MIP dataset - MIPs can be
thought of as PCR amplicons whose "primers" invert to form a circular molecule.

This dataset is a randomly sampled group of 100 samples from Tanzania, which
were originally sampled with MIPs. I selected five MIPs to serve as primers,
covering the following drug resistance loci: k13_561, mdr1_86, dhfr_164,
dhps_581, crt_76

A user can input the data here into the seekdeep repository (link above) and
then run the output table through the haplotype variant calling pipeline, here:
https://github.com/bailey-lab/haplotype_variant_calling

The final output AA table can then be visualized using the current repo.