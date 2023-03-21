# kmerDedup
Reduce genome assembly redundancy using shared mapped k-mer method.
This pipeline was inspired by khaper, but I have never seen the source code of the core program "Compress" in khaper pipeline.

Use cases
1. assembling software (flye, canu, ipa, hifiasm, masurca, wtdbg2 or any others) generate larger genome size than expected and purge_dups can not work well (no distinct depth peaks when using purge_dups or purge_haplotigs);
2. pooled sample assembly (assembling software may give you 2-4X larger genome size than expected);
3. there are several assemblies generated from different softwares, and you may want to cat then together and choose longer non-redundancy contig/scaffold set.

Notes
1. polish your assembly first if it was generated from long noisy reads;
2. don't use long noisy reads in this pipeline. You can use pacbio-hifi reads and/or ngs reads.

There are 4 steps to run this pipeline, 
1. use jellyfish2 to count and dump k-mer sequences, than filter the dump file using fa2fa.pl script; 
2. map k-mer sequences to contigs/scaffolds using bowtie2 and use BamDeal to calculate coverage;
3. use kmerDedup.pl (core program in this pipeline) to select longest contig/scaffold set;
4. (optional) use purge_dups to purge kmerDedup.pl results if the kmerDedup.pl result is still larger than expected.

TO DO: write each step clearer.
