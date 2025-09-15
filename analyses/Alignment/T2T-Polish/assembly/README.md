# Collect Assembly Summary Statistics

Collect assembly summary statistics per haplotypes.

Each sequence are assumed to have hap1 / hap2 / mat / pat as haplotypes.

Sequences with chrM or random will be ignored.


## Requirements
* Assembly: `sp_ver.fa` and `sp_ver.fa.fai`

* Merqury: Output in diploid mode, preferably on a hybrid database.

  The only files used are `prefix.asm1.qv` and `prefix.asm2.qv` - QV per sequences.

* [T2T-Ref](https://github.com/arangrhie/T2T-Ref) under `$tools`

