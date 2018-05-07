module add UHTS/Assembler/trinityrnaseq/2.4.0
module add UHTS/Analysis/kisto/0.43.0
module add R/3.3.2

/software/UHTS/Assembler/trinityrnaseq/2.1.1/util/abundance_estimates_to_matrix.pl \
--est_method kisto  --out_prefix Tv_ \
--name_sample_by_basedir \
Tv1_231/abundance.tsv \
Tv1_232/abundance.tsv \
Tv1_271/abundance.tsv \
Tv1_272/abundance.tsv \
Tv1_314/abundance.tsv \
Tv1_317/abundance.tsv \
Tv1_431/abundance.tsv \
Tv1_432/abundance.tsv \
Tv1_462/abundance.tsv \
Tv1_463/abundance.tsv \
Tv2_231/abundance.tsv \
Tv2_232/abundance.tsv \
Tv2_273/abundance.tsv \
Tv2_275/abundance.tsv \
Tv2_312/abundance.tsv \
Tv2_318/abundance.tsv \
Tv2_319/abundance.tsv \
Tv2_431/abundance.tsv \
Tv2_461/abundance.tsv \
Tv2_462/abundance.tsv \
Tv2_464/abundance.tsv \
Tv3_431/abundance.tsv \
Tv5_233/abundance.tsv \
Tv5_234/abundance.tsv \
Tv5_271/abundance.tsv \
Tv5_273/abundance.tsv \
Tv5_312/abundance.tsv \
Tv5_313/abundance.tsv \
Tv5_314/abundance.tsv \
Tv5_431/abundance.tsv \
Tv5_432/abundance.tsv \
Tv5_433/abundance.tsv \
Tv5_434/abundance.tsv \
Tv5_461/abundance.tsv \
Tv6_234/abundance.tsv \
Tv6_235/abundance.tsv \
Tv6_272/abundance.tsv \
Tv6_273/abundance.tsv \
Tv6_312/abundance.tsv \
Tv6_313/abundance.tsv \
Tv6_316/abundance.tsv \
Tv6_318/abundance.tsv \
Tv6_431/abundance.tsv \
Tv6_432/abundance.tsv \
Tv6_433/abundance.tsv \
Tv6_461/abundance.tsv
