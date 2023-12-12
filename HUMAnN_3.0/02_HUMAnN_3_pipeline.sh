######################################################

# HUMAnN 3.0 - Metabolic Pathway Prediction Pipeline

######################################################

##### Create merged FASTQ files for import to HUMAnN #####

cd $HOME/humann_uniref50/merged_reads

while IFS= read -r i || [[ -n "$i" ]]
do 
cat $HOME/amrplusplus/NonHostReads/$i.non.host.R1.fastq.gz $HOME/amrplusplus/NonHostReads/$i.non.host.R2.fastq.gz > ./$i.merged.fastq.gz
done < $HOME/ERIN_samples_IDs.txt;

##### Load the HUMAnN 3.0 Database #####

module load Anaconda/3

cd $HOME

conda activate biobakery3

humann_databases --download chocophlan full $HOME/Database/humann --update-config yes;
humann_databases --download uniref uniref90_diamond $HOME/Database/humann --update-config yes;
humann_databases --download utility_mapping full $HOME/Database/humann --update-config yes


##### Run HUMAnN 3.0 to predict metabolic pathways and infer taxonomy #####

INPUT_DIR=$SCRATCH/humann_uniref50/merged_reads
OUTPUT_DIR=$SCRATCH/humann_uniref50/humann_uniref50_output

cd $HOME

conda activate biobakery3

# Identify metabolic pathways #
while IFS= read -r i || [[ -n "$i" ]]
do 
humann --input $INPUT_DIR/${i}.merged.fastq.gz --input-format fastq.gz --search-mode uniref50 --output $OUTPUT_DIR --output-basename ${i} --verbose
done < $HOME/ERIN_samples_IDs.txt;

# Infer taxonomy associated with predicted pathways #
while IFS= read -r i || [[ -n "$i" ]]
do 
humann_infer_taxonomy --input $SCRATCH/humann/humann_results/${i}_genefamilies.tsv --output $SCRATCH/humann/humann_infer_taxa/${i}_genefamilies_genus.tsv --level Genus --database uniref90-tol-lca --mode stratified
done < $HOME/ERIN_samples_IDs.txt;

# "Renorm" tables produced by HUMAnN (obtain relative abundances) #
while IFS= read -r i || [[ -n "$i" ]]
do 
humann_renorm_table --input $SCRATCH/humann/humann_infer_taxa/${i}_genefamilies_phylum.tsv --units relab --mode community --special y --output $SCRATCH/humann/humann3_GENUS_results/humann3_renorm_table/${i}_genefamilies_genus_relabd.tsv;

humann_renorm_table --input $SCRATCH/humann/humann3_GENUS_results/humann3_output/${i}_pathabundance.tsv --units relab --mode community --special y --output $SCRATCH/humann/humann3_PHYLUM_results/humann3_renorm_table/${i}_pathabundance_genus_relabd.tsv
done < $HOME/ERIN_samples_IDs.txt;

# Join tables across samples to create a comprehensive abundance table #
humann_join_tables --input $SCRATCH/humann/humann3_GENUS_results/humann3_renorm_table/ --output $SCRATCH/humann/humann3_PHYLUM_results/ERIN_genefamilies_genus_relabd_merged.tsv --file_name genefamilies_genus_relabd;
humann_join_tables --input $SCRATCH/humann/humann3_GENUS_results/humann3_renorm_table/ --output $SCRATCH/humann/humann3_PHYLUM_results/ERIN_pathabundance_genus_relabd_merged.tsv --file_name pathabundance_genus_relabd;
humann_join_tables --input $SCRATCH/humann/humann3_GENUS_results/humann3_output/ --output $SCRATCH/humann/humann3_PHYLUM_results/ERIN_pathcoverage_genus_merged.tsv --file_name pathcoverage;

# Re-group and Re-name metabolic pathways by designated Database (MetaCyc)
humann_regroup_table --input $SCRATCH/humann/humann3_GENUS_results/ERIN_genefamilies_genus_relabd_merged.tsv --groups uniref90_rxn --function sum --ungrouped Y --protected Y --output $SCRATCH/humann/humann3_GENUS_results/humann3_regroup/ERIN_genefamilies_genus_regrouped_MetaCyc.tsv;

humann_rename_table --input $SCRATCH/humann/humann3_GENUS_results/humann3_regroup/ERIN_genefamilies_genus_regrouped_MetaCyc.tsv --names metacyc-rxn --simplify --output $SCRATCH/humann/humann3_GENUS_results/humann3_rename/ERIN_genefamilies_genus_MetaCyc_rename.tsv

# Create HUMAnN-generated relative abundance barplots for pathways of interest #
humann_barplot --input $SCRATCH/humann/humann3_GENUS_results/ERIN_pathabundance_genus_relabd_CaseFollowPairs.txt \
  --last-metadata Pathogen \
  --focal-feature PWY-5675\ # Designate pathway of intereste here
  --sort metadata \ 
  --meta-colormap $SCRATCH/humann/humann3_GENUS_results/humann3_barplot/humann3_barplot_meta_colors.clr \ # Use a pre-made CLR document with RGB color codes to designate metadata coloring
  --focal-metadata Case.status \
  --output $SCRATCH/humann/humann3_GENUS_results/humann3_barplot/PWY-5675_humann_barplot.png 
  
conda deactivate

