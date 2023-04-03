# kmer analysis

## Preparing the kmer tables
In order to speed up the analysis, precalculate for each gene what kmers it includes. 
This calculation needs to be run once, and then the kmer prediction can be run for different parameters using this. 
If you're using one of the precalculated organisms you can skip this section.

 
First, Download all sequences from ensembl biomart  
```
make download_fasta dtype=3utr ORGANISM=mmusculus_gene_ensembl FASTA_FILE=biomart_mmusculus_gene_ensembl_3utr_ensembl109.fasta
```

Next, the downloaded fasta file needs to be filtered to exclude missing sequences, sequences that are too short and keep 
only one sequence per gene.
```
make create_fastas_all dtype=3utr FASTA_FILE=biomart_mmusculus_gene_ensembl_3utr_ensembl109.fasta
```

The last step in this section includes calculating the tables themselves, and then splitting them into smaller tables so 
that the downstream calculations are faster.
```
make create_kmer_folders dtype=3utr filter=_longest directory=/path/to/current/directory
# for each of the kmer lengths you're interested in analyzing
make create_kmer_data dtype=3utr filter=_longest klen=5
make split_tables dtype=3utr
```

## Running the kmer analysis

For the analysis itself a csv/tsv file with the parameters the analysis should be run on. The file should include 4 columns:
- category
- ensembl_gene_id
- param_name
- param_val

multiple parameters can be calculated at once, just make sure to name them separately in the param_name column, with 
the respective value in the param_value column. The category can be used to run th eanalysis only on specific parameters 
from the same file.

To run analysis:

1. Calculate the k-s test p-values and the effect size:

   ```
   make run_ks_test_job path_to_kmer_matrices=3utr/kmer_matrices PARAMETERS_FILE=parameters.tsv output_path=3utr/kmer_out term=half_life 
   ```

2. Create plots:

   ```
   make ks_tests_plots ks_test_output_folder=3utr/kmer_out term=half_life output_path=3utr/kmer_out/plots
   ```

Note: instead of the ensembl_gene_id column can have gene_name instead, though often this is less accurate
