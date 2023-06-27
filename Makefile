
# ----------------------------------------------------------
# make sure to update all following parameters
# datatypes tested are - 3utr 5utr coding
# make sure to note in the file name what ensembl version was downloaded
# ----------------------------------------------------------

directory = /path/to/current/directory
ORGANISM = mmusculus_gene_ensembl
dtype = 3utr

FASTA_FILE = biomart_$(ORGANISM)_$(dtype)_ensembl109.fasta

# four columns category, gene_name/ensembl_gene_id, param_name and param_val
PARAMETERS_FILE = parameters.tsv

# ----------------------------------------------------------
# Prepare fastas fastas
# (1) make download_fasta
# (2) make create_fastas_all_job dtype=3utr
# ----------------------------------------------------------

download_fasta:
	mkdir $(dtype); \
	wget -O $(dtype)/$(FASTA_FILE) 'http://www.ensembl.org/biomart/martservice?query=<?xml version="1.0" encoding="UTF-8"?><!DOCTYPE Query><Query  virtualSchemaName = "default" formatter = "FASTA" header = "0" uniqueRows = "0" count = "" datasetConfigVersion = "0.6" ><Dataset name = "$(ORGANISM)" interface = "default" ><Attribute name = "ensembl_gene_id" /><Attribute name = "ensembl_gene_id_version" /><Attribute name = "ensembl_transcript_id" /><Attribute name = "ensembl_transcript_id_version" /><Attribute name = "$(dtype)" /><Attribute name = "external_gene_name" /></Dataset></Query>'

create_fastas_all:
	cd $(dtype); \
	mkdir fastas; \
	Rscript ../create_fastas.R $(FASTA_FILE); \
	cd ..

# ----------------------------------------------------------
# prepare kmer files
# (1) make create_kmer_folders dtype=3utr filter=_longest
# (2) make create_kmer_data_all
# (3) make split_tables_job
# ----------------------------------------------------------

create_kmer_folders:
	mkdir $(dtype)/kmer_matrices_tmp;
	mkdir $(dtype)/kmer_matrices;
	mkdir $(dtype)/kmer_out;
	mkdir $(dtype)/kmer_out/plots;
	ln -s $(directory)/$(dtype)/fastas/cluster_all_filtered$(filter).fa $(directory)/$(dtype)/kmer_matrices_tmp/cluster_all_filtered$(filter).fa;

create_kmer_data:
	Rscript create_kmer_matrices.R $(dtype)/kmer_matrices_tmp/cluster_all_filtered$(filter).fa $(klen)

create_kmer_data_all:
	rm -rf kmer_table.txt; \
	$(foreach k, $(shell seq 5 8), \
		echo "make create_kmer_data klen=$(k) dtype=3utr filter=_longest" >> kmer_table.txt; \
	) \
	sbatch.pl -a kmer_table.txt -R 64 -t 6:00:0 "module load hurcs; module load R4/4.1.3";

split_tables:
	Rscript split_tables.R $(dtype)/kmer_matrices

# ----------------------------------------------------------
# extract tables downloaded from directory
# ----------------------------------------------------------

unzip_kmer_files:
	cat *.csv.tar.gz | tar zxvf - -i; \
	rm -rf *.csv.tar.gz;

# ----------------------------------------------------------
# kmer analysis
# (1) make run_ks_test_job
# (2) make ks_tests_plots
# ----------------------------------------------------------

run_ks_test:
	#Rscript motif_ks_test.R $(dtype)/kmer_matrices $(PARAMETERS_FILE) $(dtype)/kmer_out $(term)
	Rscript motif_ks_test.R $(path_to_kmer_matrices) $(PARAMETERS_FILE) $(output_path) $(term)

ks_tests_plots:
	#Rscript motif_ks_plots.R $(dtype)/kmer_out/$(term)_ks_raw_with_stats.tsv $(dtype)/kmer_out/plots
	Rscript motif_ks_plots.R $(ks_test_output_folder)/$(term)_ks_raw_with_stats.tsv $(output_path)
