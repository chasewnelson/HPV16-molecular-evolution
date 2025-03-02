# HPV16 Molecular Evolution
Supplementary material for study of HPV16 molecular evolution

* [Custom Scripts](#custom-scripts)
	* [`BED_reduce_cov_to_seqlength.py`](#BED_reduce_cov_to_seqlength-py)
	* [`clade_assigner.py`](#clade_assigner-py)
	* [`evobioinfo.py`](#evobioinfo-py)
	* [`execute_clade_assign_replicates.py`](#execute_clade_assign_replicates-py)
	* [`extract_fasta_by_sites.pl`](#extract_fasta_by_sites-pl)
	* [`FASTA_add_metadata.py`](#FASTA_add_metadata-py)
	* [`FASTA_exclude_seqs_by_name.py`](#FASTA_exclude_seqs_by_name-py)
	* [`FASTA_extract_seqs_by_name.py`](#FASTA_extract_seqs_by_name-py)
	* [`FASTA_group_defining_sites.py`](#FASTA_group_defining_sites-py)
	* [`FASTA_inspect.py`](#FASTA_inspect-py)
	* [`FASTA_mask_sites.py`](#FASTA_mask_sites-py)
	* [`FASTA_mask_sites_by_coverage.py`](#FASTA_mask_sites_by_coverage-py)
	* [`FASTA_to_consensus.py`](#FASTA_to_consensus-py)
	* [`FASTA_to_haplotypes.py`](#FASTA_to_haplotypes-py)
	* [`generate_random_protein.py`](#generate_random_protein-py)
	* [`homoplasy_identifier.py`](#homoplasy_identifier-py)
	* [`mean_of_tables.py`](#mean_of_tables-py)
	* [`raxml-ng-log-process.py`](#raxml-ng-log-process-py)
	* [`substitute_string_set.py`](#substitute_string_set-py)
	* [`translate_FASTA.py`](#translate_FASTA-py)
	* [`VCF_splitter.py`](#VCF_splitter-py)

* [Command Line Reproducibility](#command-line)
	* [`HyPhy`](#HyPhy)
	* [`NetMHCpan-4.1`](#NetMHCpan-4.1)
	* [`OLGenie`](#OLGenie)
	* [`raxml-ng`](#raxml-ng)
	* [`SNPGenie`](#SNPGenie)

* [Acknowledgments](#acknowledgments)

* [Citation](#citation)

* [Contact](#contact)


## <a name="custom-scripts"></a>Custom Scripts

These scripts are intended to be executed from the bash command line with the specified arguments. Call with `--help` to examine all arguments and input types.

### <a name="BED_reduce_cov_to_seqlength-py"></a>BED\_reduce\_cov\_to\_seqlength.py

For BED files with ranges wrapping the circular sequence end, reduce coverage values to sequence length. Example:

	BED_reduce_cov_to_seqlength.py -b "." -s 7906 -c HPV16_Ref > BED_reduce_cov_to_seqlength.out

### <a name="clade_assigner-py"></a>clade\_assigner.py

Determine the mutually exclusive members of each clade in the tree. Example:

	clade_assigner.py -t plausible_tree_set.nw -r rep_IDs.txt -p 1000 -o output_prefix -R 5000 -T 0.95 -m 0.5

Notes: `plausible_tree_set.nw` contains the 200 final maximum likelihood trees; `rep_IDs.txt` contains exactly one sequence ID known to belong to each of the four lineages (A = PAP266425, B = PAP3508, C = IRC200686, D = PAP139245); `-p` 1000 indicates that each tree should be randomly rooted at one of the representatives and pruned 1000 times; `-R` 5000 increases the system recursion limit to accommodate large trees; `-T` 0.95 indicates that the minimum proportion of supporting permutations required to classify a sequence into a lineage is 95%; and `-m` 0.5 indicates that the maximum proportion of inconclusive permutations allowed for a sequence to be assigned to a lineage is 50%.

### <a name="evobioinfo-py"></a>evobioinfo.py

Module file containing some variables and functions imported by other scripts. Example:

	from evobioinfo import GAPS, hamming, summary_string

### <a name="execute_clade_assign_replicates-py"></a>execute\_clade\_assign\_replicates.py

Execute simulation replicates using `clade_assigner.py`. Example:

	execute_clade_assign_replicates.py -r 1000 -p 1000 -t five_clade_example.newick -R five_clade_example_reps.txt

### <a name="extract_fasta_by_sites-pl"></a>extract\_fasta\_by\_sites.pl

Output one FASTA file for each CDS record in a GTF file. Example:

	extract_fasta_by_sites.pl genomes.fasta ORFs.gtf

### <a name="FASTA_add_metadata-py"></a>FASTA\_add\_metadata.py

Add tabular metadata to FASTA sequence headers. Example:

	FASTA_add_metadata.py --seq_file=seqs.fasta --meta_file=metadata.tsv --join_key="sample_ID" --out_file=seqs_wMeta.fa > FASTA_add_metadata.out

### <a name="FASTA_exclude_seqs_by_name-py"></a>FASTA\_exclude\_seqs\_by\_name.py

Exclude a given set of sequences from a FASTA file. Example:

	FASTA_exclude_seqs_by_name.py --seq_file=seqs.fasta --seq_names=exclusions.txt --out_file=seqs_filtered.fasta > FASTA_exclude_seqs_by_name.out

### <a name="FASTA_extract_seqs_by_name-py"></a>FASTA\_extract\_seqs\_by\_name.py

Extract a given set of sequences from a FASTA file. Example:

	FASTA_extract_seqs_by_name.py -i seqs.fasta -I names.txt -o out.fasta

### <a name="FASTA_group_defining_sites-py"></a>FASTA\_group\_defining\_sites.py

Determine group-defining variants for a range of within-group variant frequencies. Example:

	FASTA_group_defining_sites.py -i A1A2_seqs.fasta A3_seqs.fasta A4_seqs.fasta B_seqs.fasta C_seqs.fasta D2D3_seqs.fasta --group_key "sublineage" --min_freq 0.9432 --max_freq 0.9432 --step_size 0 --out_dir sublineage_def --min_def_count 6 --custom_sites sites_of_interest.txt --tree "(((A1A2,A3),A4),(B,(C,D2D3)));" > sublineage_def.log

### <a name="FASTA_inspect-py"></a>FASTA\_inspect.py

Inspect a FASTA file for basic metrics of interest and lack of redundancy. Example:

	FASTA_inspect.py --seq_file=seqs.fasta --p_dist > seqs_inspect.out


### <a name="FASTA_mask_sites-py"></a>FASTA\_mask\_sites.py

Mask specified range(s) of sites in a FASTA file. Example:

	FASTA_mask_sites.py -i seqs.fasta -o seqs_masked.fasta -s 256,257,1055,1297,1936,1939,2949,3193,3313,3798,6471,6474,7895,7896 -e 256,257,1055,1297,1936,1939,2949,3193,3313,3798,6471,6474,7895,7896 > FASTA_mask_sites.out

### <a name="FASTA_mask_sites_by_coverage-py"></a>FASTA\_mask\_sites\_by\_coverage.py

Mask sites in a FASTA file based on coverage values from a BED file. Example:

	FASTA_mask_sites_by_coverage.py -s seqs.fasta -b "BED_file_dir/" -e ".bed" -o seqs_mincov10.fasta -m 10 -c HPV16REF -k N > seqs_mincov10.out

### <a name="FASTA_to_consensus-py"></a>FASTA\_to\_consensus.py

Determine the consensus sequence (majority allele at each site) for a MSA. Example:

	FASTA_to_consensus.py --aln_file alignment.fasta 2>&1 > alignment_consensus.out

### <a name="FASTA_to_haplotypes-py"></a>FASTA\_to\_haplotypes.py

Tally the unique haplotypes in a FASTA file (nucleotide or amino acid). Example:

	FASTA_to_haplotypes.py -i aligned_seqs.fasta

### <a name="generate_random_protein-py"></a>generate\_random\_protein.py

Generate randomized peptides from an input amino acid sequence. Example:

	generate_random_protein.py -i aa_aln.fasta -o aa_aln_random_l1000_n1.fasta -l 1000 -n 1 -p Lineage_ORF_ > aa_aln_random_l1000_n1.out

### <a name="homoplasy_identifier-py"></a>homoplasy\_identifier.py

Determine nucleotide or amino acid sites with homoplasies from a tree and an alignment of sequences. Example:

	homoplasy_identifier.py --fasta_file seqs.fasta --tree_file raxml.bestTree --out_file raxml.bestTree-homoplasies.txt --seq_type n > raxml.bestTree-homoplasies.log

### <a name="mean_of_tables-py"></a>mean\_of\_tables.py

Compute the mean values across all numeric cells for multiple identically formatted tables. Example:

	mean_of_tables.py -i props.tsv -k seq_name -o means.tsv

### <a name="raxml-ng-log-process-py"></a>raxml-ng-log-process.py

Process raxml-ng log files to extract ML tree information. Example:

	raxml-ng-log-process.py -i job.out -o job_metadata.txt

### <a name="substitute_string_set-py"></a>substitute\_string\_set.py

Substitute all occurrences of one string set with another in a file. Example:

	substitute_string_set.py -i seqs.fasta -m string_table.txt -f find_string -r replacement_string -o output.fasta -E ","

### <a name="translate_FASTA-py"></a>translate\_FASTA.py

Translate a nucleotide FASTA file into amino acid sequences. Example:

	translate_FASTA.py -i ORF_seqs.fasta

### <a name="VCF_splitter-py"></a>VCF\_splitter.py

Split records having multiple ALT values onto their own lines. Example:

	VCF_splitter.py -i vcf_table.txt -o vcf_table_split.txt > log.out

## <a name="command-line-reproducibility"></a>Command-Line Reproducibility

### <a name="HyPhy"></a>HyPhy

	hyphy remove-duplicates.bf --msa this_lineage_ORF_seqs.fasta --tree this_lineage.raxml.bestTree --output this_lineage_ORF_uniques.nxh 2>&1 | tee this_lineage_ORF_uniques.out 
	
	mpirun -np $SLURM_NTASKS HYPHYMPI fel --alignment this_lineage_ORF_uniques.nxh --branches Internal ENV="TOLERATE_NUMERICAL_ERRORS=1;" 2>&1 | tee this_lineage_ORF_uniques.nxh.FEL.out

	mpirun -np $SLURM_NTASKS HYPHYMPI meme --alignment this_lineage_ORF_uniques.nxh --branches Internal ENV="TOLERATE_NUMERICAL_ERRORS=1;" 2>&1 | tee this_lineage_ORF_uniques.nxh.MEME.out

### <a name="NetMHCpan-4.1"></a>NetMHCpan-4.1

	netMHCpan -f aa_seqs.fasta -a HLA-A01:01,HLA-A02:01,HLA-A03:01,HLA-A11:01,HLA-A23:01,HLA-A24:02,HLA-A26:01,HLA-B07:02,HLA-B08:01,HLA-B15:01,HLA-B27:05,HLA-B35:01,HLA-B39:01,HLA-B40:01,HLA-B44:02,HLA-B44:03,HLA-B58:01,HLA-C03:03,HLA-C04:01,HLA-C05:01,HLA-C06:02,HLA-C07:01,HLA-C07:02,HLA-C08:02,HLA-C12:03,HLA-C14:02,HLA-C15:02 -l 9 -BA

### <a name="OLGenie"></a>OLGenie

	OLGenie.pl --fasta_file ref_ORF_alignment.fasta --frame ss13 --output_file ref_ORF-OLGenie-ss13.txt --verbose

### <a name="raxml-ng"></a>raxml-ng

	raxml-ng --parse --msa seqs.fasta --model GTR+G --prefix P.seqs 

	raxml-ng --search --msa P.seqs.raxml.rba --model GTR+G --prefix lineage_name --threads auto{$SLURM_CPUS_PER_TASK} --blmin 1e-9 --blopt nr_safe

### <a name="SNPGenie"></a>SNPGenie

	snpgenie_within_group.pl --fasta_file_name=genomes.fasta --gtf_file_name=ORFs.gtf --num_bootstraps=100 --procs_per_node=4

## <a name="acknowledgments"></a>Acknowledgments

This work was supported by a Research Fellowship from the National Cancer Institute (NCI), National Institutes of Health (NIH) to C.W.N. (2021-present), Lisa Mirabello group. This product is the result of work by Laurie Burdette, Lisa Mirabello, Sambit Mishra, Chase W. Nelson, Maisa Pinheiro, and Meredith Yeager.

## <a name="citation"></a>Citation

When using this software, please refer to and cite this page:

>https://github.com/chasewnelson/HPV16-molecular-evolution


## <a name="contact"></a>Contact and troubleshooting

If you have questions about our scripts or study, please first thoroughly read the documentation and in-line comments relevant to the script of interest. If these do not answer your question, please click on the <a target="_blank" href="https://github.com/chasewnelson/HPV16-molecular-evolution/issues">Issues</a> tab at the top of this page and search to see if your question has already been answered; if not, begin a new issue, so that others might benefit from the discussion.

Other queries should be addressed to the corresponding authors: 

*  Chase W. Nelson, chase.nelson <**AT**> nih <**DOT**> gov