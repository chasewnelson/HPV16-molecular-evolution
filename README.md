# HPV16 molecular evolution
Supplementary material for study of HPV16 molecular evolution

* [Scripts](#scripts)
	* [`FASTA_add_metadata.py`](#FASTA_add_metadata-py).
	* [`FASTA_exclude_seqs_by_name.py`](#FASTA_exclude_seqs_by_name-py).
	* [`FASTA_inspect.py`](#FASTA_inspect-py).
* [Acknowledgments](#acknowledgments)
* [Citation](#citation)
* [Contact](#contact)


## <a name="scripts"></a>Scripts

Where applicable, scripts are arranged in the order they should be executed. The scripts are of two types: 

1. ***Command-line***. These scripts are intended to be executed from the bash command line with the specified arguments. 

2. ***Manual analysis***. These R scripts document the bulk of our data analyses and visualizations. They are intended to be executed manually line-by-line in R/RStudio. The user should replace path names and arguments with the appropriate values for the user's analysis and directories. Attention has been drawn to lines or variables that should be modified by the user with the flag: `CHANGE THIS`.

### <a name="FASTA_add_metadata-py"></a> `FASTA_add_metadata.py`

* `FASTA_add_metadata.py` (*command-line script*)
	* **Description**. Add tabular metadata to FASTA sequence headers.
	* **Requirements**. See import statements.
	* **Input**. 
		1. `--help`: show help message and exit
		2. `--seq_file`: FASTA file containing a multiple sequence alignment [REQUIRED]
		2. `--meta_file`: TSV file containing metadata [REQUIRED]
		3. `--join_key`: Column name from --meta_file containing key for data joining, expected at start of FASTA header. May provide multiple comma-separated columns names to be searched in order. Value cannot be "nan" or "NA" [REQUIRED]
		4. `--out_file`: Name of output FASTA file, with metadata added to headers [REQUIRED]
		5. `--meta_label`: Label to be used for the new text block in the FASTA header (default: METADATA) [OPTIONAL]
	* **Output**. 
		1. `STDOUT`: analysis log
		2. `--out_file`
	* **Examples**:

		1. IARC

				FASTA_add_metadata.py --seq_file=seq/HPV16_IARC_20200717.N-30.fasta --meta_file=seq_metadata/HPV16_lineages_phenoIARC2.tsv --join_key="NCI_repeat,NCI Sample ID" --out_file=seq/HPV16_IARC_20200717.N-30_wMeta.fa > FASTA_add_metadata_IARC.out
		
		2. PAP

				FASTA_add_metadata.py --seq_file=seq/HPV16_PAP_20200813.N-30.fasta --meta_file=seq_metadata/Lisa_hpv16meth.20190723.tsv --join_key="MERGE_ID" --out_file=seq/HPV16_PAP_20200813.N-30_wMeta.fa > FASTA_add_metadata_PAP.out

				FASTA_add_metadata.py --seq_file=seq/HPV16_PAP_20200813.N-30_wMeta.fa --meta_file=seq_metadata/allPAP_HPV16_sublineage_geneticancestry.tsv --join_key="PAP_ID" --out_file=seq/HPV16_PAP_20200813.N-30_wMeta2.fa --meta_label=METADATA2 > FASTA_add_metadata_PAP2.out 

		3. SUCCEED

				FASTA_add_metadata.py --seq_file=seq/HPV16_Succeed_20181123.N-30.fasta --meta_file=seq_metadata/HPV16_SUCCEED_pheno_lineageV3.1.tsv --join_key="ID2,ID" --out_file=seq/HPV16_Succeed_20181123.N-30_wMeta.fa > FASTA_add_metadata_SUCCEED.out

### <a name="FASTA_exclude_seqs_by_name-py"></a> `FASTA_exclude_seqs_by_name.py`

* `FASTA_exclude_seqs_by_name.py` (*command-line script*)
	* **Description**. Exclude a given set of sequences from a FASTA file.
	* **Requirements**. See import statements.
	* **Input**. 
		1. `--help`: show help message and exit
		2. `--seq_file`: FASTA file(s) [REQUIRED]
		3. `--seq_names`: A comma-separated list of sequence IDs (names; beginning of FASTA headers) [REQUIRED]
		4. `--out_file`: Output file prefix for printing sequences; otherwise STDOUT [OPTIONAL]
	* **Output**. 
		1. `STDOUT`: analysis log
		2. `--out_file`
	* **Examples**:

		1. IARC

				FASTA_exclude_seqs_by_name.py --seq_file=IARC_HPV16_n2062_meta.fasta --seq_names=../seq_metadata/IARC_exclusions_20220203.txt --out_file=IARC_HPV16_n2062_meta_exclude.fasta > IARC_HPV16_n2062_meta_exclude.out
		
		2. PAP

				FASTA_exclude_seqs_by_name.py --seq_file=PAP_HPV16_n3305_meta5.fasta --seq_names=../seq_metadata/PAP_exclusions_20220319.txt --out_file=PAP_HPV16_n3305_meta5_exclude.fasta > PAP_HPV16_n3305_meta5_exclude.out
				
		3. SUCCEED
		
				FASTA_exclude_seqs_by_name.py --seq_file=SUCCEED_HPV16_n1044_meta.fasta --seq_names=../seq_metadata/SUCCEED_exclusions_20220206.txt --out_file=SUCCEED_HPV16_n1044_meta_exclude.fasta > SUCCEED_HPV16_n1044_meta_exclude.out
				
### <a name="FASTA_inspect-py"></a> `FASTA_inspect.py`

* `FASTA_inspect.py` (*command-line script*)
	* **Description**. Inspect a FASTA file for basic metrics of interest and lack of redundancy.
	* **Requirements**. See import statements.
	* **Input**. 
		1. `--help`: show help message and exit
		2. `--seq_file`: FASTA file containing sequence(s) [REQUIRED]
		3. `--p_dist`: Activate with alignments to calculate p-distance between sequences with identical IDs [OPTIONAL]
		1. `STDOUT`: analysis log
	* **Examples**:

		1. IARC

				FASTA_inspect.py --seq_file=IARC_HPV16_n2062_meta_exclude_n1966.fasta --p_dist > IARC_HPV16_n2062_meta_exclude_n1966_inspect.out
		
		2. PAP

				FASTA_inspect.py --seq_file=PAP_HPV16_n3305_meta5_exclude_n3220.fasta --p_dist > PAP_HPV16_n3305_meta5_exclude_n3220_inspect.out
				
		3. SUCCEED
		
				FASTA_inspect.py --seq_file=IARC_HPV16_n2062_meta_exclude_n1966.fasta --p_dist > IARC_HPV16_n2062_meta_exclude_n1966_inspect.out

## <a name="acknowledgments"></a>Acknowledgments

This work was supported by a Research Fellowship from the National Cancer Institute (NCI), National Institutes of Health (NIH) to C.W.N. (2021-present), Lisa Mirabello group. This product is the result of work by Laurie Burdette, Lisa Mirabello, Sambit Mishra, Chase W. Nelson, Maisa Pinheiro, and Meredith Yeager.

## <a name="citation"></a>Citation

When using this software, please refer to and cite this page:

>https://github.com/chasewnelson/HPV16-molecular-evolution


## <a name="contact"></a>Contact and troubleshooting

If you have questions about our scripts or study, please first thoroughly read the documentation and in-line comments relevant to the script of interest. If these do not answer your question, please click on the <a target="_blank" href="https://github.com/chasewnelson/HPV16-molecular-evolution/issues">Issues</a> tab at the top of this page and search to see if your question has already been answered; if not, begin a new issue, so that others might benefit from the discussion.

Other queries should be addressed to the corresponding authors: 

*  Chase W. Nelson, chase.nelson <**AT**> nih <**DOT**> gov