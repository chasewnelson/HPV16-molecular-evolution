# HPV16 molecular evolution
Supplementary material for study of HPV16 molecular evolution

* [Supplementary data](#supplementary-data)
* [Supplementary scripts](#supplementary-scripts)
	* [`FASTA_add_metadata.py`](#FASTA_add_metadata-py).

* [Figures](#figures)
	* [**Figure X**. Informative title](#figure-X).
		* `name.bash`
		* `name.R`
* [Acknowledgments](#acknowledgments)
* [Citation](#citation)
* [Contact](#contact)

## <a name="supplementary-data"></a>Supplementary data

All supplementary are available at <a target="_blank" href="https://zenodo.org">Zenodo</a> under record ID `XXX`.

For easy access, the four most important supplementary data files are available in the `/data/` directory of this repository:

1. `name.fasta`: whole-genome multiple sequence alignment of *n*=X genomes of the species *X*. See [manuscript](#citation) for details.

2. `name.gtf`: Gene Transfer Format (GTF) file giving gene positions within **name.fasta**.

3. `name2.gtf`: Gene Transfer Format (GTF) file giving gene positions within the <a target="_blank" href="https://www.ncbi.nlm.nih.gov/nuccore/">NCBI</a> genome. Records are ordered by start site in the genome.

4. `Supplementary_Tables.xlsx`: Supplementary Tables referred to in the [manuscript](#citation).


## <a name="supplementary-scripts"></a>Supplementary scripts

Scripts are arranged by Figure, and therefore by analysis. Every effort has been made to describe all steps and input, here or in each script's comments. Where applicable, scripts are arranged in the order they should be executed. The scripts are of two types: 

1. ***Command-line***. These scripts are intended to be executed from the bash command line with the specified arguments. 

2. ***Manual analysis***. These R scripts document the bulk of our data analyses and visualizations. They are intended to be executed manually line-by-line in R/RStudio. The user should replace path names and arguments with the appropriate values for the user's analysis and directories. Attention has been drawn to lines or variables that should be modified by the user with the flag: `CHANGE THIS`.

### <a name="FASTA_add_metadata-py"></a> `FASTA_add_metadata.py`

* `FASTA_add_metadata.py` (*command-line script*)
	* **Description**. Import HPV16 sequence data and metadata to label FASTA headers.
	* **Requirements**. See import statements.
	* **Input**. 
		1. `--seq_file`: FASTA file containing multiple sequence alignment [REQUIRED]
		2. `--meta_file`: TSV file containing metadata [REQUIRED]
		3. `--join_key`: Column name from `--meta_file` containing key for data joining, expected at start of FASTA header. May provide multiple comma-separated columns names to be searched in order [REQUIRED]
		4. `--out_file`: Name of output FASTA file, with metadata added to headers [REQUIRED]
		5. `--meta_label`: Label to be used for the new text block in the FASTA header (default: METADATA) [OPTIONAL]
	* **Output**. 
		1. `STDOUT`: analysis log
		2. `--out_file`
	* **Examples**:

		1. PAP

				$ FASTA_add_metadata.py --seq_file=seq/HPV16_PAP_20200813.N-30.fasta --meta_file=seq_metadata/Lisa_hpv16meth.20190723.tsv --join_key="MERGE_ID" --out_file=seq/HPV16_PAP_20200813.N-30_wMeta.fa > FASTA_add_metadata_PAP.out

				$ FASTA_add_metadata.py --seq_file=seq/HPV16_PAP_20200813.N-30_wMeta.fa --meta_file=seq_metadata/allPAP_HPV16_sublineage_geneticancestry.tsv --join_key="PAP_ID" --out_file=seq/HPV16_PAP_20200813.N-30_wMeta2.fa --meta_label=METADATA2 > FASTA_add_metadata_PAP2.out 

		2. SUCCEED

				$ FASTA_add_metadata.py --seq_file=seq/HPV16_Succeed_20181123.N-30.fasta --meta_file=seq_metadata/HPV16_SUCCEED_pheno_lineageV3.1.tsv --join_key="ID2,ID" --out_file=seq/HPV16_Succeed_20181123.N-30_wMeta.fa > FASTA_add_metadata_SUCCEED.out

		3. IARC

				$ FASTA_add_metadata.py --seq_file=seq/HPV16_IARC_20200717.N-30.fasta --meta_file=seq_metadata/HPV16_lineages_phenoIARC2.tsv --join_key="NCI_repeat,NCI Sample ID" --out_file=seq/HPV16_IARC_20200717.N-30_wMeta.fa > FASTA_add_metadata_IARC.out



## <a name="figures"></a>Figures

### <a name="figure-X"></a>Figure X. Description

* `name.bash` (*command-line script*)
	* **Description**. Produces Figure X using PyGenomeTracks.
	* **Requirements**. PyGenomeTracks, Seqkit.
	* **Input**. The following files need to be in the working directory: 
		1. `name.fasta`: multiple alignment `.fasta` file of n=X genomes of the species X.
		2. `name.gtf`: Gene Transfer Format (GTF) file giving gene positions within `X`.
		3. `name.nw`: Newick tree for alignment
		4. `name.txt`
		5. `name2.txt`
	* **Output**. 
		1. `FigX.png`
	* **Example**:

			figX.bash

* `name.R` (*manual analysis script*)
	* **Description**. Analyze X to produce Figure X.
	* **Requirements**. tidyverse, RColorBrewer, ggrepel, patchwork.
	* **Input**. 
		1. `name.txt`, produced by `x`.
	* **Output**. 
		1. `x.tsv`, description.

## <a name="acknowledgments"></a>Acknowledgments

This work was supported by a Research Fellowship from the National Cancer Institute (NCI), National Institutes of Health (NIH) to C.W.N. (2021-present), Lisa Mirabello group. This product is the result of work by Laurie Burdette, Lisa Mirabello, Sambit Mishra, Chase W. Nelson, Maisa Pinheiro, and Meredith Yeager.

## <a name="citation"></a>Citation

When using this software, please refer to and cite this page:

>https://github.com/chasewnelson/HPV16-molecular-evolution


## <a name="contact"></a>Contact and troubleshooting

If you have questions about our scripts or study, please first thoroughly read the documentation and in-line comments relevant to the script of interest. If these do not answer your question, please click on the <a target="_blank" href="https://github.com/chasewnelson/HPV16-molecular-evolution/issues">Issues</a> tab at the top of this page and search to see if your question has already been answered; if not, begin a new issue, so that others might benefit from the discussion.

Other queries should be addressed to the corresponding authors: 

*  Chase W. Nelson, chase.nelson <**AT**> nih <**DOT**> gov