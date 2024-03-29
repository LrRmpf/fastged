
<!-- README.md is generated from README.Rmd. Please edit that file -->

# fastged

<!-- badges: start -->
<!-- badges: end -->

The goal of fastged is to provide a fast implementation for gene expression discretization. It uses gaussian mixture modeling on the estimated signal density curve for determining the tresholds.

The gene expression discretization uses either binary classification (values 0,1) corresponding to expressed and non-expressed genes respectively or 3 classes (values -1, 0, 1) corresponding to non-expressed, normally expressed and highly expressed genes. The user can choose the option to apply deseq (across sample normalization) and tpm normalization before the discretization calculation.

## How to use?

An input file containing a gene expression matrix (Rows: Gene IDs, Columns: Sample IDs) has to be specified in console command. The other arguments are optional. Run Rscript discretization.R with option -h to show all options:

    Rscript <Path to discretization.R> -h

Console output:    

Options:
        -f CHARACTER, --file=CHARACTER
                dataset file name

        -n, --normalize
                Should the expression data be TPM normalized first? [default FALSE]

        -g CHARACTER, --genome=CHARACTER
                Options 'h' for Hsapiens, 'mm' for Mmusculus genome when no gene length file is provided for TPM normalisation. [default h]

        -v CHARACTER, --vector=CHARACTER
                gene length file name

        -p, --plot
                Should the program plot the resulting densities of each sample? [default FALSE]

        -b, --binary
                Should the program's output be a binary classification? [default FALSE]

        -o CHARACTER, --out=CHARACTER
                output file name discretization [default= NA]

        -t CHARACTER, --tpmout=CHARACTER
                output file name normalization [default= NA]

        -l CHARACTER, --lout=CHARACTER
                output file name gene lengths [default= NA]

        -h, --help
                Show this help message and exit

Example commands:

        Rscript <Path to discretization.R> -f <Path to input file>

        Rscript <Path to discretization.R> -f <Path to input file> -n -g 'mm' -p

        Rscript <Path to discretization.R> -f <Path to input file> -p -b

        Rscript <Path to discretization.R> -f <Path to input file> -p -o <Path to output file + filename>


If -p or --plot is chosen the plots are saved to the folder "/.fastged/Plots".

The file 'tpm.txt' that stores the results of the normalization can be found in the same folder as the input data if no other specification is given. If no gene length file for the tpm normalization is applied, the two genome options for the automated gene length query (Biomart) are mouse (-g 'mm') and human (-g 'h'). The resulting file 'gene_lengths.txt' is also stored in the input data folder. Gene ids in count matrix (and in length vector when not specified) are converted to ENSEMBL ids. If the gene symbol to ENSEMBL id mapping is not unique or no ENSEMBL id can be found, the gene is removed from the counts file. Removed genes and the cause of removal can be found in the log file 'deleted_genes.log'.   

Specifications for the gene length file if provided by the user: gene ids (ENSEMBL ids!) as rownames, one column containing the corresponding gene lengths. The gene ids have to be ENSEMBL Ids, no gene symbols!   

The binary classification is chosen with -b --binary.

A name an location for the output file containing the discretized matrix can be chosen with -o or --out. If none is specified the file is saved to "/.fastged/GDE/gde.txt".

## Credits

This approach of gene expression discretization is based on the following research/ projects:

Li G, Ma Q, Tang H, Paterson AH, Xu Y. QUBIC: a qualitative biclustering algorithm for analyses of gene expression data. Nucleic Acids Res. 2009;37(15):e101. doi:10.1093/nar/gkp491; QUBIC project git clone: [https://git.bioconductor.org/packages/QUBIC](https://git.bioconductor.org/packages/QUBIC)

Maria Pires Pacheco, Tamara Bintener, Dominik Ternes, Dagmar Kulms, Serge Haan, Elisabeth Letellier, Thomas Sauter,
Identifying and targeting cancer-specific metabolism with network-based drug target prediction, EBioMedicine, Volume 43,
2019, Pages 98-106, ISSN 2352-3964, [https://doi.org/10.1016/j.ebiom.2019.04.046](https://doi.org/10.1016/j.ebiom.2019.04.046).
[https://www.sciencedirect.com/science/article/pii/S2352396419302853](https://www.sciencedirect.com/science/article/pii/S2352396419302853); rFASTCORMICS project: [https://wwwen.uni.lu/research/fstc/life_sciences_research_unit/research_areas/systems_biology/software/rfastcormics](https://wwwen.uni.lu/research/fstc/life_sciences_research_unit/research_areas/systems_biology/software/)
