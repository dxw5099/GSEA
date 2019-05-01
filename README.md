# GSEA

Gene Set Enrichment Analysis (GSEA) is particularly suitable and is recommended when ranks are available for all or most of the genes in the genome (e.g., for RNA-seq data).

GSEA searches for pathways whose genes are enriched at the top or bottom of the ranked gene list, more so than exempted expected by chance alone.It calculates enrichment score (ES) for a pathway, GSEA progressively examines genes from the top to the bottom of the ranked list, increasing the ES if a gene is part of pathway and decreasing the score otherwise. THE ES score is calculated as the maximum value of the running sum and normalized relative to pathway size, resulting in a normalized enrichment score (NES) that reflects the enrichment of the pathway in the list. 
