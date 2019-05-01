# GSEA

Gene Set Enrichment Analysis (GSEA) is particularly suitable and is recommended when ranks are available for all or most of the genes in the genome (e.g., for RNA-seq data). GSEA searches for pathways whose genes are enriched at the top or bottom of the ranked gene list, more so than exempted expected by chance alone. 

* To calculate enrichment score (ES) for a pathway, GSEA progressively examines genes from the top to the bottom of the ranked list, increasing the ES if a gene is part of pathway and decreasing the score otherwise. 
* These running sum values are weighted, so that enrichment in the very top- (and bottom-) ranking genes is amplified, whereas enrichment in genes with more moderate ranks are not amplified.
* THE ES score is calculated as the maximum value of the running sum and normalized relative to pathway size, resulting in a normalized enrichment score (NES) that reflects the enrichment of the pathway in the list. 


Two input options:
1. pre-ranked list
2. expression file and phenotype file

Output summary:

<img width="568" alt="Screen Shot 2019-05-01 at 9 32 39 AM" src="https://user-images.githubusercontent.com/19800554/57028494-16c96380-6bf4-11e9-8523-b79c607651f4.png">



An example enrichment plot for the top pathway in the G1 set:

![image](https://user-images.githubusercontent.com/19800554/57028992-69574f80-6bf5-11e9-9081-a944ed6ca0d5.png)

*Enrichment plot: ROSTY_CERVICAL_CANCER_PROLIFERATION_CLUSTER (Profile of the Running ES Score & Positions of GeneSet Members on the Rank Ordered List)

An example enrichment plot for the top pathway in the G2 set:

![image](https://user-images.githubusercontent.com/19800554/57029100-a28fbf80-6bf5-11e9-9a72-85c93ad67b27.png)

*Enrichment plot: LEE_NEURAL_CREST_STEM_CELL_DN(Profile of the Running ES Score & Positions of GeneSet Members on the Rank Ordered List)

![image](https://user-images.githubusercontent.com/19800554/57028930-3614c080-6bf5-11e9-9fca-61aa90ba26d4.png)

