working_dir <- file.path(getwd(), "test")
#Load Expression Data

data=read.table(file.path(working_dir, "G1_vs_G2_DEGs_all.csv"),sep=',',header=T,check=F,comment.char="")

gene_list <- gsub("\\..*","",data$Row.names) #remove all after "."
data$Row.names <- gene_list
colnames(data)[1]<-"Gene"
colnames(data)[3]<-"log2FC"
colnames(data)[5]<-"p_adj"
head(data)
# add gene symbol column
data$symbol = mapIds(org.Hs.eg.db,
                     keys=data$Gene, 
                     column="SYMBOL",
                     keytype="ENSEMBL",
                     multiVals="first")
#add Entrez ID column
data$entrez = mapIds(org.Hs.eg.db,
                     keys=data$Gene, 
                     column="ENTREZID",
                     keytype="ENSEMBL",
                     multiVals="first")
#add gene name column
data$name =   mapIds(org.Hs.eg.db,
                     keys=data$Gene, 
                     column="GENENAME",
                     keytype="ENSEMBL",
                     multiVals="first")
#creat ranks file
# calculate ranks
ranks_RNAseq = sign(data$log2FC) * -log10(data$pvalue)
#Create a two-column rank (.RNK) file of all gene IDs and corresponding scores to for GSEA pre-ranked analysis
ranks_RNAseq <- cbind(data$symbol, ranks_RNAseq)
colnames(ranks_RNAseq) <- c("GeneName","rank")
#remove rows is NA
ranks_RNAseq_wo_NA <- ranks_RNAseq[complete.cases(ranks_RNAseq[,2]),]

#sort ranks in decreasing order
ranks_RNAseq_final <- ranks_RNAseq_wo_NA[order(as.numeric(ranks_RNAseq_wo_NA[,2]),decreasing = TRUE),]
rnk_file <- file.path(working_dir,"RNASeq_ranks.rnk")
write.table(ranks_RNAseq_final, rnk_file,col.name = TRUE, sep="\t", row.names = FALSE, quote = FALSE)

#Creat input expression file for GSEA
GSEA_expression <- data[, c(12, 14,6:11)]
colnames(GSEA_expression)[1] <- "Name"
colnames(GSEA_expression)[2] <- "Description"
input_file <- file.path(working_dir,"RNAseq_expression_input.txt")
write.table(GSEA_expression, input_file,col.name = TRUE, sep="\t",row.names = FALSE, quote = FALSE)

#Creat class file for GSEA

sample_info <- read.table(file.path(working_dir, "sample_info.csv"),sep=',',header=T)
GSEA_class <- file(
  file.path(working_dir,"GSEA_classes.cls"))
writeLines(c(paste(length(sample_info$Sample_Name), length(levels(sample_info$Group)), "1"),
                         paste("# ", unique(sample_info[,"Group"])[1], " ",unique(sample_info[,"Group"])[2])), GSEA_class)
write.table(t(as.character(sample_info$Group)),
            file.path(working_dir,"GSEA_classes.cls"),
            col.name=FALSE, sep="\t",row.names=FALSE, quote=FALSE, append=TRUE)
close(GSEA_class)


########## run GSEA ##############
gsea_jar <- "/Applications/gsea-3.0.jar"
run_gsea = TRUE
gsea_directory = ""
#TODO: change this to update to the latest gmt file.
dest_gmt_file <- file.path(working_dir,"msigdb.v6.2.symbols.gmt" )
num_randomizations <- 1000
analysis_name <- "G1_vs_G2"
#Run GSEA with gene set randomization
timestamp()
start_gs_perm <- Sys.time()
### input option #1: pre-ranked gene list #######
if(run_gsea){
  command <- paste("java -Xmx8G -cp",gsea_jar, "xtools.gsea.GseaPreranked -gmx",
                   dest_gmt_file, "-rnk" ,rnk_file,
                   "-cls",file.path(working_dir,"GSEA_classes.cls"),
                   "-collapse false -nperm ",num_randomizations,
                   " -permute gene_set -scoring_scheme weighted -rpt_label ",
                   paste(analysis_name,"gsrand",sep="_"),
                   " -num 100 -plot_top_x 20 -rnd_seed 12345 -set_max 200 -set_min 15 -zip_report false -out" ,
                   file.path(working_dir,"apr30"), "-gui false > gsea_output.txt",sep=" ")
  system(command)
}

### input option #2: expression file and phenotype file #######
#if(run_gsea){
#  command <- paste("java -Xmx8G -cp",gsea_jar, "xtools.gsea.GseaPreranked -gmx",
#                   dest_gmt_file, "-rnk" ,rnk_file,
#                   "-cls",file.path(working_dir,"GSEA_classes.cls"),
#                   "-collapse false -nperm ",num_randomizations,
#                   " -permute gene_set -scoring_scheme weighted -rpt_label ",
#                   paste(analysis_name,"gsrand",sep="_"),
#                   " -num 100 -plot_top_x 20 -rnd_seed 12345 -set_max 200 -set_min 15 -zip_report false -out" ,
#                   file.path(working_dir,"apr30"), "-gui false > gsea_output.txt",sep=" ")
#  system(command)
#}

stop_gs_perm <- Sys.time()
timestamp()
