## Installing and loading netZooR ##

## install netZooR
install.packages("devtools")
library(devtools)
devtools::install_github("netZoo/netZooR", build_vignettes = FALSE)

## install other packages used in this tutorial
# if (!requireNamespace("BiocManager", quietly = TRUE))   
# install.packages("BiocManager",repos = "http://cran.us.r-project.org")  
# BiocManager::install("fgsea")   
# install.packages("ggplot2")  
# install.packages("reshape2")
# install.packages("visNetwork")

# load packages used in this tutorial
library('netZooR') # To load PANDA and LIONESS
library('fgsea') # For gene enrichment analysis
library('ggplot2') # For plotting
library('reshape2') # For data processing
library('visNetwork') # for network visualization


#### Input data pre-processing ####

# read in the input data
motif <- read.delim("./data/motif_subset.txt", stringsAsFactors=F, header=F);
ppi <- read.delim("./data/ppi_subset.txt", stringsAsFactors=F, header=F);
exp <- read.delim("./data/expression_tpm_lcl_blood_subset.txt", stringsAsFactors = F, check.names = F);

# pre-process the expression data and remove lowly expressed genes
exp <- log2(exp+1);
zero_na_counts <- apply(exp, MARGIN = 1, FUN = function(x) length(x[(!is.na(x) & x!=0) ]));
exp <- exp[zero_na_counts > 20,];

# subset the expression data to LCL and whole blood samples
lcl_samples <-read.delim("./data/LCL_samples.txt", header=FALSE, stringsAsFactors=FALSE);
lcl_exp <- exp[,colnames(exp) %in% lcl_samples[,1]];
wblood_samples <-read.delim("./data/WholeBlood_samples.txt", header=FALSE, stringsAsFactors=FALSE);
wb_exp <- exp[,colnames(exp) %in% wblood_samples[,1]];



#### Building a PANDA regulatory network ####
pandaLCL <- panda(motif, lcl_exp, ppi, mode="intersection");
pandaWB <- panda(motif, wb_exp, ppi, mode="intersection");


#### Visualize the 200 highest edge weights in LCL ####

nDiffs= 200 # top edges to plot (top edges with largest absolute value)
diffNet = pandaLCL@regNet;
nTFs  = dim(diffNet)[1];

# Edges data frame
edges = matrix(0L, nDiffs, 3);
colnames(edges) = c("from","to","value");
edges = as.data.frame(edges);
aa = order(as.matrix(abs(diffNet)), decreasing = TRUE);
bb = sort(as.matrix(abs(diffNet)), decreasing = TRUE);
edges$value = as.matrix(diffNet)[aa[1:nDiffs]];
geneIdsTop = (aa[1:nDiffs] %/% dim(diffNet)[1]) + 1;
tfIdsTop = aa[1:nDiffs] %% dim(diffNet)[1];
tfIdsTop[tfIdsTop == 0] = nTFs;
edges$to = colnames(diffNet)[geneIdsTop];
edges$from = rownames(diffNet)[tfIdsTop]  ;                                
edges$arrows = "to"; 
edges$value = edges$value;

# Nodes data frame
nodes = data.frame(id = unique(as.vector(as.matrix(edges[,c(1,2)]))), label=unique(as.vector(as.matrix(edges[,c(1,2)]))));
nodes$group = ifelse(nodes$id %in% edges$from, "TF", "gene");

# Plot network
net <- visNetwork(nodes, edges, width = "100%");
net <- visGroups(net, groupname = "TF", shape = "triangle", color = list(background = "purple", border="black"));
net <- visGroups(net, groupname = "gene", shape = "dot", color = list(background = "teal", border="black"));
visLegend(net, main="Legend", position="right", ncol=1);


#### Visualize the top differential edges ####

nDiffs= 200 # top differential edges to plot 
diffNet = pandaLCL@regNet - pandaWB@regNet

# Edges data frame
edges = matrix(0L, nDiffs, 3)
colnames(edges) = c("from","to","value")
edges = as.data.frame(edges)
aa = order(as.matrix(abs(diffNet)), decreasing = TRUE)
bb  = sort(as.matrix(abs(diffNet)), decreasing = TRUE)
edges$value = as.matrix(diffNet)[aa[1:nDiffs]]
geneIdsTop = (aa[1:nDiffs] %/% dim(diffNet)[1]) + 1
tfIdsTop = aa[1:nDiffs] %% dim(diffNet)[1]
tfIdsTop[tfIdsTop == 0] = nTFs
edges$to = colnames(diffNet)[geneIdsTop]
edges$from = rownames(diffNet)[tfIdsTop]                                  
edges$arrows = "to"   
edges$color = ifelse(edges$value > 0, "green", "red")
edges$value = abs(edges$value)

# Nodes data frame
nodes = data.frame(id = unique(as.vector(as.matrix(edges[,c(1,2)]))), label=unique(as.vector(as.matrix(edges[,c(1,2)]))))
nodes$group = ifelse(nodes$id %in% edges$from, "TF", "gene")

# Plot network
net <- visNetwork(nodes, edges, width = "100%")
net <- visGroups(net, groupname = "TF", shape = "triangle", color = list(background = "purple", border="black"))
net <- visGroups(net, groupname = "gene", shape = "dot", color = list(background = "teal", border="black"))
visLegend(net, main="Legend", position="right", ncol=1) 


#### Gene set enrichment analysis on differentially-targeted genes ####

# demonstration of how to calculate TF and gene degree in each network
lcl_outdegree <- calcDegree(pandaLCL, type="tf");
wb_outdegree <- calcDegree(pandaWB, type="tf");
lcl_indegree <- calcDegree(pandaLCL, type="gene");
wb_indegree <- calcDegree(pandaWB, type="gene");
# demonstration of how to calculate the difference in degree between two networks
degreeDiff <- calcDegreeDifference(pandaLCL, pandaWB, type="gene")

# read in pathway data
pathways <- gmtPathways("./data/c2.cp.kegg.v7.0.symbols.gmt");
# read in degree difference across all genes, as calculated in Lopes-Ramos 2017 (and format for fgsea analysis)
degreeDiff_all <- read.delim("./data/lclWB_indegreeDifference.rnk",stringsAsFactors = F,header=F);
degreeDiff_all <- setNames(degreeDiff_all[,2], degreeDiff_all[,1]);

# run fgsea
fgseaRes <- fgsea(pathways, degreeDiff_all, minSize=15, maxSize=500, nperm=1000);

# visualize the results as a bubble plot
dat <- data.frame(fgseaRes)

# bubble plot settings
fdrcut <- 0.05 # FDR cut-off to use as output for significant signatures
dencol_neg <- "blue" # bubble plot color for negative ES
dencol_pos <- "red" # bubble plot color for positive ES
signnamelength <- 4 # set to remove prefix from signature names (2 for "GO", 4 for "KEGG", 8 for "REACTOME")
asp <- 3 # aspect ratio of bubble plot
charcut <- 100 # cut signature name in heatmap to this nr of characters

# Make signature names more readable
a <- as.character(dat$pathway) # 'a' is a great variable name to substitute row names with something more readable
for (j in 1:length(a)){
  a[j] <- substr(a[j], signnamelength+2, nchar(a[j]))
}
a <- tolower(a) # convert to lower case (you may want to comment this out, it really depends on what signatures you are looking at, c6 signatures contain gene names, and converting those to lower case may be confusing)
for (j in 1:length(a)){
  if(nchar(a[j])>charcut) { a[j] <- paste(substr(a[j], 1, charcut), "...", sep=" ")}
} # cut signature names that have more characters than charcut, and add "..."
a <- gsub("_", " ", a)
dat$NAME <- a

# Determine what signatures to plot (based on FDR cut)
dat2 <- dat[dat[,"padj"]<fdrcut,]
dat2 <- dat2[order(dat2[,"padj"]),] 
dat2$signature <- factor(dat2$NAME, rev(as.character(dat2$NAME)))
# Determine what labels to color
sign_neg <- which(dat2[,"NES"]<0)
sign_pos <- which(dat2[,"NES"]>0)
# Color labels
signcol <- rep(NA, length(dat2$signature))
signcol[sign_neg] <- dencol_neg # text color of negative signatures
signcol[sign_pos] <- dencol_pos # text color of positive signatures
signcol <- rev(signcol) # need to revert vector of colors, because ggplot starts plotting these from below

# Make bubble plot
g<-ggplot(dat2, aes(x=padj,y=signature,size=size))
g+geom_point(aes(fill=NES), shape=21, colour="white")+
  theme_bw()+ # white background, needs to be placed before the "signcol" line
  xlim(0,fdrcut)+
  scale_size_area(max_size=10,guide="none")+
  scale_fill_gradient2(low=dencol_neg, high=dencol_pos)+
  theme(axis.text.y = element_text(colour=signcol))+
  theme(aspect.ratio=asp, axis.title.y=element_blank()) # test aspect.ratio
