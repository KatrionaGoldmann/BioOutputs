


#install.packages("devtools")
library("devtools")
#devtools::install_github("klutometis/roxygen")
library(roxygen2)

#Initialize the package
setwd("~/Documents/Scripts/Git/")
#create_package("BioOutputs")


# Load/update the package
package.name="BioOutputs"
setwd(paste("~/Documents/Scripts/Git/", package.name, sep=""))
document()

# Examples
?bio_frequency

metadata <- readRDS("~/V:/PEAC/PEAC_Imputed_data.rds")
data <- metadata
data <- data[, colnames(data) !="CD21"]
data$path2 <- data$Pathotype
bio_frequency(data, "Pathotype")
bio_frequency(data, columns=c("Pathotype", "path2"), freq.percent = "percent")

which(colnames(data) == "Pathotype")
bio_frequency(data, c(6, 73), include.na=FALSE)
bio_frequency(data, "Pathotype", remove.vars = c("Ungraded"))
bio_frequency(data, "Pathotype", remove.vars = c("Ungraded"), include.na=FALSE)
bio_frequency(data, columns=colnames(data)[grepl("CD", colnames(data))], include.na=TRUE, freq.percent = "both")



?bio_volcano


load("~/V:/Katriona/ShinyDocker/CurrentWebsite/PEAC HPC 1.85/Complete_web_data.rdata")
save(LVsF, file="/home/katrionagoldmann/Documents/Scripts/Git/BioOutputs/data/grav.rda")
toptable=LVsF
fc.col="log2FoldChange"
fc.cutoff=1
padj.cutoff=0.05
main="hello"
label.p.cutoff=NULL
add.lines=TRUE
label.row.indices=order(toptable$padj)[1:50] # The 50 most significant genens

bio_volcano(toptable, fc.col="log2FoldChange", padj.col=NULL, padj.method="fdr",
											 padj.cutoff=0.01, label.p.cutoff=NULL, label.row.indices=order(toptable$padj)[1:10],
											 fc.cutoff=1, add.lines=TRUE, line.colour="grey14", label.colour="black",
											 main="hello world", xlims=NULL, ylims=NULL, marker.colour=c("grey60", "violet", "salmon", "darkturquoise"), 
											 legend.labs=NULL)


?bio_corr
data <- data.frame("exp"=as.numeric(as.character(syn.rldMatrix[1, ])), "clin"=as.numeric(as.character(syn.metadata.final$age_inclusion)))
x="clin"
y="exp"



library("devtools")
#devtools::install_github("klutometis/roxygen")
library(roxygen2)

#Initialize the package
setwd("~/Documents/Scripts/Git/")
#create_package("BioOutputs")


# Load/update the package
package.name="BioOutputs"
setwd(paste("~/Documents/Scripts/Git/", package.name, sep=""))
document()

bio_corr(data, x, y)
