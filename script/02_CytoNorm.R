
library("flowCore")
library("FlowSOM")
library("flowViz")
library("flowVS")
library("flowAI")
library("PeacoQC")
library("CytoML")
library("flowWorkspace")
library("CytoNorm")
library("dplyr")
library("Seurat")
library("uwot")
library("ggplot2")
library("tidyr")
library("cowplot")
library("tidyverse")
library("plotly")
library("ggpubr")
library("knitr")

rm(list=ls())


###Transformation ----

fcs.dir<- file.path("") #file path of .fcs files that should be transformed
fcs_data <- read.flowSet(path=fcs.dir, pattern="*.fcs", transformation = FALSE, truncate_max_range = FALSE) 
panel<-read.csv("")#file path of .csv with colors and markers
markerstotransform <- panel$fcs_colname[c(7:35)] ##select the markers to transform

##calculate Cofactors for Arcsinh transformation (after downsampling, takes some time to run)

fcs_data_small <- Downsampling_FlowSet(x=fcs_data, samplesize = 2000) #samplesize is the number of cells included, you can include more cells. (will need more computation; this is just for the transformation)
cofactors <- estParamFlowVS(fcs_data_small, channels=markerstotransform)
cofactordata <- data.frame(fcs_colname=markerstotransform, cofactors)%>%left_join(panel, by="fcs_colname")

#manual adjustments: 
cofactordata$cofactors<-ifelse(cofactordata$antigen=="CD3"&!is.na(cofactordata$antigen),10000,
                               ifelse(cofactordata$antigen=="CD122"&!is.na(cofactordata$antigen),10000,
                                      ifelse(cofactordata$antigen=="CD8"&!is.na(cofactordata$antigen),10000,
                                             ifelse(cofactordata$antigen=="CD19"&!is.na(cofactordata$antigen),10000,
                                                    ifelse(cofactordata$antigen=="CD126"&!is.na(cofactordata$antigen),3000,
                                                           ifelse(cofactordata$fcs_colname=="FJComp-Zombie UV-A",40000,
                                                                  cofactordata$cofactors))))))

#transformation
fcs_transform <- transFlowVS(fcs_data, channels = markerstotransform, cofactordata$cofactors) 
filenames <- sampleNames(fcs_data)
sampleNames(fcs_transform) <- filenames
markernames(fcs_transform)<-unlist((markernames(fcs_data)[2]))

#output
setwd(file.path("")) #set working directory to save transformed .fcs files and cofactordata
write.csv(x=cofactordata, file="cofactordata.csv")

pdf(file="density_plot_before_transformation.pdf", width=20, height=20)
densityplot(~., fcs_data[[1]])
dev.off()

pdf(file="density_plot_after_transformation.pdf", width=20, height=20)
densityplot(~., fcs_transform[[1]])
dev.off()

outdir <- file.path(getwd()) 
filenames <- paste(gsub(".fcs", "",fcs_data@phenoData@data$name), "_transformed", sep="")
write.flowSet(fcs_transform, outdir = outdir, filename = filenames) 

# CytoNorm Clustering with selected markers for FlowSOM----

fcs.dir<- file.path("") #file path of transformed .fcs files
files <- list.files(fcs.dir, pattern = ".fcs")
train_files <- files[grep("control",files)] # select training files 
validation_files <- files[grep("Acute|6months|12months|vacc|Healthy", files)]
files[which(!files%in% c(train_files, validation_files))]

data <- data.frame(File = files,
                   Path = file.path(fcs.dir, files),
                   Batch = strtrim(files, 5),
                   stringsAsFactors = FALSE)

data$Type<-ifelse(grepl("control", data$File)&!grepl("new_batch", data$File), "control",
                  ifelse(grepl("Acute|6months|12months|Healthy", data$File),"Patient",
                         ifelse(grepl("vacc", data$File),"Vaccination",
                                NA)))
train_data<-data%>%filter(Type=="control")
markerstocluster <- panel%>% filter(antigen %in% c("CD3", "CD4", "CD8", "CD56", "CD45RA", "CCR7"))%>%select(fcs_colname)%>%unlist(use.names=F) #define markers used for clustering

#check 5x5 SOM 
fsom <- prepareFlowSOM(train_data$Path, colsToUse = markerstocluster, transformList = NULL, FlowSOM.params = list(xdim=5,ydim=5, nClus=22, scale=FALSE)) 
cvs <- testCV(fsom,cluster_values = 3:22) 
dev.off()


##Calculate model----
dir <- file.path("") #define input directory
setwd(dir)
files <- list.files(dir, pattern = "fcs")
paths<-file.path(dir, files)

data <- data.frame(File = files,
                   Path = file.path(dir, files),
                   Batch = strtrim(files, 5),
                   stringsAsFactors = FALSE)


data$Type<-ifelse(grepl("control", data$File)&!grepl("new_batch", data$File), "control",
                  ifelse(grepl("Acute|6months|12month|Healthy", data$File),"Patient",
                         ifelse(grepl("vacc", data$File),"Vaccination",
                                NA)))

train_data <- data%>%filter(Type=="control")
validation_data <- data%>%filter(Type %in%c("Patient","Vaccination"))

ff <- flowCore::read.FCS(data$Path[1], truncate_max_range=FALSE)
markers<-flowCore::colnames(ff)
channels <- grep("FJComp-", flowCore::colnames(ff), value = TRUE) 
panel<-data.frame(marker=flowCore::markernames(ff)%>%unlist(use.names=F), channel=grep("FJComp-", flowCore::colnames(ff), value = TRUE))
markers<-flowCore::markernames(ff)
channels_selected<-channels
channels_clustering<-panel%>%filter(marker %in% c("CD3", "CD4", "CD8", "CD56", "CD45RA", "CCR7") )%>%select(channel)%>%unlist(use.names=F)
transformList <- flowCore::transformList(channels_clustering,cytofTransform) 
transformList.reverse <- flowCore::transformList(channels_selected, cytofTransform.reverse) 

model <- CytoNorm.train(train_data$Path ,
                        labels = train_data$Batch,
                        channels = channels_selected,
                        transformList = transformList,
                        FlowSOM.params = list(nCells = 1000000,
                                              xdim = 5,
                                              ydim = 5,
                                              nClus = 7,
                                              scale = FALSE,
                                              colsToUse=channels_clustering),
                        normParams = list(nQ = 101),
                        seed = 123,
                        plot=TRUE,
                        clean=FALSE)

#Apply model---- 

outputDir_manual<-file.path("") #set output directory for normalized .fcs files

#apply model on samples:

CytoNorm.normalize(model = model,
                   files = validation_data$Path,
                   labels = validation_data$Batch,
                   transformList = transformList,
                   transformList.reverse = transformList.reverse,
                   verbose = TRUE,
                   outputDir = outputDir_manual,
                   prefix = "Norm_")

# apply the model on control samples: 

CytoNorm.normalize(model = model,
                   files = train_data$Path,
                   labels = train_data$Batch,
                   transformList = transformList,
                   transformList.reverse = transformList.reverse,
                   verbose = TRUE,
                   outputDir = outputDir_manual,
                   prefix = "Norm_")