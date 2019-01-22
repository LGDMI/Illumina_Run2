#### USER SETTINGS ####
silvadb <- 128 # set the version of the SILVA release that was used for classification
rdprel <- 16 #set the version of the RDP release that was used
fldr <- "/media/projects2/LisaM/Batch_Experiments/OTU-recovery2_AllBatch/2_Processing/" #relative location from the report project to the data.
dataname <- "Batch"
metadatafn <- "MetaData.xlsx" #relative path (from fldr) to select other Metadata files
#during the SOP workshop we'll tell people to use the pre-trimmed file in the 
#taxonomy folders on the server to have a more reproducible way of calling this.

#### load required packages ####
library(scales)   #for percent formatting
library(ggplot2)  #for adequate plotting
library(plyr)     #data wrangling (mapvalues)
suppressPackageStartupMessages(library(dplyr))    #data wrangling
library(tidyr)    #tidy data
suppressPackageStartupMessages(library(reshape2)) #for melt function
library(vegan)    #for ecological calculations
library(phyloseq) #for microbiome census data processing
library(ade4)     #for ecological calculations
library(gplots)   #for heatmap.2
library(splitstackshape) #for csplit
library(knitr)    #for simple markdown tables
library(xtable)   #for more advanced tables
library(ape)      #dealing with phylogenetic trees
library(xlsx)     #handling Excel files
library(openxlsx) #handling Excel files without Java deps
library(readxl)   #faster handling of excel files
library(data.table) #data wrangling
library(SPECIES)  #alpha diversity estimators
library(parallel) #parallel computation in R
library(purrr) #functional programming
library(DESeq2)  #normalization procedures to cope with differences in smp depth
library(NMF)
library(devtools) #nicer session info
library(viridis) # for color-blind nice color palettes
library(Phenoflow)
library(pander)


#### user defined functions ####
#need to be integrated in dedicated package at https://github.ugent.be/LabMETNGS/CMETNGS_package
preformattax <- function(taxonomy)     #reformats column with taxonomy (see excelfile) and splits it into different columns
{
  #Input: mothur taxonomy with probs=T (default)
  #Output: splitted and stripped tax file for R processing (e.g. phyloseq)
  tax.good <- cSplit(taxonomy,"Taxonomy",";")
  tax.good.probs <- do.call(cbind, lapply(tax.good[,3:ncol(tax.good)], 
                                          function(x){
                                            do.call(rbind, 
                                                    strsplit(as.character(x), 
                                                             "\\((?=\\d)", 
                                                             perl = TRUE))
                                          }
  )
  )
  # tax.good.probs.nona <- subset(tax.good.probs,
  #           select=which(colSums(apply(tax.good.probs,2,is.na))==0))
  tax.final <- as.data.frame(apply(tax.good.probs,2,
                                   function(x) sub(")","",x)))
  otunames <- tax.good$OTU
  rownames(tax.final) <- otunames
  colnames(tax.final) <- c("Regnum","Prob_R","Phylum","Prob_P",
                           "Classis","Prob_C","Ordo","Prob_O",
                           "Familia","Prob_F","Genus","Prob_G")
  if(ncol(tax.final)==12){
    #for RDP & SILVA taxonomies which stop at genus level
    tax.final <- as.data.frame(mutate(tax.final,Species=NA,Prob_S=NA))
    rownames(tax.final) <- otunames
  }else{
    #for greengenes taxonomy
    colnames(tax.final)[13:14] <- c("Species","Prob_S")
    tax.final <- as.data.frame(apply(tax.final,
                                     2,
                                     function(x)sub("[a-z]__","",x)))
    tax.final$Species <- paste(tax.final$Genus,tax.final$Species) #just the species is not very sensible
  }
  #to comply to mothur 1.38 tax format
  charvectgenus <- as.character(tax.final$Genus)
  probvectgenus <- as.numeric(as.character(tax.final$Prob_G))
  nasgenuslev <- which(is.na(charvectgenus))
  for(i in nasgenuslev)
  {
    if(grepl("unclassified",as.character(tax.final$Familia[i])))
    {
      if(!is.na(as.character(tax.final$Familia[i])))
      {
        charvectgenus[i]<-as.character(tax.final$Familia[i])
      }else{
        if(grepl("unclassified",as.character(tax.final$Ordo[i])))
        {
          charvectgenus[i]<-as.character(tax.final$Ordo[i])
        }
      }
      
    }else{
      if(!is.na(as.character(tax.final$Familia[i])))
      {
        charvectgenus[i]<-paste0(as.character(tax.final$Familia[i]),"_unclassified")
      }else{
        charvectgenus[i]<-paste0(as.character(tax.final$Ordo[i]),"_unclassified")
      }
      
    }
    probvectgenus[i] <- as.numeric(as.character(tax.final$Prob_F[i]))
  }
  tax.final$Genus <- factor(charvectgenus)
  tax.final$Prob_G <- probvectgenus
  
  charvectfamilia <- as.character(tax.final$Familia)
  probvectfamilia <- as.numeric(as.character(tax.final$Prob_F))
  nasfamlev <- which(is.na(charvectfamilia))
  for(i in nasfamlev)
  {
    if(grepl("unclassified",as.character(tax.final$Ordo[i])))
    {
      if(!is.na(as.character(tax.final$Ordo[i])))
      {
        charvectfamilia[i]<-as.character(tax.final$Ordo[i])
      }else{
        if(grepl("unclassified",as.character(tax.final$Classis[i])))
        {
          charvectfamilia[i]<-as.character(tax.final$Classis[i])
        }
      }
      
    }else{
      if(!is.na(as.character(tax.final$Ordo[i])))
      {
        charvectfamilia[i]<-paste0(as.character(tax.final$Ordo[i]),"_unclassified")
      }else{
        charvectfamilia[i]<-paste0(as.character(tax.final$Classis[i]),"_unclassified")
      }
      
    }
    probvectfamilia[i] <- as.numeric(as.character(tax.final$Prob_O[i]))
  }
  tax.final$Familia <- factor(charvectfamilia)
  tax.final$Prob_F <- probvectfamilia
  return(tax.final)
}

preformattax.new <- function(taxonomy)
{
  #Input: mothur taxonomy with probs=T (default)
  #Output: splitted and stripped tax file for R processing (e.g. phyloseq)
  tax.good <- cSplit(taxonomy,"Taxonomy",";")
  tax.good.probs <- do.call(cbind, lapply(tax.good[,3:ncol(tax.good)], 
                                          function(x){
                                            do.call(rbind, 
                                                    strsplit(as.character(x), 
                                                             "\\((?=\\d)", 
                                                             perl = TRUE))
                                          }
  )
  )
  #above lies difference with regular preformattax, checks for number after bracket (line above and below)
  tax.final <- as.data.frame(apply(tax.good.probs,2,
                                   function(x) sub(")","",x)))
  
  # tax.final.noprobs <- as.data.frame(apply(tax.good.probs,2,
  #                                 function(x) sub(")","",x)))
  # tax.final <- suppressMessages(suppressWarnings(as.data.frame(sapply(tax.final.noprobs,
  #                                     function(x){mapvalues(x,from="00",to="100")}))))
  otunames <- tax.good$OTU
  #tax.final <- subset(tax.final,select=-c(OTU,Size))
  rownames(tax.final) <- otunames
  colnames(tax.final) <- c("Regnum","Prob_R","Phylum","Prob_P",
                           "Classis","Prob_C","Ordo","Prob_O",
                           "Familia","Prob_F","Genus","Prob_G")
  if(ncol(tax.final)==12){
    #for RDP & SILVA taxonomies which stop at genus level
    tax.final <- as.data.frame(mutate(tax.final,Species=NA,Prob_S=NA))
    rownames(tax.final) <- otunames
  }else{
    #for greengenes taxonomy
    colnames(tax.final)[13:14] <- c("Species","Prob_S")
    tax.final <- as.data.frame(apply(tax.final,
                                     2,
                                     function(x)sub("[a-z]__","",x)))
    tax.final$Species <- paste(tax.final$Genus,tax.final$Species) #just the species is not very sensible
  }
  #to comply to mothur 1.38 tax format
  charvectgenus <- as.character(tax.final$Genus)
  probvectgenus <- as.numeric(as.character(tax.final$Prob_G))
  nasgenuslev <- which(is.na(charvectgenus))
  for(i in nasgenuslev)
  {
    if(grepl("unclassified",as.character(tax.final$Familia[i])))
    {
      if(!is.na(as.character(tax.final$Familia[i])))
      {
        charvectgenus[i]<-as.character(tax.final$Familia[i])
      }else{
        if(grepl("unclassified",as.character(tax.final$Ordo[i])))
        {
          charvectgenus[i]<-as.character(tax.final$Ordo[i])
        }
      }
      
    }else{
      if(!is.na(as.character(tax.final$Familia[i])))
      {
        charvectgenus[i]<-paste0(as.character(tax.final$Familia[i]),"_unclassified")
      }else{
        charvectgenus[i]<-paste0(as.character(tax.final$Ordo[i]),"_unclassified")
      }
      
    }
    probvectgenus[i] <- as.numeric(as.character(tax.final$Prob_F[i]))
  }
  tax.final$Genus <- factor(charvectgenus)
  tax.final$Prob_G <- probvectgenus
  
  charvectfamilia <- as.character(tax.final$Familia)
  probvectfamilia <- as.numeric(as.character(tax.final$Prob_F))
  nasfamlev <- which(is.na(charvectfamilia))
  for(i in nasfamlev)
  {
    if(grepl("unclassified",as.character(tax.final$Ordo[i])))
    {
      if(!is.na(as.character(tax.final$Ordo[i])))
      {
        charvectfamilia[i]<-as.character(tax.final$Ordo[i])
      }else{
        if(grepl("unclassified",as.character(tax.final$Classis[i])))
        {
          charvectfamilia[i]<-as.character(tax.final$Classis[i])
        }
      }
      
    }else{
      if(!is.na(as.character(tax.final$Ordo[i])))
      {
        charvectfamilia[i]<-paste0(as.character(tax.final$Ordo[i]),"_unclassified")
      }else{
        charvectfamilia[i]<-paste0(as.character(tax.final$Classis[i]),"_unclassified")
      }
      
    }
    probvectfamilia[i] <- as.numeric(as.character(tax.final$Prob_O[i]))
  }
  tax.final$Familia <- factor(charvectfamilia)
  tax.final$Prob_F <- probvectfamilia
  return(tax.final)
}

source("bargraphGG.R") #will be integrated into CMETNGS package

#### load data ####
if(.Platform$OS.type == "unix") {
  crfn <- system(paste("ls",fldr," | grep contigs.report"),intern=TRUE)
} else {
  filelist <- list.files(fldr)
  crfn <- grep(".*contigs.report",filelist,value=TRUE)
}
fn <- sub(".contigs.report","",crfn)
ini <- fread(paste(fldr,"/",crfn,sep=""),header=TRUE)
csum <- fread(paste(fldr,"/",fn,
                    ".trim.contigs.summary",sep=""),header=TRUE)
firsttrimsum <- fread(paste(fldr,"/",fn,
                            ".trim.contigs.good.summary",sep=""),
                      header=TRUE)
uniquesum <- fread(paste(fldr,"/",fn,
                         ".trim.contigs.good.unique.summary",sep=""),
                   header=TRUE)
postalnsum <- fread(paste(fldr,"/",fn,
                          ".trim.contigs.good.unique.good.filter.summary",sep=""),
                    header=TRUE)
preclussum <- fread(paste(fldr,"/",fn,
                          ".trim.contigs.good.unique.good.filter.unique.precluster.summary",sep=""),
                    header=TRUE)
postuchimeclasssum <- fread(paste(fldr,"/",fn,".trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.summary",sep=""),
                            header=TRUE)




# By default (if only one taxonomy, run the code below)
# otutaxonomy <- fread(paste(fldr,"/",fn,
# ".trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.an.unique_list.0.03.cons.taxonomy",
#             sep=""),header=TRUE)


# otherwise, we default to SILVA (for now)
otutaxonomy <- fread(paste(fldr,"/",fn,
                           ".trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.0.03.cons.taxonomy",
                           sep=""),header=TRUE)

taxonomy.spl <- preformattax.new(otutaxonomy)
taxonomy.np <- taxonomy.spl %>% dplyr::select(-dplyr::contains("Prob"))    #np = no probabilities

# read in reads per sample
shared <- fread(paste(fldr,"/",fn,
                      ".trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.shared",sep=""),
                header=TRUE)
shared <- as.data.frame(shared)
desgroups <- shared$Group
shared.x <- shared[,4:ncol(shared)]    # OTU-table
rownames(shared.x) <- desgroups
shared.t <- as.data.frame(t(shared.x))   # transposed: OTU's are rows, samples are columns
shared.t.ns <- shared.t[which(rowSums(shared.t)!=1),]    # ns = no singletons
sharedns <- data.frame(label=rep(0.03,ncol(shared.t.ns)),
                       Group=colnames(shared.t.ns),
                       numOtus=nrow(shared.t.ns),t(shared.t.ns))
suppressWarnings(try(write.table(x=sharedns,file=paste(fldr,"/sharedns.shared",sep=""),
                                 row.names=FALSE,quote=FALSE),silent=TRUE))


#write.xlsx2(x = shared.t, file = "Results.xlsx", sheetName = "OTU_table" )
#write.xlsx2(x = sharedns, file = "Results.xlsx",
#            sheetName = "OTU_tablenosingletons",
#            append = TRUE)

#filter taxonomy (i.e. remove singletons from taxonomy)
taxonomy.np.ns <- taxonomy.np[which(rownames(taxonomy.np) 
                                    %in% rownames(shared.t.ns)),]

metadata <- readxl::read_excel("MetaData.xlsx",sheet = "ForR")
metadata.tibble <- readxl::read_excel(paste0(fldr,metadatafn),sheet="ForR")
factdescs <- readxl::read_excel(paste0(fldr,metadatafn),sheet="FactDesc")

#TODO(fpkerckh) : generic naming & format for MD? 
metadata <- as.data.frame(metadata.tibble) #to avoid warnings/errors with rownames
if(is.numeric(metadata$SampleName)) #check for fully numeric sample names
{
  metadata$SampleName <- paste(dataname,metadata$SampleName,sep="")
}
rownames(metadata) <- metadata$SampleName
metadata.smpdat <- sample_data(metadata)
colnames(shared.t.ns) <- plyr::mapvalues(colnames(shared.t.ns),
                                         from=as.character(metadata$Code),
                                         to=as.character(metadata$SampleName)) #rename

#the part above maps sample names to codes that were used for sequencing at LGC, regardless of the order in the shared file



#### phyloseq constructors ####
otumat.ns <- as.matrix(shared.t.ns)
taxmat.ns <- as.matrix(taxonomy.np.ns)
OTU       <- otu_table(otumat.ns,taxa_are_rows = TRUE)
TAX       <- tax_table(taxmat.ns)
physeqobj <- phyloseq(OTU,TAX)
physeqobj.meta <- merge_phyloseq(physeqobj,metadata.smpdat)

#### filtering according to WNWN ####
## prevalence = (fraction of samples in which an OTU is observed minimum 1 time)
minobs=1
prevalence <- apply(as.matrix(shared.t.ns),1,function(x,minobs){sum(x>=minobs)},minobs)/ncol(shared.t.ns)
prevalencefilter <- prevalence>0.05
sharedminsingletonwnwn <- shared.t.ns[prevalencefilter,]


##Read counts should exceed 0.5 times the number of samples
sharedfilteredwnwn <- sharedminsingletonwnwn[rowSums(sharedminsingletonwnwn)>0.5*ncol(sharedminsingletonwnwn),]
# deseq normalise ==> very dependent upon design
metadata$factor1 <- factor(metadata$Factor1)
metadata$factor2 <- factor(metadata$Factor2)
deseqdata <- as.matrix(sharedfilteredwnwn +1)
# deseqdata <- DESeqDataSetFromMatrix(deseqdata,colData=metadata,design= ~ factor1 + factor2)
# deseqdata <- estimateSizeFactors(deseqdata)
# sharedfiltereddeseq <- counts(deseqdata,normalized=TRUE)
# deseqdata <- estimateDispersions(deseqdata,fitType="local")
# sharedfiltereddeseqvartransf <- varianceStabilizingTransformation(deseqdata, blind = FALSE)
# sharedfiltereddeseqvartransfmatrix <- assay(sharedfiltereddeseqvartransf)
# wnwndeseqvartransf <- getVarianceStabilizedData(deseqdata)
# jaccard_deseq <- vegdist(t(sharedfiltereddeseq),method="jaccard",binary="FALSE") # not interpretable on negative counts!

# set global chunk options: 
library(knitr)
opts_chunk$set(cache=TRUE, autodep = TRUE)
#dep_auto()

#in this way for every chunck if any previous chunck is changed, it will be recompiled
#given that only generally the first chunck should change, on first run the cache will be 
#rebuilt by default but subsequent downstream changes should be fine

hist(csum$nbases,col = "lightblue",main="",xlab="Number of nucleotides")
box()

hist(firsttrimsum$nbases,col = "lightblue",main="",xlab="Number of nucleotides")
box()

hist(postalnsum$nbases,col = "lightblue",main="",xlab="Number of nucleotides")
box()


countdfprepr <- data.frame(uniqueseqs=c(nrow(csum),nrow(firsttrimsum),
                                        nrow(uniquesum),nrow(postalnsum),
                                        nrow(preclussum),
                                        nrow(postuchimeclasssum)),
                           totalseqs=c(nrow(csum),nrow(firsttrimsum),
                                       sum(uniquesum$numSeqs),
                                       sum(postalnsum$numSeqs),
                                       sum(postalnsum$numSeqs),
                                       sum(postuchimeclasssum$numSeqs)),
                           step=c("Contigs","Initial trim",
                                  "Unique","Alignment",
                                  "Precluster","UChime"))
countdfprepr.m <- melt(countdfprepr)
p <- ggplot(aes(x=step,y=value,fill=variable),
            data=countdfprepr.m) +
  geom_bar(stat="identity",position="dodge") +
  xlab("Nucleotide") + 
  
  
  
  
  groupcounts <- data.frame(Reads = sort(colSums(shared.t.ns),
                                         decreasing=TRUE))
kable(groupcounts,
      caption = "Number of reads for all groups after removing singletons")
write.xlsx2(x = groupcounts, file = "Results.xlsx",
            sheetName = "groupcounts", append = TRUE)


groupcountdf <- data.frame(groupcounts,samplename=rownames(groupcounts))
metadatadf <- data.frame(metadata,samplename=rownames(metadata))
countsdfplot <- dplyr::full_join(groupcountdf,metadatadf,by="samplename")
countsdfplot$Factor1 <- as.factor(countsdfplot$Factor1)
countsdfplot$Factor2 <- as.factor(countsdfplot$Factor2)
gglibsize <- ggplot(countsdfplot,aes(x=Factor1,y=Reads,fill=Factor1))
gglibsize + geom_boxplot() + facet_grid(~Factor2,drop=TRUE) + 
  theme_bw() +
  theme(axis.title.x = element_blank(),axis.text.x = element_blank()) +
  scale_fill_discrete(name=factdescs$Desc[1])+
  ggtitle(paste("Library size by",factdescs$Desc[1],"and",factdescs$Desc[2]))



#TODO(fpkerckh) : make it possible to increase n 
barplot.genus <- makebargraphrawggplot2(tax=taxonomy.np.ns,shared=shared.t.ns,
                                        topn=8,taxlevel="Genus",shared.abs=TRUE,
                                        tax.prob=FALSE,samples=dataname,
                                        plot=TRUE)
  ylab("Number of nucleotides") + theme_bw() +
  scale_fill_discrete(labels=c("Unique",
                               "Total"),
                      name="") + 
  scale_x_discrete(limits=c("Contigs",
                            "Initial trim",
                            "Unique",
                            "Alignment",
                            "Precluster",
                            "UChime")) +
  theme(legend.title=element_blank())
p
write.xlsx2(x = countdfprepr.m, file = "Results.xlsx", 
            sheetName = "filteredReadcounts", append = TRUE)



#TODO(fpkerckh) : make it possible to increase n 
barplot.genus <- makebargraphrawggplot2(tax=taxonomy.np.ns,shared=shared.t.ns,
                                        topn=8,taxlevel="Genus",shared.abs=TRUE,
                                        tax.prob=FALSE,samples=dataname,
                                        plot=TRUE)