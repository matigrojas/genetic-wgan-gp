library("edgeR")
library('vctrs')
library(scales)
library(data.table)

setwd("") #CHANGE WORKING DIRECTORY TO GET THE CODE WORKING

alpha = 1e-2

classes = read.csv('./BATCH/PRAD/classes.csv',sep = ';', header = TRUE)
gene_names = read.csv('./BATCH/PRAD/gene_names.csv',sep = ';', header = TRUE)

classes = as.factor(classes$class)

disenio=model.matrix(~0+classes)
colnames(disenio)= levels(classes)
head(disenio)
apply(disenio, 2, sum)

batches = 30

freq <- data.frame(matrix(0,dim(gene_names)[1],ncol = 4))
rownames(freq)<-gene_names$gene
colnames(freq)<-c('in_real','in_synth','both', 'none')

##### ANALYSIS

for (i in 0:(batches-1)){
  #Reading real matrix
  file_path = './BATCH/PRAD/real_matrix'
  file_path  = paste(file_path,'_',i,'.csv',sep="", collapse=NULL)
  real_mat = read.csv(file_path,sep = ';', header = TRUE)
  rownames(real_mat)=gene_names$gene

  #Reading synthetic matrix
  file_path = './BATCH/PRAD/synth_matrix'
  file_path  = paste(file_path,'_',i,'.csv',sep="", collapse=NULL)
  synth_mat = read.csv(file_path,sep = ';', header = TRUE)
  rownames(synth_mat)=gene_names$gene

  #Starting Diff Exp analysis for real data
  real_dgelist <- DGEList(counts=real_mat, group=classes, genes = gene_names)
  real_dgelist <- calcNormFactors(real_dgelist)
  real_dgelist <- estimateGLMCommonDisp(real_dgelist, disenio)
  real_dgelist <- estimateGLMTagwiseDisp(real_dgelist, disenio)

  de_real <- exactTest(real_dgelist)

  tablaEDif_real = topTags(de_real,n='inf')[["table"]]
  tablaORDth_real <- tablaEDif_real[order(tablaEDif_real$FDR),]

  SIthREAL = rownames(tablaORDth_real[which(tablaORDth_real$FDR<=alpha),])
  NOthREAL = rownames(tablaORDth_real[which(tablaORDth_real$FDR>alpha),])

  #Starting Diff Exp analysis for synthetic data
  synth_dgelist <- DGEList(counts=synth_mat, group=classes, genes = gene_names)
  synth_dgelist <- calcNormFactors(synth_dgelist)
  synth_dgelist <- estimateGLMCommonDisp(synth_dgelist, disenio)
  synth_dgelist <- estimateGLMTagwiseDisp(synth_dgelist, disenio)

  de_synth <- exactTest(synth_dgelist)

  tablaEDif_synth = topTags(de_synth,n='inf')[["table"]]
  tablaORDth_synth <- tablaEDif_synth[order(tablaEDif_synth$FDR),]

  SIthSYNTH = rownames(tablaORDth_synth[which(tablaORDth_synth$FDR<=alpha),])
  NOthSYNTH = rownames(tablaORDth_synth[which(tablaORDth_synth$FDR>alpha),])

  for (gene in setdiff(SIthREAL,SIthSYNTH)){
    freq[gene,'in_real']<-freq[gene,'in_real']+1
  }

  for (gene in setdiff(SIthSYNTH,SIthREAL)){
    freq[gene,'in_synth']<-freq[gene,'in_synth']+1
  }

  for (gene in intersect(SIthREAL,SIthSYNTH)){
    freq[gene,'both']<-freq[gene,'both']+1
  }

  for (gene in intersect(NOthREAL,NOthSYNTH)){
    freq[gene,'none']<-freq[gene,'none']+1
  }
}

