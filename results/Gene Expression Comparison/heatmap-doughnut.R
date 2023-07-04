library(reshape)
library(ggplot2)
library("scales")
library("ggsci")

setwd("") #CHANGE WORKING DIRECTORY TO GET THE CODE WORKING

#Loading matrices
classes = read.csv('./GENERATED_MATRICES/THCA/classes.csv',sep = ';', header = TRUE)
gene_names = read.csv('./GENERATED_MATRICES/THCA/gene_names.csv',sep = ';', header = TRUE)
synth_mat = read.csv('./GENERATED_MATRICES/THCA/synth_exp_matrix.csv',sep = ';', header = TRUE)
real_mat = read.csv('./GENERATED_MATRICES/THCA/real_exp_matrix.csv',sep = ';', header = TRUE)

#Making classes as factor
classes = as.factor(classes$class)

#Defining matrix design 
disenio=model.matrix(~0+classes)
colnames(disenio)= levels(classes)
head(disenio)
apply(disenio,�2,�sum)

#Selecting randomly a column from synth matrix
synth_mat_col_name <- sample(colnames(synth_mat[,disenio[,'Primary Tumor'] == 1]),1)
synth_mat<-as.data.frame(synth_mat[,synth_mat_col_name])
synth_mat['gene'] <- gene_names$gene

#Selecting randomly a column from real matrix
real_mat_col_name <- sample(colnames(real_mat[,disenio[,'Primary Tumor'] == 1]),1)
real_mat<-as.data.frame(real_mat[,real_mat_col_name])
real_mat['gene'] <- gene_names$gene
colnames(real_mat) <- c('values', 'gene')

#Add numerical label to data

synth_mat['label'] <- data.frame(matrix(1,dim(gene_names)[1],ncol = 1))
colnames(synth_mat) <- c('values', 'gene', 'label')

real_mat['label'] <- data.frame(matrix(2,dim(gene_names)[1],ncol = 1))
colnames(real_mat) <- c('values', 'gene', 'label')

#Joining both datasets

merged_matrix = rbind(synth_mat,real_mat)
merged_matrix$var2 = as.numeric(merged_matrix$label) + 1.5

#### Computing angles of the gene names
y_labels = merged_matrix$label
y_breaks = seq_along(y_labels) + 1.5

merged_matrix.labs <- subset(merged_matrix, label==1)
merged_matrix.labs <- merged_matrix.labs[order(merged_matrix.labs$gene),]
merged_matrix.labs$ang <- seq(from=(360/nrow(merged_matrix.labs))/1.5, to=(1.5*(360/nrow(merged_matrix.labs)))-360, length.out=nrow(merged_matrix.labs))+80
merged_matrix.labs$hjust <- 0
merged_matrix.labs$hjust[which(merged_matrix.labs$ang < -90)] <- 1
merged_matrix.labs$ang[which(merged_matrix.labs$ang < -90)] <- (180+merged_matrix.labs$ang)[which(merged_matrix.labs$ang < -90)]

##### PLOTING

p2 = ggplot(merged_matrix, aes(x=gene, y=var2, fill=values)) +
  geom_tile(colour="white") +
  geom_text(data=merged_matrix.labs, aes(x=gene, y=var2+1.5, #1.5 permite que se aleje del centro
                               label=gene, angle=ang, hjust=hjust), size=3) +
  scale_fill_gradient(low = "white", high = "#B71B1BFF") +
  ylim(c(0, max(merged_matrix$var2)+1.5)) +
  scale_y_discrete(breaks=y_breaks, labels=y_labels, limits=factor(1,2)) +
  coord_polar(theta="x") +
  ggtitle("BREAST CANCER") + 
  theme(panel.background=element_blank(),
        axis.title=element_blank(),
        panel.grid=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks=element_blank(),
        axis.text.y=element_text(size=5))+
  annotate('label', x = 0.5, y = c(3.5,2.5), label = c('FAKE','REAL'), fill="gray", alpha=0.8)

print(p2)
ggsave(filename="./breast-hm.eps", plot=p2, height=7, width=7)