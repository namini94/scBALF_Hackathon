# Linear Discriminant Analysis (LDA) 

Here we explain how LDA can be performed to classify cells based on expression of top DE genes.

``` R 
memory.limit(size=150000)
library(dplyr)
library(ggplot2)
library(MASS)
#install.packages("PRROC")
#library(PRROC)
data<-read.csv(file = "/scratch/user/naminiyakan/hackathon/covid-selected-data.csv")
data<-as.matrix(data)
rownames(data)<-data[,1]
data<-data[,-1]
n<-rownames(data)
data<-apply(data,2,as.numeric)
rownames(data)<-n
data<-t(data)
dim(data)

true_lables<-read.csv(file="/scratch/user/naminiyakan/hackathon/covid-selected-data-labels.csv")
true_lables<-as.matrix(true_lables)
rownames(true_lables)<-true_lables[,1]
true_lables<-true_lables[,-1]
n<-rownames(data)
rownames(data)<-n
true_lables<-as.matrix(true_lables)
colnames(true_lables)<-c("cell_type")
true_lables<-as.data.frame(true_lables)
summary(as.factor(true_lables$cell_type))

true_lables<-true_lables %>% mutate(Nor_vs_rest = case_when(
  cell_type == c("Normal") ~ "Normal",
  cell_type != c("Normal") ~ "Rest"))

true_lables<-true_lables %>% mutate(Sev_vs_rest = case_when(
  cell_type == c("Severe") ~ "Severe",
  cell_type != c("Severe") ~ "Rest"))

true_lables<-true_lables %>% mutate(Mil_vs_rest = case_when(
  cell_type == c("Mild") ~ "Mild",
  cell_type != c("Mild") ~ "Rest"))

imp_genes<-c("SFTA2","RSAD2","IL21","PEG10","IGHV3.69.1","SPON1","IFNL2","SLAMF7",
             "ADGRF5","CCL8","CALHM6","RGL1","CD52","CXCL10","IL4I1","IGFL2","IDO1", 
             "GCH1","CXCL11","SFTPB","FABP4","SFTA1P","SGK1","GJA1","HKDC1","CRYBA4",
             "GBP5","MARCKS","CD3E","ANKRD22","CTSE","SFN","CCL3","S100A14","PFKFB3",
             "SDS","LGMN","KRT7","CYR61","IFIT1","CTSB","CTSL","SAA2.SAA4","ATF5","IFNL3",
             "GBP1","GBP4","TRAV38.1","CCL2", "TSHZ2","CCL4","NINJ1","CD38","MAL2","ALOX5AP",
             "CD7","SPRR1B","MAFB","TRBV7.4","ATF3","RNASE7","CXCL9","GZMA","FBXW10","KRT16",
             "IGFBP4","FGL2","KRTAP3.1","CLDN3","LAMC2","PLA2G7","MARCO","IER3","NUPR1","CD2",
             "PTPRCAP","SPHK1","ICAM1","TMPRSS11D","CD247","IGHV5.51","SPOCK2","CEACAM7","CD48",
             "SFTPA2","CCR4","SPP1","LYZ")

data_imp<-t(data)[,intersect(rownames(data),imp_genes)]
#data_imp_nor<-cbind(true_lables$Nor_vs_rest,data_imp)
data_imp_nor<-cbind(as.numeric(as.factor(true_lables$Mil_vs_rest)),data_imp)
#data_imp<-as.data.frame(data_imp)

smp_size_raw<-floor(0.75*nrow(data_imp_nor))
set.seed(1000)
train_ind_raw<-sample(nrow(data_imp_nor),size = smp_size_raw)
train_raw.df<-as.data.frame(data_imp_nor[train_ind_raw,])
test_raw.df<-as.data.frame(data_imp_nor[-train_ind_raw,])

print("pre-processing done!")
gc()

model_nor <- lda(V1~., data = train_raw.df)
print("Model fitting done!")
test.lda.predict_nor<-predict(model_nor,test_raw.df)
print("Prediction done!")


pdf(file = "/scratch/user/naminiyakan/hackathon/88genes_lda_mil_1000.pdf")
rownames(true_lables)<-colnames(data)
lda.data<-cbind(true_lables[intersect(rownames(train_raw.df),rownames(true_lables)),],predict(model_nor)$x)
#lda.data <- cbind(train_raw.df, predict(model_nor)$x)
ggplot(lda.data, aes(LD1, 0)) +
  geom_point(aes(color = Mil_vs_rest))
dev.off()

mean(test.lda.predict_nor$class==test_raw.df$V1)

save.image(file="/scratch/user/naminiyakan/hackathon/88genes_lda_mil_1000.RData")
        
```
