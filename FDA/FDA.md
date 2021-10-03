# Flexible Discriminant Analysis (LDA) 

Here we explain how FDA can be performed to classify cells based on expression of top DE genes.

``` R 
#### Flexible Discriminant Analysis (FDA)
memory.limit(size=150000)
library(dplyr)
library(ggplot2)
library(MASS)
library(mda)
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
             "GCH1","CXCL11","SFTPB"
             ,"FABP4","SFTA1P","SGK1","GJA1","HKDC1","CRYBA4",
             "GBP5","MARCKS","CD3E","ANKRD22","CTSE","SFN","CCL3","S100A14","PFKFB3",
             "SDS","LGMN","KRT7","CYR61","IFIT1","CTSB","CTSL","SAA2.SAA4","ATF5","IFNL3",
             "GBP1","GBP4","TRAV38.1","CCL2", "TSHZ2","CCL4","NINJ1","CD38","MAL2","ALOX5AP",
             "CD7","SPRR1B","MAFB","TRBV7.4","ATF3","RNASE7","CXCL9","GZMA","FBXW10","KRT16",
             "IGFBP4","FGL2","KRTAP3.1","CLDN3","LAMC2","PLA2G7","MARCO","IER3","NUPR1","CD2",
             "PTPRCAP","SPHK1","ICAM1","TMPRSS11D","CD247","IGHV5.51","SPOCK2","CEACAM7","CD48",
             "SFTPA2","CCR4","SPP1","LYZ")

imp_genes<-c("RSAD2","CXCL10","IDO1","GCH1","CXCL11","CRYBA4","CCL3","LGMN","IFIT1","CTSB",
"GBP1","CCL2")

imp_genes<-c( "SFTA2",      "RSAD2",     "IL21",       "PEG10",      "IGHV3.69.1",
              "SPON1",      "IFNL2"  ,    "SLAMF7",     "ADGRF5",     "CCL8",
              "CALHM6",     "RGL1"    ,   "CD52"   ,    "CXCL10" ,    "IL4I1",
              "IGFL2" ,     "IDO1"    ,   "GCH1"   ,    "CXCL11" ,    "SFTPB",
              "FABP4" ,     "SFTA1P"  ,   "SGK1"   ,    "GJA1"   ,    "HKDC1",
              "CRYBA4",     "GBP5"    ,   "MARCKS" ,    "CD3E"   ,    "TNFSF10",
              "ANKRD22",    "CTSE"    ,   "SFN"    ,    "CCL3"   ,    "S100A14",
              "PFKFB3" ,    "SDS"     ,   "LGMN"   ,    "KRT7"   ,    "CYR61",
              "IFIT1"  ,    "CTSB"    ,   "CTSL"   ,    "SAA2.SAA4",  "ATF5" ,
              "IFNL3"  ,    "GBP1"    ,   "GBP4"   ,    "TRAV38.1",   "CCL2"  ,
              "TSHZ2"  ,    "CCL4"    ,   "NINJ1"  ,    "CD38"   ,    "MAL2"   ,
              "ALOX5AP" ,   "CD7"     ,   "SPRR1B" ,    "MAFB"   ,    "TRBV7.4" ,
              "MX1"     ,   "IFIT3"   ,   "ATF3"   ,    "RNASE7" ,    "CXCL9"   ,
              "AC089983.1", "GZMA"    ,   "FBXW10" ,    "KRT16"  ,    "VAMP5"   ,
              "IGFBP4"  ,   "HES4"    ,   "FGL2"   ,    "KRTAP3.1",   "CLDN3"   ,
              "LAMC2"   ,   "PLA2G7"  ,   "IL1RN"  ,    "MARCO"   ,   "IER3"    ,
              "NUPR1"   ,   "CD2"     ,   "PTPRCAP",    "ISG20"   ,   "SPHK1"   ,
              "TYROBP"  ,   "ICAM1"   ,   "TMPRSS11D",  "CD247"   ,   "IGHV5.51" ,
              "KLF6"    ,   "TNFSF13B",   "SPOCK2" ,    "CEACAM7" ,   "CCL4L2" ,
              "CD48"    ,   "SFTPA2" ,    "CCR4"   ,    "SPP1"    ,   "LYZ"    ,
              "FNIP2"    ,  "TGFBI"   ,   "RGS1"    ,   "CST3"     ,  "SCIN"    ,
              "CEACAM6"  ,  "GPR183"  ,   "PALM"    ,   "TM4SF1"   ,  "CAV1"    ,
              "SSTR2"    ,  "RAB25"   ,   "S100A12" ,   "APOBEC3A" ,  "IFITM1"  ,
              "FFAR2"    ,  "C1QA"    ,   "IFITM2"  ,   "C1QB"     ,  "ACOD1"   ,
              "S100A8"   ,  "S100A9"  ,   "CLU"     ,   "C1QC"     ,  "MT.CO1"  ,
              "FCN1"     ,  "MT1M"    ,   "IFITM3"  ,   "MT.CO2"   ,  "GPNMB"   ,
              "FPR1"     ,  "CD74"    ,   "FTL"     ,   "DEFB1"    ,  "CTSD"    ,
              "LILRA5"   ,  "MT.ND3"  ,   "APOC1"   ,   "APOE"     ,  "TIMD4"   ,
              "IFIT2"    ,  "MT.ND4L" ,   "SLC25A37",   "MT.CO3"   ,  "HPSE"    ,
              "FEZ1"     ,  "HLA.DQA2",   "CCL3L1" ,    "MT1X"     ,  "IGLV3.9" ,
              "CLEC4E"   ,  "SELL"    ,   "FCGR3B" ,    "RPL21"    ,  "ISG15"   ,
              "NCF1"     ,  "RPS27"   ,   "MT.ATP8",    "IL1R2"    ,  "CCRL2"   ,
              "MT2A"     ,  "MSR1"    ,   "TRAV34",     "HLA.DPB1" ,  "CLEC4D"  ,
              "CCR1"     ,  "NR1H3"   ,   "CD177",      "PLD3"     ,  "MT1G"    ,
              "BCL2A1"   ,  "C15orf48",   "PLEK",       "RNASE1"   ,  "HLA.DRA" ,
              "OLIG1"    ,  "SOCS1"  ,    "HLA.DRB1",   "IFNA17"   ,  "RNASE2"  ,
              "GSN"       , "TAOK1"   ,   "IFNA7",      "G0S2"      , "IFNA2"    ,
              "PPP1R15A" ,  "ANXA3"  ,    "FOLR3",      "NAMPT"    ,  "MT1E"    ,
              "TNFAIP6"  ,  "MT.ATP6",    "SAT1",       "AC136475.9")



data_imp<-t(data)[,intersect(rownames(data),imp_genes)]
#data_imp_nor<-cbind(true_lables$Nor_vs_rest,data_imp)
data_imp_nor<-cbind(as.numeric(as.factor(true_lables$Sev_vs_rest)),data_imp)
#data_imp<-as.data.frame(data_imp)

smp_size_raw<-floor(0.75*nrow(data_imp_nor))
set.seed(1000)
train_ind_raw<-sample(nrow(data_imp_nor),size = smp_size_raw)
train_raw.df<-as.data.frame(data_imp_nor[train_ind_raw,])
test_raw.df<-as.data.frame(data_imp_nor[-train_ind_raw,])

print("pre-processing done!")
gc()

model_nor <- fda(V1~., data = train_raw.df, method=mars)
print("Model fitting done!")
test.lda.predict_nor<-predict(model_nor,test_raw.df,type='posterior')
print("Prediction done!")

pdf(file = "/scratch/user/naminiyakan/hackathon/88genes_fda_sev_1000.pdf")
rownames(true_lables)<-colnames(data)
lda.data<-cbind(true_lables[intersect(rownames(train_raw.df),rownames(true_lables)),],predict(model_nor)$x)
#lda.data <- cbind(train_raw.df, predict(model_nor)$x)
ggplot(lda.data, aes(LD1, 0)) +
  geom_point(aes(color = Mil_vs_rest))
dev.off()

mean(test.lda.predict_nor$class==test_raw.df$V1)

#save.image(file="/scratch/user/naminiyakan/hackathon/88genes_fda_sev_1000.RData")

#### AUCROC
#PRROC_obj <- roc.curve(scores.class0 = (test.lda.predict_nor[,1]), weights.class0=as.numeric(as.factor(test_raw.df$V1)),
#                       curve=TRUE)
#print(PRROC_obj$auc)

library(pROC)
auc(as.numeric(as.factor(test_raw.df$V1)),(test.lda.predict_nor[,1]))

'''
