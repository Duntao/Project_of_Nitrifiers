setwd("D:/Omicsdata Sequencing/SCI papers/Datas of twelve paper/Randomforest")

library(multcomp)
library(phyloseq)
attach(cholesterol)
library(tidyverse)
library(ggplot2)
library(agricolae)
library(reshape2)
library(vegan)
library(plyr)
library(dplyr)
library(randomForest)
library(rfPermute)
library(rfUtilities) 
library(ggcor)
library(relaimpo)
library(gridExtra)
library(ggpubr)
library(Hmisc)
library(RColorBrewer)
getwd()
#devtools::install_git("https://gitee.com/dr_yingli/ggcor")
packageVersion("ggcor")
#input
beta<-read.csv("PCoA_nitrifying_guilds.csv",header=T,row.names=1)
soil<-read.csv("soil_re.csv",header=T,row.names=1)
gene<-read.csv("gene.csv",header=T)
crop<-read.csv("crop_production.csv",header=TRUE)
View(crop)
View(soil)
multi_index<-read.csv("incdx_nmds_gene_pn.csv",header=TRUE,row.names=1)
View(multi_index)
df3_crop<-multi_index
df3_crop$yield<-crop$yield
View(df3_crop)

df_env_crop<-cbind(soil,df3_crop)
View(df_env_crop)
df_env_crop<-df_env_crop[,-c(14:15)]
View(beta)
View(df_env_crop)
beta_all<-cbind(beta[,-c(6,7,8,12,13)],df_env_crop[,c(14:26,35:44)])
View(beta_all)
write.csv(beta_all,file="beta_all_final.csv")

####ramdorest for PNO VS Environmental factors######
set.seed(123)
fit1.rp <- rfPermute(PN ~ .,ntree = 1000, nrep = 9999, num.cores = 4,
                     data = beta_all[,c(1:8,18)])#num.cores=1或2
rmp1<-rp.importance(fit1.rp)
View(rmp1)
plot(rp.importance(fit1.rp, scale = TRUE))

rf.perm1<-rf.significance(fit1.rp,beta_all[c("PN")],nperm=999,ntree=1000)
summary(rf.perm1)

###ggplot for df.imp1
df.imp1 <- as.data.frame(rmp1[, 1:2])
#df.imp1 <- mutate(df.imp1, type = row.names(df.imp1))
df.imp1<-df.imp1[order(rownames(df.imp1),decreasing=F),]
colnames(df.imp1) <- c("MSE", "MSE_pval")
df.imp2 <- mutate(df.imp1, sig= ifelse (df.imp1[2] < 0.05&df.imp1[2]>=0.01, "*",
                                        ifelse(df.imp1[2]>=0.001&df.imp1[2]<0.01,"**",
                                               ifelse(df.imp1[2]<=0.001,"***",""))))

row.names(df.imp2) <- row.names(df.imp1)
df.imp2$id<-rownames(df.imp2)
df.imp2<-df.imp2[order(df.imp2$MSE,decreasing=F),]
df.imp2$id<-factor(df.imp2$id,levels=df.imp2$id)
#df.imp3<-separate(df.imp2,id,c("gene","alpha"),sep="_")
View(df.imp2)
#df.imp2$type<-c(rep("aoa",2),rep("aob",2),rep("com",2),rep("nob",2))
#View(df.imp3)
#df.imp3$alpha<-factor(df.imp2$id,level=c("Chao","Shannon","Chao","Shannon",
#                                            "Chao","Shannon","Chao","Shannon"))
p1<-ggplot(df.imp2,aes(x=id,y=abs(MSE)))+
  geom_bar(stat="identity",width=0.75,fill="black")+
  geom_text(aes(y=abs(MSE)*1.05,label=sig),position="stack",stat="identity")+
  ylab("Increase in MSE (%)")+
  ggtitle("Explain Variation: 85.25%")+
  scale_y_continuous(position="right")+
  annotate("text",x=-Inf,y=Inf,hjust=2,vjust=-9,
           label="R^2=0.852, P=0.001")+coord_flip()+
  theme_bw()+
  theme(axis.text.x=element_text(color="black",hjust=0.7,vjust=0),
        axis.text.y=element_text(color="black"),axis.title.y=element_blank(),
        axis.title.x=element_text(),legend.title=element_blank())
p1
ggsave("PN_Env.pdf",p1,width=3.5,height=6)

##################################################################
##################################################################
##################################################################
####ramdorest for crop VS Environmental factors###################
set.seed(123)
fit2.rp <- rfPermute(yield ~ .,ntree = 1000, nrep = 999, num.cores = 4,
                     data = beta_all[,c(1:8,40)])#num.cores=1或2
fit2.rp#67.71_67.22_0.001
rmp2<-rp.importance(fit2.rp)
View(rmp2)
plot(rp.importance(fit2.rp, scale = TRUE))

rf.perm2<-rf.significance(fit2.rp,df_env_crop[c("yield")],nperm=999,ntree=1000)
summary(rf.perm2)

##ggplot for rmp2
###ggplot for df.imp1
df.imp2 <- as.data.frame(rmp2[, 1:2])
#df.imp2 <- mutate(df.imp2, type = row.names(df.imp2))
df.imp2<-df.imp2[order(rownames(df.imp2),decreasing=F),]
colnames(df.imp2) <- c("MSE", "MSE_pval")
df.imp3 <- mutate(df.imp2, sig= ifelse (df.imp2[2] < 0.05&df.imp2[2]>=0.01, "*",
                                        ifelse(df.imp2[2]>=0.001&df.imp2[2]<0.01,"**",
                                               ifelse(df.imp2[2]<=0.001,"***",""))))

row.names(df.imp3) <- row.names(df.imp2)
df.imp3$id<-rownames(df.imp3)
df.imp3<-df.imp3[order(df.imp3$MSE,decreasing=F),]
df.imp3$id<-factor(df.imp3$id,levels=df.imp3$id)
#df.imp3<-separate(df.imp3,id,c("gene","alpha"),sep="_")
View(df.imp3)
#df.imp3$type<-c(rep("aoa",2),rep("aob",2),rep("com",2),rep("nob",2))
#View(df.imp3)
#df.imp3$alpha<-factor(df.imp3$id,level=c("Chao","Shannon","Chao","Shannon",
#                                            "Chao","Shannon","Chao","Shannon"))
p2<-ggplot(df.imp3,aes(x=id,y=abs(MSE)))+
  geom_bar(stat="identity",width=0.75,fill="black")+
  geom_text(aes(y=abs(MSE)*1.05,label=sig),position="stack",stat="identity")+
  ylab("Increase in MSE (%)")+
  ggtitle("Explain Variation: 67.71%")+
  scale_y_continuous(position="right")+
  annotate("text",x=-Inf,y=Inf,hjust=2,vjust=-9,
           label="R^2=0.672, P=0.001")+coord_flip()+
  theme_bw()+
  theme(axis.text.x=element_text(color="black",hjust=0.7,vjust=0),
        axis.text.y=element_text(color="black"),axis.title.y=element_blank(),
        axis.title.x=element_text(),legend.title=element_blank())
p2
ggsave("Yield_Env.pdf",p2,width=3.5,height=6)



####ramdorest for AOA NMDS1 VS Environmental factors######
set.seed(123)
fit3.rp <- rfPermute(aoa_MDS1 ~ .,ntree = 1000, nrep = 9999, num.cores = 4,
                     data = beta_all[,c(1:8,9)])#num.cores=1或2
fit3.rp
#14.12
plot(rp.importance(fit3.rp, scale = TRUE))
rmp3<-rp.importance(fit3.rp)
View(rmp3)

rf.perm3<-rf.significance(fit3.rp,beta[c("aoa_MDS1")],nperm=999,ntree=1000)
summary(rf.perm3)
#R_14.04 p=0.047
##ggplot for AOA_NMDS1
##ggplot for rmp2

###ggplot for df.imp1
df.imp3 <- as.data.frame(rmp3[, 1:2])
#df.imp3 <- mutate(df.imp3, type = row.names(df.imp3))
df.imp3<-df.imp3[order(rownames(df.imp3),decreasing=F),]
colnames(df.imp3) <- c("MSE", "MSE_pval")
df.imp4 <- mutate(df.imp3, sig= ifelse (df.imp3[2] < 0.05&df.imp3[2]>=0.01, "*",
                                        ifelse(df.imp3[2]>=0.001&df.imp3[2]<0.01,"**",
                                               ifelse(df.imp3[2]<=0.001,"***",""))))

row.names(df.imp4) <- row.names(df.imp3)
df.imp4$id<-rownames(df.imp4)
df.imp4<-df.imp4[order(df.imp4$MSE,decreasing=F),]
df.imp4$id<-factor(df.imp4$id,levels=df.imp4$id)
#df.imp4<-separate(df.imp4,id,c("gene","alpha"),sep="_")
View(df.imp4)
#df.imp4$type<-c(rep("aoa",2),rep("aob",2),rep("com",2),rep("nob",2))
#View(df.imp4)
#df.imp4$alpha<-factor(df.imp4$id,level=c("Chao","Shannon","Chao","Shannon",
#                                            "Chao","Shannon","Chao","Shannon"))
p3<-ggplot(df.imp4,aes(x=id,y=abs(MSE)))+
  geom_bar(stat="identity",width=0.75,fill="black")+
  geom_text(aes(y=abs(MSE)*1.05,label=sig),position="stack",stat="identity")+
  ylab("Increase in MSE (%)")+
  ggtitle("Explain Variation: 14.12%")+
  scale_y_continuous(position="right")+
  annotate("text",x=-Inf,y=Inf,hjust=2,vjust=-9,
           label="R^2=0.140, P=0.047")+coord_flip()+
  theme_bw()+
  theme(axis.text.x=element_text(color="black",hjust=0.7,vjust=0),
        axis.text.y=element_text(color="black"),axis.title.y=element_blank(),
        axis.title.x=element_text(),legend.title=element_blank())
p3
ggsave("AOA_NMDS1_Env.pdf",p3,width=3.5,height=6)


####ramdorest for AOB NMDS1 VS Environmental factors######
set.seed(123)
fit4.rp <- rfPermute(aob_MDS1 ~ .,ntree = 1000, nrep = 9999, num.cores = 4,
                     data = beta_all[,c(1:8,11)])#num.cores=1或2
fit4.rp
#73.88
plot(rp.importance(fit4.rp, scale = TRUE))
rmp4<-rp.importance(fit4.rp)
View(rmp4)

rf.perm4<-rf.significance(fit4.rp,beta[c("aob_MDS1")],nperm=999,ntree=1000)
summary(rf.perm4)
#R2=0.738 p=0.001

##ggplot for AOB NMDS1
df.imp4 <- as.data.frame(rmp4[, 1:2])
#df.imp4 <- mutate(df.imp4, type = row.names(df.imp4))
df.imp4<-df.imp4[order(rownames(df.imp4),decreasing=F),]
colnames(df.imp4) <- c("MSE", "MSE_pval")
df.imp5 <- mutate(df.imp4, sig= ifelse (df.imp4[2] < 0.05&df.imp4[2]>=0.01, "*",
                                        ifelse(df.imp4[2]>=0.001&df.imp4[2]<0.01,"**",
                                               ifelse(df.imp4[2]<=0.001,"***",""))))

row.names(df.imp5) <- row.names(df.imp4)
df.imp5$id<-rownames(df.imp5)
df.imp5<-df.imp5[order(df.imp5$MSE,decreasing=F),]
df.imp5$id<-factor(df.imp5$id,levels=df.imp5$id)
#df.imp5<-separate(df.imp5,id,c("gene","alpha"),sep="_")
View(df.imp5)
#df.imp5$type<-c(rep("aoa",2),rep("aob",2),rep("com",2),rep("nob",2))
#View(df.imp5)
#df.imp5$alpha<-factor(df.imp5$id,level=c("Chao","Shannon","Chao","Shannon",
#                                            "Chao","Shannon","Chao","Shannon"))
p4<-ggplot(df.imp5,aes(x=id,y=abs(MSE)))+
  geom_bar(stat="identity",width=0.75,fill="black")+
  geom_text(aes(y=abs(MSE)*1.05,label=sig),position="stack",stat="identity")+
  ylab("Increase in MSE (%)")+
  ggtitle("Explain Variation: 73.9%")+
  scale_y_continuous(position="right")+
  annotate("text",x=-Inf,y=Inf,hjust=2,vjust=-9,
           label="R^2=0.738, P=0.001")+coord_flip()+
  theme_bw()+
  theme(axis.text.x=element_text(color="black",hjust=0.7,vjust=0),
        axis.text.y=element_text(color="black"),axis.title.y=element_blank(),
        axis.title.x=element_text(),legend.title=element_blank())
p4
ggsave("AOB_NMDS1_Env.pdf",p4,width=3.5,height=6)


####ramdorest for nob NMDS1 VS Environmental factors######
set.seed(123)
fit5.rp <- rfPermute(nob_MDS2 ~ .,ntree = 1000, nrep = 9999, num.cores = 4,
                     data = beta_all[,c(1:8,14)])#num.cores=1或2
fit5.rp
#31.56
plot(rp.importance(fit5.rp, scale = TRUE))
rmp5<-rp.importance(fit5.rp)
View(rmp5)

rf.perm5<-rf.significance(fit5.rp,beta[c("nob_MDS1")],nperm=999,ntree=1000)
summary(rf.perm5)
#31.4 0.002

##ggplot for NOB NMDS1
df.imp5 <- as.data.frame(rmp5[, 1:2])
#df.imp5 <- mutate(df.imp5, type = row.names(df.imp5))
df.imp5<-df.imp5[order(rownames(df.imp5),decreasing=F),]
colnames(df.imp5) <- c("MSE", "MSE_pval")
df.imp6 <- mutate(df.imp5, sig= ifelse (df.imp5[2] < 0.05&df.imp5[2]>=0.01, "*",
                                        ifelse(df.imp5[2]>=0.001&df.imp5[2]<0.01,"**",
                                               ifelse(df.imp5[2]<=0.001,"***",""))))

row.names(df.imp6) <- row.names(df.imp5)
df.imp6$id<-rownames(df.imp6)
df.imp6<-df.imp6[order(df.imp6$MSE,decreasing=F),]
df.imp6$id<-factor(df.imp6$id,levels=df.imp6$id)
#df.imp6<-separate(df.imp6,id,c("gene","alpha"),sep="_")
View(df.imp6)
#df.imp6$type<-c(rep("aoa",2),rep("aob",2),rep("com",2),rep("nob",2))
#View(df.imp6)
#df.imp6$alpha<-factor(df.imp6$id,level=c("Chao","Shannon","Chao","Shannon",
#                                            "Chao","Shannon","Chao","Shannon"))
p5<-ggplot(df.imp6,aes(x=id,y=abs(MSE)))+
  geom_bar(stat="identity",width=0.75,fill="black")+
  geom_text(aes(y=abs(MSE)*1.05,label=sig),position="stack",stat="identity")+
  ylab("Increase in MSE (%)")+
  ggtitle("Explain Variation: 73.9%")+
  scale_y_continuous(position="right")+
  annotate("text",x=-Inf,y=Inf,hjust=2,vjust=-9,
           label="R^2=0.738, P=0.001")+coord_flip()+
  theme_bw()+
  theme(axis.text.x=element_text(color="black",hjust=0.7,vjust=0),
        axis.text.y=element_text(color="black"),axis.title.y=element_blank(),
        axis.title.x=element_text(),legend.title=element_blank())
p5
ggsave("NOB_NMDS2_Env.pdf",p5,width=3.5,height=6)


####ramdorest for cmx NMDS1 VS Environmental factors######
set.seed(123)
View(beta_all)
fit6.rp <- rfPermute(com_MDS2 ~ .,ntree = 1000, nrep = 9999, num.cores = 4,
                     data = beta_all[,c(1:8,16)])#num.cores=1或2
fit6.rp
plot(rp.importance(fit6.rp, scale = TRUE))
rmp6<-rp.importance(fit6.rp)
View(rmp6)
write.csv(rmp6,file="rmp6.csv")

rf.perm6<-rf.significance(fit6.rp,beta[c("com_MDS1")],nperm=999,ntree=1000)
summary(rf.perm6)

##ggplot for comammox NMDS1
df.imp6 <- as.data.frame(rmp6[, 1:2])
#df.imp6 <- mutate(df.imp6, type = row.names(df.imp6))
df.imp6<-df.imp6[order(rownames(df.imp6),decreasing=F),]
colnames(df.imp6) <- c("MSE", "MSE_pval")
df.imp7 <- mutate(df.imp6, sig= ifelse (df.imp6[2] < 0.05&df.imp6[2]>=0.01, "*",
                                        ifelse(df.imp6[2]>=0.001&df.imp6[2]<0.01,"**",
                                               ifelse(df.imp6[2]<=0.001,"***",""))))

row.names(df.imp7) <- row.names(df.imp6)
df.imp7$id<-rownames(df.imp7)
df.imp7<-df.imp7[order(df.imp7$MSE,decreasing=F),]
df.imp7$id<-factor(df.imp7$id,levels=df.imp7$id)
#df.imp7<-separate(df.imp7,id,c("gene","alpha"),sep="_")
View(df.imp7)
#df.imp7$type<-c(rep("aoa",2),rep("aob",2),rep("com",2),rep("nob",2))
#View(df.imp7)
#df.imp7$alpha<-factor(df.imp7$id,level=c("Chao","Shannon","Chao","Shannon",
#                                            "Chao","Shannon","Chao","Shannon"))
p6<-ggplot(df.imp7,aes(x=id,y=abs(MSE)))+
  geom_bar(stat="identity",width=0.75,fill="black")+
  geom_text(aes(y=abs(MSE)*1.05,label=sig),position="stack",stat="identity")+
  ylab("Increase in MSE (%)")+
  ggtitle("Explain Variation: 73.9%")+
  scale_y_continuous(position="right")+
  annotate("text",x=-Inf,y=Inf,hjust=2,vjust=-9,
           label="R^2=0.738, P=0.001")+coord_flip()+
  theme_bw()+
  theme(axis.text.x=element_text(color="black",hjust=0.7,vjust=0),
        axis.text.y=element_text(color="black"),axis.title.y=element_blank(),
        axis.title.x=element_text(),legend.title=element_blank())
p6
ggsave("COM_NMDS2_Env.pdf",p5,width=3.5,height=6)


###########################################################
###########################################################
####ramdorest for AOA gene  VS Environmental factors######
View(beta_env)
View(beta_all)
beta2$PH<-log10(beta2$PH)+1
set.seed(123)
fit7.rp <- rfPermute(AOA_amoA ~ .,ntree = 1000, nrep = 9999, num.cores = 5,
                     data = beta_all[,c(1:8,31)])#num.cores=1或2
plot(rp.importance(fit7.rp, scale = TRUE))
rmp7<-rp.importance(fit7.rp)
View(rmp7)

rf.perm7<-rf.significance(fit7.rp,df_env_crop[c("AOA_amoA")],nperm=999,ntree=1000)
summary(rf.perm7)


####ramdorest for AOB gene VS Environmental factors######
set.seed(123)
fit8.rp <- rfPermute(AOB_amoA ~ .,ntree = 1000, nrep = 9999, num.cores = 4,
                     data = beta_all[,c(1:8,32)])#num.cores=1或2
plot(rp.importance(fit8.rp, scale = TRUE))
rmp8<-rp.importance(fit8.rp)
View(rmp8)

rf.perm8<-rf.significance(fit8.rp,df_env_crop[c("AOB_amoA")],nperm=999,ntree=1000)
summary(rf.perm8)


####ramdorest for NOB gene VS Environmental factors######
set.seed(123)
fit9.rp <- rfPermute(Nitrobacter_nxrA ~ .,ntree = 1000, nrep = 9999, num.cores = 4,
                     data = beta_all[,c(1:8,33)])#num.cores=1或2
plot(rp.importance(fit9.rp, scale = TRUE))
rmp9<-rp.importance(fit9.rp)
View(rmp9)

rf.perm9<-rf.significance(fit9.rp,df_env_crop[c("Nitrobacter_nxrA")],nperm=999,ntree=1000)
summary(rf.perm9)

####ramdorest for CMX gene VS Environmental factors######
set.seed(123)
fit10.rp <- rfPermute(Comammox ~ .,ntree = 1000, nrep = 9999, num.cores = 4,
                     data = beta_all[,c(1:8,34)])#num.cores=1或2
fit10.rp
#92.54
plot(rp.importance(fit10.rp, scale = TRUE))
rmp10<-rp.importance(fit10.rp)
View(rmp10)

rf.perm10<-rf.significance(fit10.rp,df_env_crop[c("Comammox")],nperm=999,ntree=1000)
summary(rf.perm10)

##ggplot for comammox NMDS1
df.imp6 <- as.data.frame(rmp10[, 1:2])
#df.imp6 <- mutate(df.imp6, type = row.names(df.imp6))
df.imp6<-df.imp6[order(rownames(df.imp6),decreasing=F),]
colnames(df.imp6) <- c("MSE", "MSE_pval")
df.imp7 <- mutate(df.imp6, sig= ifelse (df.imp6[2] < 0.05&df.imp6[2]>=0.01, "*",
                                        ifelse(df.imp6[2]>=0.001&df.imp6[2]<0.01,"**",
                                               ifelse(df.imp6[2]<=0.001,"***",""))))

row.names(df.imp7) <- row.names(df.imp6)
df.imp7$id<-rownames(df.imp7)
df.imp7<-df.imp7[order(df.imp7$MSE,decreasing=F),]
df.imp7$id<-factor(df.imp7$id,levels=df.imp7$id)
#df.imp7<-separate(df.imp7,id,c("gene","alpha"),sep="_")
View(df.imp7)
#df.imp7$type<-c(rep("aoa",2),rep("aob",2),rep("com",2),rep("nob",2))
#View(df.imp7)
#df.imp7$alpha<-factor(df.imp7$id,level=c("Chao","Shannon","Chao","Shannon",
#                                            "Chao","Shannon","Chao","Shannon"))
p6<-ggplot(df.imp7,aes(x=id,y=abs(MSE)))+
  geom_bar(stat="identity",width=0.75,fill="black")+
  geom_text(aes(y=abs(MSE)*1.05,label=sig),position="stack",stat="identity")+
  ylab("Increase in MSE (%)")+
  ggtitle("Explain Variation: 92.54%")+
  scale_y_continuous(position="right")+
  annotate("text",x=-Inf,y=Inf,hjust=2,vjust=-9,
           label="R^2=0.915, P=0.001")+coord_flip()+
  theme_bw()+
  theme(axis.text.x=element_text(color="black",hjust=0.7,vjust=0),
        axis.text.y=element_text(color="black"),axis.title.y=element_blank(),
        axis.title.x=element_text(),legend.title=element_blank())
p6
ggsave("COM_GENE_Env.pdf",p6,width=3.5,height=6)

p_yield_pn<-ggarrange(p1,p2,nrow=1,ncol=2,labels=c("a","b"))
p_yield_pn
ggsave("p_yield_pn.pdf",p_yield_pn,height=6,width=8.5)


p_nitrifying<-ggarrange(p3,p4,p5,p6,nrow=2,ncol=2,labels=c("a","b","c","d"))
p_nitrifying
ggsave("p_nitrifying.pdf",p_nitrifying,height=13,width=8.5)
ggsave("p_nitrifying_2.pdf",p_nitrifying,height=9,width=10)
ggsave("p_nitrifying_3.pdf",p_nitrifying,height=8,width=10)

###############################################################
#################################################################
##############################################################
##################################################################

corr_meta_core<-function(df1,df2){
  corr_1<-rcorr(as.matrix(df1),as.matrix(df2),type = "spearman")
  #str(corr_1)
  corr_1r<-corr_1$r[-c(9:20),-c(1:8)]
  corr_1p<-corr_1$P[-c(9:20),-c(1:8)]
  #corr_1r[abs(corr_1r)<0.4]<-0
  #View(corr_1r)
  corr_1p[corr_1p<0.05]<-1
  corr_1p[corr_1p<1]<-0
  #View(corr_1p)
  corr_1_rp<-corr_1r*corr_1p
  #View(corr_1_rp)
  return(corr_1_rp)
}
##dt_16s
corr_dt_16s_mp<-corr_meta_core(mp[1:15,],core_dt7_16s)
write.csv(corr_dt_16s_mp,"corr_dt_16s_mp.csv")

##Spearman correlation for alpha diversity#############
corr_1<-rcorr(as.matrix(beta_all[,1:8]),as.matrix(beta_all[,19:30]),type = "spearman")
corr_1r<-corr_1$r[-c(9:20),-c(1:8)]
corr_1p<-corr_1$P[-c(9:20),-c(1:8)]
View(corr_1p)
corr_1p[corr_1p<0.05]<-1
corr_1p[corr_1p<1]<-0
#View(corr_1p)
corr_alpha<-corr_1r*corr_1p
View(corr_alpha)
write.csv(corr_alpha,"corr_alpha.csv")
##Spearman correlation for beta diversity###############
View(beta_all)
corr_2<-rcorr(as.matrix(beta_all[,1:8]),as.matrix(beta_all[,9:16]),type = "spearman")
corr_2r<-corr_2$r[-c(9:16),-c(1:8)]
corr_2p<-corr_2$P[-c(9:16),-c(1:8)]
View(corr_2p)
corr_2p[corr_2p<0.05]<-1
corr_2p[corr_2p<1]<-0
#View(corr_1p)
corr_beta<-corr_2r*corr_2p
View(corr_beta)
write.csv(corr_beta,"corr_alpha.csv")

##Spearman correlation for gene abundance, crop and PNR##
corr_3<-rcorr(as.matrix(beta_all[,1:8]),as.matrix(beta_all[,c(18,31:34,40)]),type = "spearman")
corr_3r<-corr_3$r[-c(9:14),-c(1:8)]
corr_3p<-corr_3$P[-c(9:14),-c(1:8)]
View(corr_3$r)
corr_3p[corr_3p<0.05]<-1
corr_3p[corr_3p<1]<-0
#View(corr_1p)
corr_gene<-corr_3r*corr_3p
View(corr_gene)
write.csv(corr_beta,"corr_alpha.csv")

##mantel test for alpha diversity########################


##mantel test for beta diversity#########################


##mantel test for gene abundance, crop and PNR###########


#################ggcor for aoa,aob,com,nob###############
df_ecgd<-read.csv("df_env_crop_gene_divesity.csv",row.names=1)
View(df_ecgd)
df_env<-df_ecgd[,c(1:14,44)]
View(df_env)
df_env2<-df_env[,-c(7,8,13)]
View(df_env2)
df_aoa<-df_ecgd[,c(16,27,35)]
View(df_aoa)
df_aob<-df_ecgd[,c(19,29,36)]
df_nob<-df_ecgd[,c(25,33,37)]
df_cmx<-df_ecgd[,c(22,31,38)]


#aoa#
#packageVersion("ggcor")
corr_aoa <- fortify_cor(df_aoa, df_env2)
str(corr_aoa)
##spec.select选择最终呈现变量所在的单列或者多列
mantel <- mantel_test(df_aoa, df_env2,
                      spec.select = list(a_diversity = 1,
                                         B_diversity = 2,
                                         Gene_AOA = 3
                                         )) %>% 
  mutate(rd = cut(r, breaks = c(-Inf, 0.2, 0.4, Inf),
                  labels = c("< 0.2", "0.2 - 0.4", ">= 0.4")),
         pd = cut(p.value, breaks = c(-Inf, 0.01, 0.05, Inf),
                  labels = c("< 0.01", "0.01 - 0.05", ">= 0.05")))

aoa<-quickcor(df_env2, type = "upper") +
  geom_square() +
  anno_link(aes(colour = pd, size = rd), data = mantel) +
  scale_size_manual(values = c(0.5, 1, 2)) +
  scale_colour_manual(values = c("#D95F02", "#1B9E77", "#A2A2A288")) +
  guides(size = guide_legend(title = "Mantel's r",
                             override.aes = list(colour = "grey35"), 
                             order = 2),
         colour = guide_legend(title = "Mantel's p", 
                               override.aes = list(size = 3), 
                               order = 1),
         fill = guide_colorbar(title = "Pearson's r", 
                               order = 3))
aoa
ggsave("aoa_2.pdf",aoa,width=9,height=7)


#AOB
mantel <- mantel_test(df_aob, df_env2,
                      spec.select = list(a_diversity = 1,
                                         B_diversity = 2,
                                         Gene_AOB = 3
                      )) %>% 
  mutate(rd = cut(r, breaks = c(-Inf, 0.2, 0.4, Inf),
                  labels = c("< 0.2", "0.2 - 0.4", ">= 0.4")),
         pd = cut(p.value, breaks = c(-Inf, 0.01, 0.05, Inf),
                  labels = c("< 0.01", "0.01 - 0.05", ">= 0.05")))

aob<-quickcor(df_env2, type = "upper") +
  geom_square() +
  anno_link(aes(colour = pd, size = rd), data = mantel) +
  scale_size_manual(values = c(0.5, 1, 2)) +
  scale_colour_manual(values = c("#D95F02", "#1B9E77", "#A2A2A288")) +
  guides(size = guide_legend(title = "Mantel's r",
                             override.aes = list(colour = "grey35"), 
                             order = 2),
         colour = guide_legend(title = "Mantel's p", 
                               override.aes = list(size = 3), 
                               order = 1),
         fill = guide_colorbar(title = "Pearson's r", 
                               order = 3))
aob
ggsave("aob_2.pdf",aob,width=9,height=7)

#######nob########
mantel <- mantel_test(df_nob, df_env2,
                      spec.select = list(a_diversity = 1,
                                         B_diversity = 2,
                                         Gene_nob = 3
                      )) %>% 
  mutate(rd = cut(r, breaks = c(-Inf, 0.2, 0.4, Inf),
                  labels = c("< 0.2", "0.2 - 0.4", ">= 0.4")),
         pd = cut(p.value, breaks = c(-Inf, 0.01, 0.05, Inf),
                  labels = c("< 0.01", "0.01 - 0.05", ">= 0.05")))

nob<-quickcor(df_env2, type = "upper") +
  geom_square() +
  anno_link(aes(colour = pd, size = rd), data = mantel) +
  scale_size_manual(values = c(0.5, 1, 2)) +
  scale_colour_manual(values = c("#D95F02", "#1B9E77", "#A2A2A288")) +
  guides(size = guide_legend(title = "Mantel's r",
                             override.aes = list(colour = "grey35"), 
                             order = 2),
         colour = guide_legend(title = "Mantel's p", 
                               override.aes = list(size = 3), 
                               order = 1),
         fill = guide_colorbar(title = "Pearson's r", 
                               order = 3))
nob
ggsave("nob_2.pdf",nob,width=9,height=7)

#######cmx########
mantel <- mantel_test(df_cmx, df_env2,
                      spec.select = list(a_diversity = 1,
                                         B_diversity = 2,
                                         Gene_cmx = 3
                      )) %>% 
  mutate(rd = cut(r, breaks = c(-Inf, 0.2, 0.4, Inf),
                  labels = c("< 0.2", "0.2 - 0.4", ">= 0.4")),
         pd = cut(p.value, breaks = c(-Inf, 0.01, 0.05, Inf),
                  labels = c("< 0.01", "0.01 - 0.05", ">= 0.05")))
set_scale()
cmx<-quickcor(df_env2, type = "upper") +
  geom_square() +
  anno_link(aes(colour = pd, size = rd), data = mantel) +
  scale_size_manual(values = c(0.5, 1, 2)) +
  scale_colour_manual(values = c("#D95F02", "#1B9E77", "#A2A2A288")) +
  guides(size = guide_legend(title = "Mantel's r",
                             override.aes = list(colour = "grey35"), 
                             order = 2),
         colour = guide_legend(title = "Mantel's p", 
                               override.aes = list(size = 3), 
                               order = 1),
         fill = guide_colorbar(title = "Pearson's r", 
                               order = 3))
  scale_fill_brewer(palette="Set2")
cmx
ggsave("cmx_2.pdf",cmx,width=9,height=7)



corr <- fortify_cor(varespec[1:3], varechem)
quickcor(varechem, type = "upper", cor.test = TRUE) +
  geom_square(colour = NA) +
  geom_cross(colour = "grey60") +
  anno_link(aes(colour = r < 0, size = abs(r),
                linetype = r < 0), data = corr) +
  anno_link_label(nudge_x = 0.5) +
  scale_colour_manual(values = c("#4DAF4A", "#FF7F00")) +
  scale_size(range = c(0.25, 1.5), limits = c(0, 1), breaks = c(0, 0.5, 1))



############################################################################
############################################################################
####ggcor for alpha of aoa,aob,cmx,nob###############################
View(df_env2)
df_env2$yield<-log10(df_env2$yield)
View(beta_all)
alpha_index<-beta_all[,c(19:30)]
View(alpha_index)
beta_index<-beta_all[,c(9:16)]
View(beta_index)
gene_index<-beta_all[,c(31:34)]
View(gene_index)

##ggcor for alpha
corr_alpha<- fortify_cor(beta_index, df_env2)
##spec.select选择最终呈现变量所在的单列或者多列
mantel <- mantel_test(beta_index, df_env2,
                      spec.select = list(AOA = 1:3,
                                         AOB = 4:6,
                                         Comammox = 7:9,
                                         NOB=10:11
                      )) %>% 
  mutate(rd = cut(r, breaks = c(-Inf, 0.2, 0.4, Inf),
                  labels = c("< 0.2", "0.2 - 0.4", ">= 0.4")),
         pd = cut(p.value, breaks = c(-Inf, 0.01, 0.05, Inf),
                  labels = c("< 0.01", "0.01 - 0.05", ">= 0.05")))
set_scale()
beta<-quickcor(df_env2, type = "upper") +
  geom_square() +
  anno_link(aes(colour = pd, size = rd), data = mantel) +
  scale_size_manual(values = c(0.5, 1, 2)) +
  scale_colour_manual(values = c("#D95F02", "#1B9E77", "#A2A2A288")) +
  guides(size = guide_legend(title = "Mantel's r",
                             override.aes = list(colour = "grey35"), 
                             order = 2),
         colour = guide_legend(title = "Mantel's p", 
                               override.aes = list(size = 3), 
                               order = 1),
         fill = guide_colorbar(title = "Pearson's r", 
                               order = 3))
beta
ggsave("alpha_2.pdf",alpha,width=9,height=7)


##ggcor for beta
View(beta_index)
corr_beta<- fortify_cor(beta_index, df_env2)
##spec.select选择最终呈现变量所在的单列或者多列
mantel <- mantel_test(beta_index, df_env2,
                      spec.select = list(AOA = 1:2,
                                         AOB = 3:4,
                                         Comammox = 7:8,
                                         NOB=5:6
                      )) %>% 
  mutate(rd = cut(r, breaks = c(-Inf, 0.2, 0.4, Inf),
                  labels = c("< 0.2", "0.2 - 0.4", ">= 0.4")),
         pd = cut(p.value, breaks = c(-Inf, 0.01, 0.05, Inf),
                  labels = c("< 0.01", "0.01 - 0.05", ">= 0.05")))
set_scale()
beta<-quickcor(df_env2, type = "upper") +
  geom_square() +
  anno_link(aes(colour = pd, size = rd), data = mantel) +
  scale_size_manual(values = c(0.5, 1, 2)) +
  scale_colour_manual(values = c("#D95F02", "#1B9E77", "#A2A2A288")) +
  guides(size = guide_legend(title = "Mantel's r",
                             override.aes = list(colour = "grey35"), 
                             order = 2),
         colour = guide_legend(title = "Mantel's p", 
                               override.aes = list(size = 3), 
                               order = 1),
         fill = guide_colorbar(title = "Pearson's r", 
                               order = 3))
beta
ggsave("beta_2.pdf",beta,width=9,height=7)



##############################################
##ggcor for gene###############################
View(df_env2)
View(gene_index)
df_env3<-df_env2[,1:10]
corr_gene<- fortify_cor(gene_index, df_env3)
##spec.select选择最终呈现变量所在的单列或者多列
mantel <- mantel_test(gene_index, df_env3,
                      spec.select = list(AOA = 1,
                                         AOB = 2,
                                         Comammox = 4,
                                         NOB=3
                      )) %>% 
  mutate(rd = cut(r, breaks = c(-Inf, 0.2, 0.4, Inf),
                  labels = c("< 0.2", "0.2 - 0.4", ">= 0.4")),
         pd = cut(p.value, breaks = c(-Inf, 0.01, 0.05, Inf),
                  labels = c("< 0.01", "0.01 - 0.05", ">= 0.05")))
set_scale()
gene<-quickcor(df_env3, type = "upper") +
  geom_square() +
  anno_link(aes(colour = pd, size = rd),data = mantel) +
  scale_size_manual(values = c(0.5, 1, 2)) +
  scale_colour_manual(values = c("#D95F02", "#1B9E77", "#A2A2A288")) +
  guides(size = guide_legend(title = "Mantel's r",
                             override.aes = list(colour = "grey35"), 
                             order = 2),
         colour = guide_legend(title = "Mantel's p", 
                               override.aes = list(size = 3), 
                               order = 1),
         fill = guide_colorbar(title = "Pearson's r", 
                               order = 3))
gene
ggsave("gene_2.pdf",gene,width=9,height=7)
ggsave("gene_3.pdf",gene,width=9,height=7)

###################################################################
#####################ggcor for community##############################
aoa<-read.csv("otu_rel_t_AOA.csv",row.names=1)
dim(aoa)
aob<-read.csv("otu_rel_t_AOB.csv",row.names=1)
dim(aob)
nob<-read.csv("otu_rel_t_NOB.csv",row.names=1)
dim(nob)
com<-read.csv("otu_rel_t_COM.csv",row.names=1)
dim(com)

df_group<-read.csv("group.csv",row.names=1)
View(df_group)
colnames(aoa);colnames(aob);colnames(nob);colnames(com)
df_nf<-rbind(aoa,aob,nob,com)
View(df_nf)
df_nft<-t(df_nf)
rownames(df_nft)<-df_group$sample
View(df_nft)
df_nft<-df_nft[order(rownames(df_nft)),]
#aoa:1:2191
#aob
View(df_nft[,2192:10184])
#nob:10185:11867
View(df_nft[,10185:11867])
#com:11867:13131
View(df_nft[,11867:13131])

##ggcor for gene###############################
View(df_env2)
df_env3<-df_env2[,1:10]
corr_community<- fortify_cor(df_nft, df_env3)
##spec.select选择最终呈现变量所在的单列或者多列
mantel <- mantel_test(df_nft, df_env3,
                      spec.select = list(AOA = 1:2191,
                                         AOB = 2192:10184,
                                         Comammox = 11867:13131,
                                         NOB=10185:11867
                      )) %>% 
  mutate(rd = cut(r, breaks = c(-Inf, 0.2, 0.4, Inf),
                  labels = c("< 0.2", "0.2 - 0.4", ">= 0.4")),
         pd = cut(p.value, breaks = c(-Inf, 0.01, 0.05, Inf),
                  labels = c("< 0.01", "0.01 - 0.05", ">= 0.05")))
set_scale()
community3<-quickcor(df_env3, type = "upper") +
  geom_square() +
  anno_link(aes(colour = pd, size = rd), data = mantel) +
  scale_size_manual(values = c(0.5, 1, 2)) +
  scale_colour_manual(values = c("#D95F02", "#1B9E77", "#A2A2A288")) +
  guides(size = guide_legend(title = "Mantel's r",
                             override.aes = list(colour = "grey35"), 
                             order = 2),
         colour = guide_legend(title = "Mantel's p", 
                               override.aes = list(size = 3), 
                               order = 1),
         fill = guide_colorbar(title = "Pearson's r", 
                               order = 3))
community3
ggsave("community_3.pdf",community3,width=9,height=7)


##################ggcor for PN_YIELD with env_factor, alhpa,beta,gene############
View(df_env2)
View(beta_all)
df_env3<-df_env2[,1:10]
df_pn_yield<-df_env2[,11:12]
df_alpha<-beta_all[,c(19:30)]
df_beta<-beta_all[,c(9:16)]
df_gene<-beta_all[,c(31:39)]

###########ggcor for PN_YIELD with env_factor
corr_py<- fortify_cor(df_pn_yield, df_env3)
##spec.select选择最终呈现变量所在的单列或者多列
mantel <- mantel_test(df_pn_yield, df_env3,
                      spec.select = list(PN = 1,
                                         yield=2
                      )) %>% 
  mutate(rd = cut(r, breaks = c(-Inf, 0.2, 0.4, Inf),
                  labels = c("< 0.2", "0.2 - 0.4", ">= 0.4")),
         pd = cut(p.value, breaks = c(-Inf, 0.01, 0.05, Inf),
                  labels = c("< 0.01", "0.01 - 0.05", ">= 0.05")))
set_scale()
py<-quickcor(df_env3, type = "upper") +
  geom_square() +
  anno_link(aes(colour = pd, size = rd), data = mantel) +
  scale_size_manual(values = c(0.5, 1, 2)) +
  scale_colour_manual(values = c("#D95F02", "#1B9E77", "#A2A2A288")) +
  guides(size = guide_legend(title = "Mantel's r",
                             override.aes = list(colour = "grey35"), 
                             order = 2),
         colour = guide_legend(title = "Mantel's p", 
                               override.aes = list(size = 3), 
                               order = 1),
         fill = guide_colorbar(title = "Pearson's r", 
                               order = 3))
py
ggsave("py_2.pdf",py,width=9,height=7)

##############################################
###########ggcor for PN_YIELD with alpha#############
corr_py_alpha<- fortify_cor(df_pn_yield, df_alpha)
##spec.select选择最终呈现变量所在的单列或者多列
mantel <- mantel_test(df_pn_yield, df_alpha,
                      spec.select = list(PN = 1,
                                         yield=2
                      )) %>% 
  mutate(rd = cut(r, breaks = c(-Inf, 0.2, 0.4, Inf),
                  labels = c("< 0.2", "0.2 - 0.4", ">= 0.4")),
         pd = cut(p.value, breaks = c(-Inf, 0.01, 0.05, Inf),
                  labels = c("< 0.01", "0.01 - 0.05", ">= 0.05")))
set_scale()
py_alpha<-quickcor(df_alpha, type = "upper") +
  geom_square() +
  anno_link(aes(colour = pd, size = rd), data = mantel) +
  scale_size_manual(values = c(0.5, 1, 2)) +
  scale_colour_manual(values = c("#D95F02", "#1B9E77", "#A2A2A288")) +
  guides(size = guide_legend(title = "Mantel's r",
                             override.aes = list(colour = "grey35"), 
                             order = 2),
         colour = guide_legend(title = "Mantel's p", 
                               override.aes = list(size = 3), 
                               order = 1),
         fill = guide_colorbar(title = "Pearson's r", 
                               order = 3))
py_alpha
ggsave("py_alpha.pdf",py_alpha,width=9,height=7)

###################################################
###########ggcor for PN_YIELD with beta
corr_py_alpha<- fortify_cor(df_pn_yield, df_beta)
##spec.select选择最终呈现变量所在的单列或者多列
mantel <- mantel_test(df_pn_yield, df_beta,
                      spec.select = list(PN = 1,
                                         yield=2
                      )) %>% 
  mutate(rd = cut(r, breaks = c(-Inf, 0.2, 0.4, Inf),
                  labels = c("< 0.2", "0.2 - 0.4", ">= 0.4")),
         pd = cut(p.value, breaks = c(-Inf, 0.01, 0.05, Inf),
                  labels = c("< 0.01", "0.01 - 0.05", ">= 0.05")))
set_scale()
py_beta<-quickcor(df_beta, type = "upper") +
  geom_square() +
  anno_link(aes(colour = pd, size = rd), data = mantel) +
  scale_size_manual(values = c(0.5, 1, 2)) +
  scale_colour_manual(values = c("#D95F02", "#1B9E77", "#A2A2A288")) +
  guides(size = guide_legend(title = "Mantel's r",
                             override.aes = list(colour = "grey35"), 
                             order = 2),
         colour = guide_legend(title = "Mantel's p", 
                               override.aes = list(size = 3), 
                               order = 1),
         fill = guide_colorbar(title = "Pearson's r", 
                               order = 3))
py_beta
ggsave("py_beta.pdf",py_beta,width=9,height=7)

########################################################
###########ggcor for PN_YIELD with gene
View(df_gene)
df_gene$AOA_AOB<-log10(df_gene$AOA_AOB)+1
df_gene$COM_AOB<-log10(df_gene$COM_AOB)+1
df_gene$COM_AOB_NOB<-log10(df_gene$COM_AOB_NOB)+1
corr_py_alpha<- fortify_cor(df_pn_yield, df_gene)
##spec.select选择最终呈现变量所在的单列或者多列
mantel <- mantel_test(df_pn_yield, df_gene,
                      spec.select = list(PN = 1,
                                         yield=2
                      )) %>% 
  mutate(rd = cut(r, breaks = c(-Inf, 0.2, 0.4, Inf),
                  labels = c("< 0.2", "0.2 - 0.4", ">= 0.4")),
         pd = cut(p.value, breaks = c(-Inf, 0.01, 0.05, Inf),
                  labels = c("< 0.01", "0.01 - 0.05", ">= 0.05")))
set_scale()
py_gene<-quickcor(df_gene, type = "upper") +
  geom_square() +
  anno_link(aes(colour = pd, size = rd), data = mantel) +
  scale_size_manual(values = c(0.5, 1, 2)) +
  scale_colour_manual(values = c("#D95F02", "#1B9E77", "#A2A2A288")) +
  guides(size = guide_legend(title = "Mantel's r",
                             override.aes = list(colour = "grey35"), 
                             order = 2),
         colour = guide_legend(title = "Mantel's p", 
                               override.aes = list(size = 3), 
                               order = 1),
         fill = guide_colorbar(title = "Pearson's r", 
                               order = 3))
py_gene
ggsave("py_gene.pdf",py_gene,width=9,height=7)
