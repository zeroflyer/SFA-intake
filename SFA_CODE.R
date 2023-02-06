##case-contol study: logistic regression##
CC<-fread("cc.csv",header = T,check.names = F)
head(CC)
log_fix <- glm(pheno ~factor(prs3_sa2)+sex+enroll_age+family_history+factor(alcohol_intake_frequency)+bmi+factor(smoke)+factor(education)+factor(race)+red_meat_s+disease+factor(ukb_finance),data=CC,family = "binomial")
log_fix <- glm(pheno ~SFxPRS,data=CC,family = "binomial")
log_fix <- glm(pheno ~SFxPRS+sex+enroll_age,data=CC,family = "binomial")
summary(log_fix)
exp(coef(log_fix)[2]),ci
exp(coef(log_fix)[3])


#cohort study: Cox proportional hazard models## 
library("survival")
library("survminer")
library(data.table)
CRC<-fread("cohort.csv",header = T,check.names = F)
res.cox<-coxph(Surv(time,pheno)~PRS_iceqc2ld781w+sex+enroll_age+family_history+factor(alcohol_intake)+bmi+factor(smoke)+factor(education)+factor(race)+red_meat_s+disease+factor(ukb_finance),data=CRC)
summary(res.cox)
res.cox<-coxph(Surv(time,pheno)~factor(prs3)+sex+enroll_age+family_history+factor(alcohol_intake)+bmi+factor(smoke)+factor(education)+factor(race)+red_meat_s+disease+factor(ukb_finance),data=CRC)
res.cox<-coxph(Surv(time,pheno)~sfxprs,data=CRC)
res.cox<-coxph(Surv(time,pheno)~sfxprs+sex+enroll_age,data=CRC)
summary(res.cox)


##density plot##
library(ggplot2)
library(hrbrthemes)
library(ggsci)
BiocManager::install('hrbrthemes')
sort(CRC$pheno,decreasing = TRUE)
rank(CRC$pheno,)
library(data.table)
CRC<-fread("cohort.csv",header = T,check.names = F)
names(CRC)
CRC$LOG2SFA=log2(CRC$saturated_fat)
p<-ggplot(data=CRC, aes(x = b_society,..density..))
p + geom_density(aes(color = group,fill=group,alpha=0.2,))+
  theme_bw()+ theme(panel.grid=element_blank())
A=ggplot(data=CRC, aes(x=prsice,fill=group,color=group)) +
  geom_density(adjust=1.5, alpha=0.5,lwd=1.5,linetype = 1) +
  theme_ipsum()


##Cumulative incidence curve##
library(ggplot2) 
library("survival")
library("survminer")
library(survminer)
library(data.table)
CRC<-fread("cohort.csv",header = T,check.names = F)
ggsurvplot(survfit(Surv(time,pheno)~prs3,data=CRC),        
           fun = "event",
           conf.int = TRUE,  
           conf.int.style="ribbon",  
           conf.int.alpha=0.45,     
           surv.plot.height= 0.7, 
           surv.plot.weight= 0.7,
           main = "Survival curve",           
           legend.title = "PRS",           
           legend.labs = c("Low", "Intermediate", "High"),    
           
           xlab="Years",         
           ylab="Cumulative Event Rate (%)", 
           xlim = c(0,14), 
           ylim = c(0,0.015),      
           risk.table = TRUE,       
           tables.height = 0.15, 
           tables.theme = theme_cleantable(),     
           font.main = c(14, "bold", "darkblue"),     
           font.x = c(14, "plain", "black"),          
           break.x.by = 2, break.y.by = 0.005,         
           font.y = c(14, "plain", "black"),          
           font.tickslab = c(12, "plain", "black"),   
           risk.table.col="strata",       
           risk.table.height=0.2,         
           palette = "npg",       
           
           ggtheme=theme_survminer()
)

ggsurvplot(survfit(Surv(time,pheno)~sfa1_4,data=CRC),        
           fun = "event",
           conf.int = TRUE,  
           conf.int.style="ribbon",  
           conf.int.alpha=0.45,     
           surv.plot.height= 0.7,   surv.plot.weight= 0.7,       
           main = "Survival curve",           
           legend.title = "Saturated fat",           
           legend.labs = c("1", "4"),    
           pval = TRUE, 
           xlab="Years",         
           ylab="Cumulative Event Rate (%)", 
           xlim = c(0,14), 
           ylim = c(0,0.012), 
           break.x.by = 2, break.y.by = 0.005,
           risk.table = TRUE,       
           tables.height = 0.15, 
           tables.theme = theme_cleantable(),     
           font.main = c(14, "bold", "darkblue"),     
           font.x = c(14, "plain", "black"),          
           
           font.y = c(14, "plain", "black"),          
           font.tickslab = c(12, "plain", "black"),   
           risk.table.col="strata",       
           risk.table.height=0.2,         
           palette = "npg",       
           ggtheme=theme_survminer()
)


##Restricted cubic spline##
dt=fread("cc.csv",head=T)
install.packages("rms")
library(rms) #RCS
library(survminer)
library(ggplot2)
library(ggsci) 
dd<-datadist(dt) 
options(datadist='dd')
names(dt)
#for case-control study
fit <- lrm(pheno ~ rcs(prsice,4)+sex+enroll_age+family_history+alcohol_intake_frequency+bmi+smoke+education+race+red_meat_s+disease+ukb_finance ,data=dt)
an<-anova(fit)
plot(Predict(fit, prsice,fun=exp), anova=an, pval=T)
ggcoxzph(cox.zph(fit, "rank"))
Pre0 <-rms::Predict(fit,prsice,fun=exp,type="predictions",ref.zero=TRUE,conf.int = 0.95,digits=2)
ggplot(Pre0)
View(Pre0)
#for cohort study
S <- Surv(dt$time,dt$pheno==1)
fit <- cph(S ~ rcs(prsice,4) +sex+enroll_age+family_history+alcohol_intake+bmi+smoke+education+race+red_meat_s+disease+ukb_finance, x=TRUE, y=TRUE,data=dt)
an<-anova(fit)
plot(Predict(fit, prsice,fun=exp), anova=an, pval=T)
ggcoxzph(cox.zph(fit, "rank"))
Pre0 <-rms::Predict(fit,prsice,fun=exp,type="predictions",ref.zero=TRUE,conf.int = 0.95,digits=2)
ggplot(Pre0)
View(Pre0)
anova(fit)


#Area Under Curve, AUC
setwd("G:\\ukb")
library(data.table)
library(pROC)
CRC<-fread("cc.csv",header = T,check.names = F)
names(CRC)
library("survival")
library("survminer")
#prs
pre1 <- glm(pheno ~factor(prs3),family=binomial(link = "logit"),data = CRC)
summary(pre1)
real_T <- CRC$pheno
predict._T1 <- predict.glm(pre1,type='response',newdata=CRC)
summary(predict._T1)
modelroc_CRC1 <- roc(real_T ,predict._T1)
plot(modelroc_CRC1,xlim=c(1,0),ylim=c(0,1),col="aquamarine4",print.thres=F, print.auc=F)
#sex age
pre2 <- glm(pheno ~sex+enroll_age,family=binomial(link = "logit"),data = CRC)
predict._T2 <- predict.glm(pre2,type='response',newdata=CRC)
modelroc_CRC2 <- roc(real_T ,predict._T2)
plot(modelroc_CRC2,xlim=c(1,0),ylim=c(0,1),col="cornflowerblue",print.thres=F, print.auc=F,add=T)
#SFA
pre3 <- glm(pheno ~factor(saturatedfat4),family=binomial(link = "logit"),data = CRC)
predict._T3 <- predict.glm(pre3,type='response',newdata=CRC)
modelroc_CRC3 <- roc(real_T ,predict._T3)
plot(modelroc_CRC3,xlim=c(1,0),ylim=c(0,1),col="bisque2",print.thres=F, print.auc=F,add=T)
#11
pre4 <- glm(pheno ~sex+enroll_age+family_history+factor(alcohol_intake_frequency)+bmi+factor(smoke)+factor(education)+factor(race)+red_meat_s+disease+factor(ukb_finance),family=binomial(link = "logit"),data = CRC)
predict._T4 <- predict.glm(pre4,type='response',newdata=CRC)
modelroc_CRC4 <- roc(real_T ,predict._T4)
plot(modelroc_CRC4,xlim=c(1,0),ylim=c(0,1),col="darkorange1",print.thres=F, print.auc=F,add=T)
#prs age sex
pre5 <- glm(pheno ~factor(prs3)+sex+enroll_age,family=binomial(link = "logit"),data = CRC)
predict._T5 <- predict.glm(pre5,type='response',newdata=CRC)
modelroc_CRC5 <- roc(real_T ,predict._T5)
plot(modelroc_CRC5,xlim=c(1,0),ylim=c(0,1),col="darkolivegreen1",print.thres=F, print.auc=F,add=T)
#prs sfa
pre6 <- glm(pheno ~factor(prs3)+factor(saturatedfat4),family=binomial(link = "logit"),data = CRC)
predict._T6 <- predict.glm(pre6,type='response',newdata=CRC)
modelroc_CRC6 <- roc(real_T ,predict._T6)
plot(modelroc_CRC6,xlim=c(1,0),ylim=c(0,1),col="darkseagreen3",print.thres=F, print.auc=F,add=T)
#prs 11
pre7 <- glm(pheno ~factor(prs3)+sex+enroll_age+family_history+factor(alcohol_intake_frequency)+bmi+factor(smoke)+factor(education)+factor(race)+red_meat_s+disease+factor(ukb_finance),family=binomial(link = "logit"),data = CRC)
predict._T7 <- predict.glm(pre7,type='response',newdata=CRC)
modelroc_CRC7 <- roc(real_T ,predict._T7)
plot(modelroc_CRC7,xlim=c(1,0),ylim=c(0,1),col="brown4",print.thres=F, print.auc=F,add=T)
#prs 11 sfa
pre8 <- glm(pheno ~factor(prs3)+factor(saturatedfat4)+sex+enroll_age+family_history+factor(alcohol_intake_frequency)+bmi+factor(smoke)+factor(education)+factor(race)+red_meat_s+disease+factor(ukb_finance),family=binomial(link = "logit"),data = CRC)
predict._T8 <- predict.glm(pre8,type='response',newdata=CRC)
modelroc_CRC8 <- roc(real_T ,predict._T8)
plot(modelroc_CRC8,xlim=c(1,0),ylim=c(0,1),col="steelblue1",print.thres=F, print.auc=F,add=T)
# 11 sfa
pre9 <- glm(pheno ~factor(saturatedfat4)+sex+enroll_age+family_history+factor(alcohol_intake_frequency)+bmi+factor(smoke)+factor(education)+factor(race)+red_meat_s+disease+factor(ukb_finance),family=binomial(link = "logit"),data = CRC)
predict._T9 <- predict.glm(pre9,type='response',newdata=CRC)
modelroc_CRC9 <- roc(real_T ,predict._T9)
plot(modelroc_CRC9,xlim=c(1,0),ylim=c(0,1),col="darkviolet",print.thres=F, print.auc=F,add=T)
# prs sfa model2
pre10 <- glm(pheno ~factor(saturatedfat4)+factor(prs3)+sex+enroll_age,family=binomial(link = "logit"),data = CRC)
predict._T10 <- predict.glm(pre10,type='response',newdata=CRC)
modelroc_CRC10 <- roc(real_T ,predict._T10)
plot(modelroc_CRC10,xlim=c(1,0),ylim=c(0,1),col="lightgoldenrod2",print.thres=F, print.auc=F,add=T)
pdf("auc10_cc.pdf")
roc.test(modelroc_CRC8,modelroc_CRC9)
auc(modelroc_CRC10);ci(modelroc_CRC10)

