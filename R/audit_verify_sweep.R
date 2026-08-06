## audit_verify_sweep.R — verification sweep (diagnostic only, no fixes).
## NOT part of the numbered run order. Six of its seven blocks are superseded or
## belong to the deleted predictive section; only the MFA block is live content.
## Blocks: VIF | component-level S3/S4 | Comp7 AUCs | 10-fold CV | 18-40 subgroup
##          | Table S5 selection | SMOTE justification | MFA PSQI×Sex dims.
suppressPackageStartupMessages({
  library(mice); library(dplyr); library(caret); library(themis); library(recipes)
  library(pROC); library(car); library(randomForest); library(FactoMineR)
  library(factoextra)   # owns get_eigenvalue(); missing import broke the MFA block
})
root <- here::here()
logf <- file.path(root,"outputs/reports/audit_verify_sweep.txt")
say <- function(...){l<-sprintf(...);cat(l,"\n");cat(l,"\n",file=logf,append=TRUE)}
tryb <- function(label,expr){ say("\n===== %s =====",label); tryCatch(expr, error=function(e) say("  ERROR/cannot-determine: %s", conditionMessage(e))) }
cat("",file=logf); say("VERIFICATION SWEEP  START %s", format(Sys.time()))

comp <- paste0("Comp.",1:7)
imp <- readRDS(file.path(root,"data/processed/mice_imputation_m100.rds"))
long <- complete(imp,"long")
long$bmi  <- long$weight/((long$height/100)^2)
long$bmic <- factor(ifelse(long$bmi<25,"Normal",ifelse(long$bmi<30,"Overweight","Obesity")),levels=c("Normal","Overweight","Obesity"))
long$psqi.total <- rowSums(long[,comp])
long$psqi <- factor(ifelse(long$psqi.total>5,"p.sleep","g.sleep"),levels=c("g.sleep","p.sleep"))
long$sex  <- factor(long$sex, levels=c("female","male"))
n_imp <- max(long$.imp)

## ---- VIF (car::vif on glm per imputation, averaged) ----
tryb("VIF (target: sex 2.07 / age 1.04 / bmic 1.04 / bf 1.65 / hb 1.67)",{
  V <- matrix(NA,n_imp,5, dimnames=list(NULL,c("sex","age","bmic","bf","hb")))
  for(i in 1:n_imp){
    di <- long[long$.imp==i,]
    v <- car::vif(glm(psqi~sex+age+bmic+bf+hb, data=di, family=binomial()))
    if(is.matrix(v)){ # factor present -> GVIF^(1/(2Df)) col 3
      V[i,"bmic"] <- v["bmic",3]^2  # square the (1/(2Df)) adj to get VIF-scale? report both below
      V[i,"sex"]  <- v["sex",1]; V[i,"age"]<-v["age",1]; V[i,"bf"]<-v["bf",1]; V[i,"hb"]<-v["hb",1]
      V[i,"bmic"] <- v["bmic",1]   # store raw GVIF for bmic
    } else { V[i,names(v)] <- v }
  }
  cm <- colMeans(V,na.rm=TRUE)
  say("  mean over %d imps: sex=%.2f age=%.2f bmic(GVIF)=%.2f bf=%.2f hb=%.2f", n_imp, cm["sex"],cm["age"],cm["bmic"],cm["bf"],cm["hb"])
  # also GVIF^(1/2Df) for bmic
  b <- sapply(1:n_imp,function(i){di<-long[long$.imp==i,]; car::vif(glm(psqi~sex+age+bmic+bf+hb,data=di,family=binomial()))["bmic",3]})
  say("  bmic GVIF^(1/(2*Df)) mean = %.2f (this is the caret/car-printed comparable value)", mean(b))
})

## ---- Component-level (09): dichotomize 0 vs 1-3, glm pooled ----
tryb("Component-level OR/p (target C6 OR 0.21 p .020 ; C7 OR 2.44 p .046)",{
  for(k in 1:7){
    cc <- comp[k]
    ml <- lapply(1:n_imp,function(i){di<-long[long$.imp==i,]
      di$out <- factor(ifelse(di[[cc]]==0,"No","Yes"),levels=c("No","Yes"))
      glm(out~sex+age+bmic+bf+hb, data=di, family=binomial())})
    pl <- summary(pool(as.mira(ml)), conf.int=TRUE, exponentiate=TRUE)
    sig <- pl[pl$term!="(Intercept)" & pl$p.value<0.05,,drop=FALSE]
    line <- if(nrow(sig)) paste(sprintf("%s OR=%.3f(%.3f-%.3f)p=%.4f",sig$term,sig$estimate,sig$`2.5 %`,sig$`97.5 %`,sig$p.value),collapse=" | ") else "none<.05"
    say("  %s: %s", cc, line)
  }
})

## ---- Comp7 dedicated (13): no-SMOTE OR + no-SMOTE & SMOTE AUC ----
tryb("Comp7 (target OR Overwt 2.44 p .046 ; AUC 0.549 vs 0.546)",{
  set.seed(123456); fi <- long[long$.imp==1,]
  fi$c7 <- factor(ifelse(fi$Comp.7==0,"No","Yes"),levels=c("No","Yes"))
  idx <- as.integer(createDataPartition(fi$c7,p=0.80,list=FALSE))
  ml <- lapply(1:n_imp,function(i){di<-long[long$.imp==i,]
    di$c7 <- factor(ifelse(di$Comp.7==0,"No","Yes"),levels=c("No","Yes"))
    glm(c7~sex+age+bmic+bf+hb,data=di,family=binomial())})
  pl <- summary(pool(as.mira(ml)),conf.int=TRUE,exponentiate=TRUE)
  ow <- pl[pl$term=="bmicOverweight",]
  say("  no-SMOTE pooled Overweight OR=%.3f (%.3f-%.3f) p=%.4f", ow$estimate,ow$`2.5 %`,ow$`97.5 %`,ow$p.value)
  aN <- aS <- numeric(n_imp)
  for(i in 1:n_imp){ di<-long[long$.imp==i,]; di$c7<-factor(ifelse(di$Comp.7==0,"No","Yes"),levels=c("No","Yes"))
    tr<-di[idx,c("sex","age","bmic","bf","hb","c7")]; te<-di[-idx,c("sex","age","bmic","bf","hb","c7")]
    mN<-glm(c7~sex+age+bmic+bf+hb,data=tr,family=binomial())
    aN[i]<-as.numeric(roc(te$c7,predict(mN,te,type="response"),quiet=TRUE,levels=c("No","Yes"),direction="<")$auc)
    trb<-tryCatch({set.seed(123);bake(prep(recipe(c7~.,data=tr)%>%step_smotenc(c7,over_ratio=1,seed=123,neighbors=5),training=tr,verbose=FALSE),new_data=NULL)},error=function(e)tr)
    mS<-glm(c7~sex+age+bmic+bf+hb,data=trb,family=binomial())
    aS[i]<-as.numeric(roc(te$c7,predict(mS,te,type="response"),quiet=TRUE,levels=c("No","Yes"),direction="<")$auc)
  }
  say("  Comp7 AUC no-SMOTE=%.3f | SMOTE=%.3f  -> map to pub 0.549 / 0.546", mean(aN), mean(aS))
})

## ---- 10-fold CV global RF (07 STEP7): sampling=smote ----
tryb("10-fold CV global (target AUC 0.683 / sens 0.600 / spec 0.682)",{
  fi <- long[long$.imp==1,c("sex","age","bmic","bf","hb","psqi")]
  levels(fi$psqi) <- c("Good_Sleep","Poor_Sleep")
  set.seed(123456)
  cvm <- train(psqi~sex+age+bmic+bf+hb, data=fi, method="rf",
               trControl=trainControl(method="cv",number=10,classProbs=TRUE,summaryFunction=twoClassSummary,sampling="smote"),
               metric="ROC", ntree=500)
  b <- which.max(cvm$results$ROC)
  say("  CV-ROC=%.3f sens=%.3f spec=%.3f (best mtry=%s)", cvm$results$ROC[b], cvm$results$Sens[b], cvm$results$Spec[b], cvm$results$mtry[b])
})

## ---- 18-40 subgroup (¶214: n=239, AUC ~0.70) ----
tryb("18-40 subgroup (target n=239, AUC ~0.70)",{
  n239 <- sum(long$.imp==1 & long$age>=18 & long$age<=40)
  say("  subgroup n (age 18-40, imp1) = %d", n239)
  set.seed(123456); fi <- long[long$.imp==1,]; sub_mask <- fi$age>=18 & fi$age<=40
  idx <- as.integer(createDataPartition(fi$psqi[sub_mask],p=0.80,list=FALSE))
  a <- numeric(n_imp)
  for(i in 1:n_imp){ di<-long[long$.imp==i,]; di<-di[di$age>=18&di$age<=40,c("sex","age","bmic","bf","hb","psqi")]
    tr<-di[idx,]; te<-di[-idx,]
    trb<-tryCatch({set.seed(123);bake(prep(recipe(psqi~.,data=tr)%>%step_smotenc(psqi,over_ratio=1,seed=123,neighbors=5),training=tr,verbose=FALSE),new_data=NULL)},error=function(e)tr)
    m<-glm(psqi~sex+age+bmic+bf+hb,data=trb,family=binomial())
    a[i]<-as.numeric(roc(te$psqi,predict(m,te,type="response"),quiet=TRUE,levels=c("g.sleep","p.sleep"),direction="<")$auc)
  }
  say("  subgroup mean AUC (SMOTE LR) = %.3f", mean(a))
})


## ---- SMOTE justification global (11): LR no-SMOTE vs SMOTE ----
tryb("SMOTE justification (target AUC .690 vs .702 p.006; spec .727 vs .750 p.022)",{
  set.seed(123456); fi<-long[long$.imp==1,]; idx<-as.integer(createDataPartition(fi$psqi,p=0.80,list=FALSE))
  aN<-aS<-spN<-spS<-snN<-snS<-numeric(n_imp)
  for(i in 1:n_imp){ di<-long[long$.imp==i,c("sex","age","bmic","bf","hb","psqi")]; tr<-di[idx,]; te<-di[-idx,]
    mN<-glm(psqi~sex+age+bmic+bf+hb,data=tr,family=binomial()); pN<-predict(mN,te,type="response")
    rN<-roc(te$psqi,pN,quiet=TRUE,levels=c("g.sleep","p.sleep"),direction="<"); aN[i]<-as.numeric(rN$auc)
    coN<-coords(rN,"best",best.method="youden",ret=c("sensitivity","specificity")); snN[i]<-coN$sensitivity[1]; spN[i]<-coN$specificity[1]
    trb<-tryCatch({set.seed(123);bake(prep(recipe(psqi~.,data=tr)%>%step_smotenc(psqi,over_ratio=1,seed=123,neighbors=5),training=tr,verbose=FALSE),new_data=NULL)},error=function(e)tr)
    mS<-glm(psqi~sex+age+bmic+bf+hb,data=trb,family=binomial()); pS<-predict(mS,te,type="response")
    rS<-roc(te$psqi,pS,quiet=TRUE,levels=c("g.sleep","p.sleep"),direction="<"); aS[i]<-as.numeric(rS$auc)
    coS<-coords(rS,"best",best.method="youden",ret=c("sensitivity","specificity")); snS[i]<-coS$sensitivity[1]; spS[i]<-coS$specificity[1]
  }
  say("  AUC: noSMOTE=%.3f SMOTE=%.3f Wilcoxon p=%.4f", mean(aN),mean(aS),wilcox.test(aS,aN,paired=TRUE)$p.value)
  say("  Sens: noSMOTE=%.3f SMOTE=%.3f p=%.4f | Spec: noSMOTE=%.3f SMOTE=%.3f p=%.4f",
      mean(snN),mean(snS),wilcox.test(snS,snN,paired=TRUE)$p.value, mean(spN),mean(spS),wilcox.test(spS,spN,paired=TRUE)$p.value)
})

## ---- MFA PSQI x Sex dimensions (04) ----
tryb("MFA PSQI×Sex (target Dim1 31.6% / Dim2 17.9%)",{
  fi <- long[long$.imp==1,]
  num <- fi[,c("age","height","weight","bmi","wc","bf","mm","water","hb", comp,"psqi.total")]
  d <- data.frame(psqi=fi$psqi, sex=fi$sex, num)
  d$mean <- rowMeans(num); d$std <- apply(num,1,sd)
  res <- MFA(d, group=c(2,9,8,2), type=c("n","s","s","s"),
             name.group=c("BBF","Physio","PSQI","Illus"), num.group.sup=4, graph=FALSE)
  ev <- get_eigenvalue(res)
  say("  Dim1=%.1f%% Dim2=%.1f%% (cumulative Dim1-2 = %.1f%%)", ev[1,2], ev[2,2], ev[2,3])
})

say("\nEND %s", format(Sys.time()))
