################################################################################
# 08_prospective_diagnosis.R
# Cumulative psychosocial exposome score (cPES) study - BHRC.
# Logistic models of psychiatric diagnosis at follow-up, adjusted for baseline CBCL.
#
# Faithful, path-genericised copy of the analysis code used for the manuscript.
# CONFIGURE PATHS: set the placeholders below to your local tree before running.
#   <PROJECT_ROOT>     project directory
#   ${cPES_OUT}        output tree (e.g. <PROJECT_ROOT>/FINAL_WITH_SUBSTANCE)
#   <BHRC_INPUT_DIR>   controlled BHRC data (.sav) - not shareable
#   <COGNITION_DIR>    cognitive-score file
# R 4.4.1; seed 241094 where stochastic.
# NB: any FDR columns written here use the per-model pool; the MANUSCRIPT FDR
#     families and the developmental-years metric are applied in
#     10_final_tables_corrected.R (run last). See README.md.
################################################################################

################################################################################
# Paper-ready deliverables: (A) domain-contribution figure (architecture
# description, NOT independent prediction); (B) prospective diagnosis table
# (unadjusted + adjusted for baseline CBCL).
################################################################################
suppressWarnings(suppressMessages({library(haven);library(dplyr);library(readr);library(janitor);library(ggplot2)}))
OUT<-"${cPES_OUT}/_EXTRA_ANALYSES"

## ---- (A) Domain figure ----
joint<-read_csv(file.path(OUT,"01_domain_joint_model.csv"),show_col_types=FALSE)
uni<-read_csv(file.path(OUT,"01_domain_univariate_R2.csv"),show_col_types=FALSE)
joint$domain<-gsub("`","",joint$domain)
p<-ggplot(joint,aes(x=std_beta,y=reorder(domain,std_beta)))+
  geom_col(fill="#0072B2")+
  geom_text(aes(label=sprintf("%.2f",std_beta)),hjust=-0.15,size=3.4)+
  labs(x="Standardized contribution to CBCL total (joint model)",y=NULL,
       title="Cumulative architecture: all seven domains contribute",
       subtitle="Descriptive decomposition of the cumulative psychosocial exposome score (all p<0.001)")+
  xlim(0,max(joint$std_beta)*1.15)+ theme_bw(12)+
  theme(plot.title=element_text(face="bold"),panel.grid.major.y=element_blank())
ggsave(file.path(OUT,"FIG_domain_contribution.png"),p,width=8,height=4.5,dpi=300,bg="white")
ggsave(file.path(OUT,"FIG_domain_contribution.pdf"),p,width=8,height=4.5)
cat("Saved domain figure.\n")

## ---- (B) Prospective table (unadj + adj for baseline CBCL) ----
auc<-function(y,p){y<-y[!is.na(p)];p<-p[!is.na(p)];n1<-sum(y==1);n0<-sum(y==0);if(n1==0||n0==0)return(NA);r<-rank(p);(sum(r[y==1])-n1*(n1+1)/2)/(n1*n0)}
cPES<-read_csv("${cPES_OUT}/01_EXWAS_EN_cPES/POOLED/07_cpes_oof_predictions.csv",show_col_types=FALSE)|>transmute(ident=as.character(ident),cPES_z=cPES_score_OOF_z)
B<-read_sav("<BHRC_INPUT_DIR>/Grassi-Oliveira_167BHRC_2023_12_09.sav",encoding="latin1")|>clean_names()
wv<-function(x)c("wave0_arm_1"=0,"wave1_arm_1"=1,"wave2_arm_1"=2)[x]
D<-B|>filter(redcap_event_name%in%c("wave0_arm_1","wave1_arm_1","wave2_arm_1"))|>transmute(ident=as.character(ident),wave=wv(redcap_event_name),age=as.numeric(age),gender=factor(gender),site=factor(site_original),dcany=dplyr::case_when(as.numeric(dcany)==2~1L,as.numeric(dcany)==0~0L,TRUE~NA_integer_))
base<-B|>filter(redcap_event_name=="wave0_arm_1")|>transmute(ident=as.character(ident),cbcl0=as.numeric(cl_tot))
# Age rule: age at assessment (visit age) is used; for the few participants without a recorded
# visit age (109 at the 8-year wave), baseline age is substituted. Results are identical under
# baseline-age-for-all (OR 1.54 at 8 years either way); the fallback is used for consistency with
# the longitudinal models and the reported Ns.
cpes0<-D|>filter(wave==0)|>transmute(ident,age0=age,gender,site)|>inner_join(cPES,by="ident")|>inner_join(base,by="ident")
tab<-lapply(0:2,function(w){dd<-D|>filter(wave==w)|>transmute(ident,dcany,age_fu=age)|>inner_join(cpes0,by="ident")|>filter(!is.na(dcany),!is.na(cbcl0))|>mutate(age=ifelse(is.na(age_fu),age0,age_fu))
  mu<-glm(dcany~cPES_z+age+gender+site,dd,family=binomial)
  if(w==0){data.frame(wave="W0",n=nrow(dd),cases=sum(dd$dcany==1),OR_unadj=round(exp(coef(mu)["cPES_z"]),2),AUC=round(auc(dd$dcany,predict(mu,type="response")),3),OR_adj=NA,adj_lo=NA,adj_hi=NA,adj_p=NA)}
  else{ma<-glm(dcany~cPES_z+scale(cbcl0)+age+gender+site,dd,family=binomial);ci<-exp(confint.default(ma)["cPES_z",])
    data.frame(wave=paste0("W",w),n=nrow(dd),cases=sum(dd$dcany==1),OR_unadj=round(exp(coef(mu)["cPES_z"]),2),AUC=round(auc(dd$dcany,predict(mu,type="response")),3),OR_adj=round(exp(coef(ma)["cPES_z"]),2),adj_lo=round(ci[1],2),adj_hi=round(ci[2],2),adj_p=signif(summary(ma)$coef["cPES_z",4],3))}})|>bind_rows()
write_csv(tab,file.path(OUT,"04_prospective_dcany_adjusted.csv"))
cat("\n=== Prospective diagnosis table (OR per SD of cPES) ===\n");print(tab,row.names=FALSE)
cat("\nDONE\n")
