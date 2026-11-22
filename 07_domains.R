################################################################################
# 07_domains.R
# Cumulative psychosocial exposome score (cPES) study - BHRC.
# Descriptive in-sample decomposition of the score by the seven presentation domains.
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
# #1 Domain decomposition of the cPES: contribution to the score + predictive
#    value of each of the 7 psychosocial domains.
################################################################################
suppressWarnings(suppressMessages({library(haven);library(dplyr);library(readr);library(janitor);library(caret);library(tidyr)}))
set.seed(241094)
PATH_BHRC<-"<BHRC_INPUT_DIR>/Grassi-Oliveira_167BHRC_2023_12_09.sav"
PATH_EXPO<-"<PROJECT_ROOT>/Output/BHRCS_final_Exposome.sav"
COEF<-"${cPES_OUT}/01_EXWAS_EN_cPES/POOLED/06_enet_nonzero_coefs.csv"
DOM<-"<DOWNLOADS_DIR>/Copy_cPES_domain_map_with_descriptions.csv"
OUT<-"${cPES_OUT}/_EXTRA_ANALYSES"; if(!dir.exists(OUT))dir.create(OUT,recursive=TRUE)

coef<-read_csv(COEF,show_col_types=FALSE)|>transmute(variable,weight=beta)
dom<-read_csv(DOM,show_col_types=FALSE)|>transmute(variable,domain)
cw<-coef|>left_join(dom,by="variable")
cw$domain[is.na(cw$domain)]<-"Unassigned"
cat("=== vars per domain (with non-zero EN weight) & sum|weight| ===\n")
contrib<-cw|>group_by(domain)|>summarise(n_vars=n(),sum_abs_weight=round(sum(abs(weight)),3))|>arrange(desc(sum_abs_weight))
print(as.data.frame(contrib),row.names=FALSE); write_csv(contrib,file.path(OUT,"01_domain_contribution.csv"))

# analytic sample + standardized exposures (interpretation space, full-sample)
B<-read_sav(PATH_BHRC,encoding="latin1")|>clean_names()
w0<-B|>filter(redcap_event_name=="wave0_arm_1")|>transmute(ident=as.character(ident),cl_tot=as.numeric(cl_tot),age=as.numeric(age))
Expo<-read_sav(PATH_EXPO)|>clean_names();Expo$ident<-as.character(Expo$ident)
vars<-intersect(cw$variable,names(Expo))
dat<-Expo|>select(ident,all_of(vars))|>inner_join(w0,by="ident")|>filter(!is.na(cl_tot))
X<-dat|>select(all_of(vars))
for(v in vars) X[[v]][is.na(X[[v]])]<-median(X[[v]],na.rm=TRUE)
pp<-caret::preProcess(X,method=c("YeoJohnson","center","scale"));Xz<-predict(pp,X)

# per-domain sub-score = weighted sum of that domain's standardized exposures
doms<-unique(cw$domain)
sub<-sapply(doms,function(dm){vs<-cw$variable[cw$domain==dm & cw$variable%in%vars];if(!length(vs))return(rep(0,nrow(Xz)));as.numeric(as.matrix(Xz[,vs,drop=FALSE])%*%cw$weight[match(vs,cw$variable)])})
sub<-as.data.frame(sub);names(sub)<-doms;sub$cl_tot<-dat$cl_tot

# univariate predictive R2 per domain
uni<-lapply(doms,function(dm){r<-suppressWarnings(cor(sub[[dm]],sub$cl_tot,use="complete.obs"));data.frame(domain=dm,cor_with_CBCL=round(r,3),R2=round(r^2,3))})|>bind_rows()|>arrange(desc(R2))
cat("\n=== per-domain sub-score: univariate prediction of CBCL total ===\n");print(as.data.frame(uni),row.names=FALSE)
write_csv(uni,file.path(OUT,"01_domain_univariate_R2.csv"))

# joint model: standardized contribution of each domain
subz<-as.data.frame(scale(sub[,doms]));subz$y<-as.numeric(scale(sub$cl_tot))
mj<-lm(y~.,subz);cj<-summary(mj)$coef
joint<-data.frame(domain=rownames(cj)[-1],std_beta=round(cj[-1,1],3),p=signif(cj[-1,4],3))|>arrange(desc(abs(std_beta)))
cat(sprintf("\n=== joint model of all domains -> CBCL (total model R2=%.3f) ===\n",summary(mj)$r.squared))
print(joint,row.names=FALSE);write_csv(joint,file.path(OUT,"01_domain_joint_model.csv"))
cat("\nDONE\n")
