setwd("")
install.packages("devtools")
library(devtools)
devtools::install_github("MRCIEU/TwoSampleMR")
library(TwoSampleMR)
devtools::install_github("MRCIEU/MRInstruments")
library(MRInstruments)

ao <- available_outcomes()

#chronotype 
exposure_dat <- read_exposure_data("chronotype_female.txt", sep = "\t",
                                   snp_col = "SNP", beta_col = "beta",
                                   se_col = "se", effect_allele_col = "effect_allele",
                                   other_allele_col = "other_allele", eaf_col = "eaf", pval_col = "pval",
                                   samplesize_col = "N")

outcome_dat <- extract_outcome_data(exposure_dat$SNP, c(1126))

dat_chrono <- harmonise_data(exposure_dat, outcome_dat, action = 2)
summary(dat_chrono$proxy.outcome)

res_chrono <- mr(dat_chrono)
or_chrono <- generate_odds_ratios(res_chrono)
write.csv(or_chrono, "or_chrono_female.csv", row.names=F, quote=F)

outcome_dat <- extract_outcome_data(exposure_dat$SNP, c(1127))

dat_chrono_erp <- harmonise_data(exposure_dat, outcome_dat, action = 2)

res_chrono_erp <- mr(dat_chrono_erp)
or_chrono_erp <- generate_odds_ratios(res_chrono_erp)
write.csv(or_chrono_erp, "or_chrono_erp_female.csv", row.names=F, quote=F)

outcome_dat <- extract_outcome_data(exposure_dat$SNP, c(1128))

dat_chrono_ern <- harmonise_data(exposure_dat, outcome_dat, action = 2)

res_chrono_ern <- mr(dat_chrono_ern)
or_chrono_ern <- generate_odds_ratios(res_chrono_ern)
write.csv(or_chrono_ern, "or_chrono_ern_female.csv", row.names=F, quote=F)

mr_heterogeneity(dat_chrono)
mr_pleiotropy_test(dat_chrono)

dat_chrono$outcome <- as.character(dat_chrono$outcome)
dat_chrono$outcome <- replace(dat_chrono$outcome, dat_chrono$outcome == "Breast cancer (Combined Oncoarray; iCOGS; GWAS meta analysis) || BCAC || 2017", "Breast cancer")
dat_chrono$exposure <- as.character(dat_chrono$exposure)
dat_chrono$exposure <- replace(dat_chrono$exposure, dat_chrono$exposure == "exposure", "Chronotype")

res_chrono <- mr(dat_chrono, method_list = c("mr_ivw", "mr_egger_regression", "mr_weighted_mode", "mr_weighted_median"))
res_chrono$outcome <- as.character(res_chrono$outcome)
res_chrono$outcome <- replace(res_chrono$outcome, res_chrono$outcome == "Breast cancer (Combined Oncoarray; iCOGS; GWAS meta analysis) || BCAC || 2017", "Breast cancer")
res_chrono$exposure <- as.character(res_chrono$exposure)
res_chrono$exposure <- replace(res_chrono$exposure, res_chrono$exposure == "exposure", "Chronotype")

scatter_chrono <- mr_scatter_plot(res_chrono, dat_chrono)
scatter_chrono[[1]]

loo <- mr_leaveoneout(dat_chrono)
loo_chrono <- mr_leaveoneout_plot(loo)
loo_chrono[[1]]

s_chrono <- mr_singlesnp(dat_chrono)
forest_chrono <- mr_forest_plot(s_chrono)
forest_chrono[[1]]

funnel_chrono <- mr_funnel_plot(s_chrono)
funnel_chrono[[1]]

#insomnia
exposure_dat <- read_exposure_data("insomniasnps_ordinal_female.txt", sep = "\t",
                                   snp_col = "SNP", beta_col = "beta",
                                   se_col = "se", effect_allele_col = "effect_allele",
                                   other_allele_col = "other_allele", eaf_col = "EAF", pval_col = "pval",
                                   samplesize_col = "N")

outcome_dat <- extract_outcome_data(exposure_dat$SNP, c(1126))

dat_any <- harmonise_data(exposure_dat, outcome_dat, action = 2)
summary(dat_any$proxy.outcome)

res_any <- mr(dat_any)
or_any <- generate_odds_ratios(res_any)
write.csv(or_any, "or_any_female.csv", row.names=F, quote=F)

outcome_dat <- extract_outcome_data(exposure_dat$SNP, c(1127))

dat_any_erp <- harmonise_data(exposure_dat, outcome_dat, action = 2)

res_any_erp <- mr(dat_any_erp)
or_any_erp <- generate_odds_ratios(res_any_erp)
write.csv(or_any_erp, "or_any_erp_female.csv", row.names=F, quote=F)

outcome_dat <- extract_outcome_data(exposure_dat$SNP, c(1128))

dat_any_ern <- harmonise_data(exposure_dat, outcome_dat, action = 2)

res_any_ern <- mr(dat_any_ern)
or_any_ern <- generate_odds_ratios(res_any_ern)
write.csv(or_any_ern, "or_any_ern_female.csv", row.names=F, quote=F)

mr_heterogeneity(dat_any)
mr_pleiotropy_test(dat_any)

res_any$outcome <- as.character(res_any$outcome)
res_any$outcome <- replace(res_any$outcome, res_any$outcome == "Breast cancer (Combined Oncoarray; iCOGS; GWAS meta analysis) || BCAC || 2017", "Breast cancer")
res_any$exposure <- as.character(res_any$exposure)
res_any$exposure <- replace(res_any$exposure, res_any$exposure == "exposure", "Insomnia")

res_any <- mr(dat_any, method_list = c("mr_ivw", "mr_egger_regression", "mr_weighted_mode", "mr_weighted_median"))
dat_any$outcome <- as.character(dat_any$outcome)
dat_any$outcome <- replace(dat_any$outcome, dat_any$outcome == "Breast cancer (Combined Oncoarray; iCOGS; GWAS meta analysis) || BCAC || 2017", "Breast cancer")
dat_any$exposure <- as.character(dat_any$exposure)
dat_any$exposure <- replace(dat_any$exposure, dat_any$exposure == "exposure", "Insomnia")

scatter_any <- mr_scatter_plot(res_any, dat_any)
scatter_any[[1]]

loo <- mr_leaveoneout(dat_any)
loo_any <- mr_leaveoneout_plot(loo)
loo_any[[1]]

s_any <- mr_singlesnp(dat_any)
forest_any <- mr_forest_plot(s_any)
forest_any[[1]]

funnel_any <- mr_funnel_plot(s_any)
funnel_any[[1]]

#duration 
exposure_dat <- read_exposure_data("sleepdurationhits_snps_female.txt", sep = "\t",
                                   snp_col = "SNP", beta_col = "beta",
                                   se_col = "se", effect_allele_col = "effect_allele",
                                   other_allele_col = "other_allele", eaf_col = "eaf", pval_col = "pval",
                                   samplesize_col = "N")

outcome_dat <- extract_outcome_data(exposure_dat$SNP, c(1126))

dat_duration <- harmonise_data(exposure_dat, outcome_dat, action = 2)
res_duration <- mr(dat_duration)
or_duration <- generate_odds_ratios(res_duration)
write.csv(or_duration, "or_duration_female.csv", row.names=F, quote=F)

outcome_dat <- extract_outcome_data(exposure_dat$SNP, c(1127))

dat_duration_erp <- harmonise_data(exposure_dat, outcome_dat, action = 2)

res_duration_erp <- mr(dat_duration_erp)
or_duration_erp <- generate_odds_ratios(res_duration_erp)
write.csv(or_duration_erp, "or_duration_erp_female.csv", row.names=F, quote=F)

outcome_dat <- extract_outcome_data(exposure_dat$SNP, c(1128))

dat_duration_ern <- harmonise_data(exposure_dat, outcome_dat, action = 2)

res_duration_ern <- mr(dat_duration_ern)
or_duration_ern <- generate_odds_ratios(res_duration_ern)
write.csv(or_duration_ern, "or_duration_ern_female.csv", row.names=F, quote=F)


mr_heterogeneity(dat_duration)
mr_pleiotropy_test(dat_duration)

res_duration$outcome <- as.character(res_duration$outcome)
res_duration$outcome <- replace(res_duration$outcome, res_duration$outcome == "Breast cancer (Combined Oncoarray; iCOGS; GWAS meta analysis) || BCAC || 2017", "Breast cancer")
res_duration$exposure <- as.character(res_duration$exposure)
res_duration$exposure <- replace(res_duration$exposure, res_duration$exposure == "exposure", "Sleep duration")

res_duration <- mr(dat_duration, method_list = c("mr_ivw", "mr_egger_regression", "mr_weighted_mode", "mr_weighted_median"))
dat_duration$outcome <- as.character(dat_duration$outcome)
dat_duration$outcome <- replace(dat_duration$outcome, dat_duration$outcome == "Breast cancer (Combined Oncoarray; iCOGS; GWAS meta analysis) || BCAC || 2017", "Breast cancer")
dat_duration$exposure <- as.character(dat_duration$exposure)
dat_duration$exposure <- replace(dat_duration$exposure, dat_duration$exposure == "exposure", "Sleep duration")

scatter_duration <- mr_scatter_plot(res_duration, dat_duration)
scatter_duration[[1]]

loo <- mr_leaveoneout(dat_duration)
loo_duration <- mr_leaveoneout_plot(loo)
loo_duration[[1]]

s_duration <- mr_singlesnp(dat_duration)
forest_duration <- mr_forest_plot(s_duration)
forest_duration[[1]]

funnel_duration <- mr_funnel_plot(s_duration)
funnel_duration[[1]]

#MR-Raps 
devtools::install_github('qingyuanzhao/mr.raps')
mr_raps <- mr(dat_chrono, method_list = c("mr_raps"))
or_mr_raps <- generate_odds_ratios(mr_raps)

mr_raps <- mr(dat_duration, method_list = c("mr_raps"))
or_mr_raps <- generate_odds_ratios(mr_raps)

mr_raps <- mr(dat_any, method_list = c("mr_raps"))
or_mr_raps <- generate_odds_ratios(mr_raps)

#Radial plots 
devtools::install_github("WSpiller/RadialMR")
library(RadialMR)

#chronotype 
res_single <- mr_singlesnp(dat_chrono)

dat <- dat_chrono[dat_chrono$SNP%in%res_single$SNP,]

raddat <- format_radial(dat$beta.exposure, dat$beta.outcome, dat$se.exposure, dat$se.outcome, dat$SNP)
ivwrad <- ivw_radial(raddat, alpha=0.05/305, weights=3)
dim(ivwrad$outliers)[1] 
#36 outliers at 0.05, 6 at bonf  

eggrad <- egger_radial(raddat, alpha=0.05/305, weights=3)
eggrad$coef 
dim(eggrad$outliers)[1] 
#36 outliers at 0.05, 6 at bonf 

#plot_radial(ivwrad, TRUE, FALSE, TRUE)
plot_radial(c(ivwrad,eggrad), TRUE, FALSE, TRUE)

ivwrad$qstatistic 
ivwrad$sortoutliers <- ivwrad$outliers[order(ivwrad$outliers$p.value),]
ivwrad$sortoutliers$Qsum <- cumsum(ivwrad$sortoutliers$Q_statistic)
ivwrad$sortoutliers$Qdif <- ivwrad$sortoutliers$Qsum - ivwrad$qstatistic
write.csv(ivwrad$sortoutliers, "chrono_outliers.csv", row.names=F, quote=F)

eggrad$qstatistic 
eggrad$sortoutliers <- eggrad$outliers[order(eggrad$outliers$p.value),]
eggrad$sortoutliers$Qsum <- cumsum(eggrad$sortoutliers$Q_statistic)
eggrad$sortoutliers$Qdif <- eggrad$sortoutliers$Qsum - eggrad$qstatistic
write.csv(eggrad$sortoutliers, "chrono_outliers_egger.csv", row.names=F, quote=F)

#Remove top outliers 
dat2 <- dat_chrono[!dat_chrono$SNP %in% ivwrad$outliers$SNP,]
mr_chrono2 <- mr(dat2)
or_chrono2 <- generate_odds_ratios(mr_chrono2)
mr_heterogeneity(dat2)
mr_pleiotropy_test(dat2)
write.csv(or_chrono2, "or_chrono_female2.csv", row.names=F, quote=F)

#duration
res_single <- mr_singlesnp(dat_duration)
dat <- dat_duration[dat_duration$SNP%in%res_single$SNP,]

raddat <- format_radial(dat$beta.exposure, dat$beta.outcome, dat$se.exposure, dat$se.outcome, dat$SNP)
ivwrad <- ivw_radial(raddat, alpha=0.05/82, weights=3)
dim(ivwrad$outliers)[1] 
#3 outliers at bonf   

eggrad <- egger_radial(raddat, alpha=0.05/82, weights=3)
eggrad$coef 
dim(eggrad$outliers)[1] 
#3 outliers at bonf 

#plot_radial(ivwrad, TRUE, FALSE, TRUE)
plot_radial(c(ivwrad,eggrad), TRUE, FALSE, TRUE)

ivwrad$qstatistic 
ivwrad$sortoutliers <- ivwrad$outliers[order(ivwrad$outliers$p.value),]
ivwrad$sortoutliers$Qsum <- cumsum(ivwrad$sortoutliers$Q_statistic)
ivwrad$sortoutliers$Qdif <- ivwrad$sortoutliers$Qsum - ivwrad$qstatistic
write.csv(ivwrad$sortoutliers, "duration_outliers.csv", row.names=F, quote=F)

eggrad$qstatistic 
eggrad$sortoutliers <- eggrad$outliers[order(eggrad$outliers$p.value),]
eggrad$sortoutliers$Qsum <- cumsum(eggrad$sortoutliers$Q_statistic)
eggrad$sortoutliers$Qdif <- eggrad$sortoutliers$Qsum - eggrad$qstatistic
write.csv(eggrad$sortoutliers, "duration_outliers_egger.csv", row.names=F, quote=F)

#Remove top outliers 
dat2 <- dat_duration[!dat_duration$SNP %in% ivwrad$outliers$SNP,]
mr_duration2 <- mr(dat2)
or_duration2 <- generate_odds_ratios(mr_duration2)
mr_heterogeneity(dat2)
mr_pleiotropy_test(dat2)
write.csv(or_duration2, "or_duration_female2.csv", row.names=F, quote=F)

#insomnia
res_single <- mr_singlesnp(dat_any)
dat <- dat_any[dat_any$SNP%in%res_single$SNP,]

raddat <- format_radial(dat$beta.exposure, dat$beta.outcome, dat$se.exposure, dat$se.outcome, dat$SNP)
ivwrad <- ivw_radial(raddat, alpha=0.05/50, weights=3)
dim(ivwrad$outliers)[1] 
#2 outliers at bonf   

eggrad <- egger_radial(raddat, alpha=0.05/50, weights=3)
eggrad$coef 
dim(eggrad$outliers)[1] 
#2 outliers at bonf 

#plot_radial(ivwrad, TRUE, FALSE, TRUE)
plot_radial(c(ivwrad,eggrad), TRUE, FALSE, TRUE)

ivwrad$qstatistic 
ivwrad$sortoutliers <- ivwrad$outliers[order(ivwrad$outliers$p.value),]
ivwrad$sortoutliers$Qsum <- cumsum(ivwrad$sortoutliers$Q_statistic)
ivwrad$sortoutliers$Qdif <- ivwrad$sortoutliers$Qsum - ivwrad$qstatistic
write.csv(ivwrad$sortoutliers, "any_outliers.csv", row.names=F, quote=F)

eggrad$qstatistic 
eggrad$sortoutliers <- eggrad$outliers[order(eggrad$outliers$p.value),]
eggrad$sortoutliers$Qsum <- cumsum(eggrad$sortoutliers$Q_statistic)
eggrad$sortoutliers$Qdif <- eggrad$sortoutliers$Qsum - eggrad$qstatistic
write.csv(eggrad$sortoutliers, "any_outliers_egger.csv", row.names=F, quote=F)

#Remove top outliers 
dat2 <- dat_any[!dat_any$SNP %in% ivwrad$outliers$SNP,]
mr_any2 <- mr(dat2)
or_any2 <- generate_odds_ratios(mr_any2)
mr_heterogeneity(any2)
mr_pleiotropy_test(any2)
write.csv(or_any2, "or_insomnia_female2.csv", row.names=F, quote=F)

#MR presso 
devtools::install_github("rondolab/MR-PRESSO")
library(MRPRESSO)

mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = dat_chrono, NbDistribution = 10000,  SignifThreshold = 0.05)

mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = dat_duration, NbDistribution = 10000,  SignifThreshold = 0.05)

#winner's curse analysis 
#chronotype 
exposure_dat <- read_exposure_data("chronotype_female2.txt", sep = "\t",
                                   snp_col = "SNP", beta_col = "beta",
                                   se_col = "se", effect_allele_col = "effect_allele",
                                   other_allele_col = "other_allele", eaf_col = "eaf", pval_col = "pval",
                                   samplesize_col = "N")

outcome_dat <- extract_outcome_data(exposure_dat$SNP, c(1126))

dat_chrono3 <- harmonise_data(exposure_dat, outcome_dat, action = 2)
#dat <- harmonise_data(exposure_dat, outcome_dat, action = 1)
#dat_chrono <- clump_data(dat_chrono)

res_chrono3 <- mr(dat_chrono3)
or_chrono3 <- generate_odds_ratios(res_chrono3)
write.csv(or_chrono3, "or_chrono_female_wc.csv", row.names=F, quote=F)


#binary sleep 
exposure_dat <- read_exposure_data("shortduration_female.txt", sep = "\t",
                                   snp_col = "SNP", beta_col = "logOR",
                                   se_col = "se", effect_allele_col = "effect_allele",
                                   other_allele_col = "other_allele", eaf_col = "eaf", pval_col = "pval",
                                   samplesize_col = "N")

outcome_dat <- extract_outcome_data(exposure_dat$SNP, c(1126))

dat_short <- harmonise_data(exposure_dat, outcome_dat, action = 2)
res_short <- mr(dat_short)
or_short <- generate_odds_ratios(res_short)
write.csv(or_short, "or_short_female.csv", row.names=F, quote=F)

exposure_dat <- read_exposure_data("longduration_female.txt", sep = "\t",
                                   snp_col = "SNP", beta_col = "logOR",
                                   se_col = "se", effect_allele_col = "effect_allele",
                                   other_allele_col = "other_allele", eaf_col = "eaf", pval_col = "pval",
                                   samplesize_col = "N")

outcome_dat <- extract_outcome_data(exposure_dat$SNP, c(1126))

dat_long <- harmonise_data(exposure_dat, outcome_dat, action = 2)
res_long <- mr(dat_long)
or_long <- generate_odds_ratios(res_long)
write.csv(or_long, "or_long_female.csv", row.names=F, quote=F)


#Objective measures 
#L5 timing
exposure_dat <- read_exposure_data("L5_timing.csv", sep = ",", phenotype_col ="Exposure",
                                     snp_col = "SNP", beta_col = "beta",
                                     se_col = "se", effect_allele_col = "effect_allel",
                                     other_allele_col = "other_allele", eaf_col = "eaf", pval_col = "p-val")

exposure_dat$beta.exposure <- exposure_dat$beta.exposure*-1 
outcome_dat <- extract_outcome_data(exposure_dat$SNP, c(1126))

dat_L5 <- harmonise_data(exposure_dat, outcome_dat, action = 2)
#dat_L5 <- harmonise_data(exposure_dat, outcome_dat, action = 1)
res_L5 <- mr(dat_L5)
or_L5 <- generate_odds_ratios(res_L5)
write.csv(or_L5, "or_L5.csv", row.names=F, quote=F)
mr_heterogeneity(dat_L5)
mr_pleiotropy_test(dat_L5)

res_L5$outcome <- as.character(res_L5$outcome)
res_L5$outcome <- replace(res_L5$outcome, res_L5$outcome == "Breast cancer (Combined Oncoarray; iCOGS; GWAS meta analysis) || id:1126", "Breast cancer")

res_L5 <- mr(dat_L5, method_list = c("mr_ivw", "mr_egger_regression", "mr_weighted_mode", "mr_weighted_median"))
dat_L5$outcome <- as.character(dat_L5$outcome)
dat_L5$outcome <- replace(dat_L5$outcome, dat_L5$outcome == "Breast cancer (Combined Oncoarray; iCOGS; GWAS meta analysis) || id:1126", "Breast cancer")

scatter_L5 <- mr_scatter_plot(res_L5, dat_L5)
scatter_L5[[1]]

loo <- mr_leaveoneout(dat_L5)
loo_L5 <- mr_leaveoneout_plot(loo)
loo_L5[[1]]

s_L5 <- mr_singlesnp(dat_L5)
forest_L5 <- mr_forest_plot(s_L5)
forest_L5[[1]]

funnel_L5 <- mr_funnel_plot(s_L5)
funnel_L5[[1]]

#MR-Raps 
devtools::install_github('qingyuanzhao/mr.raps')
mr_raps <- mr(dat_L5, method_list = c("mr_raps"))
or_mr_raps <- generate_odds_ratios(mr_raps)

#Sleep duration 
  exposure_dat <- read_exposure_data("sleep_duration.csv", sep = ",", phenotype_col ="Exposure", snp_col = "SNP", beta_col = "beta", se_col = "se", effect_allele_col = "effect_allel", other_allele_col = "other_allele", eaf_col = "eaf", pval_col = "p-val")

outcome_dat <- extract_outcome_data(exposure_dat$SNP, c(1126))

dat_sd <- harmonise_data(exposure_dat, outcome_dat, action = 2)
#dat_sd <- harmonise_data(exposure_dat, outcome_dat, action = 1)
res_sd <- mr(dat_sd)
or_sd <- generate_odds_ratios(res_sd)
write.csv(or_sd, "or_sd.csv", row.names=F, quote=F)
mr_heterogeneity(dat_sd)
mr_pleiotropy_test(dat_sd)

res_sd$outcome <- as.character(res_sd$outcome)
res_sd$outcome <- replace(res_sd$outcome, res_sd$outcome == "Breast cancer (Combined Oncoarray; iCOGS; GWAS meta analysis) || id:1126", "Breast cancer")

res_sd <- mr(dat_sd, method_list = c("mr_ivw", "mr_egger_regression", "mr_weighted_mode", "mr_weighted_median"))
dat_sd$outcome <- as.character(dat_sd$outcome)
dat_sd$outcome <- replace(dat_sd$outcome, dat_sd$outcome == "Breast cancer (Combined Oncoarray; iCOGS; GWAS meta analysis) || id:1126", "Breast cancer")


scatter_sd <- mr_scatter_plot(res_sd, dat_sd)
scatter_sd[[1]]

loo <- mr_leaveoneout(dat_sd)
loo_sd <- mr_leaveoneout_plot(loo)
loo_sd[[1]]

s_sd <- mr_singlesnp(dat_sd)
forest_sd <- mr_forest_plot(s_sd)
forest_sd[[1]]

funnel_sd <- mr_funnel_plot(s_sd)
funnel_sd[[1]]

mr_raps2 <- mr(dat_sd, method_list = c("mr_raps"))
or_mr_raps2 <- generate_odds_ratios(mr_raps2)

#Sleep fragmentation
  exposure_dat <- read_exposure_data("sleep_fragmentation.csv", sep = ",", phenotype_col ="Exposure",
                                     snp_col = "SNP", beta_col = "beta",
                                     se_col = "se", effect_allele_col = "effect_allel",
                                     other_allele_col = "other_allele", eaf_col = "eaf", pval_col = "p-val")

outcome_dat <- extract_outcome_data(exposure_dat$SNP, c(1126))

dat_sf <- harmonise_data(exposure_dat, outcome_dat, action = 2)
#dat_sf <- harmonise_data(exposure_dat, outcome_dat, action = 1)
res_sf <- mr(dat_sf)
or_sf <- generate_odds_ratios(res_sf)
write.csv(or_sf, "or_sf.csv", row.names=F, quote=F)
mr_heterogeneity(dat_sf)
mr_pleiotropy_test(dat_sf)

res_sf$outcome <- as.character(res_sf$outcome)
res_sf$outcome <- replace(res_sf$outcome, res_sf$outcome == "Breast cancer (Combined Oncoarray; iCOGS; GWAS meta analysis) || id:1126", "Breast cancer")

res_sf <- mr(dat_sf, method_list = c("mr_ivw", "mr_egger_regression", "mr_weighted_mode", "mr_weighted_median"))
dat_sf$outcome <- as.character(dat_sf$outcome)
dat_sf$outcome <- replace(dat_sf$outcome, dat_sf$outcome == "Breast cancer (Combined Oncoarray; iCOGS; GWAS meta analysis) || id:1126", "Breast cancer")

scatter_sf <- mr_scatter_plot(res_sf, dat_sf)
scatter_sf[[1]]

loo <- mr_leaveoneout(dat_sf)
loo_sf <- mr_leaveoneout_plot(loo)
loo_sf[[1]]

s_sf <- mr_singlesnp(dat_sf)
forest_sf <- mr_forest_plot(s_sf)
forest_sf[[1]]

funnel_sf <- mr_funnel_plot(s_sf)
funnel_sf[[1]]

mr_raps3 <- mr(dat_sf, method_list = c("mr_raps"))
or_mr_raps3 <- generate_odds_ratios(mr_raps3)


