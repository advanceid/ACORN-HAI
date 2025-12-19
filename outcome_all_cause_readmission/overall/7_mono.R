# Clear
rm(list = ls())

# Load packages
suppressPackageStartupMessages({
  require(pacman)
  pacman::p_load(dplyr, 
                 magrittr, 
                 crrSC, 
                 cmprsk,
                 ggsurvfit,
                 patchwork,
                 rms,
                 tidyverse,
                 extrafont,
                 purrr,
                 forcat,
                 Cairo)
})

# Load data
df <- readRDS("data/readmission90_data.RData")

# Load fonts
loadfonts()

df$hpd_admreason <- fct_recode(
  df$hpd_admreason,
  "INF" = "Infectious disease",
  "GIT" = "Gastrointestinal disorder",
  "PMD" = "Pulmonary disease",
  "OTH" = "Others"
)

# 
df <- df %>% 
  select("age_new", "sex", 
         "country_region", "country_income",
         "hpd_admreason", "comorbidities_CCI",
         "severity_score_scale",
         "icu_hd_ap", "infection_types",
         "aci_car", "ent_thir", "ent_car",
         "readm_death_time", "readm_event")

x.factors <- model.matrix(~ df$sex + 
                            df$country_region + df$country_income +
                            df$hpd_admreason + df$icu_hd_ap +
                            df$infection_types +
                            df$aci_car + df$ent_thir + df$ent_car)[,-1]


x.factors <- as.matrix(data.frame(x.factors, df$age_new, 
                                  df$comorbidities_CCI,
                                  df$severity_score_scale))

# Set labels
# Revise labels
set_labels <- function(df, labels) {
  for (col in names(labels)) {
    df <- labelled::set_variable_labels(df, !!sym(col) := labels[[col]])
  }
  return(df)
}

labels <- list(
  sex = "Sex",
  age_new = "Age (years)",
  country_region = "Region",
  country_income = "World Bank income status",
  
  hpd_admreason = "Primary admission reason",
  comorbidities_CCI = "Charlson comorbidity index",
  
  severity_score_scale = "Severity score of disease",
  icu_hd_ap = "Admission to ICU/HD at enrollment",
  infection_types = "Infection syndromes",

  aci_car = "CRA", 
  ent_thir = "3GCRE", 
  ent_car = "CRE"
)

df <- set_labels(df, labels)

#
dd <- datadist(df)
options(datadist = "dd")

cox_mod <- cph(Surv(readm_death_time, readm_event == 1) ~ 
                 age_new + sex + country_region + country_income + hpd_admreason +
                 comorbidities_CCI + severity_score_scale +icu_hd_ap + infection_types +
                 aci_car + ent_thir + ent_car,
               data = df, x = TRUE, y = TRUE, 
               surv = TRUE, time.inc = 90)

nom <- nomogram(cox_mod,
                fun=function(x)1/(1+exp(-x)),
                fun.at=c(.001,.01,.05,seq(0,1,by=.1),.95,.99,.999),
                funlabel="Risk of 90-day readmission",
                conf.int=F,
                abbrev=F,
                minlength=1,
                lp=F) 

###
# For death
cox_mod_death <- cph(Surv(readm_death_time, readm_event == 2) ~ 
                       age_new + sex + country_region + country_income + hpd_admreason +
                       comorbidities_CCI + severity_score_scale +icu_hd_ap + infection_types +
                       aci_car + ent_thir + ent_car,
                     data = df, x = TRUE, y = TRUE, 
                     surv = TRUE, time.inc = 90)

nom_death <- nomogram(cox_mod_death,
                      fun=function(x)1/(1+exp(-x)),
                      fun.at=c(.001,.01,.05,seq(0,1,by=.1),.95,.99,.999),
                      funlabel="Risk of 90-day death",
                      conf.int=F,
                      abbrev=F,
                      minlength=1,
                      lp=F) 

#family = "Times"
CairoPDF("output/pdf/nomo_readmission.pdf", width = 23, height = 9)
par(mar = c(1, 0, 3, 4), family = "Times New Roman")
plot(nom, xfrac=.25, 
     total.points.label="Sum of all points", 
     cex.axis = 1.3,
     cex = 1.4,
     force.label = TRUE,
     tcl = -0.3,
     lmgp = 0.2,
     vnames="labels",
     col.grid=gray(c(0.85,0.95)))
dev.off()

CairoPDF("output/pdf/nomo_death.pdf", width = 23, height = 9)
par(mar = c(2.5, 0, 3, 4), family = "Times New Roman")
plot(nom_death, xfrac=.25, 
     total.points.label="Sum of all points", 
     cex.axis = 1.3,
     cex = 1.4,
     force.label = TRUE,
     tcl = -0.3,
     lmgp = 0.2,
     vnames="labels",
     col.grid=gray(c(0.85,0.95)))
dev.off()
