#'---
#' title: "Time to BV Diagnosis Innate Immune SNP Genotype"
#' author: "Erin dela Cruz"

#' **Program**: Cohort1_TimeToFirstBV_TLR5Deficiency.R
#' 
#' **Purpose**: analyze relationship between innate immune SNP genotype and time to first
#' BV Diagnosis in Cohort 1
#' 
#' -----------------------------------------------------------------------------------------------
#' 
#+ setup-libs, include=FALSE
#' **Prep**: Clear data from memory
rm(list = ls())
# Tidyverse Libs
library(tidyverse)
# For survival analysis
library(survival)
# For epidemiologic analysis
library(epitools)
# ggplot add-ons
library(GGally)
library(ggthemes)
library(ggpubr)
library(cowplot)

#' For use on Windows
fp <- "H:/code/tlr5-deficiency-BV_PUBLICATION/cohort1/"

source(paste0(fp,"src/cohort1-graphfunctions.R"))

#' Time to first diagnosis of BV by Amsel's criteria (clinical)
amsel.dx.tlr5 <- read.csv(paste0(fp, "indata/time-to-first-amsel-tlr5.csv"), header=TRUE, check.names=TRUE, as.is=TRUE) 
table(amsel.dx.tlr5$TLR5_deficient, useNA="ifany")

tlr5.surv.amsel <- surv.analysis(DtaFrame=amsel.dx.tlr5, TollVar="TLR5_deficient", DayVar="days", 
                           OutcomeVar="final.amsel", "TLR5 rs5744168", "CC (Sufficient)\nn=77", "CT/TT (Deficient)\nn=37") +
  theme(legend.position="bottom")

#' Time to first diagnosis of BV by Nugent score (microbiological)
nugent.dx.tlr5 <- read.csv(paste0(fp, "indata/time-to-first-nugent-tlr5.csv"), header=TRUE, check.names=TRUE, as.is=TRUE)
table(nugent.dx.tlr5$TLR5_deficient, useNA="ifany")
tlr5.surv.nugent <- surv.analysis.n(DtaFrame=nugent.dx.tlr5, TollVar="TLR5_deficient", DayVar="days", 
                             OutcomeVar="final.dx",  "TLR5 rs5744168", "CC (Sufficient)\nn=71", "CT/TT (Deficient)\nn=40") +
  theme(legend.position="bottom")

#' Time to first BV with Nugent score 9 or 10 (presence of Mobiluncus morphotypes)
nugent910.dx.tlr5 <- read.csv(paste0(fp, "indata/time-to-first-nugent9_10-tlr5.csv"), header=TRUE, check.names=TRUE, as.is=TRUE)
table(nugent910.dx.tlr5$TLR5_deficient, useNA="ifany")
tlr5.surv.nug910 <- surv.analysis.n(DtaFrame=nugent910.dx.tlr5, TollVar="TLR5_deficient", DayVar="days", 
                             OutcomeVar="final.dx",  "TLR5 rs5744168", "CC (Sufficient)\nn=113", "CT/TT (Deficient)\nn=47") +
  theme(legend.position="bottom") + scale_y_continuous("Proportion without Mobiluncus morphotypes\n(Nugent Score 9-10)", limits=c(0,1)) +
  labs(subtitle="p=0.02 (log-rank test)")

#' Additional statistics for reporting
summary(survfit(formula = Surv(nugent910.dx.tlr5$days, nugent910.dx.tlr5$final.dx) ~ nugent910.dx.tlr5$TLR5_deficient, conf.type = "none"))
fut <- nugent910.dx.tlr5 %>% group_by(TLR5_deficient) %>%
  summarise(tot.fu.time = sum(days))
fut
