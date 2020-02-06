#'---
#' title: "Cohort 2: TLR5 Deficiency and BV"
#' author: "Erin dela Cruz"
#' -----------------------------------------------------------------------------------------------
#' 
#' **Program**: Cohort2_TLR5_deficiency_analysis.R
#' 
#' **Purpose**: Explore relationship between TLR5 deficiency and BV/BVAB colonization in Cohort #2
#' 
#' -----------------------------------------------------------------------------------------------
#' 
#+ setup-libs, include=FALSE
#' **Prep**: Load necessary libraries and clear data from memory
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
# for statistical analysis
library(lme4)
library(MatchIt)

#' For use on Windows
fp <- "H:/code/tlr5-deficiency-BV_PUBLICATION/cohort2/"

source(paste0(fp,"src/cohort2_graphfunctions.R"))

#' -----------------------------------------------------------------------------------------------
rsnumber <- "TLR5 rs5744168"
tlrvariable <- "TLR5_deficient"
wildtype.bl <- "CC (Sufficient)\nn=103"
wildtype.m <- "CC (Sufficient)\nn=32"
deficient <- "CT/TT (Deficient)\nn=8"
#' -----------------------------------------------------------------------------------------------

#' ***********************************************************************************************
#' Section 1: Load Data
#' ***********************************************************************************************

snp.bl <- read.csv(paste0(fp, "indata/cohort2-bl-data.csv"), header=TRUE, check.names=TRUE, as.is=TRUE) %>%
  mutate(TLR5_deficient=as.factor(TLR5_deficient))


#' ***********************************************************************************************
#' Section 2: Match TLR5 deficient with sufficient controls in 4:1 ratio
#' Used for colonization analyses
#' ***********************************************************************************************

#' Match TLR5 deficient individuals with TLR5 sufficient controls
tlr5.match <- snp.bl %>%
  filter(is.na(nugentdx)==FALSE &
           is.na(TLR5_rs5744168)==FALSE &
           is.na(age)==FALSE &
           is.na(contraception_horm)==FALSE &
           is.na(white)==FALSE &
           is.na(douch_ever)==FALSE) %>% 
  mutate(TLR5_deficient = (TLR5_rs5744168=="C/T")) %>%
  select(pubID, nugentdx, TLR5_deficient, age, contraception_horm, white, douch_ever)
set.seed(918)
mi <- matchit(TLR5_deficient~age + contraception_horm + white + douch_ever, 
                   data=tlr5.match,
                   ratio=4,
              method = 'nearest')
summary(mi)
tlr5.matched <- match.data(mi)

#' Subset genotype, BL data to only matched controls
snp.matched <- semi_join(snp.bl, tlr5.matched, by="pubID")




#' ***********************************************************************************************
#' Section 3.1: Import qPCR data from first set (BVAB1, Gvag, Lcrisp)
#' ***********************************************************************************************
#' Load qPCR data - BVAB1, Gvag, Lcrisp
qpcr.snp <- read.csv(paste0(fp, "indata/cohort2-qPCR.csv"), header=TRUE, check.names=TRUE, as.is=TRUE) %>%
  right_join(., snp.matched, by="pubID")  %>%
  filter(time=="B")

#' Replace with full assay name
qpcr.snp$assay <- as.character(qpcr.snp$assay)
qpcr.snp$assay[qpcr.snp$assay=="bvab1"] <- "BVAB1"
qpcr.snp$assay[qpcr.snp$assay=="gvag"] <- "G. vaginalis"
qpcr.snp$assay[qpcr.snp$assay=="lcrisp"] <- "L. crispatus"

#' Create separate datasets for each qPCR assay
bvab1.snp <- as.data.frame(qpcr.snp[which(qpcr.snp$assay == "BVAB1"),])
gvag.snp <- as.data.frame(qpcr.snp[which(qpcr.snp$assay == "G. vaginalis"),])
lcrisp.snp <- as.data.frame(qpcr.snp[which(qpcr.snp$assay == "L. crispatus"),])

#' ***********************************************************************************************
#' Section 3.2: Import Mobiluncus qPCR data (assay developed at later date)
#' ***********************************************************************************************
#' Import Mobiluncus qPCR data
mobis <- read.csv(paste0(fp,"indata/2018-02-11_McMmCombo.csv"), 
                  header=TRUE, stringsAsFactors=FALSE, strip.white=TRUE) %>%
  mutate(copies_per_rxn = ifelse(copies_per_rxn<lod | is.na(copies_per_rxn), (lod/2), copies_per_rxn)) %>%
  mutate(copies.per.uL = copies_per_rxn/2) %>%
  mutate(copies_per_swab = copies.per.uL*150) %>% 
  select(pubID, copies_per_swab, assay) %>%
  mutate(log.copies.per.swab = log10(copies_per_swab))

mobis.bvab1 <- mobis %>% select(pubID, log.copies.per.swab, assay) %>%
  spread(assay, log.copies.per.swab) %>%
  inner_join(., bvab1.snp, by="pubID") %>%
  rename(BVAB1 = log.copies.per.swab) %>%
  select(BVAB1, M.curtisii, M.mulieris, TLR5_deficient)
flag.pairs <- ggpairs(mobis.bvab1, aes(colour=TLR5_deficient))

#' Create seprate datasets for Mmulieris, Mcurtisii
m.mulie <- mobis %>% filter(assay=="M.mulieris") %>% left_join(., snp.bl, by="pubID")
m.curt <- mobis %>% filter(assay=="M.curtisii") %>% left_join(., snp.bl, by="pubID")





#' ****************************************************************************************************
#' Section 4.1: Graph Relationships between TLR5 Deficiency and BV
#' ****************************************************************************************************
snp.bl <- as.data.frame(snp.bl)
diagnosis <- bvriskgraph(snp.bl, tlrvariable, "nugentdx", "amsel", rsnumber, wildtype.bl, deficient)

hi.nug.summary <- snp.bl %>% group_by(TLR5_deficient) %>%
  summarise(pct = mean(hi.nugent, na.rm=TRUE),
             sd = sd(hi.nugent),
             n= n()) %>%
  mutate(se = sd/sqrt(n)) %>%
  na.omit() %>%
  mutate(dx="Nugent 9-10")
model <- glm(hi.nugent~TLR5_deficient,
             family=binomial(link='logit'),
             data=snp.bl)
summary(model)
label <- substitute(paste("OR=", estimate, ", CI=[", lo, ", ", hi, "] (p=", pvalue, ")"),
                    list(estimate = signif(exp(summary(model)$coefficient[2,1]), 3),
                         lo = signif(exp(summary(model)$coefficient[2,1] - 1.96*summary(model)$coefficient[2,2]), 3),
                         hi = signif(exp(summary(model)$coefficient[2,1] + 1.96*summary(model)$coefficient[2,2]), 3),
                         pvalue = signif(summary(model)$coefficient[2,4], 3)))
nugent.9.10 <- ggplot(hi.nug.summary, aes(x=TLR5_deficient, y=pct, color=factor(TLR5_deficient))) +
  geom_pointrange(aes(ymin=pct-1.96*se, ymax=pct+1.96*se),
                  size=.3, 
                  position=position_dodge(.9)) +
  labs(subtitle = label, y="Proportion with Mobiluncus morphotypes\n(Nugent Score 9-10)", x="") + 
  scale_color_manual(labels=c(paste(wildtype.bl), paste(deficient)), values=c("#f1a340", "#998ec3")) +
  scale_x_discrete(labels=c("TLR5 rs5744168\nCC (Sufficient)\nn=103", "TLR5 rs5744168\nCT/TT (Deficient)\nn=8")) +
  theme(legend.position = "none")

t <- table(snp.bl$TLR5_deficient, snp.bl$hi.nugent)
riskratio(t)

#' ****************************************************************************************************
#' Section 4.2: Graph Relationships between BVAB Colonization and TLR5 Deficiency
#' ****************************************************************************************************
bvab1 <- mini.smear(bvab1.snp, tlrvariable, HypothesisDirection = "two.sided",
                    rsnumber, wildtype.m, deficient, "BVAB1")
gvag <- mini.smear(gvag.snp, tlrvariable, HypothesisDirection = "two.sided",
                   rsnumber, wildtype.m, deficient, "Gardnerella vaginalis")
lcrisp <- mini.smear(lcrisp.snp, tlrvariable, HypothesisDirection = "two.sided", 
                     rsnumber, wildtype.m, deficient, "Lactobacillus crispatus")
mc <- mini.smear(m.curt, tlrvariable, HypothesisDirection = "two.sided", 
                 rsnumber, wildtype.m, deficient, "Mobiluncus curtisii")
mm <- mini.smear(m.mulie, tlrvariable, HypothesisDirection = "two.sided", 
                 rsnumber, wildtype.m, deficient, "Mobiluncus mulieris")

#' ****************************************************************************************************
#' Section 4.3: Correlation between IL-8 [flagellated BVAB]
#' TLR5 sufficient subset with IL-8 measurement by ELISA
#' ****************************************************************************************************
#' Plot cytokine data, test for statistically significant differences
m.mulie.merge <- m.mulie %>% select(pubID, copies_per_swab) %>%
  rename(m.mulie = copies_per_swab)
m.curt.merge <- m.curt %>% select(pubID, copies_per_swab) %>%
  rename(m.curt = copies_per_swab)
FlaAs <- left_join(bvab1.snp, m.mulie.merge, by = "pubID") %>%
  left_join(., m.curt.merge, by="pubID")

il8.cohort <- read.csv(paste0(fp,"indata/2018-05-02_R56-TLR5Cohort_IL8.csv"), 
                header=TRUE, stringsAsFactors=FALSE, strip.white=TRUE) %>%
  group_by(pubID) %>%
  summarise(mean.il8 = mean(il8.final))

#' Create separate/wide format datasets for graphing
FlaAs.il8 <- left_join(FlaAs, il8.cohort, by = "pubID")
gv.snp.il8 <- left_join(gvag.snp, il8.cohort, by = "pubID")

#' Analysis restricted to TLR5 sufficient individuals, as TLR5 deficient would not be mounting TLR5-dependent
#' IL-8 inflammatory response to flagellated BVAB
tlr5.suf <- FlaAs.il8 %>% filter(TLR5_deficient==0)
tlr5.suf.il8 <- tlr5.suf %>% filter(is.na(mean.il8)==FALSE)
gv.il8 <- gv.snp.il8 %>% filter(TLR5_deficient==0 & is.na(mean.il8)==FALSE)

#' Graph: correlation between IL8 and [Mobiluncus mulieris]
label <- substitute(paste("Correlation (Spearman) = ", corr, ", p = ", p),
                    list(corr = signif(cor.test(tlr5.suf.il8$m.mulie, tlr5.suf.il8$mean.il8, method = c("spearman"))$estimate, 3),
                         p = signif(cor.test(tlr5.suf.il8$m.mulie, tlr5.suf.il8$mean.il8, method = c("spearman"))$p.value, 3)))
mm.il8 <- ggplot(tlr5.suf.il8, aes(x=m.mulie, y=mean.il8, color=factor(TLR5_deficient), shape = factor(nugentdx))) +
  guides(color=FALSE) +
  geom_point(size=4) + scale_x_log10("Mobiluncus mulieris\n16S rRNA gene copies per swab") + scale_y_log10("IL-8 (pg/mL)") +
  scale_color_manual("TLR5 rs5744168", labels = c("CC (Sufficient)\nn=31", "CT/TT (Deficient)\nn=8"), values=c("#f1a340", "#998ec3")) +
  scale_shape_manual("BV (Nugent)", labels = c("Negative", "Positive"), values= c(5, 16)) + 
  labs(subtitle = label) +
  theme(legend.position = "bottom", legend.box = "vertical", plot.subtitle=element_text(color="gray36"))

#' Graph: correlation between IL8 and [Mobiluncus curtisii]
label <- substitute(paste("Correlation (Spearman) = ", corr, ", p = ", p),
                    list(corr = signif(cor.test(tlr5.suf.il8$m.curt, tlr5.suf.il8$mean.il8, method = c("spearman"))$estimate, 3),
                         p = signif(cor.test(tlr5.suf.il8$m.curt, tlr5.suf.il8$mean.il8, method = c("spearman"))$p.value, 3)))
mc.il8 <- ggplot(tlr5.suf.il8, aes(x=m.curt, y=mean.il8, color=factor(TLR5_deficient), shape = factor(nugentdx))) +
  guides(color=FALSE) +
  geom_point(size=4) + scale_x_log10("Mobiluncus curtisii\n16S rRNA gene copies per swab") + scale_y_log10("IL-8 (pg/mL)") +
  scale_color_manual("TLR5 rs5744168", labels = c("CC (Sufficient)\nn=31", "CT/TT (Deficient)\nn=8"), values=c("#f1a340", "#998ec3")) +
  scale_shape_manual("BV (Nugent)", labels = c("Negative", "Positive"), values= c(5, 16)) + 
  labs(subtitle = label) +
  theme(legend.position = "bottom", legend.box = "vertical", plot.subtitle=element_text(color="gray36"))

#' Graph: correlation between IL8 and [BVAB1]
label <- substitute(paste("Correlation (Spearman) = ", corr, ", p = ", p),
                    list(corr = signif(cor.test(tlr5.suf.il8$copies_per_swab, tlr5.suf.il8$mean.il8, method = c("spearman"))$estimate, 3),
                         p = signif(cor.test(tlr5.suf.il8$copies_per_swab, tlr5.suf.il8$mean.il8, method = c("spearman"))$p.value, 3)))
bvab1.il8 <- ggplot(tlr5.suf.il8, aes(x=copies_per_swab, y=mean.il8, color=factor(TLR5_deficient), shape = factor(nugentdx))) +
  guides(color=FALSE) +
  geom_point(size=4) + scale_x_log10("BVAB1\n16S rRNA gene copies per swab") + scale_y_log10("IL-8 (pg/mL)") +
  scale_color_manual("TLR5 rs5744168", labels = c("CC (Sufficient)\nn=31", "CT/TT (Deficient)\nn=8"), values=c("#f1a340", "#998ec3")) +
  scale_shape_manual("BV (Nugent)", labels = c("Negative", "Positive"), values= c(5, 16)) + 
  labs(subtitle = label) +
  theme(legend.position = "bottom", legend.box = "vertical", plot.subtitle=element_text(color="gray36"))

#' Graph: correlation betwee IL8 and Gvag
label <- substitute(paste("Correlation (Spearman) = ", corr, ", p = ", p),
                    list(corr = signif(cor.test(gv.il8$copies_per_swab, gv.il8$mean.il8, method = c("spearman"))$estimate, 3),
                         p = signif(cor.test(gv.il8$copies_per_swab, gv.il8$mean.il8, method = c("spearman"))$p.value, 3)))
gv.il8 <- ggplot(gv.il8, aes(x=copies_per_swab, y=mean.il8, color=factor(TLR5_deficient), shape = factor(nugentdx))) +
  guides(color=FALSE) +
  geom_point(size=4) + scale_x_log10("Gardnerella vaginalis\n16S rRNA gene copies per swab") + scale_y_log10("IL-8 (pg/mL)") +
  scale_color_manual("TLR5 rs5744168", labels = c("CC (Sufficient)\nn=31", "CT/TT (Deficient)\nn=8"), values=c("#f1a340", "#998ec3")) +
  scale_shape_manual("BV (Nugent)", labels = c("Negative", "Positive"), values= c(5, 16)) + 
  labs(subtitle = label) +
  theme(legend.position = "bottom", legend.box = "vertical", plot.subtitle=element_text(color="gray36"))

#' ****************************************************************************************************
#' Section 4.3: Correlations between BVAB1, Mm, Mc
#' TLR5 sufficient subset
#' ****************************************************************************************************
label <- substitute(paste("Correlation (Spearman) = ", corr, ", p = ", p),
                    list(corr = signif(cor.test(tlr5.suf$m.curt, tlr5.suf$copies_per_swab, method = c("spearman"))$estimate, 3),
                         p = signif(cor.test(tlr5.suf$m.curt, tlr5.suf$copies_per_swab, method = c("spearman"))$p.value, 3)))
bvab1.mc.corr <- ggplot(subset(tlr5.suf, TLR5_deficient==0), aes(x=m.curt, y=copies_per_swab, color=factor(TLR5_deficient), shape = factor(nugentdx))) +
  guides(color=FALSE) +
  geom_point(size=4) + scale_x_log10("Mobiluncus curtisii\n16S rRNA gene copies per swab") + scale_y_log10("BVAB1\n16S rRNA gene copies per swab") +
  scale_color_manual("BV (Nugent)", labels = c("Negative", "Positive"), values=c("#f1a340", "#998ec3")) +
  scale_shape_manual("BV (Nugent)", labels = c("Negative", "Positive"), values= c(5, 16)) + 
  labs(subtitle = label) +
  theme(legend.position = "bottom", legend.box = "vertical", plot.subtitle=element_text(color="gray36"))


label <- substitute(paste("Correlation (Spearman) = ", corr, ", p = ", p),
                    list(corr = signif(cor.test(tlr5.suf$m.mulie, tlr5.suf$copies_per_swab, method = c("spearman"))$estimate, 3),
                         p = signif(cor.test(tlr5.suf$m.mulie, tlr5.suf$copies_per_swab, method = c("spearman"))$p.value, 3)))
bvab1.mm.corr <- ggplot(subset(tlr5.suf, TLR5_deficient==0), aes(x=m.mulie, y=copies_per_swab, color=factor(TLR5_deficient), shape = factor(nugentdx))) +
  guides(color=FALSE) +
  geom_point(size=4) + scale_x_log10("Mobiluncus mulieris\n16S rRNA gene copies per swab") + scale_y_log10("BVAB1\n16S rRNA gene copies per swab") +
  scale_color_manual("BV (Nugent)", labels = c("Negative", "Positive"), values=c("#f1a340", "#998ec3")) +
  scale_shape_manual("BV (Nugent)", labels = c("Negative", "Positive"), values= c(5, 16)) + 
  labs(subtitle = label) +
  theme(legend.position = "bottom", legend.box = "vertical", plot.subtitle=element_text(color="gray36"))


label <- substitute(paste("Correlation (Spearman) = ", corr, ", p = ", p),
                    list(corr = signif(cor.test(tlr5.suf$m.mulie, tlr5.suf$m.curt, method = c("spearman"))$estimate, 3),
                         p = signif(cor.test(tlr5.suf$m.mulie, tlr5.suf$m.curt, method = c("spearman"))$p.value, 3)))
mc.mm.corr <- ggplot(subset(tlr5.suf, TLR5_deficient==0), aes(x=m.mulie, y=m.curt, color=factor(TLR5_deficient), shape = factor(nugentdx))) +
  guides(color=FALSE) +
  geom_point(size=4) + scale_x_log10("Mobiluncus mulieris\n16S rRNA gene copies per swab") + scale_y_log10("Mobiluncus curtisii\n16S rRNA gene copies per swab") +
  scale_color_manual("BV (Nugent)", labels = c("Negative", "Positive"), values=c("#f1a340", "#998ec3")) +
  scale_shape_manual("BV (Nugent)", labels = c("Negative", "Positive"), values= c(5, 16)) + 
  labs(subtitle = label) +
  theme(legend.position = "bottom", legend.box = "vertical", plot.subtitle=element_text(color="gray36"))
