#'---
#' title: "Cohort 1: TLR5 deficiency and Colonization with BVAB"
#' author: "Erin dela Cruz"
#' ----
#' 
#' **Program**: Cohort1_TLR5_BVAB1_analysis.R
#' 
#' **Purpose**: Analysis of TLR5 deficiency and colonization with BVAB1 and other vaginal bacteria in Cohort 1
#' 
#+ setup-libs, include=FALSE
#' **Prep**: Clear data from memory
rm(list = ls())
# Tidyverse Libs
library(tidyverse)
# ggplot add-ons
library(GGally)
library(ggthemes)
library(ggpubr)
library(cowplot)
# For matching genotype groups
library(MatchIt)
# For Survival Analysis
library(survival)
# For multilevel modeling
library(lme4)
library(lmerTest)

#' For use on Windows
fp <- "H:/code/tlr5-deficiency-BV_PUBLICATION/cohort1/"

source(paste0(fp,"src/cohort1-graphfunctions.R"))

#' -----------------------------------------------------------------------------------------------
#' Set TLR labels for all graphs
tlrvariable <- "TLR5_deficient"
rsnumber <- "TLR5 rs5744168"
wildtype <- "CC (Sufficient)\nn=32"
deficient <- "CT/TT (Deficient)\nn=16"
wildtype.wo.n <- "CC (Sufficient)"
deficient.wo.n <- "CT/TT (Deficient)"
#' -----------------------------------------------------------------------------------------------

#' Load dataset - all observations at start of episode of BV
start.bv <- read.csv(paste0(fp, "indata/tlr5-colonization.csv")) %>%
  group_by(pubID)

#' -----------------------------------------------------------------------------------------------
#' Section #1: Match TLR5 deficient subjects to sufficient subjects
#' 
#' Output dataset: gene.matched
#' -----------------------------------------------------------------------------------------------

#' Create % Days BV+ variable in dataset with one line per pubID
#' Exclude any observations where TLR5 genotype did not resolve
start.bv.obs <- start.bv %>% group_by(pubID, TLR5_deficient, nugent.dx, age, age_cat, sexual_pref, wsw, wsm, con.any.hormonal, black, white) %>% 
  summarise(days.bv = sum(hc.nug.bv.dx, na.rm=TRUE), 
            person.time=n(), 
            pct.bv.days = 100*(days.bv/person.time)) %>%
  filter(is.na(TLR5_deficient)==FALSE) %>%
  ungroup()

#' Find matched controls for TLR5 at risk study participants
tlr5.matching.df <- start.bv.obs %>% 
  filter(is.na(nugent.dx)==FALSE &
           is.na(TLR5_deficient)==FALSE &
           is.na(age_cat)==FALSE &
           is.na(sexual_pref)==FALSE &
           is.na(con.any.hormonal)==FALSE &
           is.na(white)==FALSE) %>%
  select(pubID, nugent.dx, TLR5_deficient, age_cat, sexual_pref, con.any.hormonal, white)
tlr5.matching.df <- as.data.frame(tlr5.matching.df)
set.seed(1809)
mi <- matchit(TLR5_deficient~nugent.dx + age_cat + white + 
                sexual_pref + con.any.hormonal, 
              data=tlr5.matching.df,
              ratio=2)
summary(mi)
tlr5.matched <- match.data(mi, group = "all") %>% select(pubID, distance, weights)
gene.matched <- semi_join(start.bv, tlr5.matched, by="pubID")


#' -----------------------------------------------------------------------------------------------
#' Section #4: Check balance between groups
#' -----------------------------------------------------------------------------------------------

#' Check balance of baseline characteristics - used in table included in supplement
start.bv.day0 <- gene.matched %>% top_n(1, post.day) 
tlr5 <- start.bv.day0 %>% select(pubID, TLR5_deflabel) %>% distinct(.)
#' Describe Study Population
table(tlr5$TLR5_deflabel, useNA="ifany")
study.pop.chars <- left_join(tlr5, start.bv.obs, by="pubID") %>%
  group_by(TLR5_deflabel) %>%
  summarise(pct.nugent.pos = mean(nugent.dx)*100,
            mean.age = mean(age),
            pct.black = mean(black)*100,
            pct.wsm = mean(wsm)*100,
            pct.wsw = mean(wsw)*100,
            pct.hc = mean(con.any.hormonal)*100, 
            n=n())

#' Survival analysis of longitudinal risk of BV: TLR5 deficiency vs sufficiency
ts.data <- read.csv(paste0(fp, "indata/time-to-next-BV_colonization-analysis.csv"), header=TRUE, check.names=TRUE, as.is=TRUE)
clin.amsel <- left_join(tlr5, ts.data, by="pubID")
clin.amsel.df <- as.data.frame(clin.amsel)

survival <- surv.analysis.nocens(DtaFrame=clin.amsel.df, TollVar=tlrvariable, DayVar="days", 
                                 OutcomeVar="final.amsel", PrrLabel=rsnumber, 
                                 PrrWTLabel=wildtype, PrrDefLabel=deficient)
survival

#' Check # of observations included at each day post-analysis
num.pts <- gene.matched %>% group_by(post.day, TLR5_deficient) %>%
  summarise(n=n())
ggplot(num.pts, aes(x=as.numeric(post.day), y=n, color=factor(TLR5_deficient))) +
  geom_line() + facet_wrap(~TLR5_deficient, scales="free_y")

#' Make a dataset specifically for ploting kinetics of colonization
start.bv.longitudinal <- gene.matched %>% filter(post.day<=90)
sbl <- as.data.frame(start.bv.longitudinal)

#' % of each group that is BV+ by nugent at any given timepoint
nugent <- nugent.plot(sbl, tlrvariable, rsnumber, wildtype, deficient)
nugent

#' Create table describing balance of baseline characteristics as well as % days with BV by AA versus AG/GG
days.bv.cohort <- gene.matched %>% group_by(TLR5_deflabel) %>%
  summarise(days.bv = sum(hc.nug.bv.dx, na.rm=TRUE), person.time=n(), pct.bv.days = 100*(days.bv/person.time)) %>%
  full_join(., study.pop.chars, by="TLR5_deflabel") %>%
  select(TLR5_deflabel, n, pct.nugent.pos, person.time, pct.bv.days, mean.age, pct.black, pct.wsm, pct.wsw, pct.hc)
cohort.print <- as.data.frame(t(days.bv.cohort)) %>%
  mutate(V1=as.numeric(as.character(V1))) %>% mutate(V2=as.numeric(as.character(V2))) %>%
  filter(is.na(V1)==FALSE)
cohort.table <- ggtexttable(signif(cohort.print,3), rows=c("N", "% Nugent+ at Enrollment", "Person-Days Observed", "% Days Nugent+",
                                                           "Age (Mean)", "% Black" ,"% WSM", "% WSW", "% Hormonal Contraception"),
                            cols=c(paste0(rsnumber,"\n", wildtype.wo.n), paste0(rsnumber,"\n", deficient.wo.n)),
                            theme=ttheme("minimal"))
cohort.table



#' -----------------------------------------------------------------------------------------------
#' Section #5.1: Modeling/Statistical Testing: Relationship between colonization/levels of BVAB and TLR5
#' -----------------------------------------------------------------------------------------------


#' ***
#' **Modeling/Statistical Testing**
#' Relationship between colonization/levels of BVAB and TLR5
#' Create a matrix to hold betas and SEs
beta.se <- matrix(nrow=0,ncol=3)
anova <- matrix(nrow=0, ncol=1)
names <- c("BVAB1", "G. vaginalis", "L. crispatus")
gene.matched <- as.data.frame(gene.matched)
for (i in seq(25,27)) {
  model <- lmer(log10(gene.matched[,i])~gene.matched$TLR5_deficient + (1 | gene.matched$pubID),
                control=lmerControl(
                  optimizer="Nelder_Mead"
                ))
  model.full <- lmer(
    log10(gene.matched[,i])~gene.matched$TLR5_deficient + (1 | gene.matched$pubID), 
    REML=FALSE,
    control=lmerControl(
      optimizer="Nelder_Mead"
    )
  )
  model.null <- lmer(
    log10(gene.matched[,i])~(1 | gene.matched$pubID), 
    REML = FALSE,
    control=lmerControl(
      optimizer="Nelder_Mead"
    )
  )
  temp <- c(summary(model)$coef[2,1], summary(model)$coef[2,2], anova(model.null, model.full)$"Pr(>Chisq)"[[2]])
  beta.se <- rbind(beta.se, temp)
}
beta.se <- cbind(beta.se, names)
b.se.graph <- as.data.frame(beta.se, row.names = TRUE)
b.se.graph <- b.se.graph %>%
  mutate(b = (as.numeric(as.character(V1)))) %>%
  mutate(se = (as.numeric(as.character(V2)))) %>%
  mutate(p.anova = (as.numeric(as.character(V3)))) %>%
  rename(bacteria=names) %>%
  mutate(bacteria=as.character(bacteria)) %>%
  select(b, se, p.anova, bacteria)
b.se.graph$bacteria <- factor(b.se.graph$bacteria, levels = b.se.graph$bacteria[order(-b.se.graph$b)])
effect.sizes <- ggplot(b.se.graph, aes(x=b, y=bacteria)) +
  geom_errorbarh(aes(xmin=b-1.96*se, xmax=b+1.96*se), height=0) +
  geom_point() +
  geom_vline(xintercept=0, size=0.25) +
  labs(x="TLR5 deficiency: Mean log10 change in bacterial colonization", y="")


#' Add statistical modeling results to smear plots
#' smearplot <- function(DtaFrame, TollVar, BugVar, TollLabel, SufLabel, DefLabel, BugLabel, Ests, BugIndexString)
bvab1 <- smearplot(gene.matched, "TLR5_deficient","BVAB.1", rsnumber, wildtype, deficient,
                   "BVAB1", b.se.graph, "BVAB1")
bvab1.avg.kinetics <- average.kinetics(sbl, "TLR5_deficient", "BVAB.1", rsnumber, wildtype, deficient,
                                       "BVAB1")
gvag <- smearplot(gene.matched, "TLR5_deficient","G.vag", rsnumber, wildtype, deficient,
                  "Gardnerella vaginalis", b.se.graph, "G. vaginalis")
gv.avg.kinetics <- average.kinetics(sbl, "TLR5_deficient", "G.vag", rsnumber, wildtype, deficient,
                                    "Gardnerella vaginalis")
lcrisp <- smearplot(gene.matched, "TLR5_deficient","L.crisp", rsnumber, wildtype, deficient,
                    "Lactobacillus crispatus", b.se.graph, "L. crispatus")
lcrisp.avg.kinetics <- average.kinetics(sbl, "TLR5_deficient", "L.crisp", rsnumber, wildtype, deficient,
                                        "L. crispatus")

#' -----------------------------------------------------------------------------------------------
#' Section #6.1: Make figures for publication
#' -----------------------------------------------------------------------------------------------

#' Balance between matched cases and controls
bal <- plot_grid(survival + theme(legend.position="none"), nugent + theme(legend.position="none"),
                 align = 'vh',
                 labels = c("A", "B"),
                 hjust = -1,
                 nrow = 1)
gp.balance.legend <- get_legend(nugent+ theme(legend.position="bottom"))
bal1 <- plot_grid(bal, gp.balance.legend, ncol=1, rel_heights = c(1, .2))
bal.table <- plot_grid(NULL,cohort.table,NULL, nrow=1, rel_widths = c(1, 1,1), labels=c("","C",""))
bal2 <- plot_grid(bal1, bal.table, ncol=1, rel_heights = c(1,0.75))

#' BVAB1 figure
title1 <- ggdraw() + draw_label("Boxplots of concentration,\nordered by study participant\nand mean concentration", fontface='bold', size=12)
title2 <- ggdraw() + draw_label("Post-BV diagnosis, daily average\nconcentration across study\nparticipants, aggregated by\nTLR5 genotype", fontface='bold', size=12)
titles <- plot_grid(title1, title2, ncol=2)
b1 <- plot_grid(bvab1 + theme(legend.position="bottom"), bvab1.avg.kinetics + theme(legend.position="bottom"),
                align="vh", labels = c("A", "B"), hjust=-1, nrow=1)
qPCR <- plot_grid(titles, b1, ncol=1, rel_heights = c(0.3, 1))