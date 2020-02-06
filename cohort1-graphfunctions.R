#'---
#' author: "Erin dela Cruz"
#' ---
#' 
#' **Program**: cohort1-graphfunctions.R
#' 
#' **Purpose**: Houses a suite of custom R functions
#' 
#' **Last Edited**: January 15, 2019
#' 
#' **Dependencies:  None
#' 
theme_set(theme_tufte())

# Function: Smear Plots of Colonization with BugVar by pubID, color-coded by TLR deficiency
smearplot <- function(DtaFrame, TollVar, BugVar, TollLabel, SufLabel, DefLabel, BugLabel, Ests, BugIndexString) {
  mdf <- data.frame(pubID=as.factor(DtaFrame$pubID),
                    yvar=DtaFrame[[BugVar]],
                    tlrdef=as.factor(DtaFrame[[TollVar]]))
  
  label <- substitute(paste("Mean log10 ", Delta, "=", estimate, "; 95% CI: [", lo, ", ", hi, "]; p=", pvalue),
                      list(estimate = signif(Ests[which(Ests$bacteria==BugIndexString),]$b, 3),
                           lo = signif(Ests[which(Ests$bacteria==BugIndexString),]$b - (1.96*Ests[which(Ests$bacteria==BugIndexString),]$se), 3),
                           hi = signif(Ests[which(Ests$bacteria==BugIndexString),]$b + (1.96*Ests[which(Ests$bacteria==BugIndexString),]$se), 3),
                           pvalue = signif(Ests[which(Ests$bacteria==BugIndexString),]$p.anova, 3)))
  
  smear.d <- ggplot(mdf,aes(x=reorder(pubID,yvar,FUN=mean),
                            y=yvar, 
                            color=tlrdef,
                            shape=tlrdef)) +
    guides(fill = guide_legend(title=paste(TollLabel)), 
           linetype=guide_legend(title=paste(TollLabel)),
           color = guide_legend(title=paste(TollLabel)),
           shape = guide_legend(title=paste(TollLabel))) +
    geom_boxplot(width=0.5, lwd=0.33, outlier.shape=NA) +
    geom_point(position = position_jitter(height = 0, width = .15), size=1) +
    scale_fill_discrete(labels=c(paste(SufLabel), paste(DefLabel)), na.value=NA) +
    scale_shape_discrete(labels=c(paste(SufLabel), paste(DefLabel)), solid=FALSE) +
    scale_color_manual(labels=c(paste(SufLabel), paste(DefLabel)), values=c("#f1a340", "#998ec3")) +
    labs(y=paste0(BugLabel, "\n16S rRNA gene copies per swab"), 
         x=paste0("Study Participants\nRanked by mean copies per swab"),
         subtitle=label) +
    scale_y_log10(breaks=c(1e2, 1e4, 1e6, 1e8, 1e10)) +
    theme(axis.text.x = element_blank(), 
          axis.ticks.x = element_blank(),
          plot.subtitle=element_text(color="gray36"))
  return(smear.d)
}

# Function: Average Kinetics of Colonization with BugVar
average.kinetics <- function(DtaFrame, TollVar, BugVar, TollLabel, SufLabel, DefLabel, BugLabel) {
  mdf <- data.frame(day=as.numeric(DtaFrame$post.day),
                    yvar=DtaFrame[[BugVar]],
                    tlrdef=as.factor(DtaFrame[[TollVar]]))
  #' geom_point(position = position_jitter(height = 0, width = .15), size=1, alpha=0.5) +
  plot <- ggplot(mdf, aes(x=day, y=yvar,
                          linetype=tlrdef, color=tlrdef, fill=tlrdef, shape=tlrdef))+
    guides(fill = guide_legend(title=paste(TollLabel)), 
           linetype=guide_legend(title=paste(TollLabel)),
           color = guide_legend(title=paste(TollLabel)),
           shape = guide_legend(title=paste(TollLabel))) +
    stat_summary(fun.data=mean_cl_boot, geom = "ribbon", fun.args=list(conf.int=0.95), linetype=0) +
    stat_summary(fun.data=mean_cl_boot, geom = "smooth") +
    scale_y_log10(breaks=c(1e2, 1e4, 1e6, 1e8, 1e10)) +
    scale_x_continuous(breaks=c(0,30,60,90,120)) +
    scale_linetype_discrete(labels=c(paste(SufLabel), paste(DefLabel))) +
    scale_fill_manual(labels=c(paste(SufLabel), paste(DefLabel)), values=alpha(c("#f1a340", "#998ec3"), 0.2)) +
    scale_color_manual(labels=c(paste(SufLabel), paste(DefLabel)), values=c("#f1a340", "#998ec3")) +
    scale_shape_manual(labels=c(paste(SufLabel), paste(DefLabel)), values=c("#f1a340", "#998ec3")) +
    labs(x="Days after BV diagnosis\n(by Nugent Score)", 
         y=paste0(BugLabel, "\n16S rRNA gene copies per swab")) +
    theme(plot.subtitle=element_text(color="gray36"))
  
  return(plot)
}

# Function: nugent.plot
# Purpose: Shows moving daily average proportion Nugent+ by TLR5 deficient/sufficient
nugent.plot <- function(DtaFrame, TollVar, TollLabel, SufLabel, DefLabel) {
  mdf <- data.frame(day=as.numeric(DtaFrame$post.day),
                    yvar=as.numeric(DtaFrame$hc.nug.bv.dx),
                    tlrdef=as.factor(DtaFrame[[TollVar]]))
  
  ggplot(subset(mdf, is.na(yvar)==FALSE), aes(x=day, y=yvar,
                                              linetype=tlrdef, color=tlrdef)) +
    guides(color = guide_legend(title=paste(TollLabel)), linetype=guide_legend(title=paste(TollLabel))) +
    stat_summary(fun.y="mean", geom = "line") +
    scale_x_continuous(breaks=c(0,30,60,90,120)) +
    scale_linetype_discrete(labels=c(paste(SufLabel), paste(DefLabel))) +
    scale_color_manual(labels=c(paste(SufLabel), paste(DefLabel)), values=c("#f1a340", "#998ec3")) +
    labs(x="Days after BV diagnosis\n(by Nugent Score)", 
         y="Proportion BV+ by Nugent Score") +
    theme(plot.subtitle=element_text(color="gray36"))
}

#' Function: "surv.analysis"
# Purpose: displays survival curves - time to first episode of BV after enrollment
#          uses only observations that are BV-negative at enrollment
surv.analysis <- function(DtaFrame, TollVar, DayVar, OutcomeVar, PrrLabel, PrrWTLabel, PrrDefLabel) {
  mdf <- data.frame(day=as.numeric(DtaFrame[[DayVar]]),
                    outcome=DtaFrame[[OutcomeVar]],
                    tlrdef=DtaFrame[[TollVar]])
  surv.all <- survfit(Surv(mdf$day, mdf$outcome)~mdf$tlrdef, conf.type="none")
  print(summary(surv.all))
  print(survdiff(Surv(mdf$day, mdf$outcome)~mdf$tlrdef))
  hr <- summary(coxph(Surv(mdf$day, mdf$outcome)~mdf$tlrdef))$coef[1,2]
  p <- summary(coxph(Surv(mdf$day, mdf$outcome)~mdf$tlrdef))$coef[1,5]
  ci.low <- summary(coxph(Surv(mdf$day, mdf$outcome)~mdf$tlrdef))$conf.int[1,3]
  ci.hi <- summary(coxph(Surv(mdf$day, mdf$outcome)~mdf$tlrdef))$conf.int[1,4]
  label <- substitute(paste("Hazard Ratio=", haz, "; 95% CI: [", lo, ", ", high, "]; p=", pvalue),
                      list(haz = signif(hr, 3), lo = signif(ci.low, 3), high = signif(ci.hi, 3), pvalue = signif(p, 3)))
  
  graph <- ggsurv(surv.all, plot.cens = TRUE, order.legend = FALSE, 
                  lty.est=c(1,2), surv.col = 1) +
    guides(colour=guide_legend(title=paste(PrrLabel),nrow=1,byrow=TRUE), 
           linetype=guide_legend(title=paste(PrrLabel),nrow=1,byrow=TRUE)) +
    theme_tufte() + 
    scale_linetype_discrete(labels = c(paste(PrrWTLabel), paste(PrrDefLabel))) +
    scale_color_manual(labels = c(paste(PrrWTLabel), paste(PrrDefLabel)), values=c("#f1a340", "#998ec3")) +
    scale_x_continuous(name = "Days post-enrollment", breaks = seq(0, 500, by=50)) +
    scale_y_continuous(name = "Proportion without clinical\ndiagnosis of BV", breaks = seq(0, 1, by=0.2), limits = c(0,1)) +
    labs(subtitle=label) +
    theme(plot.subtitle=element_text(color="gray36"))
  
  
  return(graph)
}

# Function: "surv.analysis.n"
# Purpose: As above, labels as "microbiological" diagnosis (for use with time to first Nugent BV+)
surv.analysis.n <- function(DtaFrame, TollVar, DayVar, OutcomeVar, PrrLabel, PrrWTLabel, PrrDefLabel) {
  mdf <- data.frame(day=as.numeric(DtaFrame[[DayVar]]),
                    outcome=DtaFrame[[OutcomeVar]],
                    tlrdef=DtaFrame[[TollVar]])
  surv.all <- survfit(Surv(mdf$day, mdf$outcome)~mdf$tlrdef, conf.type="none")
  print(summary(surv.all))
  print(survdiff(Surv(mdf$day, mdf$outcome)~mdf$tlrdef))
  hr <- summary(coxph(Surv(mdf$day, mdf$outcome)~mdf$tlrdef))$coef[1,2]
  p <- summary(coxph(Surv(mdf$day, mdf$outcome)~mdf$tlrdef))$coef[1,5]
  ci.low <- summary(coxph(Surv(mdf$day, mdf$outcome)~mdf$tlrdef))$conf.int[1,3]
  ci.hi <- summary(coxph(Surv(mdf$day, mdf$outcome)~mdf$tlrdef))$conf.int[1,4]
  label <- substitute(paste("Hazard Ratio=", haz, "; 95% CI: [", lo, ", ", high, "]; p=", pvalue),
                      list(haz = signif(hr, 3), lo = signif(ci.low, 3), high = signif(ci.hi, 3), pvalue = signif(p, 3)))
  
  graph <- ggsurv(surv.all, plot.cens = TRUE, order.legend = FALSE, 
                  lty.est=c(1,2), surv.col = 1) +
    guides(colour=guide_legend(title=paste(PrrLabel),nrow=1,byrow=TRUE), 
           linetype=guide_legend(title=paste(PrrLabel),nrow=1,byrow=TRUE)) +
    theme_tufte() + 
    scale_linetype_discrete(labels = c(paste(PrrWTLabel), paste(PrrDefLabel))) +
    scale_color_manual(labels = c(paste(PrrWTLabel), paste(PrrDefLabel)), values=c("#f1a340", "#998ec3")) +
    scale_x_continuous(name = "Days post-enrollment", breaks = seq(0, 500, by=50)) +
    scale_y_continuous(name = "Proportion without microbiological\ndiagnosis of BV", breaks = seq(0, 1, by=0.2), limits = c(0,1)) +
    labs(subtitle=label) +
    theme(plot.subtitle=element_text(color="gray36"))
  
  
  return(graph)
}


surv.analysis.nocens <- function(DtaFrame, TollVar, DayVar, OutcomeVar, PrrLabel, PrrWTLabel, PrrDefLabel) {
  mdf <- data.frame(day=as.numeric(DtaFrame[[DayVar]]),
                    outcome=DtaFrame[[OutcomeVar]],
                    tlrdef=DtaFrame[[TollVar]])
  surv.all <- survfit(Surv(mdf$day, mdf$outcome)~mdf$tlrdef, conf.type="none")
  print(summary(surv.all))
  print(survdiff(Surv(mdf$day, mdf$outcome)~mdf$tlrdef))
  hr <- summary(coxph(Surv(mdf$day, mdf$outcome)~mdf$tlrdef))$coef[2]
  p <- summary(coxph(Surv(mdf$day, mdf$outcome)~mdf$tlrdef))$coef[5]
  ci.low <- summary(coxph(Surv(mdf$day, mdf$outcome)~mdf$tlrdef))$conf.int[3]
  ci.hi <- summary(coxph(Surv(mdf$day, mdf$outcome)~mdf$tlrdef))$conf.int[4]
  label <- substitute(paste("Hazard Ratio=", haz, "; 95% CI: [", lo, ", ", high, "]; p=", pvalue),
                      list(haz = signif(hr, 3), lo = signif(ci.low, 3), high = signif(ci.hi, 3), pvalue = signif(p, 3)))
  
  graph <- ggsurv(surv.all, plot.cens = FALSE, order.legend = FALSE, 
                  lty.est=c(1,2), surv.col = 1) +
    guides(colour=guide_legend(title=paste(PrrLabel),nrow=1,byrow=TRUE), 
           linetype=guide_legend(title=paste(PrrLabel),nrow=1,byrow=TRUE)) +
    theme_tufte() + 
    scale_linetype_discrete(labels = c(paste(PrrWTLabel), paste(PrrDefLabel))) +
    scale_color_manual(labels = c(paste(PrrWTLabel), paste(PrrDefLabel)), values=c("#f1a340", "#998ec3")) +
    scale_x_continuous(name = "Days post-enrollment", breaks = seq(0, 500, by=50)) +
    scale_y_continuous(name = "Proportion without clinical\ndiagnosis of BV", breaks = seq(0, 1, by=0.2), limits = c(0,1)) +
    labs(subtitle=label) +
    theme(plot.subtitle=element_text(color="gray36"))
  
  
  return(graph)
}