#'---
#' author: "Erin dela Cruz"
#' ---
#' 
#' **Program**: cohort2-graphfunctions.R
#' 
#' **Purpose**: Houses a suite of custom R functions for analyzing Cohort 2
#' 
#' **Last Edited**: January 16, 2020
#' ***
library(lazyeval)
theme_set(theme_tufte())


#' Risk of BV among Sufficient and Deficient
bvriskgraph <- function(BaselineDF, TLRVariable, NugentVariable, AmselVariable, RSNumber, WTString, DefString) {
  nug <- BaselineDF %>% group_by_(TLRVariable) %>%
    summarise_(pct = interp(~mean(v, na.rm=TRUE), v=as.name(NugentVariable)),
               sd=interp(~sd(v, na.rm=TRUE), v=as.name(NugentVariable)),
               n=interp(~n())) %>%
    mutate(se = sd/sqrt(n)) %>%
    na.omit() %>%
    mutate(dx="Nugent")
  print(nug)
  ams <- BaselineDF %>% group_by_(TLRVariable) %>%
    summarise_(pct = interp(~mean(v, na.rm=TRUE), v=as.name(AmselVariable)),
               sd=interp(~sd(v, na.rm=TRUE), v=as.name(AmselVariable)),
               n=interp(~n())) %>%
    mutate(se = sd/sqrt(n)) %>%
    na.omit() %>%
    mutate(dx="Amsel")
  print(ams)
  bvdxes <- bind_rows(nug, ams)
  bvdxes <- as.data.frame(bvdxes)
  
  diagnosis <- ggplot(bvdxes, aes_string(x="dx", y="pct", color=TLRVariable)) +
    guides(color = guide_legend(title=paste(RSNumber))) +
    geom_pointrange(aes(ymin=pct-1.96*se, ymax=pct+1.96*se),
                    size=.3, 
                    position=position_dodge(.9)) +
    labs(y="Proportion Diagnosed with BV", x="") + 
    scale_color_manual(labels=c(paste(WTString), paste(DefString)), values=c("#f1a340", "#998ec3")) +
    scale_x_discrete(labels=c("Amsel\n(Clinical)", "Nugent\n(Microbiological)")) +
    theme(legend.position = "bottom", plot.subtitle=element_text(color="gray36"))
  
  
  return(diagnosis)
}

#' Colonization with BVAB Plots for publication
mini.smear <- function(DataFrame, TLRVariable, HypothesisDirection, RSNumber, WTString, DefString, BugString) {
  label <- substitute(paste("Mean log10 ", Delta, "=", estimate, " (p=", pvalue, ")"),
                      list(estimate = signif(t.test(DataFrame$log.copies.per.swab~DataFrame[[TLRVariable]], alternative = c(HypothesisDirection))$estimate[2]-
                                               t.test(DataFrame$log.copies.per.swab~DataFrame[[TLRVariable]], alternative = c(HypothesisDirection))$estimate[1], 3),
                           pvalue = signif(t.test(DataFrame$log.copies.per.swab~DataFrame[[TLRVariable]], alternative = c(HypothesisDirection))$p.value, 3)))
  # label <- substitute(paste("Mean log10 ", Delta, "=", estimate, " (p=", pvalue, ")"),
  #                     list(estimate = signif(t.test(DataFrame$log.copies.per.swab~DataFrame[[TLRVariable]], alternative = c(HypothesisDirection))$estimate[2]-
  #                                              t.test(DataFrame$log.copies.per.swab~DataFrame[[TLRVariable]], alternative = c(HypothesisDirection))$estimate[1], 3),
  #                          pvalue = signif(wilcox.test(DataFrame$log.copies.per.swab~DataFrame[[TLRVariable]])$p.value, 3)))
  plot <- ggplot(DataFrame, aes_string(x=TLRVariable, y="copies_per_swab", 
                                       color=TLRVariable,
                                       shape=TLRVariable)) +
    guides(color="none", shape="none") +
    geom_boxplot(width=0.25, outlier.shape=NA) +
    geom_point(position = position_jitter(height = 0, width = .15), size=2) +
    scale_y_log10(paste0(BugString, "\n16S rRNA gene copies per swab"),
                  breaks=c(1e2,1e4,1e6,1e8,1e10)) +
    scale_fill_discrete(na.value=NA) +
    scale_shape_discrete(solid=FALSE) +
    scale_color_manual(values=c("#f1a340", "#998ec3")) +
    scale_x_discrete(paste(RSNumber), labels=c(paste(WTString), paste(DefString))) +
    labs(subtitle=label) +
    theme(plot.subtitle=element_text(color="gray36"))
  
  
  return(plot)
}