##------------------------------------------------------------
## RANDOM SURVIVAL FOREST
##------------------------------------------------------------
library(randomForestSRC)
library(ggplot2)
library(survival)
library(survminer)
library(dplyr)
library(foreign)
library(tidyr)
library(survivalROC)
library(survAUC)
library(patchwork)

proj_path = "[your project path]"

birth1 = read.dta(paste0(proj_path, "[your data]"))

process_data = function(df) {
  df = df %>%
    mutate(educ = factor(educ, labels= c("1","2","3"))) %>%
    select(duration, event, cohort, religion_important, sample, t_frsoe, educ, sexee,
           motive, lnaism, educ_p, educ_m, lnaisp, t_elev_a, t_elev_b, t_elev_d, t_elev_e,
           t_argent, t_malad, t_dispar, t_alcool, t_violen) %>%
    drop_na()
  return(df)
}

birth1 = process_data(birth1)

set.seed(131)
set.seed(123)
rfsrc_birth1 <- rfsrc(Surv(duration, event) ~ ., data = birth1, ntree= 1000, block.size = 1, nsplit = 10 , na.action = "na.impute", tree.err = TRUE,importance = TRUE)

# Out of Bag Errors (OOB)
plot(rfsrc_birth1) + theme_bw()


# ROC

roc_age20 = data.frame(survivalROC.C(Stime=birth1$duration, status=birth1$event, marker=rfsrc_birth1$predicted.oob, predict.time = 60)) %>%
  mutate(name = "Age 20") %>%
  select(name, FP, TP, AUC)

roc_age30 = data.frame(survivalROC.C(Stime=birth1$duration, status=birth1$event, marker=rfsrc_birth1$predicted.oob, predict.time = 180)) %>%
  mutate(name = "Age 30") %>%
  select(name, FP, TP, AUC)

roc_age40 = data.frame(survivalROC.C(Stime=birth1$duration, status=birth1$event, marker=rfsrc_birth1$predicted.oob, predict.time = 300)) %>%
  mutate(name = "Age 40") %>%
  select(name, FP, TP, AUC)

roc_age50 = data.frame(survivalROC.C(Stime=birth1$duration, status=birth1$event, marker=rfsrc_birth1$predicted.oob, predict.time = 420)) %>%
  mutate(name = "Age 50") %>%
  select(name, FP, TP, AUC)

roc_combined = rbind(roc_age20, roc_age30, roc_age40, roc_age50)

aucs = roc_combined %>%
  select(name, AUC) %>%
  distinct() %>%
  mutate(FP=0.6,
         TP=c(0.45,0.4,0.35,0.3))

roc_plot = ggplot(roc_combined, aes(x=FP*100, y=TP*100, color=name))+
  geom_line()+
  theme_bw()+
  geom_abline()+
  scale_color_discrete("Age (years)", labels=c("20", "30", "40", "50"))+
  xlab("False Positive (%)")+
  ylab("True Positive (%)")+
  geom_text(data = aucs,
            aes(x = FP*100, y = TP*100, label = paste0("AUC = ", round(AUC,2)), color = name),
            hjust = 0,
            show.legend = FALSE)+
  coord_fixed(expand=F)


ggsave(paste(proj_path,"/figures/roc_first_birth.png",sep=""),
       plot=roc_plot,dpi=400,height=5,width=6)


# VIMP plots

labels_list = c("age_months"="Age","sexee"="Gender","educ"="Education","t_frsoe"="Family size",
                "cohort"="Birth cohort","religion_important"="Religion is important","origind"="Origin group",
                "sample"="Generation (1,2 or natives)","m_avisit"="Visited France before arrival",
                "inat"="French nationality","motive"="Reason for migration","lnaism"="Mother born abroad",
                "educ_p"="Father's education","educ_m"="Mother's education","m_vise"="Visited country of origin",
                "lnaisp"="Father born abroad","t_elev_a"="Raised by both parents, in union",
                "t_elev_b"="Raised by both parents, separated","t_elev_c"="Raised by polygamous family",
                "t_elev_d"="Raised by mother alone","t_elev_e"="Raised by father alone",
                "t_elev_f"="Raised by father and new partner","t_elev_g"="Raised by mother and new partner",
                "t_elev_h"="Raised by grandparents","t_elev_i"="Raised by other family member",
                "t_elev_j"="Raised in institution","t_elev_k"="Raised in host family",
                "t_elev_l"="Raised by employer","t_elev_m"="Raised by other",
                "t_argent"="Money issues during childhood", "t_malad"="Parents had health issues during childhood",
                "t_dispar"="Conflicts between parents during childhood",
                "t_alcool"="Parents had alcoholism problems during childhood",
                "t_violen"="Violence during childhood")


# Variable selection: variable importance
VIMP_plot = plot(gg_vimp(rfsrc_birth1)) +
  theme(legend.position = c(0.8, 0.2)) +
  labs(fill = "VIMP > 0") + theme_bw() +
  scale_x_discrete(labels=labels_list) +
  ggtitle("A")

# Minimal depth
varsel_birth1 <- var.select(rfsrc_birth1)
gg_md <- gg_minimal_depth(varsel_birth1)
minimal_depth = plot(gg_md) +
  theme_bw() +
  scale_x_discrete(labels=labels_list) +
  ggtitle("B")

# Variable selection comparison
comparison = plot(gg_minimal_vimp(gg_md), lbls) +
  theme(legend.position=c(0.8, 0.2)) + theme_bw() +
  scale_x_discrete(labels=labels_list) +
  ggtitle("C")

combined = VIMP_plot/minimal_depth/comparison

ggsave(paste(proj_path,"/figures/first_birth_combined_v2.png",sep=""),
       plot=combined,dpi=400,height=10,width=6)


# PREDICTED PROBABILITIES

ggvar <- gg_variable(rfsrc_birth1,time=180,time.labels=c("Age 30"))
partial_coplot <- gg_partial_coplot(rfsrc_birth1, xvar = "sample",groups = ggvar$cohort,surv_type = "surv",time=180,show.plots = FALSE) %>%
  mutate(time = "Age 30")

ggvar <- gg_variable(rfsrc_birth1,time=360,time.labels=c("Age 45"))
partial_coplot2 <- gg_partial_coplot(rfsrc_birth1, xvar = "sample",groups = ggvar$cohort,surv_type = "surv",time=360,show.plots = FALSE) %>%
  mutate(time = "Age 45")

total <- rbind(partial_coplot,partial_coplot2) %>%
  mutate(sample = factor(sample, levels = c("Immigrants", "Descendants", "Natives")),
         sample_time = interaction(sample, time)) %>%
  filter(group =="1948-1959" | group=="1960-1969" | (group=="1970-1979" & time=="Age 30"))

pred_probs_first_birth = ggplot(total, aes(x=sample, y=1-yhat)) +
    geom_boxplot(fill="steelblue", size=0,outlier.shape = NA,outlier.size = 0) +labs(y = "Predicted probability of having a first birth (%)",x="Population group") + theme_bw() + theme(legend.position = "none") + scale_colour_discrete("Time") + coord_cartesian(ylim = c(0.4,1)) + facet_grid(group~time)

ggsave(paste(proj_path,"/figures/pred_probs_first_birth.png",sep=""),
       plot=pred_probs_first_birth,dpi=400,height=6,width=6)