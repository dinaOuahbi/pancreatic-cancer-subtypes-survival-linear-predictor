
#Download the IDH status and 1p/19q codeltetion status data of TCGA-GBM and TCGA-LGG==============
#================================================================
#======== Method 1:  TCGAquery_subtype ========
#https://www.bioconductor.org/packages/release/bioc/vignettes/TCGAbiolinks/inst/doc/subtypes.html
install.packages("openxlsx", dependencies=TRUE)

library(TCGAbiolinks)
library(openxlsx)

paad_subtype <- TCGAquery_subtype(tumor = "paad")
paad_subtype$`mRNA Bailey Clusters (All 150 Samples) 1squamous 2immunogenic 3progenitor 4ADEX`

write.xlsx(paad_subtype, file = "/work/shared/ptbc/CNN_Pancreas_V2/RNAseq_TCGA_analysis/data/paad_subtype.xlsx")



# SCRIPT PYTHON : /work/shared/ptbc/CNN_Pancreas_V2/RNAseq_TCGA_analysis/stats_analysis/generate_subtype_data.py


df<- read.csv('/work/shared/ptbc/CNN_Pancreas_V2/RNAseq_TCGA_analysis/data/model8/paad_subtype.csv')[-1]
df$condition=as.factor(df$condition)
df$cancer_subtype = as.factor(df$cancer_subtype)

########################## CONDITION - SUBTYPES 

subtype_groups <- table(df$condition,
                   df$cancer_subtype)

# Mosaic plots
mosaicplot(subtype_groups)


# Proportions of the Contingency Tables
prop.table(subtype_groups)


# Rows and Columns Totals
addmargins(subtype_groups)

# Statistical Tests
# in order to test if the relationship of these two variables is independent or not
chisq.test(subtype_groups)

# FISHER
fisher.test(subtype_groups)


#  Log Likelihood Ratio
library(MASS)
loglm( ~ 1 + 2, data = subtype_groups)


########################## LINEAR PREDICTION - SUBTYPES 

library(ggplot2)

# p <- ggplot(df, aes(x=cancer_subtype, y=pred, fill=cancer_subtype)) + 
#   geom_boxplot(notch=T,outlier.colour="red", outlier.shape=8,outlier.size=4)+
#   stat_summary(fun.y=mean, geom="point", shape=23, size=4)+
#   geom_jitter(shape=16, position=position_jitter(0.2))+
#   scale_color_brewer(palette="Dark2")+
#   theme(legend.position="bottom")
# p


#wilcox.test(df$pred~df$cancer_subtype)

######## anova 
library(tidyverse) #VISUALISATION
library(ggpubr) # GRAPHIQUE PRETS A LA PUBLICATION
library(rstatix) #ANALYSE STAT FACILE
set.seed(1234)

# ligne aleatoire par groupe
df %>% sample_n_by(cancer_subtype, size = 1)

# cancer subtype leverls
levels(df$cancer_subtype)

#stat descriptive

df %>%
  group_by(cancer_subtype) %>%
  get_summary_stats(pred, type = 'mean_sd')


# creer un box de pred / cancer subtype
#ggboxplot(df, x='cancer_subtype', y='pred')


#### verifier les hypotheses

#-----------> Identifier les valeurs aberante
df %>%
  group_by(cancer_subtype) %>%
  identify_outliers(pred)
# (ccl : le sous type 2 ?? des outliers) ==> test ANOVA robuste dans le package WRS2
#-----------> hypothese de normalit?? : pred de chaque sous type de cancer doit etre normal 
# Le graphique QQ plot dessine la corr??lation entre une donn??e d??finie et la distribution normale
#### en analysant les r??sidus du mod??le (QQ plot + test de normalit?? de shapiro-wilk)

# Construire le mod??le lin??aire
model  <- lm(pred ~ cancer_subtype, data = df)

# Cr??er un QQ plot des r??sidus
ggqqplot(residuals(model))
# (ccl : tous les points ne suivent pas la ligne de reference ==> on suppose une non normalit??)

# calculer le test de normalit?? de shapiro wilk
shapiro_test(residuals(model))
# (ccl : p val est significative ==> on suppose que pred n'est pas normal)

# verifier l'hypotese de la normalit?? par groupe : shapiro / groupe
df %>%
  group_by(cancer_subtype) %>%
  shapiro_test(pred)
#(ccl : le sous cancer 2 et 4 ne sont pas normal) ==>  Kruskal-Wallis

# qq plot par sous type : recommander quand l'ech > 50 car le test de shapiro devient tres sensible meme ?? un ecart mineur
ggqqplot(df, "pred", facet.by = "cancer_subtype")


#-----------> L???hypoth??se d???homog??n??it?? des variances : les variance des groupes doivent etre ??gales
# on va utiliser tjrs le graphique residuals versus prediction 
plot(model, 1)
df %>% levene_test(pred ~ cancer_subtype)
# (ccl : p = 0.04 | il ya une difference sign entre les variances des groupes | variances pas tr??s homogene) ==> test anova de welch welch_anova_test() 


#-----------> ANOVA
# * : effet d'interaction 
res.aov <- df %>% anova_test(pred ~ cancer_subtype)
res.aov
# gse = taille de la variabilit?? due au facteur ; p=pvalue ; F = valeur statistique ; DFn = degr?? de libert?? du numerateur ; DFd = DF du d??nominateur
# (ccl : 2.8% de la variation de pred pourrais etre li?? au sous type de cancer)
# (ccl : pval 0.27 : difference de pred entre les subtypes NS)


#-----------> POST HOC : faire de multiple comparaison par paires entre les groupes

pwc <- df %>% tukey_hsd(pred ~ cancer_subtype)
pwc
# (ccl : les diff des moyennes par paire ne sont pas significative )


#-----------> Rapporter

#L???ANOVA ?? un facteur a ??t?? r??alis??e pour ??valuer si le predicteur lineaire ??tait diff??rent pour les 4 
#groupes de cancer de pancreas : 1 (n = 26), 2 (n = 28) et 3 (n = 50) et 4 (n = 36)

#Les donn??es sont pr??sent??es sous forme de moyenne +/- ??cart type.
#Le predicteur lineaire n'??tait pas statistiquement significativement diff??rent entre les diff??rents groupes de tcancer, F(3, 136) = 1.301, p = 0.277, eta-carr?? g??n??ralis?? = 0,028.

# Les analyses post-hoc de Tukey ont r??v??l?? qu???aucune diff??rences inter-groupes n?????tait statistiquement significative.

# Visualisation : Boxplots avec p-values
pwc <- pwc %>% add_xy_position(x = "cancer_subtype")
ggboxplot(df, x = "cancer_subtype", y = "pred", fill = 'cancer_subtype', notch = T) +
  geom_jitter(shape=16, position=position_jitter(0.2))+
  stat_summary(fun.y=mean, geom="point", shape=23, size=4)+
  scale_color_brewer(palette="Dark2")+
  stat_pvalue_manual(pwc, hide.ns = TRUE) +
  labs(
    subtitle = get_test_label(res.aov, detailed = TRUE),
    caption = get_pwc_label(pwc)
  )


#-----------> Relaxer l???hypoth??se d???homog??n??it?? de la variance
# dans une situation ou le l'hypothese d'homogeinit?? de la variance est viol?? ==>  test ANOVA de Welch ?? 1 facteur
# Test ANOVA de Welch ?? un facteur
res.aov2 <- df %>% welch_anova_test(pred ~ cancer_subtype)
res.aov2

# Comparaisons par paires (Games-Howell)
pwc2 <- df %>% games_howell_test(pred ~ cancer_subtype)
pwc2

pwc2 <- pwc2 %>% add_xy_position(x = "cancer_subtype", step.increase = 1)
ggboxplot(df, x = "cancer_subtype", y = "pred", fill = 'cancer_subtype', notch = T) +
  geom_jitter(shape=16, position=position_jitter(0.2))+
  stat_summary(fun.y=mean, geom="point", shape=23, size=4)+
  scale_color_brewer(palette="Dark2")+
  stat_pvalue_manual(pwc2, hide.ns = TRUE) +
  labs(
    subtitle = get_test_label(res.aov, detailed = TRUE),
    caption = get_pwc_label(pwc2)
  )


# test t sans hypoth??se d?????galit?? des variances:
pwc3 <- df %>% 
  pairwise_t_test(
    pred ~ cancer_subtype, pool.sd = FALSE,
    p.adjust.method = "bonferroni"
  )
pwc3







