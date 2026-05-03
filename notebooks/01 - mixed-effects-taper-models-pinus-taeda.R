############################################################################
############################# Pacotes ######################################
############################################################################
library(ggplot2)
library(dplyr)
library(lme4)
library(nlme)
library(lmerTest)
library(cAIC4)
library(splines)
library(gamlss)
library(gridExtra)
library(RColorBrewer)
library(hnp)
library(varTestnlme)
library(car)
library(MASS) 
library(tidyverse)
library(qqplotr)
############################################################################
################################### PTA ####################################
############################################################################
DadosPTA=subset(dados,dados$Esp == "PTA")
DadosPTA2=DadosPTA[DadosPTA$DSC1_cm!=0,]
############################################################################
######################### Arranging Variables ##############################
############################################################################
DadosPTA2$DCC3=rowMeans(DadosPTA2[,c(5:6)])
DadosPTA2$DSC3=rowMeans(DadosPTA2[,c(7:8)])
DadosPTA2$Talhao=as.factor(DadosPTA2$Talhao)
DadosPTA2$RF=as.factor(DadosPTA2$RF)
DadosPTA2=DadosPTA2[DadosPTA2$RF!="Extras",]
DadosPTA2=DadosPTA2[DadosPTA2$Talhao!="X2A",]
DadosPTA2$Arvore=as.factor(DadosPTA2$Arvore)
DadosPTA2$Arvore_key = as.factor(paste(DadosPTA2$RF, DadosPTA2$Talhao, DadosPTA2$Arvore, SEP="_"))
############################################################################
#################### Model Variables #######################################
############################################################################
DadosPTA2 <- DadosPTA2 %>%
  mutate(hi_m = ifelse(hi_m == 0, 0.01, hi_m))
DadosPTA2 <- DadosPTA2 %>% mutate(DSC_DAP = DadosPTA2$DSC3/DadosPTA2$DAP_cm)
DadosPTA2 <- DadosPTA2 %>% mutate(hi_ht = DadosPTA2$hi_m/DadosPTA2$h_m)
############################################################################
############################### Removing Tree 2336 #########################
############################################################################
# When you skip these lines, you will observe the atypical profile of the trees
# removed below
DadosPTA2=DadosPTA2[DadosPTA2$Arvore_key!="BMR I7A 2336 _",]
DadosPTA2=DadosPTA2[DadosPTA2$Arvore_key!="BMO B2A 897 _",]
DadosPTA2=DadosPTA2[DadosPTA2$Arvore_key!="BMO B2A 907 _",]
DadosPTA2=DadosPTA2[DadosPTA2$Arvore_key!="BRS K9A 2082 _",]
############################################################################
################# Creating the new field variable ##########################
############################################################################
DadosPTA2$stands_new <- paste(DadosPTA2$Talhao, DadosPTA2$RF, sep = "_")
############################################################################
############# checking if there are shared stands ##########################
############################################################################
talhoes_compartilhados <- DadosPTA2 %>%
  group_by(Talhao) %>%
  filter(n_distinct(RF) > 1) %>%  
  summarise(RFs = paste(unique(RF), collapse = ", ")) %>%  
  ungroup()
print(talhoes_compartilhados)
############################################################################
############# Descriptive analysis: number of stands per FR ################
############################################################################
talhoes_por_RF <- DadosPTA2 %>%
  group_by(RF) %>%
  summarize(n_talhoes = n_distinct(stands_new)) %>%
  mutate(RF = paste0("R", row_number()))
print(talhoes_por_RF, n = 44)
table(talhoes_por_RF$n_talhoes)
talhoes_por_RF <- data.frame(
  RF = paste0("R", 1:44),
  n_talhoes = talhoes_por_RF$n_talhoes
)
talhoes_por_RF$RF <- factor(talhoes_por_RF$RF, levels = paste0("R", 1:44))
ggplot(talhoes_por_RF, aes(x = RF, y = n_talhoes)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = n_talhoes), vjust = -0.5, size = 8) + 
  scale_fill_manual(values = 'gray') +
  labs(x = "Forest Regions", y = "Number of Stands") +
  theme(
    axis.title.x = element_text(size = 30, color = "black"),
    axis.title.y = element_text(size = 30, color = "black"),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 30, color = "black"),
    plot.title = element_text(size = 45),
    legend.position = "none"
  )
result <- talhoes_por_RF %>%
  summarize(
    min_talhao = min(n_talhoes),
    max_talhao = max(n_talhoes),
    mean_talhao = mean(n_talhoes),
    var_talhao = var(n_talhoes)
  ) 
print(result)
sum(talhoes_por_RF$n_talhoes)
############################################################################
########### Descriptive analysis: number of trees per stands ###############
############################################################################
arvores_por_talhoes <- DadosPTA2 %>%
  group_by(stands_new) %>%
  summarize(n_arvores = n_distinct(Arvore_key)) 
print(arvores_por_talhoes, n = 92)
sum(arvores_por_talhoes$n_arvores)
table(arvores_por_talhoes$n_arvores)
arvores_por_talhoes <- data.frame(
  S = paste0("S", 1:92),
  n_arvores = arvores_por_talhoes$n_arvores
)
arvores_por_talhoes$S<- factor(arvores_por_talhoes$S, levels = paste0("S", 1:92))
ggplot(arvores_por_talhoes, aes(x = S, y = n_arvores)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = n_arvores), vjust = -0.5, size = 5) + 
  scale_fill_manual(values = 'gray') +
  labs(x = "Stands", y = "Number of Trees") +
  theme(
    axis.title.x = element_text(size = 30, color = "black"),
    axis.title.y = element_text(size = 30, color = "black"),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 30, color = "black"),
    plot.title = element_text(size = 45),
    legend.position = "none"
  )

result_arv <- arvores_por_talhoes %>%
  summarize(
    min_talhao = min(n_arvores),
    max_talhao = max(n_arvores),
    mean_talhao = mean(n_arvores),
    var_talhao = var(n_arvores)
  ) 
print(result_arv)
############################################################################
##################### Table: descriptive analysis ##########################
############################################################################
classes_d <- c(4, 6.6, 9.2, 11.8, 14.4, 17, 19.6, 22.2, 24.8, 27.4, 30)
classes_h <- c(4, 6.6, 9.2, 11.8, 14.4, 17, 19.6, 22.2, 24.8, 27.4, 30)
DadosPTA2 <- DadosPTA2 %>%
  mutate(
    d_class = cut(DAP_cm, breaks = classes_d, right = FALSE, include.lowest = TRUE),
    h_class = cut(h_m, breaks = classes_h, right = FALSE, include.lowest = TRUE)
  )
todas_classes <- expand.grid(
  d_class = levels(cut(DadosPTA2$DAP_cm, breaks = classes_d, right = FALSE, include.lowest = TRUE)),
  h_class = levels(cut(DadosPTA2$h_m, breaks = classes_h, right = FALSE, include.lowest = TRUE))
)
tabela_frequencia <- DadosPTA2 %>%
  distinct(Arvore_key, d_class, h_class) %>%  
  count(d_class, h_class, name = "num_arvores") %>%
  right_join(todas_classes, by = c("d_class", "h_class")) %>%
  replace_na(list(num_arvores = 0)) %>%
  pivot_wider(names_from = h_class, values_from = num_arvores, values_fill = 0) %>%
  arrange(d_class)
tabela_frequencia <- tabela_frequencia %>%
  mutate(Total = rowSums(across(-d_class))) %>%
  bind_rows(
    summarise(., d_class = "Total", across(-d_class, ~ sum(.x, na.rm = TRUE)))
  )
tabela_frequencia
############################################################################
##################### Table: descriptive analysis ##########################
############################################################################
Table <- DadosPTA2 %>%
  group_by(Classe_Idade = cut(Idade, breaks = seq(4, 20, by = 2), right = FALSE, include.lowest = TRUE)) %>%
  summarise(
    n = n_distinct(Arvore_key),          
    d_mean = mean(DAP_cm, na.rm = TRUE), 
    d_max = max(DAP_cm, na.rm = TRUE),   
    d_min = min(DAP_cm, na.rm = TRUE),   
    h_mean = mean(h_m, na.rm = TRUE),    
    h_max = max(h_m, na.rm = TRUE),     
    h_min = min(h_m, na.rm = TRUE)       
  )
print(Table)
############################################################################
##################### Graphs: descriptive analysis #########################
############################################################################
DadosPTA2 <- DadosPTA2 %>%
  mutate(Classe_Idade = cut(Idade, breaks = seq(4, 20, by = 2), right = FALSE, include.lowest = TRUE))
boxplot_d <- ggplot(DadosPTA2, aes(x = Classe_Idade, y = DAP_cm, fill = Classe_Idade)) +
  geom_boxplot() +
  labs(x = "Age Classes", y = "Diameter (cm) at breast height") +
  theme(
    axis.title.x = element_text(size = 18, color = "black"),
    axis.title.y = element_text(size = 18, color = "black"),
    axis.text.x = element_text(size = 20, angle = 45, color = "black", hjust = 1),
    axis.text.y = element_text(size = 20, color = "black"),
    plot.title = element_text(size = 45),
    legend.position = "none"
  )
violin_d <- ggplot(DadosPTA2, aes(x = Classe_Idade, y = DAP_cm, fill = Classe_Idade)) +
  geom_violin() +
  geom_boxplot(width = 0.1, color = "black", alpha = 0.7) + 
  labs(x = "Age Classes", y = "Diameter (cm) at breast height") +
  theme(
    axis.title.x = element_text(size = 18, color = "black"),
    axis.title.y = element_text(size = 18, color = "black"),
    axis.text.x = element_text(size = 20, angle = 45, color = "black", hjust = 1),
    axis.text.y = element_text(size = 20, color = "black"),
    plot.title = element_text(size = 45),
    legend.position = "none"
  )
boxplot_h <- ggplot(DadosPTA2, aes(x = Classe_Idade, y = h_m, fill = Classe_Idade)) +
  geom_boxplot() +
  labs(x = "Age Classes", y = "Total height (m)") +
  theme(
    axis.title.x = element_text(size = 18, color = "black"),
    axis.title.y = element_text(size = 18, color = "black"),
    axis.text.x = element_text(size = 20, angle = 45, color = "black", hjust = 1),
    axis.text.y = element_text(size = 20, color = "black"),
    plot.title = element_text(size = 45),
    legend.position = "none"
  )
violin_h <- ggplot(DadosPTA2, aes(x = Classe_Idade, y = h_m, fill = Classe_Idade)) +
  geom_violin() +
  geom_boxplot(width = 0.1, color = "black", alpha = 0.7) + 
  labs(x = "Age Classes", y = "Total height (m)") +
  theme(
    axis.title.x = element_text(size = 18, color = "black"),
    axis.title.y = element_text(size = 18, color = "black"),
    axis.text.x = element_text(size = 20, angle = 45, color = "black", hjust = 1),
    axis.text.y = element_text(size = 20, color = "black"),
    plot.title = element_text(size = 45),
    legend.position = "none"
  )
grid.arrange(boxplot_d, violin_d, boxplot_h, violin_h, ncol = 2)
############################################################################
##################### Graphs: descriptive analysis #########################
############################################################################
nomes_personalizados <- paste0("R", 1:44)
p11 = ggplot(DadosPTA2,
       aes(x = hi_ht,
           y = DSC_DAP, group = Arvore_key, colour = RF)) +
  geom_point(size = 0.9) +
  scale_colour_manual(values = rainbow(44), labels = nomes_personalizados,
                      guide = guide_legend(nrow = 5, byrow = TRUE, 
                                           override.aes = list(size = 3))) + 
  labs(x = expression(frac(h[ij], h[i])),
       y = expression(frac(d[ij], d[i])),
       colour = "Forest \n Regions") + 
  theme(axis.title.x = element_text(size = 30, color = "black"),
        axis.title.y = element_text(size = 30, color = "black", angle = 0, vjust = 0.5),
        axis.text.y = element_text(size = 30, color = "black"),
        axis.text.x = element_text(size = 30, color = "black"),
        plot.title = element_text(size = 45),
        strip.text.x = element_text(size = 12),
        legend.title = element_text(size = 12, hjust = 0.5),
        legend.text = element_text(size = 12),
        legend.position = "none")

DadosPTA2$DSC_DAP_2 = (DadosPTA2$DSC_DAP)^2

nomes_personalizados <- paste0("R", 1:44)
p12 = ggplot(DadosPTA2,
       aes(x = hi_ht,
           y = DSC_DAP_2, group = Arvore_key, colour = RF)) +
  geom_point(size = 0.9) +
  scale_colour_manual(values = rainbow(44), labels = nomes_personalizados,
                      guide = guide_legend(nrow = 5, byrow = TRUE, 
                                           override.aes = list(size = 3))) + 
  labs(x = 'x',
       y = 'Y^2',
       colour = "Forest \n Regions") + 
  theme(axis.title.x = element_text(size = 30, color = "black"),
        axis.title.y = element_text(size = 30, color = "black", angle = 0, vjust = 0.5),
        axis.text.y = element_text(size = 30, color = "black"),
        axis.text.x = element_text(size = 30, color = "black"),
        plot.title = element_text(size = 45),
        strip.text.x = element_text(size = 12),
        legend.title = element_text(size = 12, hjust = 0.5),
        legend.text = element_text(size = 12),
        legend.position = "none")
gridExtra::grid.arrange(p11,p12, ncol =2)

nomes_personalizados <- paste0("R", 1:44)
DadosPTA2$RF_custom <- factor(paste0("R", as.numeric(factor(DadosPTA2$RF))),
                              levels = paste0("R", 1:44))
ggplot(DadosPTA2,
       aes(x = hi_ht,
           y = DSC_DAP, group = Arvore_key, colour = RF)) +
  geom_line(size = 0.05) +
  scale_colour_manual(values = rainbow(44), labels = nomes_personalizados,
                      guide = guide_legend(nrow = 3, byrow = TRUE,
                                           override.aes = list(size = 2, linetype = 1))) +  
  labs(x = 'x',
       y = 'Y',
       colour = "Forest \n Regions") +  
  theme(axis.title.x = element_text(size = 30, color = "black"),
        axis.title.y = element_text(size = 30, color = "black", angle = 0, vjust = 0.5),
        axis.text = element_blank(),
        plot.title = element_text(size = 45),
        strip.text = element_blank(),  
        legend.title = element_text(size = 15, hjust = 0.5),
        legend.text = element_text(size = 15),
        legend.position = "none")+
  facet_wrap(~RF_custom, nrow = 11, ncol = 4)

DadosPTA2$new_variable <- DadosPTA2$DSC3 * DadosPTA2$hi_ht

nomes_personalizados <- paste0("R", 1:44)
ggplot(DadosPTA2,
       aes(x = new_variable,
           y = DSC3, group = Arvore_key, colour = RF)) +
  geom_point(size = 0.9) +
  scale_colour_manual(values = rainbow(44), labels = nomes_personalizados,
                      guide = guide_legend(nrow = 5, byrow = TRUE, 
                                           override.aes = list(size = 3))) + 
  labs(x = expression(d[i] %*% frac(h[ij], h[i])),
       y = expression(d[ij]),
       colour = "Forest \n Regions") + 
  theme(axis.title.x = element_text(size = 30, color = "black"),
        axis.title.y = element_text(size = 30, color = "black", angle = 0, vjust = 0.5),
        axis.text.y = element_text(size = 30, color = "black"),
        axis.text.x = element_text(size = 30, color = "black"),
        plot.title = element_text(size = 15),
        strip.text.x = element_text(size = 12),
        legend.title = element_text(size = 12, hjust = 0.5),
        legend.text = element_text(size = 12),
        legend.position = "none")


nomes_personalizados <- paste0("R", 1:44)
ggplot(DadosPTA2,
       aes(x = hi_ht,
           y = DSC3, group = Arvore_key, colour = RF)) +
  geom_point(size = 0.9) +
  scale_colour_manual(values = rainbow(44), labels = nomes_personalizados,
                      guide = guide_legend(nrow = 5, byrow = TRUE, 
                                           override.aes = list(size = 3))) + 
  labs(x = expression(frac(hi, h)),
       y = expression(d[j]),
       colour = "Forest \n Regions") + 
  theme(axis.title.x = element_text(size = 30, color = "black"),
        axis.title.y = element_text(size = 30, color = "black", angle = 0, vjust = 0.5),
        axis.text.y = element_text(size = 30, color = "black"),
        axis.text.x = element_text(size = 30, color = "black"),
        plot.title = element_text(size = 45),
        strip.text.x = element_text(size = 12),
        legend.title = element_text(size = 12, hjust = 0.5),
        legend.text = element_text(size = 12),
        legend.position = "none")
############################################################################
#################### Model - Fixed part based on Kozak #####################
########################## fitted mixed model and ##########################
############################################################################
# mod.without = lm(DSC3~DAP_cm*(I(hi_ht)+I(hi_ht^2)+I(hi_ht^3)+I(hi_ht^4)),
#                   DadosPTA2)
# mod.stands = lme(DSC3~DAP_cm*(I(hi_ht)+I(hi_ht^2)+I(hi_ht^3)+I(hi_ht^4)),random = ~1|stands_new, 
#                    DadosPTA2, method="REML")
# mod.stands.tree = lme(DSC3~DAP_cm*(I(hi_ht)+I(hi_ht^2)+I(hi_ht^3)+I(hi_ht^4)),random = ~1|stands_new/Arvore_key,
#                    DadosPTA2, method="REML")
# mod.tree = lme(DSC3~DAP_cm*(I(hi_ht)+I(hi_ht^2)+I(hi_ht^3)+I(hi_ht^4)),random = ~1|Arvore_key, 
#                  DadosPTA2, method="REML")
# mod.rf = lme(DSC3~DAP_cm*(I(hi_ht)+I(hi_ht^2)+I(hi_ht^3)+I(hi_ht^4)),random = ~1|RF, 
#                    DadosPTA2, method="REML")
# mod.rf.tree = lme(DSC3~DAP_cm*(I(hi_ht)+I(hi_ht^2)+I(hi_ht^3)+I(hi_ht^4)),random = ~1|RF/Arvore_key,
#                    DadosPTA2, method="REML")
# mod.rf.stands.tree = lme(DSC3~DAP_cm*(I(hi_ht)+I(hi_ht^2)+I(hi_ht^3)+I(hi_ht^4)),random = ~1|RF/stands_new/Arvore_key,
#                    DadosPTA2, method="REML")

mod.without = lm(DSC3 ~ -1 + DAP_cm*(I(hi_ht)+I(hi_ht^2)),
                 DadosPTA2)

mod.stands0 = lme(DSC3~ -1 +DAP_cm*(I(hi_ht)+I(hi_ht^2)),
                 random = ~1|stands_new, 
                 DadosPTA2, method="REML")

mod.stands.tree = lme(DSC3~ -1 +DAP_cm*(I(hi_ht)+I(hi_ht^2)),
                      random = ~1|stands_new/Arvore_key,
                      DadosPTA2, method="REML")

mod.tree = lme(DSC3~ -1 +DAP_cm*(I(hi_ht)+I(hi_ht^2)),random = ~1|Arvore_key, 
               DadosPTA2, method="REML")
mod.rf = lme(DSC3~ -1 +DAP_cm*(I(hi_ht)+I(hi_ht^2)),random = ~1|RF, 
             DadosPTA2, method="REML")
mod.rf.tree = lme(DSC3~ -1 +DAP_cm*(I(hi_ht)+I(hi_ht^2)),random = ~1|RF/Arvore_key,
                  DadosPTA2, method="REML")
mod.rf.stands.tree = lme(DSC3~DAP_cm*(I(hi_ht)+I(hi_ht^2)),random = ~1|RF/stands_new/Arvore_key,
                         DadosPTA2, method="REML")
#############################################################################
####################### create new variable #################################
############################################################################
#DadosPTA2 <- DadosPTA2 %>% dplyr::mutate(dsc_dapNew <- ifelse(DSC_DAP >= 1, 1,0)); head(DadosPTA2)
#DadosPTA2 <-DadosPTA2 %>% rename(dsc_dapNew = `dsc_dapNew <- ifelse(DSC_DAP >= 1, 1, 0)`); head(DadosPTA2)
#str(DadosPTA2$dsc_dapNew)
#DadosPTA2$dsc_dapNew = as.factor(DadosPTA2$dsc_dapNew)
mod.without = lm(DSC3 ~ -1 + DAP_cm + DAP_cm:I(hi_ht) + DAP_cm:I(hi_ht^2),
                 DadosPTA2)

mod.rf.stands.tree = lme(DSC3 ~ -1 + DAP_cm + DAP_cm:I(hi_ht) + DAP_cm:I(hi_ht^2),
                         random = ~1|RF/stands_new/Arvore_key,
                         DadosPTA2, method="REML")

mod.stands.tree = lme(DSC3 ~ -1 + DAP_cm + DAP_cm:I(hi_ht) + DAP_cm:I(hi_ht^2),
                      random = ~1|stands_new/Arvore_key,
                      DadosPTA2, method="REML")
mod.rf.tree = lme(DSC3 ~ -1 + DAP_cm + DAP_cm:I(hi_ht) + DAP_cm:I(hi_ht^2),
                  random = ~1|RF/Arvore_key,
                  DadosPTA2, method="REML")
mod.rf = lme(DSC3 ~ -1 + DAP_cm + DAP_cm:I(hi_ht) + DAP_cm:I(hi_ht^2),
             random = ~1|RF, 
             DadosPTA2, method="REML")

mod.stands = lme(DSC3 ~ -1 + DAP_cm + DAP_cm:I(hi_ht) + DAP_cm:I(hi_ht^2),
                 random = ~1|stands_new, 
                 DadosPTA2, method="REML")
mod.tree = lme(DSC3 ~ -1 + DAP_cm + DAP_cm:I(hi_ht) + DAP_cm:I(hi_ht^2),
               random = ~1|Arvore_key, 
               DadosPTA2, method="REML")


mod.het.rf.dummy <- lme(DSC3 ~ -1 + DAP_cm + DAP_cm:I(hi_ht) + DAP_cm:I(hi_ht^2),
                        random = list(RF=pdDiag(~dsc_dapNew)), 
                        DadosPTA2, method="REML")
control <- lmeControl(opt = "optim", msMaxIter = 500, msVerbose = TRUE)

mod.het.stand.dummy <- lme(DSC3 ~ -1 + DAP_cm + DAP_cm:I(hi_ht) + DAP_cm:I(hi_ht^2),
                           random = list(stands_new=pdDiag(~dsc_dapNew)), 
                           DadosPTA2, method="REML", control = control)

mod.het.tree.dummy <- lme(DSC3 ~ -1 + DAP_cm + DAP_cm:I(hi_ht) + DAP_cm:I(hi_ht^2),
                           random = list(Arvore_key=pdDiag(~dsc_dapNew)), 
                           DadosPTA2, method="REML", control = control)
####### fixed effects test #########################
mod.het.tree.dummy.0 <- lme(DSC3 ~ -1,
                          random = list(Arvore_key=pdDiag(~dsc_dapNew)), 
                          DadosPTA2, method="ML")
mod.het.tree.dummy.1 <- lme(DSC3 ~ -1 + DAP_cm,
                            random = list(Arvore_key=pdDiag(~dsc_dapNew)), 
                            DadosPTA2, method="ML")
mod.het.tree.dummy.2 <- lme(DSC3 ~ -1 + DAP_cm + DAP_cm:I(hi_ht),
                            random = list(Arvore_key=pdDiag(~dsc_dapNew)), 
                            DadosPTA2, method="ML")
mod.het.tree.dummy.3 <- lme(DSC3 ~ -1 + DAP_cm + DAP_cm:I(hi_ht) + DAP_cm:I(hi_ht^2),
                          random = list(Arvore_key=pdDiag(~dsc_dapNew)), 
                          DadosPTA2, method="ML")

anova(mod.het.tree.dummy.0, mod.het.tree.dummy.1)
anova(mod.het.tree.dummy.1, mod.het.tree.dummy.2)
anova(mod.het.tree.dummy.2, mod.het.tree.dummy.3)

logLik(mod.without)
logLik(mod.rf.stands.tree)
logLik(mod.stands.tree)
logLik(mod.rf.tree)
logLik(mod.rf)
logLik(mod.stands)
logLik(mod.tree)
logLik(mod.het.rf.dummy)
logLik(mod.het.stand.dummy)
logLik(mod.het.tree.dummy)

AIC(mod.without)
AIC(mod.rf.stands.tree)
AIC(mod.stands.tree)
AIC(mod.rf.tree)
AIC(mod.rf)
AIC(mod.stands)
AIC(mod.tree)
AIC(mod.het.rf.dummy)
AIC(mod.het.stand.dummy)
AIC(mod.het.tree.dummy)

BIC(mod.without)
BIC(mod.rf.stands.tree)
BIC(mod.stands.tree)
BIC(mod.rf.tree)
BIC(mod.rf)
BIC(mod.stands)
BIC(mod.tree)
BIC(mod.het.rf.dummy)
BIC(mod.het.stand.dummy)
BIC(mod.het.tree.dummy)

VarCorr(mod.rf)
VarCorr(mod.stands)
############################################################################
##############################  Test #######################################
############################################################################
########################## without vs. stands ##############################
############################################################################
logLik(mod.without) # -22223.77
logLik(mod.stands) # -21915.84
2*(-21915.84-(-22223.77)) # 615.86
1-0.5*pchisq(615.86,0)-0.5*pchisq(615.86,1) # p-value: 0
############################################################################
##############################  Test #######################################
############################################################################
########################## stands vs. stands/tree ##########################
############################################################################
logLik(mod.stands) # -21915.84
logLik(mod.stands.tree) # -21702.86
2*(-21702.86-(-21915.84)) # 425.96
1-0.5*pchisq(425.96,0)-0.5*pchisq(425.96,1) # p-value: 0
############################################################################
##############################  Test #######################################
############################################################################
################### stands/tree vs. rf/stands/tree #########################
############################################################################
logLik(mod.stands.tree) # -21702.86
logLik(mod.rf.stands.tree) # -21701.15
2*(-21701.15-(-21702.86)) # 3.42
1-0.5*pchisq(3.42,0)-0.5*pchisq(3.42,1) # p-value: 0.03
############################################################################
##############################  Test #######################################
############################################################################
########################## without vs. rf ##################################
############################################################################
logLik(mod.without) # -22223.77
logLik(mod.rf) # -21992.04
2*(-21992.04-(-22223.77)) # 463.46
1-0.5*pchisq(463.46,0)-0.5*pchisq(463.46,1) # p-value: 0
############################################################################
##############################  Test #######################################
############################################################################
########################## rf vs. rf/tree ##################################
############################################################################
logLik(mod.rf) # -21992.04
logLik(mod.rf.tree) # -21725.38
2*(-21725.38-(-21992.04)) # 533.32
1-0.5*pchisq(533.32,0)-0.5*pchisq(533.32,1) # p-value: 0
############################################################################
##############################  Test #######################################
############################################################################
################### rf/tree vs. rf/stands/tree #############################
############################################################################
logLik(mod.rf.tree) # -21725.38
logLik(mod.rf.stands.tree) # -21701.15
2*(-21701.15-(-21725.38)) # 48.46
1-0.5*pchisq(48.46,0)-0.5*pchisq(48.46,1) # p-value: 0
############################################################################
##############################  Test #######################################
############################################################################
########################## without vs. tree ################################
############################################################################
logLik(mod.without) # -22223.77
logLik(mod.tree) # -21798.09
2*(-21798.09-(-22223.77)) # 851.36
1-0.5*pchisq(851.36,0)-0.5*pchisq(851.36,1) # p-value: 0
#############################################################################
####################### create new variable #################################
############################################################################
DadosPTA2 <- DadosPTA2 %>% dplyr::mutate(dsc_dapNew <- ifelse(DSC_DAP >= 1, 1,0)); head(DadosPTA2)
DadosPTA2 <-DadosPTA2 %>% rename(dsc_dapNew = `dsc_dapNew <- ifelse(DSC_DAP >= 1, 1, 0)`); head(DadosPTA2)
str(DadosPTA2$dsc_dapNew)
DadosPTA2$dsc_dapNew = as.factor(DadosPTA2$dsc_dapNew)
############################################################################
################## model with heterogeneity  ###############################
############################################################################
mod.het.rf.dummy <- lme(DSC3~DAP_cm*(I(hi_ht)+I(hi_ht^2)),
            random = list(RF=pdDiag(~dsc_dapNew)), 
            DadosPTA2, method="REML")
control <- lmeControl(opt = "optim", msMaxIter = 500, msVerbose = TRUE)
mod.het.stand.dummy <- lme(DSC3~DAP_cm*(I(hi_ht)+I(hi_ht^2)),
              random = list(stands_new=pdDiag(~dsc_dapNew)), 
              DadosPTA2, method="REML", control = control)

mod.het.tree.dummy <- lme(DSC3~DAP_cm*(I(hi_ht)+I(hi_ht^2)),
              random = list(Arvore_key=pdDiag(~dsc_dapNew)), 
              DadosPTA2, method="REML")
VarCorr(mod.het.tree.dummy)
ranef(mod.het.rf.dummy)


p00 = DadosPTA2 %>%
  mutate(fitted = fitted(mod.het.rf.dummy),
         resid = rstudent(mod.het.rf.dummy)) %>%
  ggplot(aes(fitted, resid)) +
  coord_cartesian(xlim=c(0,50), ylim=c(-7, 12)) +  
  ggtitle("M1.2") +
  geom_hline(yintercept = 0, colour="red", size=1.5) +
  geom_point(size=1.2) +
  labs(x = "Fitted Values", y = "Least Confounded Residuals") + 
  theme(legend.position = "none", legend.title = element_text(size = 15),
        legend.text = element_text(size = 25),
        axis.title = element_text(size = 25),
        axis.text.x = element_text(color = "black", hjust=1),
        axis.text.y = element_text(color = "black", hjust=1),
        axis.text = element_text(size = 25),
        plot.title = element_text(size = 25, hjust = 0.5))   

p11 = DadosPTA2 %>%
  mutate(fitted = fitted(mod.stands),
         resid = residuals(mod.stands, type = "pearson")) %>%
  ggplot(aes(fitted, resid)) +
  coord_cartesian(xlim=c(0,50), ylim=c(-7, 12)) + 
  ggtitle("M6.2") +
  geom_hline(yintercept = 0, colour="red", size=1.5) +
  geom_point(size=1.2) +
  labs(x = "Fitted Values", y = "Least Confounded Residuals") + 
  theme(legend.position = "none", legend.title = element_text(size = 15),
        legend.text = element_text(size = 25),
        axis.title = element_text(size = 25),
        axis.text.x = element_text(color = "black", hjust=1),
        axis.text.y = element_text(color = "black", hjust=1),
        axis.text = element_text(size = 25),
        plot.title = element_text(size = 25, hjust = 0.5))      

p22 = DadosPTA2 %>%
  mutate(fitted = fitted(mod.stands.tree),
         resid = residuals(mod.stands.tree, type = "pearson")) %>%
  ggplot(aes(fitted, resid)) +
  coord_cartesian(xlim=c(0,50), ylim=c(-7, 12)) + 
  ggtitle("M3.2") +
  geom_hline(yintercept = 0, colour="red", size=1.5) +
  geom_point(size=1.2) +
  labs(x = "Fitted Values", y = "Least Confounded Residuals") + 
  theme(legend.position = "none", legend.title = element_text(size = 15),
        legend.text = element_text(size = 25),
        axis.title = element_text(size = 25),
        axis.text.x = element_text(color = "black", hjust=1),
        axis.text.y = element_text(color = "black", hjust=1),
        axis.text = element_text(size = 25),
        plot.title = element_text(size = 25, hjust = 0.5))    

############################################################################
##############################  Test #######################################
############################################################################
################# Tree vs. tree with heterogeneity #########################
############################################################################
logLik(mod.tree)
logLik(mod.het.tree.dummy)
2*(-19367.53-(-21798.09)) # 4861.12
1-0.5*pchisq(4861.12,0)-0.5*pchisq(4861.12,1) # 0


dados_arvore <- DadosPTA2[DadosPTA2$Arvore_key == "CPI F8A 394 _", ]
pred_het <- predict(mod.het.tree.dummy, newdata = dados_arvore)
pred_without <- predict(mod.without, newdata = dados_arvore)

ggplot(dados_arvore, aes(x = hi_m, y = DSC3, group = Arvore_key)) +
  geom_line(aes(colour = "Observed profile", linetype = "Observed profile"), size = 1.5) +
  geom_line(aes(y = pred_het, colour = "M10.2", linetype = "M10.2"), size = 1.5) +
  geom_line(aes(y = -1 * pred_het, colour = "M10.2", linetype = "M10.2"), size = 1.5) +
  geom_line(aes(y = pred_without, colour = "M7.2", linetype = "M7.2"), size = 1.5) +
  geom_line(aes(y = -1 * pred_without, colour = "M7.2", linetype = "M7.2"), size = 1.5) +
  geom_line(aes(y = -1 * DSC3, colour = "Observed profile", linetype = "Observed profile"), size = 1.5) +
  geom_hline(yintercept = 0, colour = "black", size = 1, lty = 2) + 
  geom_vline(xintercept = 0, colour = "black", size = 1, lty = 2) +
  labs(
    x = "Height (m)", 
    y = "Diameter (cm)",         
    colour = "Taper curves",
    linetype = "Taper curves"
  ) +
  coord_flip(xlim = c(0, 30), ylim = c(-50, 50)) +
  scale_y_reverse() +
  scale_y_continuous(
    breaks = c(-50, 50),        
    labels = c("0", "50")        
  ) +
  ggtitle("CPI F8A 394 _") + 
  theme(
    legend.position = c(0, 1), 
    legend.justification = c(0, 1),
    legend.title = element_text(size = 15),
    legend.text = element_text(size = 15),
    axis.title = element_text(size = 20, color = "black"),
    axis.text = element_text(size = 20, color = "black"),
    axis.title.y = element_text(size = 20, color = "black", angle = 90),
    plot.title = element_text(size = 20, hjust = 0.5),
    strip.text.x = element_text(size = 15)
  ) +
  scale_colour_manual(name = "Taper curves", 
                      values = c("Observed profile" = "blue", "M10.2" = "red", "M7.2" = "black")) +
  scale_linetype_manual(name = "Taper curves", 
                        values = c("Observed profile" = "dashed", "M10.2" = "solid", "M7.2" = "solid"))


dados_arvore <- DadosPTA2[DadosPTA2$Arvore_key == "CVL N6A 1422 _", ]
pred_het <- predict(mod.het.tree.dummy, newdata = dados_arvore)
pred_without <- predict(mod.without, newdata = dados_arvore)


ggplot(dados_arvore, aes(x = hi_m, y = DSC3, group = Arvore_key)) +
  geom_line(aes(colour = "Observed profile", linetype = "Observed profile"), size = 1.5) +
  geom_line(aes(y = pred_het, colour = "M10.2", linetype = "M10.2"), size = 1.5) +
  geom_line(aes(y = -1 * pred_het, colour = "M10.2", linetype = "M10.2"), size = 1.5) +
  geom_line(aes(y = pred_without, colour = "M7.2", linetype = "M7.2"), size = 1.5) +
  geom_line(aes(y = -1 * pred_without, colour = "M7.2", linetype = "M7.2"), size = 1.5) +
  geom_line(aes(y = -1 * DSC3, colour = "Observed profile", linetype = "Observed profile"), size = 1.5) +
  geom_hline(yintercept = 0, colour = "black", size = 1, lty = 2) + 
  geom_vline(xintercept = 0, colour = "black", size = 1, lty = 2) +
  labs(
    x = "Height (m)", 
    y = "Diameter (cm)",         
    colour = "Taper curves",
    linetype = "Taper curves"
  ) +
  coord_flip(xlim = c(0, 30), ylim = c(-50, 50)) +
  scale_y_reverse() +
  scale_y_continuous(
    breaks = c(-50, 50),          
    labels = c("0", "50")         
  ) +
  ggtitle("CVL N6A 1422 _") + 
  theme(
    legend.position = c(0, 1), 
    legend.justification = c(0, 1),
    legend.title = element_text(size = 15),
    legend.text = element_text(size = 15),
    axis.title = element_text(size = 20, color = "black"),
    axis.text = element_text(size = 20, color = "black"),
    axis.title.y = element_text(size = 20, color = "black", angle = 90),
    plot.title = element_text(size = 20, hjust = 0.5),
    strip.text.x = element_text(size = 15)
  ) +
  scale_colour_manual(name = "Taper curves", 
                      values = c("Observed profile" = "blue", "M10.2" = "red", "M7.2" = "black")) +
  scale_linetype_manual(name = "Taper curves", 
                        values = c("Observed profile" = "dashed", "M10.2" = "solid", "M7.2" = "solid"))


dados_arvore <- DadosPTA2[DadosPTA2$Arvore_key == "CMN R2A 1771 _", ]
pred_het <- predict(mod.het.tree.dummy, newdata = dados_arvore)
pred_without <- predict(mod.without, newdata = dados_arvore)

ggplot(dados_arvore, aes(x = hi_m, y = DSC3, group = Arvore_key)) +
  geom_line(aes(colour = "Observed profile", linetype = "Observed profile"), size = 1.5) +
  geom_line(aes(y = pred_het, colour = "M10.2", linetype = "M10.2"), size = 1.5) +
  geom_line(aes(y = -1 * pred_het, colour = "M10.2", linetype = "M10.2"), size = 1.5) +
  geom_line(aes(y = pred_without, colour = "M7.2", linetype = "M7.2"), size = 1.5) +
  geom_line(aes(y = -1 * pred_without, colour = "M7.2", linetype = "M7.2"), size = 1.5) +
  geom_line(aes(y = -1 * DSC3, colour = "Observed profile", linetype = "Observed profile"), size = 1.5) +
  geom_hline(yintercept = 0, colour = "black", size = 1, lty = 2) + 
  geom_vline(xintercept = 0, colour = "black", size = 1, lty = 2) +
  labs(
    x = "Height (m)", 
    y = "Diameter (cm)",         
    colour = "Taper curves",
    linetype = "Taper curves"
  ) +
  coord_flip(xlim = c(0, 30), ylim = c(-50, 50)) +
  scale_y_reverse() +
  scale_y_continuous(
    breaks = c(-50, 50),          
    labels = c("0", "50")         
  ) +
  ggtitle("CMN R2A 1771 _") + 
  theme(
    legend.position = c(0, 1), 
    legend.justification = c(0, 1),
    legend.title = element_text(size = 15),
    legend.text = element_text(size = 15),
    axis.title = element_text(size = 20, color = "black"),
    axis.text = element_text(size = 20, color = "black"),
    axis.title.y = element_text(size = 20, color = "black", angle = 90),
    plot.title = element_text(size = 20, hjust = 0.5),
    strip.text.x = element_text(size = 15)
  ) +
  scale_colour_manual(name = "Taper curves", 
                      values = c("Observed profile" = "blue", "M10.2" = "red", "M7.2" = "black")) +
  scale_linetype_manual(name = "Taper curves", 
                        values = c("Observed profile" = "dashed", "M10.2" = "solid", "M7.2" = "solid"))



dados_arvore <- DadosPTA2[DadosPTA2$Arvore_key == "DAO A4A 2573 _", ]
pred_het <- predict(mod.het.tree.dummy, newdata = dados_arvore)
pred_without <- predict(mod.without, newdata = dados_arvore)


ggplot(dados_arvore, aes(x = hi_m, y = DSC3, group = Arvore_key)) +
  geom_line(aes(colour = "Observed profile", linetype = "Observed profile"), size = 1.5) +
  geom_line(aes(y = pred_het, colour = "M10.2", linetype = "M10.2"), size = 1.5) +
  geom_line(aes(y = -1 * pred_het, colour = "M10.2", linetype = "M10.2"), size = 1.5) +
  geom_line(aes(y = pred_without, colour = "M7.2", linetype = "M7.2"), size = 1.5) +
  geom_line(aes(y = -1 * pred_without, colour = "M7.2", linetype = "M7.2"), size = 1.5) +
  geom_line(aes(y = -1 * DSC3, colour = "Observed profile", linetype = "Observed profile"), size = 1.5) +
  geom_hline(yintercept = 0, colour = "black", size = 1, lty = 2) + 
  geom_vline(xintercept = 0, colour = "black", size = 1, lty = 2) +
  labs(
    x = "Height (m)", 
    y = "Diameter (cm)",         
    colour = "Taper curves",
    linetype = "Taper curves"
  ) +
  coord_flip(xlim = c(0, 30), ylim = c(-50, 50)) +
  scale_y_reverse() +
  scale_y_continuous(
    breaks = c(-50, 50),          
    labels = c("0", "50")        
  ) +
  ggtitle("DAO A4A 2573 _") + 
  theme(
    legend.position = c(0, 1), 
    legend.justification = c(0, 1),
    legend.title = element_text(size = 15),
    legend.text = element_text(size = 15),
    axis.title = element_text(size = 20, color = "black"),
    axis.text = element_text(size = 20, color = "black"),
    axis.title.y = element_text(size = 20, color = "black", angle = 90),
    plot.title = element_text(size = 20, hjust = 0.5),
    strip.text.x = element_text(size = 15)
  ) +
  scale_colour_manual(name = "Taper curves", 
                      values = c("Observed profile" = "blue", "M10.2" = "red", "M7.2" = "black")) +
  scale_linetype_manual(name = "Taper curves", 
                        values = c("Observed profile" = "dashed", "M10.2" = "solid", "M7.2" = "solid"))


dados_arvore <- DadosPTA2[DadosPTA2$Arvore_key == "CLA A8B 439 _", ]
pred_het <- predict(mod.het.tree.dummy, newdata = dados_arvore)
pred_without <- predict(mod.without, newdata = dados_arvore)


ggplot(dados_arvore, aes(x = hi_m, y = DSC3, group = Arvore_key)) +
  geom_line(aes(colour = "Observed profile", linetype = "Observed profile"), size = 1.5) +
  geom_line(aes(y = pred_het, colour = "M10.2", linetype = "M10.2"), size = 1.5) +
  geom_line(aes(y = -1 * pred_het, colour = "M10.2", linetype = "M10.2"), size = 1.5) +
  geom_line(aes(y = pred_without, colour = "M7.2", linetype = "M7.2"), size = 1.5) +
  geom_line(aes(y = -1 * pred_without, colour = "M7.2", linetype = "M7.2"), size = 1.5) +
  geom_line(aes(y = -1 * DSC3, colour = "Observed profile", linetype = "Observed profile"), size = 1.5) +
  geom_hline(yintercept = 0, colour = "black", size = 1, lty = 2) + 
  geom_vline(xintercept = 0, colour = "black", size = 1, lty = 2) +
  labs(
    x = "Height (m)", 
    y = "Diameter (cm)",         
    colour = "Taper curves",
    linetype = "Taper curves"
  ) +
  coord_flip(xlim = c(0, 30), ylim = c(-50, 50)) +
  scale_y_reverse() +
  scale_y_continuous(
    breaks = c(-50, 50),          
    labels = c("0", "50")        
  ) +
  ggtitle("CLA A8B 439 _") + 
  theme(
    legend.position = c(0, 1), 
    legend.justification = c(0, 1),
    legend.title = element_text(size = 15),
    legend.text = element_text(size = 15),
    axis.title = element_text(size = 20, color = "black"),
    axis.text = element_text(size = 20, color = "black"),
    axis.title.y = element_text(size = 20, color = "black", angle = 90),
    plot.title = element_text(size = 20, hjust = 0.5),
    strip.text.x = element_text(size = 15)
  ) +
  scale_colour_manual(name = "Taper curves", 
                      values = c("Observed profile" = "blue", "M10.2" = "red", "M7.2" = "black")) +
  scale_linetype_manual(name = "Taper curves", 
                        values = c("Observed profile" = "dashed", "M10.2" = "solid", "M7.2" = "solid"))


dados_arvore <- DadosPTA2[DadosPTA2$Arvore_key == "BAG A6A 154 _", ]
pred_het <- predict(mod.het.tree.dummy, newdata = dados_arvore)
pred_without <- predict(mod.without, newdata = dados_arvore)

ggplot(dados_arvore, aes(x = hi_m, y = DSC3, group = Arvore_key)) +
  geom_line(aes(colour = "Observed profile", linetype = "Observed profile"), size = 1.5) +
  geom_line(aes(y = pred_het, colour = "M10.2", linetype = "M10.2"), size = 1.5) +
  geom_line(aes(y = -1 * pred_het, colour = "M10.2", linetype = "M10.2"), size = 1.5) +
  geom_line(aes(y = pred_without, colour = "M7.2", linetype = "M7.2"), size = 1.5) +
  geom_line(aes(y = -1 * pred_without, colour = "M7.2", linetype = "M7.2"), size = 1.5) +
  geom_line(aes(y = -1 * DSC3, colour = "Observed profile", linetype = "Observed profile"), size = 1.5) +
  geom_hline(yintercept = 0, colour = "black", size = 1, lty = 2) + 
  geom_vline(xintercept = 0, colour = "black", size = 1, lty = 2) +
  labs(
    x = "Height (m)", 
    y = "Diameter (cm)",         
    colour = "Taper curves",
    linetype = "Taper curves"
  ) +
  coord_flip(xlim = c(0, 30), ylim = c(-50, 50)) +
  scale_y_reverse() +
  scale_y_continuous(
    breaks = c(-50, 50),          
    labels = c("0", "50")        
  ) +
  ggtitle("BAG A6A 154 _") + 
  theme(
    legend.position = c(0, 1), 
    legend.justification = c(0, 1),
    legend.title = element_text(size = 15),
    legend.text = element_text(size = 15),
    axis.title = element_text(size = 20, color = "black"),
    axis.text = element_text(size = 20, color = "black"),
    axis.title.y = element_text(size = 20, color = "black", angle = 90),
    plot.title = element_text(size = 20, hjust = 0.5),
    strip.text.x = element_text(size = 15)
  ) +
  scale_colour_manual(name = "Taper curves", 
                      values = c("Observed profile" = "blue", "M10.2" = "red", "M7.2" = "black")) +
  scale_linetype_manual(name = "Taper curves", 
                        values = c("Observed profile" = "dashed", "M10.2" = "solid", "M7.2" = "solid"))



dados_arvore <- DadosPTA2[DadosPTA2$Arvore_key == "CCR A1A 2048 _", ]
pred_het <- predict(mod.het.tree.dummy, newdata = dados_arvore)
pred_without <- predict(mod.without, newdata = dados_arvore)


ggplot(dados_arvore, aes(x = hi_m, y = DSC3, group = Arvore_key)) +
  geom_line(aes(colour = "Observed profile", linetype = "Observed profile"), size = 1.5) +
  geom_line(aes(y = pred_het, colour = "M10.2", linetype = "M10.2"), size = 1.5) +
  geom_line(aes(y = -1 * pred_het, colour = "M10.2", linetype = "M10.2"), size = 1.5) +
  geom_line(aes(y = pred_without, colour = "M7.2", linetype = "M7.2"), size = 1.5) +
  geom_line(aes(y = -1 * pred_without, colour = "M7.2", linetype = "M7.2"), size = 1.5) +
  geom_line(aes(y = -1 * DSC3, colour = "Observed profile", linetype = "Observed profile"), size = 1.5) +
  geom_hline(yintercept = 0, colour = "black", size = 1, lty = 2) + 
  geom_vline(xintercept = 0, colour = "black", size = 1, lty = 2) +
  labs(
    x = "Height (m)", 
    y = "Diameter (cm)",         
    colour = "Taper curves",
    linetype = "Taper curves"
  ) +
  coord_flip(xlim = c(0, 30), ylim = c(-50, 50)) +
  scale_y_reverse() +
  scale_y_continuous(
    breaks = c(-50, 50),          
    labels = c("0", "50")        
  ) +
  ggtitle("CCR A1A 2048 _") + 
  theme(
    legend.position = c(0, 1), 
    legend.justification = c(0, 1),
    legend.title = element_text(size = 15),
    legend.text = element_text(size = 15),
    axis.title = element_text(size = 20, color = "black"),
    axis.text = element_text(size = 20, color = "black"),
    axis.title.y = element_text(size = 20, color = "black", angle = 90),
    plot.title = element_text(size = 20, hjust = 0.5),
    strip.text.x = element_text(size = 15)
  ) +
  scale_colour_manual(name = "Taper curves", 
                      values = c("Observed profile" = "blue", "M10.2" = "red", "M7.2" = "black")) +
  scale_linetype_manual(name = "Taper curves", 
                        values = c("Observed profile" = "dashed", "M10.2" = "solid", "M7.2" = "solid"))



dados_arvore <- DadosPTA2[DadosPTA2$Arvore_key == "ARE N3A 862 _", ]
pred_het <- predict(mod.het.tree.dummy, newdata = dados_arvore)
pred_without <- predict(mod.without, newdata = dados_arvore)

ggplot(dados_arvore, aes(x = hi_m, y = DSC3, group = Arvore_key)) +
  geom_line(aes(colour = "Observed profile", linetype = "Observed profile"), size = 1.5) +
  geom_line(aes(y = pred_het, colour = "M10.2", linetype = "M10.2"), size = 1.5) +
  geom_line(aes(y = -1 * pred_het, colour = "M10.2", linetype = "M10.2"), size = 1.5) +
  geom_line(aes(y = pred_without, colour = "M7.2", linetype = "M7.2"), size = 1.5) +
  geom_line(aes(y = -1 * pred_without, colour = "M7.2", linetype = "M7.2"), size = 1.5) +
  geom_line(aes(y = -1 * DSC3, colour = "Observed profile", linetype = "Observed profile"), size = 1.5) +
  geom_hline(yintercept = 0, colour = "black", size = 1, lty = 2) + 
  geom_vline(xintercept = 0, colour = "black", size = 1, lty = 2) +
  labs(
    x = "Height (m)", 
    y = "Diameter (cm)",         
    colour = "Taper curves",
    linetype = "Taper curves"
  ) +
  coord_flip(xlim = c(0, 30), ylim = c(-50, 50)) +
  scale_y_reverse() +
  scale_y_continuous(
    breaks = c(-50, 50),          
    labels = c("0", "50")        
  ) +
  ggtitle("ARE N3A 862 _") + 
  theme(
    legend.position = c(0, 1), 
    legend.justification = c(0, 1),
    legend.title = element_text(size = 15),
    legend.text = element_text(size = 15),
    axis.title = element_text(size = 20, color = "black"),
    axis.text = element_text(size = 20, color = "black"),
    axis.title.y = element_text(size = 20, color = "black", angle = 90),
    plot.title = element_text(size = 20, hjust = 0.5),
    strip.text.x = element_text(size = 15)
  ) +
  scale_colour_manual(name = "Taper curves", 
                      values = c("Observed profile" = "blue", "M10.2" = "red", "M7.2" = "black")) +
  scale_linetype_manual(name = "Taper curves", 
                        values = c("Observed profile" = "dashed", "M10.2" = "solid", "M7.2" = "solid"))



dados_arvore <- DadosPTA2[DadosPTA2$Arvore_key == "BFZ H5A 1821 _", ]
pred_het <- predict(mod.het.tree.dummy, newdata = dados_arvore)
pred_without <- predict(mod.without, newdata = dados_arvore)

ggplot(dados_arvore, aes(x = hi_m, y = DSC3, group = Arvore_key)) +
  geom_line(aes(colour = "Observed profile", linetype = "Observed profile"), size = 1.5) +
  geom_line(aes(y = pred_het, colour = "M10.2", linetype = "M10.2"), size = 1.5) +
  geom_line(aes(y = -1 * pred_het, colour = "M10.2", linetype = "M10.2"), size = 1.5) +
  geom_line(aes(y = pred_without, colour = "M7.2", linetype = "M7.2"), size = 1.5) +
  geom_line(aes(y = -1 * pred_without, colour = "M7.2", linetype = "M7.2"), size = 1.5) +
  geom_line(aes(y = -1 * DSC3, colour = "Observed profile", linetype = "Observed profile"), size = 1.5) +
  geom_hline(yintercept = 0, colour = "black", size = 1, lty = 2) + 
  geom_vline(xintercept = 0, colour = "black", size = 1, lty = 2) +
  labs(
    x = "Height (m)", 
    y = "Diameter (cm)",         
    colour = "Taper curves",
    linetype = "Taper curves"
  ) +
  coord_flip(xlim = c(0, 30), ylim = c(-50, 50)) +
  scale_y_reverse() +
  scale_y_continuous(
    breaks = c(-50, 50),          
    labels = c("0", "50")        
  ) +
  ggtitle("BFZ H5A 1821 _") + 
  theme(
    legend.position = c(0, 1), 
    legend.justification = c(0, 1),
    legend.title = element_text(size = 15),
    legend.text = element_text(size = 15),
    axis.title = element_text(size = 20, color = "black"),
    axis.text = element_text(size = 20, color = "black"),
    axis.title.y = element_text(size = 20, color = "black", angle = 90),
    plot.title = element_text(size = 20, hjust = 0.5),
    strip.text.x = element_text(size = 15)
  ) +
  scale_colour_manual(name = "Taper curves", 
                      values = c("Observed profile" = "blue", "M10.2" = "red", "M7.2" = "black")) +
  scale_linetype_manual(name = "Taper curves", 
                        values = c("Observed profile" = "dashed", "M10.2" = "solid", "M7.2" = "solid"))


############################################################################
############################## Plot ########################################
############################################################################
p7 = DadosPTA2 %>%
  mutate(fitted = fitted(mod.het.rf.dummy),
         resid = residuals(mod.het.rf.dummy, type = "pearson")) %>%
  ggplot(aes(fitted, resid)) +
  coord_cartesian(xlim=c(0,50), ylim=c(-7, 12)) +  
  ggtitle("M8.2") +
  geom_hline(yintercept = 0, colour="red", size=1.5) +
  geom_point(size=1.2) +
  labs(x = "Fitted Values", y = "Least Confounded Residuals") + 
  theme(legend.position = "none", legend.title = element_text(size = 15),
        legend.text = element_text(size = 25),
        axis.title = element_text(size = 25),
        axis.text.x = element_text(color = "black", hjust=1),
        axis.text.y = element_text(color = "black", hjust=1),
        axis.text = element_text(size = 25),
        plot.title = element_text(size = 25, hjust = 0.5))    
p8 = DadosPTA2 %>%
  mutate(fitted = fitted(mod.het.stand.dummy),
         resid = residuals(mod.het.stand.dummy, type = "pearson")) %>%
  ggplot(aes(fitted, resid)) +
  coord_cartesian(xlim=c(0,50), ylim=c(-7, 12)) +   
  ggtitle("M9.2") +
  geom_hline(yintercept = 0, colour="red", size=1.5) +
  geom_point(size=1.2) +
  labs(x = "Fitted Values", y = "Least Confounded Residuals") + 
  theme(legend.position = "none", legend.title = element_text(size = 15),
        legend.text = element_text(size = 25),
        axis.title = element_text(size = 25),
        axis.text.x = element_text(color = "black", hjust=1),
        axis.text.y = element_text(color = "black", hjust=1),
        axis.text = element_text(size = 25),
        plot.title = element_text(size = 25, hjust = 0.5))   
p9 = DadosPTA2 %>%
  mutate(fitted = fitted(mod.het.tree.dummy),
         resid = residuals(mod.het.tree.dummy, type = "pearson")) %>%
  ggplot(aes(fitted, resid)) +
  coord_cartesian(xlim=c(0,50), ylim=c(-7, 12)) + 
  ggtitle("M10.2") +
  geom_hline(yintercept = 0, colour="red", size=1.5) +
  geom_point(size=1.2) +
  labs(x = "Fitted Values", y = "Least Confounded Residuals") + 
  theme(legend.position = "none", legend.title = element_text(size = 15),
        legend.text = element_text(size = 25),
        axis.title = element_text(size = 25),
        axis.text.x = element_text(color = "black", hjust=1),
        axis.text.y = element_text(color = "black", hjust=1),
        axis.text = element_text(size = 25),
        plot.title = element_text(size = 25, hjust = 0.5))      
x11()
gridExtra::grid.arrange(p7,p8,p9,
                        ncol = 4)
############################################################################
############################## AIC e BIC ###################################
############################################################################
AIC(mod.without);BIC(mod.without) # 44461.54
AIC(mod.stands);BIC(mod.stands) # 43847.68
AIC(mod.stands.tree);BIC(mod.stands.tree) # 43423.72
AIC(mod.rf);BIC(mod.rf) # 44000.08
AIC(mod.tree);BIC(mod.tree) # 44000.08
AIC(mod.rf.tree);BIC(mod.rf.tree) # 43468.75
AIC(mod.rf.stands.tree);BIC(mod.rf.stands.tree) # 43422.31
AIC(mod.het.rf.dummy);BIC(mod.het.rf.dummy) # 39538.82
AIC(mod.het.stand.dummy);BIC(mod.het.stand.dummy) # 39203.92
AIC(mod.het.tree.dummy);BIC(mod.het.tree.dummy) # 38753.05

logLik(mod.het.rf.dummy)
logLik(mod.het.stand.dummy)
logLik(mod.het.tree.dummy)
