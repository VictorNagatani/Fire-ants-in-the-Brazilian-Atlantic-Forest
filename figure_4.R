########## glm por décadas


setwd("C:/Users/vhnag/OneDrive/Desktop/Doutorado Naga/Cap 1 arquivos/Exclusão replicas")




library(tidyr)
library(dplyr)
library(lme4)
  library(car)
library(multcomp)
library(MASS)
library(ggplot2)
library(dplyr)


### INVICTA

dados_invicta <- read.table("decada_invicta.txt", header = TRUE, sep = "\t")
dados_invicta$decada <- as.factor(dados_invicta$decada)

# Analysis
modelo1 <- lmer(ocorrência ~ decada + (1 | id_hex), data = dados_invicta, REML = FALSE)

summary(modelo1)
Anova(modelo1)

tukey_test <- glht(modelo1, linfct = mcp(decada = "Tukey"))
summary(tukey_test)


# GRAPH
dados_invicta$decada <- as.factor(dados_invicta$decada)
dados_filtrados <- subset(dados_invicta, ocorrência <= 500)
dados_filtrados$decada <- as.factor(dados_filtrados$decada)
estatisticas <- dados_filtrados %>% group_by(decada) %>% summarise(media = mean(ocorrência),desvio = sd(ocorrência),.groups = "drop")

letras <- c("A", "A", "A", "AB", "B")

p <- ggplot() + geom_bar(data = estatisticas, aes(x = decada, y = media), 
           stat = "identity", fill = "#FFD700", color = "#FFD700", width = 0.6, alpha = 0.4) +
  geom_errorbar(data = estatisticas, aes(x = decada, ymin = media - desvio, ymax = media + desvio), width = 0.2, color = "black", size = 1, alpha = 0.7) +
  geom_point(data = dados_filtrados, aes(x = decada, y = ocorrência), size = 2, alpha = 0.7, color = "#FFD700") +
  geom_line(data = dados_filtrados, aes(x = decada, y = ocorrência, group = id_hex), color = "#FFD700", alpha = 0.2, size = 0.8) +  
  geom_text(data = estatisticas, aes(x = decada, y = media + desvio + 2, label = letras), color = "black", size = 5) +  
  labs(x = "Década", y = "Ocorrência") + theme_minimal() +
  theme(legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(color = "black")) + 
  coord_cartesian(ylim = c(0, 20), xlim = c(0.5, length(unique(dados_filtrados$decada)) + 1)) 

print(p)
ggsave("Fig_4A.pdf", width = 8, height = 9)

#### SAEVISSIMA

dados_saevissima <- read.table("decada_saevissima.txt", header = TRUE, sep = "\t")
dados_saevissima$decada <- as.factor(dados_saevissima$decada)

# Analysis
modelo <- lmer(ocorrência ~ decada + (1 | id_hex), data = dados_saevissima, REML = FALSE)
summary(modelo)
Anova(modelo)

tukey_test <- glht(modelo, linfct = mcp(decada = "Tukey"))
summary(tukey_test)


###GRAPH

dados_filtrados <- subset(dados_saevissima, ocorrência <= 700)
estatisticas <- dados_filtrados %>% group_by(decada) %>% summarise(media = mean(ocorrência), desvio = sd(ocorrência))

p <- ggplot() +geom_bar(data = estatisticas, aes(x = decada, y = media), 
           stat = "identity", fill = "#483d8b", color = "#483d8b", width = 0.6, alpha = 0.4) + geom_errorbar(data = estatisticas, aes(x = decada, ymin = media - desvio, ymax = media + desvio),
                width = 0.2, color = "black", size = 1, alpha = 0.7) + geom_point(data = dados_filtrados, aes(x = decada, y = ocorrência), size = 2, alpha = 0.7, color = "#483d8b") +
 geom_line(data = dados_filtrados, aes(x = decada, y = ocorrência, group = id_hex), color = "#483d8b", alpha = 0.2, size = 0.8) +  
  geom_text(data = data.frame( decada = factor(c(1985, 1995, 2005, 2015, 2024), levels = c(1985, 1995, 2005, 2015, 2024)),letra = c("A", "A", "A", "A", "B")), aes(x = decada, y = 65, label = letra), 
    color = "black", size = 5, fontface = "plain", vjust = 0.5) +labs(x = "Década", y = "Ocorrência") +theme_minimal() +
  theme(legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(color = "black")) +  
  coord_cartesian(ylim = c(0, 70), xlim = c(0.5, length(unique(dados_filtrados$decada)) + 1))  

print(p)

ggsave("Fig_4B.pdf", width = 8, height = 9)

