#### ANALISES MODELOS PAISAGEM 

setwd("C:/Users/vhnag/OneDrive/Desktop/Doutorado Naga/Cap 1 arquivos/Exclusão replicas")

library (scales)
library(visdat)
library(tidyverse)
library(lattice)
library(RVAideMemoire)
library(DHARMa)
library(performance)
library(MuMIn)
library(piecewiseSEM)
library(MASS)
library(ggExtra)
library(Rmisc)
library(emmeans) 
library(sjPlot)
library(bbmle)
library(glmmTMB)
library(ordinal)
library(car)
library(ecolottery)
library(naniar)
library(vcd)
library(generalhoslem)
library(GGally)
library(ggthemes)
library(stringi)
library(piecewiseSEM)
library(AICcmodavg)
library(MuMIn)


#### Data
caminho_para_arquivo<- "dados_invicta_atual.txt"
caminho_para_arquivo1 <- "dados_saevissima_atual.txt"

dados_df <- read.table(caminho_para_arquivo, header = TRUE, sep = "\t")
dados_df1 <- read.table(caminho_para_arquivo1, header = TRUE, sep = "\t")


##INVICTA######


### S. invicta x Number of fragments of native vegetation 1985

m0 <- glm(dados_df[,"invicta_1985"]~1)
m1 <- glm(dados_df[,"invicta_1985"]~dados_df[,"nf_veg_1985"])
m2 <- glm(dados_df[,"invicta_1985"]~dados_df[,"nf_veg_1985"] + I(dados_df[,"nf_veg_1985"]^2))

call <- list(); call[[1]]=m0; call[[2]]=m1; call[[3]]=m2

(aic <- aictab(cand.set=call, modnames=c("null", "linear", "quadratic"), sort=F))


ggplot(dados_df, aes(x = nf_veg_1985, y = invicta_1985)) +
  geom_point(color = "#FFD700", size = 3, alpha = 0.5) +  
  scale_y_continuous(breaks = seq(0, 6, by = 2),  
                     limits = c(0, 6)) +   
  theme_classic() +
  theme(axis.text.x = element_text(size = 12, color = "black"),
        axis.text.y = element_text(size = 12, color = "black"),
        axis.title.x = element_text(color = "transparent"),
        axis.title.y = element_text(color = "transparent"))

ggsave("Fig_5A1985.pdf", width = 8, height = 9)

## S. invicta x Area of fragments of native vegetation 1985


m0 <- glm(dados_df[,"invicta_1985"]~1)
m1 <- glm(dados_df[,"invicta_1985"]~dados_df[,"area_veg_1985"])
m2 <- glm(dados_df[,"invicta_1985"]~dados_df[,"area_veg_1985"] + I(dados_df[,"area_veg_1985"]^2))

call <- list(); call[[1]]=m0; call[[2]]=m1; call[[3]]=m2

(aic <- aictab(cand.set=call, modnames=c("null", "linear", "quadratic"), sort=F))


ggplot(dados_df, aes(x = area_veg_1985, y = invicta_1985)) +
  geom_point(color = "#FFD700", size = 3, alpha= 0.5) + 
  geom_smooth(method = "lm", se = TRUE, fill = "#FFD700", color = "#FFD700",alpha = 0.2, size = 1.5, 
              aes(ymin = pmax(..ymin.., 0),  
                  ymax = pmin(..ymax.., 5))) +  
  scale_y_continuous(breaks = seq(0, 4, by = 1),  
                     limits = c(0, 4))   + 
  scale_x_continuous(labels = label_number(accuracy = 1), breaks = pretty(dados_df$area_urb_2015, n = 2)) + 
  theme_classic() +
  theme(axis.text.x = element_text(size = 12, color = "black"),
        axis.text.y = element_text(size = 12, color = "black"),
        axis.title.x = element_text(color = "transparent"),
        axis.title.y = element_text(color = "transparent"))

ggsave("Fig_5B1985.pdf", width = 8, height = 9)

##### S. invicta x Number of anthropic fragments 1985

m0 <- glm(dados_df[,"invicta_1985"]~1)
m1 <- glm(dados_df[,"invicta_1985"]~dados_df[,"nf_urb_1985"])
m2 <- glm(dados_df[,"invicta_1985"]~dados_df[,"nf_urb_1985"] + I(dados_df[,"nf_urb_1985"]^2))

call <- list(); call[[1]]=m0; call[[2]]=m1; call[[3]]=m2

(aic <- aictab(cand.set=call, modnames=c("null", "linear", "quadratic"), sort=F))


ggplot(dados_df, aes(x = nf_urb_1985, y = invicta_1985)) +
  geom_point(color = "#FFD700", size = 3, alpha= 0.5) +  
  geom_smooth(aes(ymin = ifelse(..ymin.. < 0, 0, ..ymin..)), method = "lm", formula = y ~ poly(x, 2), se = TRUE, 
              fill = "#FFD700", color = "#FFD700",alpha = 0.2, size = 1.5) +
  scale_y_continuous(breaks = seq(0, 5, by = 1),  
                     limits = c(0, 5)) +
  scale_x_continuous(limits = c(0, 200)) +
  theme_classic() +
  theme(axis.text.x = element_text(size = 12, color = "black"),
        axis.text.y = element_text(size = 12, color = "black"),
        axis.title.x = element_text(color = "transparent"),
        axis.title.y = element_text(color = "transparent"))

ggsave("Fig_5C1985.pdf", width = 8, height = 9)

### S. invicta x Area of anthropic fragments 1985


m0 <- glm(dados_df[,"invicta_1985"]~1)
m1 <- glm(dados_df[,"invicta_1985"]~dados_df[,"area_urb_1985"])
m2 <- glm(dados_df[,"invicta_1985"]~dados_df[,"area_urb_1985"] + I(dados_df[,"area_urb_1985"]^2))

call <- list(); call[[1]]=m0; call[[2]]=m1; call[[3]]=m2

(aic <- aictab(cand.set=call, modnames=c("null", "linear", "quadratic"), sort=F))


ggplot(dados_df, aes(x = area_urb_1985, y = invicta_1985)) +
  geom_point(color = "#FFD700", size =3 ,alpha = 0.8) +  
  scale_y_continuous(breaks = seq(0, 6, by = 1),  
                     limits = c(0, 6)) + 
  scale_x_continuous(labels = label_number(accuracy = 1), breaks = pretty(dados_df$area_urb_2015, n = 2)) +
  theme_classic() +
  theme(axis.text.x = element_text(size = 12, color = "black"),
        axis.text.y = element_text(size = 12, color = "black"),
        axis.title.x = element_text(color = "transparent"),
        axis.title.y = element_text(color = "transparent"))

ggsave("Fig_5D1985.pdf", width = 8, height = 9)

#### S. invicta x Number of fragments of native vegetation 1995

m0 <- glm(dados_df[,"invicta_1995"]~1)
m1 <- glm(dados_df[,"invicta_1995"]~dados_df[,"nf_veg_1995"])
m2 <- glm(dados_df[,"invicta_1995"]~dados_df[,"nf_veg_1995"] + I(dados_df[,"nf_veg_1985"]^2))

call <- list(); call[[1]]=m0; call[[2]]=m1; call[[3]]=m2

(aic <- aictab(cand.set=call, modnames=c("null", "linear", "quadratic"), sort=F))


ggplot(dados_df, aes(x = nf_veg_1995, y = invicta_1995)) +
  geom_point(color = "#FFD700", size = 3, alpha = 0.5) +  
  scale_y_continuous(breaks = seq(0, 6, by = 2),  
                     limits = c(0, 6)) +theme_classic() +
  theme(axis.text.x = element_text(size = 12, color = "black"),
        axis.text.y = element_text(size = 12, color = "black"),
        axis.title.x = element_text(color = "transparent"),
        axis.title.y = element_text(color = "transparent"))

ggsave("Fig_5A1995.pdf", width = 8, height = 9)


#### S. invicta x Area of fragments of native vegetation 1995

m0 <- glm(dados_df[,"invicta_1995"]~1)
m1 <- glm(dados_df[,"invicta_1995"]~dados_df[,"area_veg_1995"])
m2 <- glm(dados_df[,"invicta_1995"]~dados_df[,"area_veg_1995"] + I(dados_df[,"area_veg_1985"]^2))

call <- list(); call[[1]]=m0; call[[2]]=m1; call[[3]]=m2

(aic <- aictab(cand.set=call, modnames=c("null", "linear", "quadratic"), sort=F))


ggplot(dados_df, aes(x = area_veg_1995, y = invicta_1995)) +
  geom_point(color = "#FFD700", size = 3, alpha = 0.5) +  
  scale_y_continuous(breaks = seq(0, 6, by = 2),  
                     limits = c(0, 6)) + 
  scale_x_continuous(labels = label_number(accuracy = 1), breaks = pretty(dados_df$area_urb_2015, n = 2)) +
  theme_classic() +
  theme(axis.text.x = element_text(size = 12, color = "black"),
        axis.text.y = element_text(size = 12, color = "black"),
        axis.title.x = element_text(color = "transparent"),
        axis.title.y = element_text(color = "transparent"))

ggsave("Fig_5B1995.pdf", width = 8, height = 9)


#### S. invicta x Number of anthropic fragments 1995

m0 <- glm(dados_df[,"invicta_1995"]~1)
m1 <- glm(dados_df[,"invicta_1995"]~dados_df[,"nf_urb_1995"])
m2 <- glm(dados_df[,"invicta_1995"]~dados_df[,"nf_urb_1995"] + I(dados_df[,"nf_urb_1985"]^2))

call <- list(); call[[1]]=m0; call[[2]]=m1; call[[3]]=m2

(aic <- aictab(cand.set=call, modnames=c("null", "linear", "quadratic"), sort=F))


ggplot(dados_df, aes(x = nf_urb_1995, y = invicta_1995)) +
  geom_point(color = "#FFD700", size = 3, alpha = 0.5) +
  scale_y_continuous(breaks = seq(0, 6, by = 2),  
                     limits = c(0, 6)) +   
  theme_classic() +
  theme(axis.text.x = element_text(size = 12, color = "black"),
        axis.text.y = element_text(size = 12, color = "black"),
        axis.title.x = element_text(color = "transparent"),
        axis.title.y = element_text(color = "transparent"))

ggsave("Fig_5C1995.pdf", width = 8, height = 9)

### S. invicta x Area of anthropic fragments 1995

m0 <- glm(dados_df[,"invicta_1995"]~1)
m1 <- glm(dados_df[,"invicta_1995"]~dados_df[,"area_urb_1995"])
m2 <- glm(dados_df[,"invicta_1995"]~dados_df[,"area_urb_1995"] + I(dados_df[,"area_urb_1985"]^2))

call <- list(); call[[1]]=m0; call[[2]]=m1; call[[3]]=m2

(aic <- aictab(cand.set=call, modnames=c("null", "linear", "quadratic"), sort=F))


ggplot(dados_df, aes(x = area_urb_1995, y = invicta_1995)) +
  geom_point(color = "#FFD700", size = 3, alpha = 0.5) +  
  scale_y_continuous(breaks = seq(0, 6, by = 2),  
                     limits = c(0, 6))  + 
  scale_x_continuous(labels = label_number(accuracy = 2), breaks = pretty(dados_df$area_urb_2015, n = 4)) + 
  theme_classic() +
  theme(axis.text.x = element_text(size = 12, color = "black"),
        axis.text.y = element_text(size = 12, color = "black"),
        axis.title.x = element_text(color = "transparent"),
        axis.title.y = element_text(color = "transparent"))

ggsave("Fig_5D1995.pdf", width = 8, height = 9)

#### S. invicta x Number of fragments of native vegetation 2005

m0 <- glm(dados_df[,"invicta_2005"]~1)
m1 <- glm(dados_df[,"invicta_2005"]~dados_df[,"nf_veg_2005"])
m2 <- glm(dados_df[,"invicta_2005"]~dados_df[,"nf_veg_2005"] + I(dados_df[,"nf_veg_2005"]^2))

call <- list(); call[[1]]=m0; call[[2]]=m1; call[[3]]=m2

(aic <- aictab(cand.set=call, modnames=c("null", "linear", "quadratic"), sort=F))


ggplot(dados_df, aes(x = nf_veg_2005, y = invicta_2005)) +
  geom_point(color = "#FFD700", size = 3, alpha= 0.5) +  
  scale_y_continuous(breaks = seq(0, 18, by = 3),  
                     limits = c(0, 18)) +   
  theme_classic() +
  theme(axis.text.x = element_text(size = 12, color = "black"),
        axis.text.y = element_text(size = 12, color = "black"),
        axis.title.x = element_text(color = "transparent"),
        axis.title.y = element_text(color = "transparent"))

ggsave("Fig_5A2005.pdf", width = 8, height = 9)

#### S. invicta x Area of fragments of native vegetation 2005

m0 <- glm(dados_df[,"invicta_2005"]~1)
m1 <- glm(dados_df[,"invicta_2005"]~dados_df[,"area_veg_2005"])
m2 <- glm(dados_df[,"invicta_2005"]~dados_df[,"area_veg_2005"] + I(dados_df[,"area_veg_2005"]^2))

call <- list(); call[[1]]=m0; call[[2]]=m1; call[[3]]=m2

(aic <- aictab(cand.set=call, modnames=c("null", "linear", "quadratic"), sort=F))


ggplot(dados_df, aes(x = area_veg_1995, y = invicta_1995)) +
  geom_point(color = "#FFD700", size = 3) +  
  scale_y_continuous(breaks = seq(0, 18, by = 3),  
                     limits = c(0, 18)) +   
scale_x_continuous(labels = label_number(accuracy = 1), breaks = pretty(dados_df$area_urb_2015, n = 2)) + 
  theme_classic() +
  theme(axis.text.x = element_text(size = 12, color = "black"),
        axis.text.y = element_text(size = 12, color = "black"),
        axis.title.x = element_text(color = "transparent"),
        axis.title.y = element_text(color = "transparent"))

ggsave("Fig_5B2005.pdf", width = 8, height = 9)


### S. invicta x Number of anthropic fragments 2005

###TESTE DE MODELOS###
m0 <- glm(dados_df[,"invicta_2005"]~1)
m1 <- glm(dados_df[,"invicta_2005"]~dados_df[,"nf_urb_2005"])
m2 <- glm(dados_df[,"invicta_2005"]~dados_df[,"nf_urb_2005"] + I(dados_df[,"nf_urb_2005"]^2))

call <- list(); call[[1]]=m0; call[[2]]=m1; call[[3]]=m2

(aic <- aictab(cand.set=call, modnames=c("null", "linear", "quadratic"), sort=F))


ggplot(dados_df, aes(x = nf_urb_2005, y = invicta_2005)) +
  geom_point(color = "#FFD700", size = 3) +  
  geom_smooth(aes(ymin = ifelse(..ymin.. < 0, 0, ..ymin..)), 
              method = "lm", formula = y ~ poly(x, 2), se = TRUE, fill = "#FFD700",  color = "#FFD700", alpha = 0.2,  size = 1.5,fullrange = TRUE) +  
  scale_y_continuous(breaks = seq(0, 20, by = 5),  
                     limits = c(0, 20)) + 
  scale_x_continuous(limits = range(dados_df$nf_urb_2005)) +  
  theme_classic() +
  theme(axis.text.x = element_text(size = 12, color = "black"),
        axis.text.y = element_text(size = 12, color = "black"),
        axis.title.x = element_text(color = "transparent"),
        axis.title.y = element_text(color = "transparent"))

ggsave("Fig_5C2005.pdf", width = 8, height = 9)


####  S. invicta x Area of anthropic fragments 2005

m0 <- glm(dados_df[,"invicta_2005"]~1)
m1 <- glm(dados_df[,"invicta_2005"]~dados_df[,"area_urb_2005"])
m2 <- glm(dados_df[,"invicta_2005"]~dados_df[,"area_urb_2005"] + I(dados_df[,"area_urb_2005"]^2))

call <- list(); call[[1]]=m0; call[[2]]=m1; call[[3]]=m2

(aic <- aictab(cand.set=call, modnames=c("null", "linear", "quadratic"), sort=F))



ggplot(dados_df, aes(x = area_urb_2005, y = invicta_2005)) +
  geom_point(color = "#FFD700", size = 3, alpha= 0.5) +  
  geom_smooth(method = "lm", se = TRUE, fill = "#FFD700", color = "#FFD700",alpha = 0.2, size = 1.5, 
              aes(ymin = pmax(..ymin.., 0),  
                  ymax = pmin(..ymax.., 5))) +  
  scale_y_continuous(breaks = seq(0, 15, by = 5),  
                     limits = c(0, 15)) +
  scale_x_continuous(labels = label_number(accuracy = 1), breaks = pretty(dados_df$area_urb_2015, n = 2)) + 
  theme_classic() +
  theme(axis.text.x = element_text(size = 12, color = "black"),
        axis.text.y = element_text(size = 12, color = "black"),
        axis.title.x = element_text(color = "transparent"),
        axis.title.y = element_text(color = "transparent"))

ggsave("Fig_5D2005.pdf", width = 8, height = 9)

#### S. invicta x Number of fragments of native vegetation 2015

m0 <- glm(dados_df[,"invicta_2015"]~1)
m1 <- glm(dados_df[,"invicta_2015"]~dados_df[,"nf_veg_2015"])
m2 <- glm(dados_df[,"invicta_2015"]~dados_df[,"nf_veg_2015"] + I(dados_df[,"nf_veg_2015"]^2))

call <- list(); call[[1]]=m0; call[[2]]=m1; call[[3]]=m2

(aic <- aictab(cand.set=call, modnames=c("null", "linear", "quadratic"), sort=F))



ggplot(dados_df, aes(x = nf_veg_2015, y = invicta_2015)) +
  geom_point(color = "#FFD700", size = 3, alpha = 0.5) + 
  scale_y_continuous(breaks = seq(0, 20, by = 5),  
                     limits = c(0, 20)) +   
  theme_classic() +
  theme(axis.text.x = element_text(size = 12, color = "black"),
        axis.text.y = element_text(size = 12, color = "black"),
        axis.title.x = element_text(color = "transparent"),
        axis.title.y = element_text(color = "transparent"))

ggsave("Fig_5A2015.pdf", width = 8, height = 9)

#### S. invicta x Area of fragments of native vegetation 2015

m0 <- glm(dados_df[,"invicta_2015"]~1)
m1 <- glm(dados_df[,"invicta_2015"]~dados_df[,"area_veg_2015"])
m2 <- glm(dados_df[,"invicta_2015"]~dados_df[,"area_veg_2015"] + I(dados_df[,"area_veg_2015"]^2))

call <- list(); call[[1]]=m0; call[[2]]=m1; call[[3]]=m2

(aic <- aictab(cand.set=call, modnames=c("null", "linear", "quadratic"), sort=F))


ggplot(dados_df, aes(x = area_veg_2015, y = invicta_2015)) +
  geom_point(color = "#FFD700", size = 3, alpha = 0.5) +  
  scale_y_continuous(breaks = seq(0, 30, by = 5),  
                     limits = c(0, 30)) +   
  scale_x_continuous(labels = label_number(accuracy = 1), breaks = pretty(dados_df$area_urb_2015, n = 2)) + 
  theme_classic() +
  theme(axis.text.x = element_text(size = 12, color = "black"),
        axis.text.y = element_text(size = 12, color = "black"),
        axis.title.x = element_text(color = "transparent"),
        axis.title.y = element_text(color = "transparent"))

ggsave("Fig_5B2015.pdf", width = 8, height = 9)


#####  S. invicta x Number of anthropic fragments 2015

m0 <- glm(dados_df[,"invicta_2015"]~1)
m1 <- glm(dados_df[,"invicta_2015"]~dados_df[,"nf_urb_2015"])
m2 <- glm(dados_df[,"invicta_2015"]~dados_df[,"nf_urb_2015"] + I(dados_df[,"nf_urb_2015"]^2))

call <- list(); call[[1]]=m0; call[[2]]=m1; call[[3]]=m2

(aic <- aictab(cand.set=call, modnames=c("null", "linear", "quadratic"), sort=F))


ggplot(dados_df, aes(x = nf_urb_2015, y = invicta_2015)) +
  geom_point(color = "#FFD700", size = 3) +  # Pontos em dourado
  geom_smooth(aes(ymin = ifelse(..ymin.. < 0, 0, ..ymin..)),method = "lm", formula = y ~ poly(x, 2), se = TRUE, fill = "#FFD700", color = "#FFD700", alpha = 0.2,  size = 1.5) +
  scale_y_continuous(breaks = seq(0, 30, by = 5),  
                     limits = c(0, 30)) +  
  theme_classic() +
  theme(axis.text.x = element_text(size = 12, color = "black"),
        axis.text.y = element_text(size = 12, color = "black"),
        axis.title.x = element_text(color = "transparent"),
        axis.title.y = element_text(color = "transparent"))

ggsave("Fig_5C2015.pdf", width = 8, height = 9)


#####  S. invicta x Area of anthropic fragments 2015

m0 <- glm(dados_df[,"invicta_2015"]~1)
m1 <- glm(dados_df[,"invicta_2015"]~dados_df[,"area_urb_2015"])
m2 <- glm(dados_df[,"invicta_2015"]~dados_df[,"area_urb_2015"] + I(dados_df[,"area_urb_2015"]^2))

call <- list(); call[[1]]=m0; call[[2]]=m1; call[[3]]=m2

(aic <- aictab(cand.set=call, modnames=c("null", "linear", "quadratic"), sort=F))


ggplot(dados_df, aes(x = area_urb_2015, y = invicta_2015)) +
  geom_point(color = "#FFD700", size = 3) +  
  geom_smooth(aes(ymin = ifelse(..ymin.. < 0, 0, ..ymin..)), method = "lm", formula = y ~ poly(x, 2), se = TRUE, fill = "#FFD700",  color = "#FFD700",  alpha = 0.2, size = 1.5) +
  scale_y_continuous(breaks = seq(0, 30, by = 5),  limits = c(0, 30)) + 
  scale_x_continuous(labels = label_number(accuracy = 1), breaks = pretty(dados_df$area_urb_2015, n = 3)) + 
  theme_classic() +
  theme(axis.text.x = element_text(size = 12, color = "black"),
        axis.text.y = element_text(size = 12, color = "black"),
        axis.title.x = element_text(color = "transparent"),
        axis.title.y = element_text(color = "transparent"))


ggsave("Fig_5D2015.pdf", width = 8, height = 9)

### S. invicta x Number of fragments of native vegetation 2024

m0 <- glm(dados_df[,"invicta_2023"]~1)
m1 <- glm(dados_df[,"invicta_2023"]~dados_df[,"nf_veg_2022"])
m2 <- glm(dados_df[,"invicta_2023"]~dados_df[,"nf_veg_2022"] + I(dados_df[,"nf_veg_2022"]^2))

call <- list(); call[[1]]=m0; call[[2]]=m1; call[[3]]=m2

(aic <- aictab(cand.set=call, modnames=c("null", "linear", "quadratic"), sort=F))


ggplot(dados_df, aes(x = nf_veg_2022, y = invicta_2023)) +
  geom_point(color = "#FFD700", size = 3, alpha = 0.5) +  
  scale_y_continuous(breaks = seq(0, 30, by = 5),  
                     limits = c(0, 30)) +   
  theme_classic() +
  theme(axis.text.x = element_text(size = 12, color = "black"),
        axis.text.y = element_text(size = 12, color = "black"),
        axis.title.x = element_text(color = "transparent"),
        axis.title.y = element_text(color = "transparent"))

ggsave("Fig_5A2023.pdf", width = 8, height = 9)

#### S. invicta x Area of fragments of native vegetation 2024

m0 <- glm(dados_df[,"invicta_2023"]~1)
m1 <- glm(dados_df[,"invicta_2023"]~dados_df[,"area_veg_2022"])
m2 <- glm(dados_df[,"invicta_2023"]~dados_df[,"area_veg_2022"] + I(dados_df[,"area_veg_2022"]^2))


call <- list(); call[[1]]=m0; call[[2]]=m1; call[[3]]=m2

(aic <- aictab(cand.set=call, modnames=c("null", "linear", "quadratic"), sort=F))


ggplot(dados_df, aes(x = area_veg_2022, y = invicta_2023)) +
  geom_point(color = "#FFD700", size = 3, alpha = 0.5) +  
  scale_y_continuous(breaks = seq(0, 30, by = 5),  
                     limits = c(0, 30)) +  
  scale_x_continuous(labels = label_number(accuracy = 1), breaks = pretty(dados_df$area_urb_2015, n = 2)) + 
  theme_classic() +
  theme(axis.text.x = element_text(size = 12, color = "black"),
        axis.text.y = element_text(size = 12, color = "black"),
        axis.title.x = element_text(color = "transparent"),
        axis.title.y = element_text(color = "transparent"))

ggsave("Fig_5B2023.pdf", width = 8, height = 9)

#####  S. invicta x Number of anthropic fragments 2024

m0 <- glm(dados_df[,"invicta_2023"]~1)
m1 <- glm(dados_df[,"invicta_2023"]~dados_df[,"nf_urb_2022"])
m2 <- glm(dados_df[,"invicta_2023"]~dados_df[,"nf_urb_2022"] + I(dados_df[,"nf_urb_2022"]^2))

call <- list(); call[[1]]=m0; call[[2]]=m1; call[[3]]=m2

(aic <- aictab(cand.set=call, modnames=c("null", "linear", "quadratic"), sort=F))


ggplot(dados_df, aes(x = nf_urb_2022, y = invicta_2023)) +
  geom_point(color = "#FFD700", size = 3) +  
  geom_smooth(aes(ymin = ifelse(..ymin.. < 0, 0, ..ymin..)),method = "lm", formula = y ~ poly(x, 2), se = TRUE, fill = "#FFD700", color = "#FFD700", alpha = 0.2, size = 1.5) +
  scale_y_continuous(breaks = seq(0, 30, by = 5),  
                     limits = c(0, 30)) + theme_classic() +
  theme(axis.text.x = element_text(size = 12, color = "black"),
        axis.text.y = element_text(size = 12, color = "black"),
        axis.title.x = element_text(color = "transparent"),
        axis.title.y = element_text(color = "transparent"))

ggsave("Fig_5C2023.pdf", width = 8, height = 9)

### ###  S. invicta x Area of anthropic fragments 2024

m0 <- glm(dados_df[,"invicta_2023"]~1)
m1 <- glm(dados_df[,"invicta_2023"]~dados_df[,"area_urb_2022"])
m2 <- glm(dados_df[,"invicta_2023"]~dados_df[,"area_urb_2022"] + I(dados_df[,"area_urb_2022"]^2))

call <- list(); call[[1]]=m0; call[[2]]=m1; call[[3]]=m2

(aic <- aictab(cand.set=call, modnames=c("null", "linear", "quadratic"), sort=F))

ggplot(dados_df, aes(x = area_urb_2022, y = invicta_2023)) +
  geom_point(color = "#FFD700", size = 3) + 
  geom_smooth(aes(ymin = ifelse(..ymin.. < 0, 0, ..ymin..)),method = "lm", formula = y ~ poly(x, 2), se = TRUE, fill = "#FFD700",  color = "#FFD700",  alpha = 0.2,  size = 1.5) +
  scale_y_continuous(breaks = seq(0, 30, by = 5),  limits = c(0, 30)) +
  scale_x_continuous(labels = label_number(accuracy = 1), breaks = pretty(dados_df$area_urb_2015, n = 2)) + # Definir limites do eixo y
  theme_classic() +
  theme(axis.text.x = element_text(size = 12, color = "black"),
        axis.text.y = element_text(size = 12, color = "black"),
        axis.title.x = element_text(color = "transparent"),
        axis.title.y = element_text(color = "transparent"))

ggsave("Fig_5D2023.pdf", width = 8, height = 9)


###########  SAEVISSIMA 




####S. saevissima x Number of fragments of native vegetation 1985


m0 <- glm(dados_df1[,"saevissima_1985"]~1)
m1 <- glm(dados_df1[,"saevissima_1985"]~dados_df1[,"nf_veg_1985"])
m2 <- glm(dados_df1[,"saevissima_1985"]~dados_df1[,"nf_veg_1985"] + I(dados_df1[,"nf_veg_1985"]^2))

call <- list(); call[[1]]=m0; call[[2]]=m1; call[[3]]=m2
(aic <- aictab(cand.set=call, modnames=c("null", "linear", "quadratic"), sort=F))


ggplot(dados_df1, aes(x = nf_veg_1985, y = saevissima_1985)) +
  geom_point(color = "#483d8d", size = 3, alpha= 0.5) + 
  scale_y_continuous(breaks = seq(0, 90, by = 20),  
                     limits = c(0, 90)) +   
  theme_classic() +
  theme(axis.text.x = element_text(size = 12, color = "black"),
        axis.text.y = element_text(size = 12, color = "black"),
        axis.title.x = element_text(color = "transparent"),
        axis.title.y = element_text(color = "transparent"))

ggsave("Fig_6A1985.pdf", width = 8, height = 9)

### S. saevissima x Area of fragments of native vegetation 1985 

m0 <- glm(dados_df1[,"saevissima_1985"]~1)
m1 <- glm(dados_df1[,"saevissima_1985"]~dados_df1[,"area_veg_1985"])
m2 <- glm(dados_df1[,"saevissima_1985"]~dados_df1[,"area_veg_1985"] + I(dados_df1[,"area_veg_1985"]^2))

call <- list(); call[[1]]=m0; call[[2]]=m1; call[[3]]=m2
(aic <- aictab(cand.set=call, modnames=c("null", "linear", "quadratic"), sort=F))


ggplot(dados_df1, aes(x = area_veg_1985, y = saevissima_1985)) +
  geom_point(color = "#483d8d", size = 3, alpha= 0.5) + 
  scale_y_continuous(breaks = seq(0, 80, by = 20),  
                     limits = c(0, 80)) +    scale_x_continuous(labels = label_number(accuracy = 1), breaks = pretty(dados_df$area_urb_2015, n = 2)) + 
  theme_classic() +
  theme(axis.text.x = element_text(size = 12, color = "black"),
        axis.text.y = element_text(size = 12, color = "black"),
        axis.title.x = element_text(color = "transparent"),
        axis.title.y = element_text(color = "transparent"))

ggsave("Fig_6B1985.pdf", width = 8, height = 9)

### S. saevissima x Number of anthropic fragments 1985

m0 <- glm(dados_df1[,"saevissima_1985"]~1)
m1 <- glm(dados_df1[,"saevissima_1985"]~dados_df1[,"nf_urb_1985"])
m2 <- glm(dados_df1[,"saevissima_1985"]~dados_df1[,"nf_urb_1985"] + I(dados_df1[,"nf_urb_1985"]^2))

call <- list(); call[[1]]=m0; call[[2]]=m1; call[[3]]=m2
(aic <- aictab(cand.set=call, modnames=c("null", "linear", "quadratic"), sort=F))


ggplot(dados_df1, aes(x = nf_urb_1985, y = saevissima_1985)) +
  geom_point(color = "#483d8d", size = 3, alpha= 0.5) +  
  geom_smooth(aes(ymin = ifelse(..ymin.. < 0, 0, ..ymin..)),method = "lm", formula = y ~ poly(x, 2), se = TRUE, fill = "#483d8d", color = "#483d8d",alpha = 0.2, size = 1.5) +
  scale_y_continuous(breaks = seq(0, 80, by = 20),  
                     limits = c(0, 80)) + 
  theme_classic() +
  theme(axis.text.x = element_text(size = 12, color = "black"),
        axis.text.y = element_text(size = 12, color = "black"),
        axis.title.x = element_text(color = "transparent"),
        axis.title.y = element_text(color = "transparent"))

ggsave("Fig_6C1985.pdf", width = 8, height = 9)

### S. saevissima x Area of anthropic fragments 1985

m0 <- glm(dados_df1[,"saevissima_1985"]~1)
m1 <- glm(dados_df1[,"saevissima_1985"]~dados_df1[,"area_urb_1985"])
m2 <- glm(dados_df1[,"saevissima_1985"]~dados_df1[,"area_urb_1985"] + I(dados_df1[,"area_urb_1985"]^2))

call <- list(); call[[1]]=m0; call[[2]]=m1; call[[3]]=m2
(aic <- aictab(cand.set=call, modnames=c("null", "linear", "quadratic"), sort=F))


ggplot(dados_df1, aes(x = area_urb_1985, y = saevissima_1985)) +
  geom_point(color = "#483d8d", size = 3,alpha= 0.5) +  
  geom_smooth(aes(ymin = ifelse(..ymin.. < 0, 0, ..ymin..)),method = "lm", formula = y ~ poly(x, 2), se = TRUE, fill = "#483d8d",  color = "#483d8d", alpha = 0.2,  size = 1.5) +
  scale_y_continuous(breaks = seq(0, 90, by = 20),  
                     limits = c(0, 90)) +  scale_x_continuous(labels = label_number(accuracy = 1), breaks = pretty(dados_df$area_urb_2015, n = 2)) +  # Definir limites do eixo y
  theme_classic() +
  theme(axis.text.x = element_text(size = 12, color = "black"),
        axis.text.y = element_text(size = 12, color = "black"),
        axis.title.x = element_text(color = "transparent"),
        axis.title.y = element_text(color = "transparent"))

ggsave("Fig_6D1985.pdf", width = 8, height = 9)

#### S. saevissima x Number of fragments of native vegetation 1995

m0 <- glm(dados_df1[,"saevissima_1995"]~1)
m1 <- glm(dados_df1[,"saevissima_1995"]~dados_df1[,"nf_veg_1995"])
m2 <- glm(dados_df1[,"saevissima_1995"]~dados_df1[,"nf_veg_1995"] + I(dados_df1[,"nf_veg_1995"]^2))

call <- list(); call[[1]]=m0; call[[2]]=m1; call[[3]]=m2
(aic <- aictab(cand.set=call, modnames=c("null", "linear", "quadratic"), sort=F))

ggplot(dados_df1, aes(x = nf_veg_1995, y = saevissima_1995)) +
  geom_point(color = "#483d8d", size = 3, alpha = 0.5) +  
  scale_y_continuous(breaks = seq(0, 80, by = 20),  
                     limits = c(0, 80)) +   
  theme_classic() +
  theme(axis.text.x = element_text(size = 12, color = "black"),
        axis.text.y = element_text(size = 12, color = "black"),
        axis.title.x = element_text(color = "transparent"),
        axis.title.y = element_text(color = "transparent"))

ggsave("Fig_6A1995.pdf", width = 8, height = 9)


#### S. saevissima x Area of fragments of native vegetation 1995

m0 <- glm(dados_df1[,"saevissima_1995"]~1)
m1 <- glm(dados_df1[,"saevissima_1995"]~dados_df1[,"area_veg_1995"])
m2 <- glm(dados_df1[,"saevissima_1995"]~dados_df1[,"area_veg_1995"] + I(dados_df1[,"area_veg_1995"]^2))

call <- list(); call[[1]]=m0; call[[2]]=m1; call[[3]]=m2
(aic <- aictab(cand.set=call, modnames=c("null", "linear", "quadratic"), sort=F))


ggplot(dados_df1, aes(x = area_veg_1995, y = saevissima_1995)) +
  geom_point(color = "#483d8d", size = 3, alpha = 0.5) + 
  scale_y_continuous(breaks = seq(0, 80, by = 20),  
                     limits = c(0, 80)) +    scale_x_continuous(labels = label_number(accuracy = 1), breaks = pretty(dados_df$area_veg_1995, n = 2)) +  # Definir limites do eixo y
  theme_classic() +
  theme(axis.text.x = element_text(size = 12, color = "black"),
        axis.text.y = element_text(size = 12, color = "black"),
        axis.title.x = element_text(color = "transparent"),
        axis.title.y = element_text(color = "transparent"))

ggsave("Fig_6B1995.pdf", width = 8, height = 9)


#### S. saevissima x Number of anthropic fragments 1995

m0 <- glm(dados_df1[,"saevissima_1995"]~1)
m1 <- glm(dados_df1[,"saevissima_1995"]~dados_df1[,"nf_urb_1995"])
m2 <- glm(dados_df1[,"saevissima_1995"]~dados_df1[,"nf_urb_1995"] + I(dados_df1[,"nf_urb_1995"]^2))

call <- list(); call[[1]]=m0; call[[2]]=m1; call[[3]]=m2
(aic <- aictab(cand.set=call, modnames=c("null", "linear", "quadratic"), sort=F))


ggplot(dados_df1, aes(x = nf_urb_1995, y = saevissima_1995)) +
  geom_point(color = "#483d8d", size = 3,alpha= 0.5) +  
  geom_smooth(aes(ymin = ifelse(..ymin.. < 0, 0, ..ymin..)),method = "lm", formula = y ~ poly(x, 2), se = TRUE, fill = "#483d8d",  color = "#483d8d", alpha = 0.2, size = 1.5) +
  scale_y_continuous(breaks = seq(0, 90, by = 20),  
                     limits = c(0, 90)) +  
  theme_classic() +
  theme(axis.text.x = element_text(size = 12, color = "black"),
        axis.text.y = element_text(size = 12, color = "black"),
        axis.title.x = element_text(color = "transparent"),
        axis.title.y = element_text(color = "transparent"))

ggsave("Fig_6C1995.pdf", width = 8, height = 9)

### S. saevissima x Area of anthropic fragments 1995


m0 <- glm(dados_df1[,"saevissima_1995"]~1)
m1 <- glm(dados_df1[,"saevissima_1995"]~dados_df1[,"area_urb_1995"])
m2 <- glm(dados_df1[,"saevissima_1995"]~dados_df1[,"area_urb_1995"] + I(dados_df1[,"area_urb_1995"]^2))

call <- list(); call[[1]]=m0; call[[2]]=m1; call[[3]]=m2
(aic <- aictab(cand.set=call, modnames=c("null", "linear", "quadratic"), sort=F))


ggplot(dados_df1, aes(x = area_urb_1995, y = saevissima_1995)) +
  geom_point(color = "#483d8d", size = 3,alpha= 0.5) + 
  geom_smooth(aes(ymin = ifelse(..ymin.. < 0, 0, ..ymin..)),method = "lm", formula = y ~ poly(x, 2), se = TRUE, fill = "#483d8d",  color = "#483d8d", alpha = 0.2, size = 1.5) +
  scale_y_continuous(breaks = seq(0, 80, by = 20),  
                     limits = c(0, 80)) +    scale_x_continuous(labels = label_number(accuracy = 1), breaks = pretty(dados_df$area_urb_1995, n = 2)) +  
  theme_classic() +
  theme(axis.text.x = element_text(size = 12, color = "black"),
        axis.text.y = element_text(size = 12, color = "black"),
        axis.title.x = element_text(color = "transparent"),
        axis.title.y = element_text(color = "transparent"))

ggsave("Fig_6D1995.pdf", width = 8, height = 9)

### S. saevissima x Number of fragments of native vegetation 2005

m0 <- glm(dados_df1[,"saevissima_2005"]~1)
m1 <- glm(dados_df1[,"saevissima_2005"]~dados_df1[,"nf_veg_2005"])
m2 <- glm(dados_df1[,"saevissima_2005"]~dados_df1[,"nf_veg_2005"] + I(dados_df1[,"nf_veg_2005"]^2))

call <- list(); call[[1]]=m0; call[[2]]=m1; call[[3]]=m2
(aic <- aictab(cand.set=call, modnames=c("null", "linear", "quadratic"), sort=F))


ggplot(dados_df1, aes(x = nf_veg_2005, y = saevissima_2005)) +
  geom_point(color = "#483d8d", size = 3, alpha= 0.5) + 
  scale_y_continuous(breaks = seq(0, 80, by = 20),  
                     limits = c(0, 80)) +   
  theme_classic() +
  theme(axis.text.x = element_text(size = 12, color = "black"),
        axis.text.y = element_text(size = 12, color = "black"),
        axis.title.x = element_text(color = "transparent"),
        axis.title.y = element_text(color = "transparent"))

ggsave("Fig_6A2005.pdf", width = 8, height = 9)

#### S. saevissima x Area of fragments of native vegetation 2005

m0 <- glm(dados_df1[,"saevissima_2005"]~1)
m1 <- glm(dados_df1[,"saevissima_2005"]~dados_df1[,"area_veg_2005"])
m2 <- glm(dados_df1[,"saevissima_2005"]~dados_df1[,"area_veg_2005"] + I(dados_df1[,"area_veg_2005"]^2))

call <- list(); call[[1]]=m0; call[[2]]=m1; call[[3]]=m2
(aic <- aictab(cand.set=call, modnames=c("null", "linear", "quadratic"), sort=F))


ggplot(dados_df1, aes(x = area_veg_2005, y = saevissima_2005)) +
  geom_point(color = "#483d8d", size = 3, alpha= 0.5) + 
  geom_smooth(method = "lm", se = TRUE, fill = "#483d8d",  color = "#483d8d", alpha = 0.2, size = 1.5) +  
  scale_y_continuous(breaks = seq(0, 90, by = 20),  
                     limits = c(0, 90)) +    scale_x_continuous(labels = label_number(accuracy = 1), breaks = pretty(dados_df$area_veg_2005, n = 3)) +   
  theme_classic() +
  theme(axis.text.x = element_text(size = 12, color = "black"),
        axis.text.y = element_text(size = 12, color = "black"),
        axis.title.x = element_text(color = "transparent"),
        axis.title.y = element_text(color = "transparent"))

ggsave("Fig_6B2005.pdf", width = 8, height = 9)

#### S. saevissima x Number of anthropic fragments 2005 


m0 <- glm(dados_df1[,"saevissima_2005"]~1)
m1 <- glm(dados_df1[,"saevissima_2005"]~dados_df1[,"nf_urb_2005"])
m2 <- glm(dados_df1[,"saevissima_2005"]~dados_df1[,"nf_urb_2005"] + I(dados_df1[,"nf_urb_2005"]^2))

call <- list(); call[[1]]=m0; call[[2]]=m1; call[[3]]=m2
(aic <- aictab(cand.set=call, modnames=c("null", "linear", "quadratic"), sort=F))


ggplot(dados_df1, aes(x = nf_urb_2005, y = saevissima_2005)) +
  geom_point(color = "#483d8d", size = 3,alpha= 0.5) + 
  geom_smooth(aes(ymin = ifelse(..ymin.. < 0, 0, ..ymin..)),method = "lm", formula = y ~ poly(x, 2), se = TRUE, fill = "#483d8d",  color = "#483d8d",  alpha = 0.2,  size = 1.5) +
  scale_y_continuous(breaks = seq(0, 80, by = 20),  
                     limits = c(00, 80)) +  # Definir limites do eixo y
  theme_classic() +
  theme(axis.text.x = element_text(size = 12, color = "black"),
        axis.text.y = element_text(size = 12, color = "black"),
        axis.title.x = element_text(color = "transparent"),
        axis.title.y = element_text(color = "transparent"))

ggsave("Fig_6C2005.pdf", width = 8, height = 9)

### S. saevissima x Area of anthropic fragments 2005  

m0 <- glm(dados_df1[,"saevissima_2005"]~1)
m1 <- glm(dados_df1[,"saevissima_2005"]~dados_df1[,"area_urb_2005"])
m2 <- glm(dados_df1[,"saevissima_2005"]~dados_df1[,"area_urb_2005"] + I(dados_df1[,"area_urb_2005"]^2))

call <- list(); call[[1]]=m0; call[[2]]=m1; call[[3]]=m2
(aic <- aictab(cand.set=call, modnames=c("null", "linear", "quadratic"), sort=F))


ggplot(dados_df1, aes(x = area_urb_2005, y = saevissima_2005)) +
  geom_point(color = "#483d8d", size = 3,alpha= 0.5) + 
  geom_smooth(aes(ymin = ifelse(..ymin.. < 0, 0, ..ymin..)),method = "lm", formula = y ~ poly(x, 2), se = TRUE, fill = "#483d8d", color = "#483d8d", alpha = 0.2, size = 1.5) +
  scale_y_continuous(breaks = seq(0, 80, by = 20),  
                     limits = c(0, 80)) +    scale_x_continuous(labels = label_number(accuracy = 1), breaks = pretty(dados_df$area_urb_2005, n = 2)) + 
  theme_classic() +
  theme(axis.text.x = element_text(size = 12, color = "black"),
        axis.text.y = element_text(size = 12, color = "black"),
        axis.title.x = element_text(color = "transparent"),
        axis.title.y = element_text(color = "transparent"))

ggsave("Fig_6D2005.pdf", width = 8, height = 9)

### S. saevissima x Number of fragments of native vegetation 2015

m0 <- glm(dados_df1[,"saevissima_2015"]~1)
m1 <- glm(dados_df1[,"saevissima_2015"]~dados_df1[,"nf_veg_2015"])
m2 <- glm(dados_df1[,"saevissima_2015"]~dados_df1[,"nf_veg_2015"] + I(dados_df1[,"nf_veg_2015"]^2))

call <- list(); call[[1]]=m0; call[[2]]=m1; call[[3]]=m2
(aic <- aictab(cand.set=call, modnames=c("null", "linear", "quadratic"), sort=F))


ggplot(dados_df1, aes(x = nf_veg_2015, y = saevissima_2015)) +
  geom_point(color = "#483d8d", size = 3, alpha= 0.5) +  
  scale_y_continuous(breaks = seq(0, 80, by = 20),  
                     limits = c(0, 80)) +   
  theme_classic() +
  theme(axis.text.x = element_text(size = 12, color = "black"),
        axis.text.y = element_text(size = 12, color = "black"),
        axis.title.x = element_text(color = "transparent"),
        axis.title.y = element_text(color = "transparent"))

ggsave("Fig_6A2015.pdf", width = 8, height = 9)

### S. saevissima x Area of fragments of native vegetation 2015

m0 <- glm(dados_df1[,"saevissima_2015"]~1)
m1 <- glm(dados_df1[,"saevissima_2015"]~dados_df1[,"area_veg_2015"])
m2 <- glm(dados_df1[,"saevissima_2015"]~dados_df1[,"area_veg_2015"] + I(dados_df1[,"area_veg_2015"]^2))

call <- list(); call[[1]]=m0; call[[2]]=m1; call[[3]]=m2
(aic <- aictab(cand.set=call, modnames=c("null", "linear", "quadratic"), sort=F))

ggplot(dados_df1, aes(x = area_veg_2015, y = saevissima_2015)) +
  geom_point(color = "#483d8d", size = 3,alpha= 0.5) +  
  geom_smooth(aes(ymin = ifelse(..ymin.. < 0, 0, ..ymin..)),method = "lm", formula = y ~ poly(x, 2), se = TRUE, fill = "#483d8d", color = "#483d8d", alpha = 0.2, size = 1.5) +
  scale_y_continuous(breaks = seq(0, 80, by = 20),  
                     limits = c(0, 80)) +  scale_x_continuous(labels = label_number(accuracy = 1), breaks = pretty(dados_df$area_veg_2015, n = 2)) + 
  theme_classic() +
  theme(axis.text.x = element_text(size = 12, color = "black"),
        axis.text.y = element_text(size = 12, color = "black"),
        axis.title.x = element_text(color = "transparent"),
        axis.title.y = element_text(color = "transparent"))

ggsave("Fig_6B2015.pdf", width = 8, height = 9)

#### S. saevissima x Number of anthropic fragments 2015  

m0 <- glm(dados_df1[,"saevissima_2015"]~1)
m1 <- glm(dados_df1[,"saevissima_2015"]~dados_df1[,"nf_urb_2015"])
m2 <- glm(dados_df1[,"saevissima_2015"]~dados_df1[,"nf_urb_2015"] + I(dados_df1[,"nf_urb_2015"]^2))

call <- list(); call[[1]]=m0; call[[2]]=m1; call[[3]]=m2
(aic <- aictab(cand.set=call, modnames=c("null", "linear", "quadratic"), sort=F))

ggplot(dados_df1, aes(x = nf_urb_2015, y = saevissima_2015)) +
  geom_point(color = "#483d8d", size = 3,alpha= 0.5) + 
  geom_smooth(aes(ymin = ifelse(..ymin.. < 0, 0, ..ymin..)),method = "lm", formula = y ~ poly(x, 2), se = TRUE, fill = "#483d8d", color = "#483d8d", alpha = 0.2,  size = 1.5) +
  scale_y_continuous(breaks = seq(0, 80, by = 20),  
                     limits = c(0, 80)) +  
  theme_classic() +
  theme(axis.text.x = element_text(size = 12, color = "black"),
        axis.text.y = element_text(size = 12, color = "black"),
        axis.title.x = element_text(color = "transparent"),
        axis.title.y = element_text(color = "transparent"))

ggsave("Fig_6C2015.pdf", width = 8, height = 9)

#### S. saevissima x Area of anthropic fragments 2015  

m0 <- glm(dados_df1[,"saevissima_2015"]~1)
m1 <- glm(dados_df1[,"saevissima_2015"]~dados_df1[,"area_urb_2015"])
m2 <- glm(dados_df1[,"saevissima_2015"]~dados_df1[,"area_urb_2015"] + I(dados_df1[,"area_urb_2015"]^2))

call <- list(); call[[1]]=m0; call[[2]]=m1; call[[3]]=m2
(aic <- aictab(cand.set=call, modnames=c("null", "linear", "quadratic"), sort=F))

ggplot(dados_df1, aes(x = area_urb_2015, y = saevissima_2015)) +
  geom_point(color = "#483d8d", size = 3,alpha= 0.5) +  
  geom_smooth(method = "lm", formula = y ~ poly(x, 2), se = TRUE,  fill = "#483d8d", color = "#483d8d",alpha = 0.2, size = 1.5) +
  scale_y_continuous(breaks = seq(0, 70, by = 10),  
                     limits = c(0, 70)) + scale_x_continuous(labels = label_number(accuracy = 1), breaks = pretty(dados_df$area_urb_2015, n = 2)) + 
  theme_classic() +
  theme(axis.text.x = element_text(size = 12, color = "black"),
        axis.text.y = element_text(size = 12, color = "black"),
        axis.title.x = element_text(color = "transparent"),
        axis.title.y = element_text(color = "transparent"))

ggsave("Fig_6D2015.pdf", width = 8, height = 9)

##### S. saevissima x Number of fragments of native vegetation 2024

m0 <- glm(dados_df1[,"saevissima_2023"]~1)
m1 <- glm(dados_df1[,"saevissima_2023"]~dados_df1[,"nf_veg_2022"])
m2 <- glm(dados_df1[,"saevissima_2023"]~dados_df1[,"nf_veg_2022"] + I(dados_df1[,"nf_veg_2022"]^2))

call <- list(); call[[1]]=m0; call[[2]]=m1; call[[3]]=m2
(aic <- aictab(cand.set=call, modnames=c("null", "linear", "quadratic"), sort=F))

ggplot(dados_df1, aes(x = nf_veg_2022, y = saevissima_2023)) +
  geom_point(color = "#483d8d", size = 3, alpha= 0.5) +  
  scale_y_continuous(breaks = seq(0, 120, by = 20), limits = c(0, 120)) + theme_classic() +
  theme(axis.text.x = element_text(size = 12, color = "black"),
        axis.text.y = element_text(size = 12, color = "black"),
        axis.title.x = element_text(color = "transparent"),
        axis.title.y = element_text(color = "transparent"))

ggsave("Fig_6A2023.pdf", width = 8, height = 9)

### S. saevissima x Area of fragments of native vegetation 2024

m0 <- glm(dados_df1[,"saevissima_2023"]~1)
m1 <- glm(dados_df1[,"saevissima_2023"]~dados_df1[,"area_veg_2022"])
m2 <- glm(dados_df1[,"saevissima_2023"]~dados_df1[,"area_veg_2022"] + I(dados_df1[,"area_veg_2022"]^2))

call <- list(); call[[1]]=m0; call[[2]]=m1; call[[3]]=m2
(aic <- aictab(cand.set=call, modnames=c("null", "linear", "quadratic"), sort=F))


ggplot(dados_df1, aes(x = area_veg_2022, y = saevissima_2023)) +
  geom_point(color = "#483d8d", size = 3, alpha= 0.5) +  
  geom_smooth(method = "lm", se = TRUE, fill = "#483d8d", color = "#483d8d", alpha = 0.2, size = 1.5) +  
  scale_y_continuous(breaks = seq(0, 80, by = 20),  
                     limits = c(0, 80)) +    scale_x_continuous(labels = label_number(accuracy = 1), breaks = pretty(dados_df$area_veg_2022, n = 3)) +   
  theme_classic() +
  theme(axis.text.x = element_text(size = 12, color = "black"),
        axis.text.y = element_text(size = 12, color = "black"),
        axis.title.x = element_text(color = "transparent"),
        axis.title.y = element_text(color = "transparent"))

ggsave("Fig_6B2023.pdf", width = 8, height = 9)

### S. saevissima x Number of anthropic fragments 2024  

m0 <- glm(dados_df1[,"saevissima_2023"]~1)
m1 <- glm(dados_df1[,"saevissima_2023"]~dados_df1[,"nf_urb_2022"])
m2 <- glm(dados_df1[,"saevissima_2023"]~dados_df1[,"nf_urb_2022"] + I(dados_df1[,"nf_urb_2022"]^2))

call <- list(); call[[1]]=m0; call[[2]]=m1; call[[3]]=m2
(aic <- aictab(cand.set=call, modnames=c("null", "linear", "quadratic"), sort=F))

ggplot(dados_df1, aes(x = nf_urb_2022, y = saevissima_2023)) +
  geom_point(color = "#483d8d", size = 3,alpha= 0.5) + 
  geom_smooth(aes(ymin = ifelse(..ymin.. < 0, 0, ..ymin..)),method = "lm", formula = y ~ poly(x, 2), se = TRUE, fill = "#483d8d", color = "#483d8d", alpha = 0.2,  size = 1.5) +
  scale_y_continuous(breaks = seq(0, 70, by = 10),  
                     limits = c(0, 70)) + 
  theme_classic() +
  theme(axis.text.x = element_text(size = 12, color = "black"),
        axis.text.y = element_text(size = 12, color = "black"),
        axis.title.x = element_text(color = "transparent"),
        axis.title.y = element_text(color = "transparent"))

ggsave("Fig_6C2023.pdf", width = 8, height = 9)

####  S. saevissima x Area of anthropic fragments 2024  

m0 <- glm(dados_df1[,"saevissima_2023"]~1)
m1 <- glm(dados_df1[,"saevissima_2023"]~dados_df1[,"area_urb_2022"])
m2 <- glm(dados_df1[,"saevissima_2023"]~dados_df1[,"area_urb_2022"] + I(dados_df1[,"area_urb_2022"]^2))

call <- list(); call[[1]]=m0; call[[2]]=m1; call[[3]]=m2
(aic <- aictab(cand.set=call, modnames=c("null", "linear", "quadratic"), sort=F))

ggplot(dados_df1, aes(x = area_urb_2022, y = saevissima_2023)) +
  geom_point(color = "#483d8d", size = 3,alpha= 0.5) + 
  geom_smooth(method = "lm", formula = y ~ poly(x, 2), se = TRUE, fill = "#483d8d", color = "#483d8d", alpha = 0.2,  size = 1.5) +
  scale_y_continuous(breaks = seq(0, 70, by = 10),  
                     limits = c(0, 70)) +  scale_x_continuous(labels = label_number(accuracy = 1), breaks = pretty(dados_df$area_urb_2022, n = 3)) +   
  theme_classic() +
  theme(axis.text.x = element_text(size = 12, color = "black"),
        axis.text.y = element_text(size = 12, color = "black"),
        axis.title.x = element_text(color = "transparent"),
        axis.title.y = element_text(color = "transparent"))


ggsave("Fig_6D2023.pdf", width = 8, height = 9)
