########################################
# Script for generating results for manuscript:
#        Modelling extreme desiccation tolerance in a marine tardigrade
#
# Thomas L. Hygum, Robyn Margaret Stuart, Aslak Jørgensen, Nadja Møbjerg
# Date last modified: Feb 28, 2018
########################################

########################################
# Load libraries
########################################
library (nlme)
library (drc)
library (readr)
library (ggplot2)
library(readxl)
library (lsmeans)

########################################
# Set global variables and read in data
########################################
setwd("~/Documents/git/tardigrades")
desiccation <- read_excel("desiccation-data.xlsx", sheet="5 - All data")
controls <- read_excel("desiccation-data.xlsx", sheet = "Controls")

desiccation <- desiccation[-165,1:8]
desiccation <- desiccation[desiccation$Surface!="Filter paper dH20",]
desiccation$Surface = as.factor(desiccation$Surface)

st.err <- function(x) { sd(x)/sqrt(length(x)) }

########################################
## Strand 1 - CONTROLS
########################################
means <- aggregate(controls$Prop, list(controls$Day), mean)
stdevs <- aggregate(controls$Prop, list(controls$Day), st.err)
(cbind(means,stdevs[,2]))

# A Wilcoxon rank sum test...
t.test(controls$Prop[controls$Day==1],controls$Prop[controls$Day==3])
wilcox.test(controls$Prop[controls$Day==1],controls$Prop[controls$Day==3])

t.test(controls$Prop[controls$Day==1],controls$Prop[controls$Day==2])
wilcox.test(controls$Prop[controls$Day==1],controls$Prop[controls$Day==2])

(g<- ggplot(aes(y = Prop, x = as.factor(Day)), data = controls)+ 
    geom_boxplot() + 
    geom_jitter(aes(y = Prop, x = as.factor(Day)), shape=16, position=position_jitter(0.2)) + 
    theme_bw() + 
    ylab("Proportion active") +
    xlab("Day") +
    theme(axis.text    = element_text(face = "bold"),
          legend.title = element_blank()) +
    coord_cartesian(ylim = c(0, 1)) )
fn <- "Figs/Fig2.png"
ggsave(fn, g, device="png", dpi=1200)


########################################
## Strand 2
########################################
desshort <- desiccation[desiccation$DesiccationLength<49&desiccation$Replicate!="TH-FP-PP-2-001",]
meansdes <- aggregate(desshort$Prop48, list(desshort$DesiccationLength, desshort$Surface), mean)
stdevsdes <- aggregate(desshort$Prop48, list(desshort$DesiccationLength, desshort$Surface), st.err)
cbind(meansdes,stdevsdes[,3])

desshort$DesLenFactor = as.factor(desshort$DesiccationLength)
glm48 <- glm(cbind(Active48,Inactive48) ~ DesLenFactor+Surface, family=binomial, data=desshort)
plot(glm48)
drop1(glm48, test="Chisq")

# Test for significance of surface
desshortglass <- subset(desshort,Surface=="Glass")
glm48surf <- glm(cbind(Active48,Inactive48) ~ DesLenFactor+WaterVolume, family=binomial, data=desshortglass)
drop1(glm48, test="Chisq")

(sur = cld(lsmeans(glm48,~Surface:DesLenFactor, type="response")))

g <- ggplot(data = desshort) + 
  geom_boxplot(aes(y = Prop48, x = as.factor(DesiccationLength))) + 
  geom_jitter(aes(y = Prop48, x = as.factor(DesiccationLength)), shape=16, position=position_jitter(0.2)) + 
  facet_wrap(~Surface) + 
  theme_bw() + 
  ylab("Proportion active after 48h") +
  xlab("Time spent desiccated") +
  theme(axis.text    = element_text(face = "bold"),
        legend.title = element_blank()) +
  coord_cartesian(ylim = c(0, 1)) 
g
fn <- "Figs/Fig3.png"
ggsave(fn, g, device="png", dpi=1200)


########################
# Strand 3
########################
desmed <- desiccation[desiccation$DesiccationLength<337&desiccation$Surface!="Hedgehog",]
meansdes <- aggregate(desmed$Prop48, list(desmed$DesiccationLength, desmed$Surface), mean)
stdevsdes <- aggregate(desmed$Prop48, list(desmed$DesiccationLength, desmed$Surface), st.err)
cbind(meansdes,stdevsdes[,3])

# Figure showing the results by water volume
g <- ggplot(data = subset(desmed,Surface=="Glass")) + 
  geom_boxplot(aes(y = Prop48, x = as.factor(DesiccationLength))) + 
  geom_jitter(aes(y = Prop48, x = as.factor(DesiccationLength)), shape=16, position=position_jitter(0.2)) + 
  facet_wrap(~WaterVolume) + 
  theme_bw() + 
  ylab("Proportion active after 48h") +
  xlab("Hours spent desiccated") +
  theme(axis.text    = element_text(face = "bold"),
        legend.title = element_blank()) +
  coord_cartesian(ylim = c(0, 1)) 
g
fn <- "Figs/Fig4.png"
ggsave(fn, g, device="png", dpi=1200)

desmed$DesLenFactor = as.factor(desmed$DesiccationLength)
glm2a <- glm(cbind(Active48,Inactive48) ~ DesLenFactor+WaterVolume, family=binomial, data=subset(desmed,Surface=="Glass"))
drop1(glm2a, test="Chisq")
(sur = cld(lsmeans(glm2a,~DesLenFactor:WaterVolume, type="response")))

glm2 <- glm(cbind(Active48,Inactive48) ~ DesLenFactor+Surface, family=binomial, data=desmed)
drop1(glm2, test="Chisq")
(sur = cld(lsmeans(glm2,~Surface:DesLenFactor, type="response")))
contrasts(lsmeans(glm2,"Surface", type="response"))

g <- ggplot(data = desmed) + 
  geom_boxplot(aes(y = Prop48, x = as.factor(DesiccationLength))) + 
  geom_jitter(aes(y = Prop48, x = as.factor(DesiccationLength)), shape=16, position=position_jitter(0.2)) + 
  facet_wrap(~Surface) + 
  theme_bw() + 
  ylab("Proportion active after 48h") +
  xlab("Hours spent desiccated") +
  theme(axis.text    = element_text(face = "bold"),
        legend.title = element_blank()) +
  coord_cartesian(ylim = c(0, 1)) 
g
fn <- "Figs/Fig5.png"
ggsave(fn, g, device="png", dpi=1200)




########################
# Strand 4-1
########################
deslong <- desiccation[desiccation$Surface=="Filter paper",]
meansdes <- aggregate(deslong$Prop48, list(deslong$DesiccationLength), mean)
mediansdes <- aggregate(deslong$Prop48, list(deslong$DesiccationLength), median)
stdevsdes <- aggregate(deslong$Prop48, list(deslong$DesiccationLength), st.err)
cbind(meansdes,stdevsdes[,2])

glm3 <- glm(cbind(Active48,Inactive48) ~ factor(DesiccationLength), family=binomial, data=deslong)
drop1(glm3, test="Chisq")

(sur = cld(lsmeans(glm3,~DesiccationLength, type="response")))

g <- ggplot(data = deslong) + 
  geom_boxplot(aes(y = Prop48, x = as.factor(DesiccationLength))) + 
  geom_jitter(aes(y = Prop48, x = as.factor(DesiccationLength)), shape=16, position=position_jitter(0.2)) + 
  theme_bw() + 
  ylab("Proportion active after 48h") +
  xlab("Hours spent desiccated") +
  theme(axis.text    = element_text(face = "bold"),
        legend.title = element_blank()) +
  coord_cartesian(ylim = c(0, 1)) 
g
fn <- "Figs/Fig6.png"
ggsave(fn, g, device="png", dpi=1200)

########################
# Strand 4-2
########################
## Make a dose-reponse model
drm1 <- drm(data=deslong, Prop48~DesiccationLength, fct=LL.3(names = c("b","d","ed50")), na.action = na.omit)

## Model summary
summary(drm1)
print(drm1)
plot(fitted(drm1), residuals(drm1))
hist(residuals(drm1))
cbind(estimate=coef(drm1),confint(drm1))
modelFit(drm1)

## Make a logit model for comparison
glm3 <- glm(cbind(Active48,Inactive48) ~ DesiccationLength, family=binomial, data=deslong)
mselect(drm1, sorted = c("IC", "Res var", "Lack of fit", "no"), linreg = FALSE, icfct = AIC)
summary(glm3)

## Make predictions
new.dat <- data.frame(DesiccationLength=seq(0,8760,10))
pm <- predict(drm1, newdata=new.dat, interval="confidence") 
new.dat$p <- pm[,1]
new.dat$pmin <- pm[,2]
new.dat$pmax <- pm[,3]

## Plot
ggplot(deslong, aes(x = DesiccationLength, y = Prop48)) +
  geom_point() +
  geom_line(data=new.dat, aes(x=DesiccationLength, y=p)) +
  geom_ribbon(data=new.dat, aes(x=DesiccationLength, y=p, ymin=pmin, ymax=pmax), alpha=0.2)+
  labs(x="Time spent desiccated", y="Proportion active after 48 hours")  + 
  theme_bw()

fn <- "Figs/Fig7.png"
ggsave(fn, device="png", dpi=1200)
