########################################
# Script for analysing tardigrade heat tolerance
#
# Ricardo C. Neves, Thomas L. Sørensen-Hygum, Robyn M Stuart, Nadja Møbjerg
# Date last modified: Apr 17, 2020
########################################

########################################
# Load libraries
########################################
library (nlme)
library (ggplot2)
library (emmeans)
library (dplyr)
library (ggpubr)

########################################
# Set global variables and read in data
########################################
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # Set working directory
source("emmeans_LD.R") # Read in emmeans file

# Read in data
heatingdess_168h <- read.csv2("Data/Desiccation & heating 1 week.csv", comment.char="#", stringsAsFactors=FALSE)
heatingdess_24h  <- read.csv2("Data/heatingdess_24h.csv", comment.char="#", stringsAsFactors=FALSE)


heatingdess_24h$exposure = "24 hours"
heatingdess_168h$exposure = "1 week"
heatingdess_168h$temp[heatingdess_168h$temp==5] = 23

# Combine into one dataset and just take the 48hr observations
alldata <- bind_rows(heatingdess_24h, heatingdess_168h)
alldata$time = factor(alldata$time)
alldata$exposure = factor(alldata$exposure)
alldata$propact = alldata$sumact/alldata$N
data48 <- filter(alldata, time == "48h")
data48$tempf <- factor(data48$temp)


########################################
## Modelling effects of exposure time on tardigrade activity
########################################
# Data visualization
(g<- ggplot(aes(y = sumact/N, x = tempf), data = data48)+ 
    geom_boxplot() + 
    geom_point(aes(y = sumact/N, x = as.factor(temp)), shape=16)+
    theme_pubclean() + 
   theme(text=element_text(size=16,  family="Optima")) + 
   ylab("Activity") +
    xlab("Temperature (ºC)") +
    facet_wrap(~exposure) +
    theme(axis.text    = element_text(face = "bold"),
          legend.title = element_blank()) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) #+
#    coord_cartesian(ylim = c(0, 1)) 
)
fn <- "Figs/activity_by_exposure.eps"
ggsave(fn, g, device="eps", dpi=300)

# First, make a model of the specimens at 5, 37, 40 degrees with temp as a factor
m1 <- glm(cbind(sumact,in.) ~ exposure + tempf, data=subset(data48,tempf%in%c("23","40","50","60","70")), family=binomial)
drop1(m1, test="Chisq") # Everything significant

# Second, make a model of the effects of acclimatization and temperature on proportion active with temp as a continuous variable
m2 <- glm(cbind(sumact,in.) ~ exposure + temp, data=data48, family=binomial)
drop1(m2, test="Chisq") # Everything significant

# Get table of parameter estimates
m2pars <- coef(summary(m2))[,1:2]
m2pars

LD50 <- emmeans_LD(m2,    ~exposure,  list(list(temp=0),list(temp=1)),  p=0.5)
myLD50 <- as.data.frame(LD50)
myLD50

pairs(emmeans(LD50,~exposure))

# Do a simple proportion test on the 60deg specimens
surviveto60 = sum(alldata$sumact[alldata$exposure=="24 hours" & alldata$time=="2h" & alldata$temp==60])
prop.test(x=c(0,surviveto60),n=c(100,100))

# Make predictions 
new.dat1 <- data.frame(temp=seq(5,100,.1), exposure="24 hours")
new.dat2 <- data.frame(temp=seq(5,100,.1), exposure="1 week")
new.dat <- rbind(new.dat1, new.dat2)

preds <- predict(m2, newdata = new.dat, type = "link", se.fit = TRUE)
critval <- 1.96 ## approx 95% CI
new.dat$pmax <- m2$family$linkinv(preds$fit + (critval * preds$se.fit))
new.dat$pmin <- m2$family$linkinv(preds$fit - (critval * preds$se.fit))
new.dat$p <- m2$family$linkinv(preds$fit)

## Plot the model predictions
(g <- ggplot() +
  geom_point(data=data48, aes(x=temp, y=sumact/N, col=exposure)) +
  geom_line(data=new.dat, aes(x=temp, y=p, col=exposure)) +
  geom_ribbon(data=new.dat, aes(x=temp, ymin=pmin, ymax=pmax, col=exposure), alpha=0.2)+
    theme_pubclean() + 
    theme(text=element_text(size=16,  family="Optima")) + 
    labs(x="Temperature (ºC)", y="Activity")  + 
  theme(legend.title = element_blank()) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) #+
)
(g <- g + theme(legend.title = element_blank()))
fn <- "Figs/activity_by_exposure.png_dr.pdf"
ggsave(fn, g, device="pdf", dpi=300)


