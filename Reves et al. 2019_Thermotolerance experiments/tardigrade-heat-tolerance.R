########################################
# Script for running statistical analyses for:
#     Neves RC, Hvidepil LKB, Sørensen-Hygum TL, Stuart RM & Møbjerg N. 
#     Thermotolerance experiments on active and desiccated states of 
#     Ramazzottius varieornatus emphasize that tardigrades are sensitive to high temperatures. 
#
# Date last modified: Dec 9, 2019
########################################

########################################
# Load libraries
########################################
library (nlme)
library (ggplot2)
library (emmeans)
library (dplyr)


########################################
# Set global variables and read in data
########################################
setwd("~/Documents/git/tardigrades")
print("Loading your data... ")
active <- read.csv2("Data/active_specimens.csv", comment.char="#")
acclimatize <- read.csv("Data/acclimatize.csv", sep=";")
heatingdess <- read.csv2("Data/heatingdess.csv", comment.char="#")
heatingdessshort <- read.csv2("Data/heatingshort_time_dess.csv", comment.char="#")
print("... data loaded.")

# Rename activity in the final dataset so it's the same as the others
names(heatingdessshort)[names(heatingdessshort) == "A."] <- "act."
names(heatingdessshort)[names(heatingdessshort) == "A"] <- "act"

# Add a variable specifying treatment type
active$treatment <- "No acclimatization"
acclimatize$treatment <- "Acclimatization"
heatingdess$treatment <- "Desiccation (24h)"
heatingdessshort$treatment <- "Desiccation (1h)"

# Combine into one dataset and just take the 48hr observations
alldata <- bind_rows(active, acclimatize, heatingdess, heatingdessshort)
data48 <- filter(alldata, time == "48h")
data48$tempf <- factor(data48$temp)


########################################
## Modelling effects of high temperatures and desiccation
########################################
dataacc <- filter(data48, treatment %in% c("No acclimatization", "Acclimatization"))

neworder <- c("No acclimatization", "Acclimatization")
dataacc <- arrange(mutate(dataacc,
                           treatment=factor(treatment,levels=neworder)),treatment)

# Data visualization
(g<- ggplot(aes(y = sumact/N, x = tempf), data = dataacc)+ 
    geom_boxplot() + 
    geom_point(aes(y = sumact/N, x = as.factor(temp)), shape=16)+
    theme_bw() + 
    ylab("Proportion active") +
    xlab("Temperature") +
    facet_wrap(~treatment) +
    theme(axis.text    = element_text(face = "bold"),
          legend.title = element_blank()) +
    coord_cartesian(ylim = c(0, 1)) )
fn <- "Figs/acclimatization_by_temp.png"
ggsave(fn, g, device="png", dpi=300)

# First, make a model of the specimens at 5, 37, 40 degrees with temp as a factor
modelacc <- glm(cbind(sumact,in.) ~ treatment + tempf, data=subset(dataacc,tempf%in%c("5","37","40")), family=binomial)
drop1(modelacc, test="Chisq") # Everything significant

# Second, make a model of the effects of acclimatization and temperature on proportion active with temp as a continuous variable
modelacc2 <- glm(cbind(sumact,in.) ~ treatment + temp, data=dataacc, family=binomial)
drop1(modelacc2, test="Chisq") # Everything significant

# Get table of parameter estimates
accparameters <- coef(summary(modelacc2))[,1:2]
accparameters

# Get the LD50s
-accparameters[1,1]/accparameters[3,1]
-(accparameters[1,1]+accparameters[2,1])/accparameters[3,1]

# Make predictions 
new.dat1 <- data.frame(temp=seq(5,60,.1), treatment="No acclimatization")
new.dat2 <- data.frame(temp=seq(5,60,.1), treatment="Acclimatization")
new.dat <- rbind(new.dat1, new.dat2)

preds <- predict(modelacc2, newdata = new.dat, type = "link", se.fit = TRUE)
critval <- 1.96 ## approx 95% CI
new.dat$pmax <- modelacc2$family$linkinv(preds$fit + (critval * preds$se.fit))
new.dat$pmin <- modelacc2$family$linkinv(preds$fit - (critval * preds$se.fit))
new.dat$p <- modelacc2$family$linkinv(preds$fit)

## Plot the model predictions
(g <- ggplot() +
  geom_point(data=dataacc, aes(x=temp, y=sumact/N, col=treatment)) +
  geom_line(data=new.dat, aes(x=temp, y=p, col=treatment)) +
  geom_ribbon(data=new.dat, aes(x=temp, ymin=pmin, ymax=pmax, col=treatment), alpha=0.2)+
  labs(x="Temperature", y="Proportion active")  + 
  theme_bw())
fn <- "Figs/acclimatization_by_temp_dr.jpg"
ggsave(fn, g, device="png", dpi=300)


########################################
## Modelling effects of high temperatures and desiccation
########################################
datades <- filter(data48, treatment %in% c("Desiccation (24h)", "Desiccation (1h)"))

# Boxplot according to different treatment types
(g<- ggplot(aes(y = sumact/N, x = tempf), data = datades)+ 
    geom_boxplot() + 
    geom_point(aes(y = sumact/N, x = as.factor(temp)), shape=16)+
    theme_bw() + 
    ylab("Proportion active") +
    xlab("Temperature") +
    facet_wrap(~treatment) +
    theme(axis.text    = element_text(face = "bold"),
          legend.title = element_blank()) +
    coord_cartesian(ylim = c(0, 1)) )
fn <- "Figs/desiccation_by_temp.png"
ggsave(fn, g, device="png", dpi=300)

# Do a simple proportion test on the 70deg specimens
prop.test(x=c(0,(19+19+20+19+18)),n=c(100,100))

# Make a model of the effects of desiccation time and temperature (factor) on proportion active
datadescut <- datades
datadescut <- filter(datades, (treatment=="Desiccation (1h)" & temp>23) | treatment=="Desiccation (24h)")

# Make a model of the effects of desiccation time and temperature (continuous) on proportion active
modeldes <- glm(cbind(sumact, in.) ~ temp + treatment, data=datadescut, family=binomial(link="logit"))
drop1(modeldes, test="Chisq") # Everything significant

# Get table of parameter estimates
desparameters <- coef(summary(modeldes))[,1:2]
desparameters

# Get the LD50s
-desparameters[1,1]/desparameters[2,1]
-(desparameters[1,1]+desparameters[3,1])/desparameters[2,1]

# Make predictions
new.dat1 <- data.frame(temp=seq(5,120,.1), treatment="Desiccation (24h)")
new.dat2 <- data.frame(temp=seq(5,120,.1), treatment="Desiccation (1h)")
new.dat <- rbind(new.dat1, new.dat2)
preds <- predict(modeldes, newdata = new.dat, se.fit = TRUE)
critval <- 1.96 ## approx 95% CI
new.dat$pmax <- modeldes$family$linkinv(preds$fit + (critval * preds$se.fit))
new.dat$pmin <- modeldes$family$linkinv(preds$fit - (critval * preds$se.fit))
new.dat$p <- modeldes$family$linkinv(preds$fit)

(g<- ggplot() +
    geom_jitter(data=datades, aes(x=temp, y=sumact/N, col=treatment)) +
    geom_line(data=new.dat, aes(x=temp, y=p, col=treatment)) +
    geom_ribbon(data=new.dat, aes(x=temp, ymin=pmin, ymax=pmax, col=treatment), alpha=0.2)+
    labs(x="temperature", y="Proportion active")  + 
    theme_bw())
fn <- "Figs/desiccation_by_temp_dr.png"
ggsave(fn, g, device="jpg", dpi=300)


