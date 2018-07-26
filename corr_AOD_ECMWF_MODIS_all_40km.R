
library(dplyr)
library(readr)
library(lubridate)
library(ggplot2)
library(tidyr)
library(splines)
library(plyr)


# bind all regridded satellite data from ECMWF
# load ECMWF + AOD(MODIS) data for year 2013; 2014; 2015; 2016

ECMWF_MODIS_2013 <- read_csv("Z:/_SHARED_FOLDERS/Air Quality/Phase 2/ECMWF_CAMS/2013/AOD_2013_MODIS_all_40km.csv")
ECMWF_MODIS_2014 <- read_csv("Z:/_SHARED_FOLDERS/Air Quality/Phase 2/ECMWF_CAMS/2014/AOD_2014_MODIS_all_40km.csv")
ECMWF_MODIS_2015 <- read_csv("Z:/_SHARED_FOLDERS/Air Quality/Phase 2/ECMWF_CAMS/2015/AOD_2015_MODIS_all_40km.csv")
ECMWF_MODIS_2016 <- read_csv("Z:/_SHARED_FOLDERS/Air Quality/Phase 2/ECMWF_CAMS/2016/AOD_2016_MODIS_all_40km.csv")

# bind all data together---only AOD
ECMWF_MODIS_AOD <- rbind(ECMWF_MODIS_2013, ECMWF_MODIS_2014, 
                     cbind(ECMWF_MODIS_2015[,1:5],ECMWF_MODIS_2015[,8: 12]),
                     cbind(ECMWF_MODIS_2016[1:5], ECMWF_MODIS_2016[,8:12]))


# bind all data together---also includes PM10 and PM2.5
ECMWF_MODIS_PM <- rbind(ECMWF_MODIS_2015, ECMWF_MODIS_2016)


#################################################################################
#################################################################################

### make correlation between ECMWF and MODIS Terra-AQUA data~~~~~~~~~~~~~~~~~~~~~~~~


## get the months of observations
ECMWF_MODIS_AOD$month <- factor(format(ECMWF_MODIS_AOD$Date, format = "%b"), levels = month.abb)

## Define seasons
ECMWF_MODIS_AOD$season <- character(length = nrow(ECMWF_MODIS_AOD))
ECMWF_MODIS_AOD$season[ECMWF_MODIS_AOD$month %in% month.abb[c(1:2)]] <- "winter"
ECMWF_MODIS_AOD$season[ECMWF_MODIS_AOD$month %in% month.abb[c(12)]] <- "winter"
ECMWF_MODIS_AOD$season[ECMWF_MODIS_AOD$month %in% month.abb[c(3:5)]] <- "spring"
ECMWF_MODIS_AOD$season[ECMWF_MODIS_AOD$month %in% month.abb[c(6:8)]] <- "summer"
ECMWF_MODIS_AOD$season[ECMWF_MODIS_AOD$month %in% month.abb[c(9:11)]] <- "fall"
ECMWF_MODIS_AOD$season <- factor(ECMWF_MODIS_AOD$season, levels = c("winter","spring","summer","fall"))


## get the months of observations
ECMWF_MODIS_PM$month <- factor(format(ECMWF_MODIS_PM$Date, format = "%b"), levels = month.abb)

## Define seasons
ECMWF_MODIS_PM$season <- character(length = nrow(ECMWF_MODIS_PM))
ECMWF_MODIS_PM$season[ECMWF_MODIS_PM$month %in% month.abb[c(1:2)]] <- "winter"
ECMWF_MODIS_PM$season[ECMWF_MODIS_PM$month %in% month.abb[c(12)]] <- "winter"
ECMWF_MODIS_PM$season[ECMWF_MODIS_PM$month %in% month.abb[c(3:5)]] <- "spring"
ECMWF_MODIS_PM$season[ECMWF_MODIS_PM$month %in% month.abb[c(6:8)]] <- "summer"
ECMWF_MODIS_PM$season[ECMWF_MODIS_PM$month %in% month.abb[c(9:11)]] <- "fall"
ECMWF_MODIS_PM$season <- factor(ECMWF_MODIS_PM$season, levels = c("winter","spring","summer","fall"))



####### ECMWF(AOD) vs MODIS(AOD) ###------------------------------------------

plot(ECMWF_MODIS_AOD$AOD_ECMWF, ECMWF_MODIS_AOD$AOD_MODIS)


# check your data~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
plot <- ggplot(ECMWF_MODIS_AOD, aes(season, AOD_ECMWF)) +
  theme_bw() +
  geom_boxplot() + 
  guides(fill=FALSE) +   # no legend
  ylim(0, 2) 
plot


# check your data~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
plot <- ggplot(ECMWF_MODIS_AOD, aes(season, AOD_MODIS)) +
  theme_bw() +
  geom_boxplot() + 
  guides(fill=FALSE)  +
  ylim(0, 2) 
plot


plot(ECMWF_MODIS_AOD$AOD_ECMWF, ECMWF_MODIS_AOD$AOD_MODIS)

ECMWF_MODIS_AOD <- ECMWF_MODIS_AOD %>%
  filter(AOD_ECMWF < 2) %>%
  filter(AOD_MODIS < 3)

plot(ECMWF_MODIS_AOD$AOD_ECMWF, ECMWF_MODIS_AOD$AOD_MODIS)



#### fit function and label for AOD ----------------------------------------------
#### this funtion FORCE regression to pass through the origin #######################

lm_eqn <- function(df){
  m <- lm(AOD_MODIS ~ -1 + AOD_ECMWF, df);
  eq <- substitute(italic(y) == b %.% italic(x)*","~~italic(r)^2~"="~r2, 
                   list(b = format(coef(m)[1], digits = 2), 
                        r2 = format(summary(m)$r.squared, digits = 3)))
  as.character(as.expression(eq));                 
}


## this function includes the intercept~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# lm_eqn <- function(df){
#   m <- lm(AOD_MODIS ~  AOD_ECMWF, df);
#   eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2,
#                    list(a = format(coef(m)[2], digits = 2),
#                         b = format(coef(m)[1], digits = 2),
#                         r2 = format(summary(m)$r.squared, digits = 3)))
#   as.character(as.expression(eq));
# }
#####################################################################################

# plot with regression line-----


jpeg('Z:/_SHARED_FOLDERS/Air Quality/Phase 2/ECMWF_CAMS/MODIS_vs_AOD_2013_2016.jpg',    
     quality = 100, bg = "white", res = 200, width = 7, height = 7, units = "in")
par(mar=c(4, 10, 9, 2) + 0.3)
oldpar <- par(las=1)


# define regression equation for each season
eq <- ddply(ECMWF_MODIS_AOD, .(season),lm_eqn)


# ggplot(AOD_data, aes(x=AOD_ECMWF, y=AOD_MODIS, color = season)) +
  ggplot(ECMWF_MODIS_AOD, aes(y=AOD_MODIS, x=AOD_ECMWF)) +
  theme_bw() +
 # geom_point(size = 2) +
# geom_( bins = 30) +
  geom_jitter(colour=alpha("black",0.15) ) +
  facet_grid(season ~ .) +
  theme( strip.text = element_text(size = 18)) + 
  scale_color_manual(values = c("#ff0000", "#0000ff", "#000000", "#ffb732")) + 
  # geom_smooth(method="lm") +  # Add linear regression line
  geom_smooth(method = "lm", formula = y ~ -1 + x) +  # force fit through the origin
  ylab(expression("AOD (ECMWF)")) +
  xlab(expression("AOD (MODIS)")) +
  ylim(c(0,2)) + 
  xlim(c(0,2)) +
  theme(axis.title.y = element_text(face="bold", colour="black", size=12),
        axis.text.y  = element_text(angle=0, vjust=0.5, size=13)) +
  theme(axis.title.x = element_text(face="bold", colour="black", size=12),
        axis.text.x  = element_text(angle=0, vjust=0.5, size=13)) +
  
  geom_text(data = eq, aes(x = 1.7, y = 0.25, label = V1),
            parse = TRUE, inherit.aes=FALSE, size = 4, color = "black" ) +
  facet_grid(season ~.)


par(oldpar)
dev.off()

#############################################################################


####### ECMWF(PM2.5) vs ECMWF(AOD) ###------------------------------------------

plot(ECMWF_MODIS_PM$AOD_ECMWF, ECMWF_MODIS_PM$PM25_ECMWF)

# subtract DUST from AOD
# ECMWF_MODIS_PM$AOD_ECMWF <- ECMWF_MODIS_PM$AOD_ECMWF - ECMWF_MODIS_PM$DUST_ECMWF

# check your data~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
plot <- ggplot(ECMWF_MODIS_PM, aes(season, PM25_ECMWF)) +
  theme_bw() +
  geom_boxplot() + 
  guides(fill=FALSE) +   # no legend
  ylim(0, 800) 
plot


# check your data~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
plot <- ggplot(ECMWF_MODIS_PM, aes(season, AOD_ECMWF)) +
  theme_bw() +
  geom_boxplot() + 
  guides(fill=FALSE)  +
  ylim(0, 2.5) 
plot



ECMWF_MODIS_PM <- ECMWF_MODIS_PM %>%
  filter(PM25_ECMWF < 800) 

plot(ECMWF_MODIS_PM$AOD_ECMWF, ECMWF_MODIS_PM$PM25_ECMWF)


#####################################################################################

# plot with regression line-----PM2.5


lm_eqn <- function(df){
  m <- lm(PM25_ECMWF ~ -1 + AOD_ECMWF, df);
  eq <- substitute(italic(y) == b %.% italic(x)*","~~italic(r)^2~"="~r2, 
                   list(b = format(coef(m)[1], digits = 2), 
                        r2 = format(summary(m)$r.squared, digits = 3)))
  as.character(as.expression(eq));                 
}



jpeg('Z:/_SHARED_FOLDERS/Air Quality/Phase 2/ECMWF_CAMS/PM25_vs_AOD_2013_2016_ECMWF.jpg',    
     quality = 100, bg = "white", res = 200, width = 7, height = 7, units = "in")
par(mar=c(4, 10, 9, 2) + 0.3)
oldpar <- par(las=1)


# define regression equation for each season
eq <- ddply(ECMWF_MODIS_PM, .(season),lm_eqn)


# ggplot(AOD_data, aes(x=AOD_ECMWF, y=AOD_MODIS, color = season)) +
ggplot(ECMWF_MODIS_PM, aes(y=PM25_ECMWF, x=AOD_ECMWF)) +
  theme_bw() +
  # geom_point(size = 2) +
  # geom_( bins = 30) +
  geom_jitter(colour=alpha("black",0.15) ) +
  facet_grid(season ~ .) +
  theme( strip.text = element_text(size = 18)) + 
  scale_color_manual(values = c("#ff0000", "#0000ff", "#000000", "#ffb732")) + 
  # geom_smooth(method="lm") +  # Add linear regression line
  geom_smooth(method = "lm", formula = y ~ -1 + x) +  # force fit through the origin
  ylab(expression("PM2.5 (ECMWF)")) +
  xlab(expression("AOD (ECMWF)")) +
  ylim(c(0,800)) + 
  xlim(c(0,2)) +
  theme(axis.title.y = element_text(face="bold", colour="black", size=12),
        axis.text.y  = element_text(angle=0, vjust=0.5, size=13)) +
  theme(axis.title.x = element_text(face="bold", colour="black", size=12),
        axis.text.x  = element_text(angle=0, vjust=0.5, size=13)) +
  
  geom_text(data = eq, aes(x = 1.7, y = 30, label = V1),
            parse = TRUE, inherit.aes=FALSE, size = 4, color = "black" ) +
  facet_grid(season ~.)


par(oldpar)
dev.off()

####################################################################################
#####################################################################################

# plot with regression line-----PM10

plot(ECMWF_MODIS_PM$AOD_ECMWF, ECMWF_MODIS_PM$PM10_ECMWF)

ECMWF_MODIS_PM <- ECMWF_MODIS_PM %>%
  filter(PM10_ECMWF < 800) 

plot(ECMWF_MODIS_PM$AOD_ECMWF, ECMWF_MODIS_PM$PM10_ECMWF)



lm_eqn <- function(df){
  m <- lm(PM10_ECMWF ~ -1 + AOD_ECMWF, df);
  eq <- substitute(italic(y) == b %.% italic(x)*","~~italic(r)^2~"="~r2, 
                   list(b = format(coef(m)[1], digits = 2), 
                        r2 = format(summary(m)$r.squared, digits = 3)))
  as.character(as.expression(eq));                 
}



jpeg('Z:/_SHARED_FOLDERS/Air Quality/Phase 2/ECMWF_CAMS/PM10_vs_AOD_2013_2016_ECMWF.jpg',    
     quality = 100, bg = "white", res = 200, width = 7, height = 7, units = "in")
par(mar=c(4, 10, 9, 2) + 0.3)
oldpar <- par(las=1)


# define regression equation for each season
eq <- ddply(ECMWF_MODIS_PM, .(season),lm_eqn)


# ggplot(AOD_data, aes(x=AOD_ECMWF, y=AOD_MODIS, color = season)) +
ggplot(ECMWF_MODIS_PM, aes(y=PM10_ECMWF, x=AOD_ECMWF)) +
  theme_bw() +
  # geom_point(size = 2) +
  # geom_( bins = 30) +
  geom_jitter(colour=alpha("black",0.15) ) +
  facet_grid(season ~ .) +
  theme( strip.text = element_text(size = 18)) + 
  scale_color_manual(values = c("#ff0000", "#0000ff", "#000000", "#ffb732")) + 
  # geom_smooth(method="lm") +  # Add linear regression line
  geom_smooth(method = "lm", formula = y ~ -1 + x) +  # force fit through the origin
  ylab(expression("PM10 (ECMWF)")) +
  xlab(expression("AOD (ECMWF)")) +
  ylim(c(0,800)) + 
  xlim(c(0,2)) +
  theme(axis.title.y = element_text(face="bold", colour="black", size=12),
        axis.text.y  = element_text(angle=0, vjust=0.5, size=13)) +
  theme(axis.title.x = element_text(face="bold", colour="black", size=12),
        axis.text.x  = element_text(angle=0, vjust=0.5, size=13)) +
  
  geom_text(data = eq, aes(x = 1.7, y = 30, label = V1),
            parse = TRUE, inherit.aes=FALSE, size = 4, color = "black" ) +
  facet_grid(season ~.)


par(oldpar)
dev.off()


####################################################################################
### check dust component


# plot with regression line-----DUST AOD

plot(ECMWF_MODIS_AOD$DUST_ECMWF, ECMWF_MODIS_AOD$AOD_ECMWF)


lm_eqn <- function(df){
  m <- lm(AOD_ECMWF ~ -1 + DUST_ECMWF, df);
  eq <- substitute(italic(y) == b %.% italic(x)*","~~italic(r)^2~"="~r2, 
                   list(b = format(coef(m)[1], digits = 2), 
                        r2 = format(summary(m)$r.squared, digits = 3)))
  as.character(as.expression(eq));                 
}


# # this function includes the intercept~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# lm_eqn <- function(df){
#   m <- lm(AOD_ECMWF ~  DUST_ECMWF, df);
#   eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2,
#                    list(a = format(coef(m)[2], digits = 2),
#                         b = format(coef(m)[1], digits = 2),
#                         r2 = format(summary(m)$r.squared, digits = 3)))
#   as.character(as.expression(eq));
# }


jpeg('Z:/_SHARED_FOLDERS/Air Quality/Phase 2/ECMWF_CAMS/AOD_vs_DUST_2013_2016_ECMWF.jpg',    
     quality = 100, bg = "white", res = 200, width = 7, height = 7, units = "in")
par(mar=c(4, 10, 9, 2) + 0.3)
oldpar <- par(las=1)


# define regression equation for each season
eq <- ddply(ECMWF_MODIS_AOD, .(season),lm_eqn)


# ggplot(AOD_data, aes(x=AOD_ECMWF, y=AOD_MODIS, color = season)) +
ggplot(ECMWF_MODIS_AOD, aes(y=AOD_ECMWF, x=DUST_ECMWF)) +
  theme_bw() +
  # geom_point(size = 2) +
  # geom_( bins = 30) +
  geom_jitter(colour=alpha("black",0.15) ) +
  facet_grid(season ~ .) +
  theme( strip.text = element_text(size = 18)) + 
  scale_color_manual(values = c("#ff0000", "#0000ff", "#000000", "#ffb732")) + 
  # geom_smooth(method="lm") +  # Add linear regression line
  geom_smooth(method = "lm", formula = y ~ -1 + x) +  # force fit through the origin
  ylab(expression("AOD (ECMWF)")) +
  xlab(expression("DUST (ECMWF)")) +
  ylim(c(0,2)) + 
  xlim(c(0,2)) +
  theme(axis.title.y = element_text(face="bold", colour="black", size=12),
        axis.text.y  = element_text(angle=0, vjust=0.5, size=13)) +
  theme(axis.title.x = element_text(face="bold", colour="black", size=12),
        axis.text.x  = element_text(angle=0, vjust=0.5, size=13)) +
  
  geom_text(data = eq, aes(x = 1.7, y = 0.5, label = V1),
            parse = TRUE, inherit.aes=FALSE, size = 4, color = "black" ) +
  facet_grid(season ~.)


par(oldpar)
dev.off()


####################################################################################
### check downscale PM2.5 according to the DUST component


plot(ECMWF_MODIS_PM$AOD_ECMWF, ECMWF_MODIS_PM$PM25_ECMWF)
ECMWF_MODIS_PM$PM25_ECMWF_no_dust <- ECMWF_MODIS_PM$PM25_ECMWF * (1- 1/1.4)
ECMWF_MODIS_PM$PM10_ECMWF_no_dust <- ECMWF_MODIS_PM$PM10_ECMWF * (1- 1/1.4)

lm_eqn <- function(df){
  m <- lm(PM25_ECMWF_no_dust ~ -1 + AOD_ECMWF, df);
  eq <- substitute(italic(y) == b %.% italic(x)*","~~italic(r)^2~"="~r2, 
                   list(b = format(coef(m)[1], digits = 2), 
                        r2 = format(summary(m)$r.squared, digits = 3)))
  as.character(as.expression(eq));                 
}



jpeg('Z:/_SHARED_FOLDERS/Air Quality/Phase 2/ECMWF_CAMS/PM25_no_dust_vs_AOD_2015_2016_ECMWF.jpg',    
     quality = 100, bg = "white", res = 200, width = 8, height = 7, units = "in")
par(mar=c(4, 10, 9, 2) + 0.3)
oldpar <- par(las=1)


# define regression equation for each season
eq <- ddply(ECMWF_MODIS_PM, .(season),lm_eqn)


# ggplot(AOD_data, aes(x=AOD_ECMWF, y=AOD_MODIS, color = season)) +
ggplot(ECMWF_MODIS_PM, aes(y=PM25_ECMWF_no_dust, x=AOD_ECMWF)) +
  theme_bw() +
  # geom_point(size = 2) +
  # geom_( bins = 30) +
  geom_jitter(colour=alpha("black",0.15) ) +
  facet_grid(season ~ .) +
  theme( strip.text = element_text(size = 18)) + 
  scale_color_manual(values = c("#ff0000", "#0000ff", "#000000", "#ffb732")) + 
  # geom_smooth(method="lm") +  # Add linear regression line
  geom_smooth(method = "lm", formula = y ~ -1 + x) +  # force fit through the origin
  ylab(expression("PM25 no DUST (ECMWF)")) +
  xlab(expression("AOD (ECMWF)")) +
  ylim(c(0,200)) + 
  xlim(c(0,1.5)) +
  theme(axis.title.y = element_text(face="bold", colour="black", size=12),
        axis.text.y  = element_text(angle=0, vjust=0.5, size=13)) +
  theme(axis.title.x = element_text(face="bold", colour="black", size=12),
        axis.text.x  = element_text(angle=0, vjust=0.5, size=13)) +
  
  geom_text(data = eq, aes(x = 1.3, y = 10, label = V1),
            parse = TRUE, inherit.aes=FALSE, size = 4, color = "black" ) +
  facet_grid(season ~.)


par(oldpar)
dev.off()


############ ADDITIONAL CHECKS #########################################
########################################################################

library(gclus)
library(ggplot2)
library(lubridate)
library(plyr)
library(dplyr)
library(tidyr)
library(reshape2)
library(gtools)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# functions ######

regrline = function(x,y) {
  points(x,y,pch=".")
  # abline(line(x,y),col="blue")  ### add a smooth line
  abline(lsfit(x,y, intercept = FALSE),col="red")   ### add a regression line
}


panel.cor <- function(x, y, digits=2, cex.cor)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(1, 0, 1, 0))
#  r <- abs(cor(x,y))
  m <- lm(y ~ -1 + x)
  r2 = format(summary(m)$r.squared, digits = 3)
  txt <- format(c(r2, 0.12), digits=digits)[1]
 text(0.5, 0.5,cex = 2,  paste(txt))
}




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# all AOD-------------------------------------------------------------
data <- ECMWF_MODIS_AOD %>%
  select(-Date,
         -Lat,
         -Lon,
         -month,
         -season)


jpeg('Z:/_SHARED_FOLDERS/Air Quality/Phase 2/ECMWF_CAMS/all_correlations_AOD.jpg',    
     quality = 100, bg = "white", res = 200, width = 7, height = 7, units = "in")
par(mar=c(4, 10, 9, 2) + 0.3)
oldpar <- par(las=1)


# Basic Scatterplot Matrix
pairs(data, 
      lower.panel = regrline,
      upper.panel = panel.cor,
      main="ECMWF & MODIS AOD correlations (all seasons 2013-2016)",
      cex.labels = 1, font.labels = 1)

par(oldpar)
dev.off()



# all AOD + PM -------------------------------------------------------------
data <- ECMWF_MODIS_PM %>%
  select(-Date,
         -Lat,
         -Lon,
         -month,
         -season)


jpeg('Z:/_SHARED_FOLDERS/Air Quality/Phase 2/ECMWF_CAMS/all_correlations_AOD_&_PM.jpg',    
     quality = 100, bg = "white", res = 200, width = 10, height = 7, units = "in")
par(mar=c(4, 10, 9, 2) + 0.3)
oldpar <- par(las=1)


# Basic Scatterplot Matrix
pairs(data, 
      lower.panel = regrline,
      upper.panel = panel.cor,
      main="ECMWF & MODIS AOD & PM correlations (all seasons 2015-2016)",
      cex.labels = 0.7, font.labels = 1)

par(oldpar)
dev.off()








# # check Sulphate vs AOD and PM2.5 ######################################
# 
# plot(ECMWF_MODIS_PM$AOD_ECMWF, ECMWF_MODIS_PM$SO4_ECMWF)
# plot(ECMWF_MODIS_PM$AOD_MODIS, ECMWF_MODIS_PM$SO4_ECMWF)
# plot(ECMWF_MODIS_PM$PM25_ECMWF, ECMWF_MODIS_PM$SO4_ECMWF)
# plot(ECMWF_MODIS_PM$PM25_ECMWF_no_dust, ECMWF_MODIS_PM$SO4_ECMWF)
# plot(ECMWF_MODIS_PM$AOD_ECMWF, ECMWF_MODIS_PM$SALT_ECMWF)
# 
# 
# lm_eqn <- function(df){
#   m <- lm(SO4_ECMWF ~ -1 + AOD_ECMWF, df);
#   eq <- substitute(italic(y) == b %.% italic(x)*","~~italic(r)^2~"="~r2, 
#                    list(b = format(coef(m)[1], digits = 2), 
#                         r2 = format(summary(m)$r.squared, digits = 3)))
#   as.character(as.expression(eq));                 
# }
# 
# 
# 
# # define regression equation for each season
# eq <- ddply(ECMWF_MODIS_PM, .(season),lm_eqn)
# 
# 
# ggplot(ECMWF_MODIS_PM, aes(y=SO4_ECMWF, x=AOD_ECMWF)) +
#   theme_bw() +
#   geom_jitter(colour=alpha("black",0.15) ) +
#   facet_grid(season ~ .) +
#   theme( strip.text = element_text(size = 18)) + 
#   scale_color_manual(values = c("#ff0000", "#0000ff", "#000000", "#ffb732")) + 
#   geom_smooth(method = "lm", formula = y ~ -1 + x) +  # force fit through the origin
#   ylab(expression("SO4 (ECMWF)")) +
#   xlab(expression("AOD (ECMWF)")) +
#   ylim(c(0,0.5)) + 
#   xlim(c(0,2.5)) +
#   theme(axis.title.y = element_text(face="bold", colour="black", size=12),
#         axis.text.y  = element_text(angle=0, vjust=0.5, size=13)) +
#   theme(axis.title.x = element_text(face="bold", colour="black", size=12),
#         axis.text.x  = element_text(angle=0, vjust=0.5, size=13)) +
#   
#   geom_text(data = eq, aes(x = 1.9, y = 0.1, label = V1),
#             parse = TRUE, inherit.aes=FALSE, size = 4, color = "black" ) +
#   facet_grid(season ~.)
# 
# 
# 
# 
# lm_eqn <- function(df){
#   m <- lm(PM25_ECMWF_no_dust ~ -1 + SO4_ECMWF, df);
#   eq <- substitute(italic(y) == b %.% italic(x)*","~~italic(r)^2~"="~r2, 
#                    list(b = format(coef(m)[1], digits = 2), 
#                         r2 = format(summary(m)$r.squared, digits = 3)))
#   as.character(as.expression(eq));                 
# }
# 
# 
# 
# 
# # define regression equation for each season
# eq <- ddply(ECMWF_MODIS_PM, .(season),lm_eqn)
# 
# 
# ggplot(ECMWF_MODIS_PM, aes(y=PM25_ECMWF_no_dust, x=SO4_ECMWF)) +
#   theme_bw() +
#   geom_jitter(colour=alpha("black",0.15) ) +
#   facet_grid(season ~ .) +
#   theme( strip.text = element_text(size = 18)) + 
#   scale_color_manual(values = c("#ff0000", "#0000ff", "#000000", "#ffb732")) + 
#   geom_smooth(method = "lm", formula = y ~ -1 + x) +  # force fit through the origin
#   ylab(expression("PM25 no DUST (ECMWF)")) +
#   xlab(expression("SO4 (ECMWF)")) +
#   ylim(c(0,200)) + 
#   xlim(c(0,0.5)) +
#   theme(axis.title.y = element_text(face="bold", colour="black", size=12),
#         axis.text.y  = element_text(angle=0, vjust=0.5, size=13)) +
#   theme(axis.title.x = element_text(face="bold", colour="black", size=12),
#         axis.text.x  = element_text(angle=0, vjust=0.5, size=13)) +
#   
#   geom_text(data = eq, aes(x = 0.45, y = 20, label = V1),
#             parse = TRUE, inherit.aes=FALSE, size = 4, color = "black" ) +
#   facet_grid(season ~.)

