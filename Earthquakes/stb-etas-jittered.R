library(etasFLP)
library(Rcpp)
library(scales)
library(ggplot2)
library(bayesianETAS)
library(cowplot)
library(gridExtra)
library(shotGroups)
library(SpatialEpi)
library(lubridate)
library(mapdata)
library(maps)
library(dplyr)
library(reshape)

data(world2HiresMapEnv)
w2hr <- map_data("world2Hires")
can.df <- data.frame("long"=(w2hr$long-360), "lat"=w2hr$lat, "group"=w2hr$group)
lat <- c(47.6062, 47.0379, 48.4284, 46.8523, 47.6588, 44.0582, 45.5051, 46.2087, 42.2249)
lon <- c(-122.3321, -122.9007, -123.3656, -121.7603, -117.4260, -121.3133, -122.6750, -119.1199, -121.7817)
city <- c("Seattle", "Olypmia", "Victoria", "Mt. Rainier", "Spokane", "Bend", "Portland", "Kennewick", "Klamath Falls")
cities_map <- data.frame(city, lon, lat, stringsAsFactors = FALSE)
states <- map_data("state")
usa <- map_data("usa")
canada <- map_data("worldHires", "Canada")
# setwd("C:/Users/Hank/Dropbox/PNW-Data-Processing/D. identify-swarms/Swarm Analysis/")
setwd("~/Google Drive/My Drive/M9/")
source("src/functions/functions-bayesian-ETAS.R")
source("src/functions/functions-bETAS-summaries.R")
source("src/functions/functions-to-calculate-lambda.R")
source("src/functions/functions-to-simulate-ETAS.R")
source("src/functions/functions-spatiotemporal-ETAS.R")
source("src/functions/functions-MLE-ETAS.R")
source("src/functions/functions-msmt-error.R")
# source("src/modeling-files/bayesianETAS/maxLikelihoodETAS.R")
# source("src/modeling-files/bayesianETAS/my-funcs/sampleETASPosterior-MS.R")
# sourceCpp("src/modeling-files/b-ETAS/bayesianETAS/my-funcs/bETAS-MS.cpp")
# sourceCpp("src/modeling-files/b-ETAS/bayesianETAS/my-funcs/st-bETAS-MS-v2.cpp")
sourceCpp("src/modeling-files/b-ETAS/bayesianETAS/my-funcs/st-bETAS-MS-v2-20210815.cpp")

source("~/Dropbox/PNW-Data-Processing/G. EDA/catalog-visualizations/src/lowM/functions-plot-catalog.R")

# Read in data and filter out different catalogs
# intl <- read.csv("data/processed-data/pnw-intl-aux-clusters-20210909.csv")
intl <- read.csv("data/processed-data/pnw-intl-ETAS-slab-lowM.csv")
intl <- intl[order(intl$time),]
intl.t <- intl %>% filter(lat >= 42, lat <= 49, lon >= -125, lon <= -116.5, mag >= 2, is.na(remove),
                          potential_duplicates_20210412 == F)
intl.aux <- intl %>% filter(lat >= 41, lat <= 50, lon >= -126, lon <= -115.5, mag >= 2, is.na(remove),
                            potential_duplicates_20210412 == F)

# Read in the spatial integrals Gevals (needed for the Bayesian computation)
pnw.Gevals.long <- read.csv("data/bayesianETAS-output/2020-experiments/st-model/pnw/Gevals-long-pnw-target.csv")
pnw.Gevals.mat <- as.matrix(pnw.Gevals.long)
pnw.Gevals.mat <- unname(pnw.Gevals.mat)
pnw.Nt.cat.noCswarms <- intl %>% filter(lat >= 45, lat <= 49, lon >= -125, lon <= -116.5, mag >= 2, 
                                        is.na(remove), date >= as.Date("1985-01-01"),
                                        !scheme2020 %in% c(1),
                                        potential_duplicates_20210412 == F)
# pnw.Nt.cat.noCswarms <- pnw.Nt.cat.noCswarms[order(pnw.Nt.cat.noCswarms$time),]
pnw.Nt.cat.noCswarms$time <- pnw.Nt.cat.noCswarms$time - min(pnw.Nt.cat.noCswarms$time)
pnw.Nt.cat.noCswarms$x <- pnw.Nt.cat.noCswarms$x - min(pnw.Nt.cat.noCswarms$x) + 1
pnw.Nt.cat.noCswarms$y <- pnw.Nt.cat.noCswarms$y - min(pnw.Nt.cat.noCswarms$y) + 1
which.pnw.Nt.cat.noCswarms <- which(intl.t$id %in% pnw.Nt.cat.noCswarms$id)
pnw.Nt.cat.noCswarms.Gevals.mat <- pnw.Gevals.mat[, c(1:2, which.pnw.Nt.cat.noCswarms+2)]

pnw.Gevals.aux.long <- read.csv("data/bayesianETAS-output/2020-experiments/st-model/pnw/Gevals-long-pnw-aux.csv")
pnw.Gevals.aux.mat <- as.matrix(pnw.Gevals.aux.long)
pnw.Gevals.aux.mat <- unname(pnw.Gevals.aux.mat)

# Analysis paraemters
cat3D.true <- c(0.1, 0.006, 2.2999, 0.05, 1.08, 0.1, 1.5)
sims <- 5000
numMCMCSamples <- 500 # 500=What Gordon uses
M0 <- 2.0; approx=FALSE
burnin <- 500
sims <- sims+burnin
maxT.pnw <- round(as.numeric(max(intl$time)-min(intl$time) + 1), 0)
maxT.pnw.target <- round(as.numeric(max(intl.t$time)-min(intl.t$time) + 1), 0)
cat2inits.ss <- c(0.1, 0.02, 0.935*log(10), 0.05, 1.08, 0.32, 1.5)
# cat2inits.ss <- c(0.1, 0.02, 0.935*log(10), 0.05, 1.08, 0.92, 1.1)

maxT.pnw.Nt <- as.numeric(as.Date("2019-01-01")-as.Date("1985-01-01") + 1)


intl.Nt.box <- latlong2grid(rbind(c(-125, 49), c(-116.5, 45)))
intl.Nt.box$x <- intl.Nt.box$x - min(intl.Nt.box$x)
intl.Nt.box$y <- intl.Nt.box$y - min(intl.Nt.box$y)
intl.Nt.area <- intl.Nt.box$x[2]*intl.Nt.box$y[1]


head(jLMs1.aux[, c("x", "y", "lon", "lat")])
head(intl.aux[, c("x", "y", "lon", "lat")])

summary(jLMs1.aux[jLMs1.aux$out.of.win,])
plot(jLMs1.aux[jLMs1.aux$out.of.win & jLMs1.aux$mag >=1.8, c("x", "y")])
plotAnyBubble(data = jLMs1.aux %>% filter(out.of.win, mag >=1.8), 
              var="mag", var.name="Magnitude", cat.name="Target PNW, No MC Swarms, Herr AND Merr is NA")

 
   
jLMs2.pnw.Nt.noCswarms <- simCatSTBETASPNW(init.cat=intl.aux, region.label="north", swarms="noCswarms", target=T, 
                                           random.seed=2, file.label="pnw-Nt-noCswarms-Nt-cat2initsSS-fixalpha")
runSTBETASPlotsPNW(cat=jLMs2.pnw.Nt.noCswarms, 
                   cat.label = "jLM2pnw-Nt-noCswarms-Nt-cat2initsSS-fixalpha", 
                   nsims=5000, this.M0=2, these.init.vals=cat2inits.ss, 
                   this.spat.area = intl.Nt.area, this.isAlphaFixed=T,
                   cat.Gevals.mat = pnw.Nt.cat.noCswarms.Gevals.mat, this.maxT=maxT.pnw.Nt, this.plot.loc="pnw/j-cats-pnw/north",
                   plot.title="PNW North Target (45-49N) JITTERED Events M2+, 1985-2018, Low Proposal Sds (=0.05), Self-Similar Alpha") 

jLMs3.pnw.Nt.noCswarms <- simCatSTBETASPNW(init.cat=intl.aux, region.label="north", swarms="noCswarms", target=T, 
                                           random.seed=3, file.label="pnw-Nt-noCswarms-Nt-cat2initsSS-fixalpha")
runSTBETASPlotsPNW(cat=jLMs3.pnw.Nt.noCswarms, 
                   cat.label = "jLM3pnw-Nt-noCswarms-Nt-cat2initsSS-fixalpha", 
                   nsims=5000, this.M0=2, these.init.vals=cat2inits.ss, 
                   this.spat.area = intl.Nt.area, this.isAlphaFixed=T,
                   cat.Gevals.mat = pnw.Nt.cat.noCswarms.Gevals.mat, this.maxT=maxT.pnw.Nt, this.plot.loc="pnw/j-cats-pnw/north",
                   plot.title="PNW North Target (45-49N) JITTERED Events M2+, 1985-2018, Low Proposal Sds (=0.05), Self-Similar Alpha") 

jLMs4.pnw.Nt.noCswarms <- simCatSTBETASPNW(init.cat=intl.aux, region.label="north", swarms="noCswarms", target=T, 
                                           random.seed=4, file.label="pnw-Nt-noCswarms-Nt-cat2initsSS-fixalpha")
runSTBETASPlotsPNW(cat=jLMs4.pnw.Nt.noCswarms, 
                   cat.label = "jLM4pnw-Nt-noCswarms-Nt-cat2initsSS-fixalpha", 
                   nsims=5000, this.M0=2, these.init.vals=cat2inits.ss, 
                   this.spat.area = intl.Nt.area, this.isAlphaFixed=T,
                   cat.Gevals.mat = pnw.Nt.cat.noCswarms.Gevals.mat, this.maxT=maxT.pnw.Nt, this.plot.loc="pnw/j-cats-pnw/north",
                   plot.title="PNW North Target (45-49N) JITTERED Events M2+, 1985-2018, Low Proposal Sds (=0.05), Self-Similar Alpha") 

jLMs5.pnw.Nt.noCswarms <- simCatSTBETASPNW(init.cat=intl.aux, region.label="north", swarms="noCswarms", target=T, 
                                           random.seed=5, file.label="pnw-Nt-noCswarms-Nt-cat2initsSS-fixalpha")
runSTBETASPlotsPNW(cat=jLMs5.pnw.Nt.noCswarms, 
                   cat.label = "jLM5pnw-Nt-noCswarms-Nt-cat2initsSS-fixalpha", 
                   nsims=5000, this.M0=2, these.init.vals=cat2inits.ss, 
                   this.spat.area = intl.Nt.area, this.isAlphaFixed=T,
                   cat.Gevals.mat = pnw.Nt.cat.noCswarms.Gevals.mat, this.maxT=maxT.pnw.Nt, this.plot.loc="pnw/j-cats-pnw/north",
                   plot.title="PNW North Target (45-49N) JITTERED Events M2+, 1985-2018, Low Proposal Sds (=0.05), Self-Similar Alpha") 

jLMs6.pnw.Nt.noCswarms <- simCatSTBETASPNW(init.cat=intl.aux, region.label="north", swarms="noCswarms", target=T, 
                                           random.seed=6, file.label="pnw-Nt-noCswarms-Nt-cat2initsSS-fixalpha")
runSTBETASPlotsPNW(cat=jLMs6.pnw.Nt.noCswarms, 
                   cat.label = "jLM6pnw-Nt-noCswarms-Nt-cat2initsSS-fixalpha", 
                   nsims=5000, this.M0=2, these.init.vals=cat2inits.ss, 
                   this.spat.area = intl.Nt.area, this.isAlphaFixed=T,
                   cat.Gevals.mat = pnw.Nt.cat.noCswarms.Gevals.mat, this.maxT=maxT.pnw.Nt, this.plot.loc="pnw/j-cats-pnw/north",
                   plot.title="PNW North Target (45-49N) JITTERED Events M2+, 1985-2018, Low Proposal Sds (=0.05), Self-Similar Alpha") 

jLMs7.pnw.Nt.noCswarms <- simCatSTBETASPNW(init.cat=intl.aux, region.label="north", swarms="noCswarms", target=T, 
                                           random.seed=7, file.label="pnw-Nt-noCswarms-Nt-cat2initsSS-fixalpha")
runSTBETASPlotsPNW(cat=jLMs7.pnw.Nt.noCswarms, 
                   cat.label = "jLM7pnw-Nt-noCswarms-Nt-cat2initsSS-fixalpha", 
                   nsims=5000, this.M0=2, these.init.vals=cat2inits.ss, 
                   this.spat.area = intl.Nt.area, this.isAlphaFixed=T,
                   cat.Gevals.mat = pnw.Nt.cat.noCswarms.Gevals.mat, this.maxT=maxT.pnw.Nt, this.plot.loc="pnw/j-cats-pnw/north",
                   plot.title="PNW North Target (45-49N) JITTERED Events M2+, 1985-2018, Low Proposal Sds (=0.05), Self-Similar Alpha") 

jLMs8.pnw.Nt.noCswarms <- simCatSTBETASPNW(init.cat=intl.aux, region.label="north", swarms="noCswarms", target=T, 
                                           random.seed=8, file.label="pnw-Nt-noCswarms-Nt-cat2initsSS-fixalpha")
runSTBETASPlotsPNW(cat=jLMs8.pnw.Nt.noCswarms, 
                   cat.label = "jLM8pnw-Nt-noCswarms-Nt-cat2initsSS-fixalpha", 
                   nsims=5000, this.M0=2, these.init.vals=cat2inits.ss, 
                   this.spat.area = intl.Nt.area, this.isAlphaFixed=T,
                   cat.Gevals.mat = pnw.Nt.cat.noCswarms.Gevals.mat, this.maxT=maxT.pnw.Nt, this.plot.loc="pnw/j-cats-pnw/north",
                   plot.title="PNW North Target (45-49N) JITTERED Events M2+, 1985-2018, Low Proposal Sds (=0.05), Self-Similar Alpha") 

jLMs9.pnw.Nt.noCswarms <- simCatSTBETASPNW(init.cat=intl.aux, region.label="north", swarms="noCswarms", target=T, 
                                           random.seed=9, file.label="pnw-Nt-noCswarms-Nt-cat2initsSS-fixalpha")
runSTBETASPlotsPNW(cat=jLMs9.pnw.Nt.noCswarms, 
                   cat.label = "jLM9pnw-Nt-noCswarms-Nt-cat2initsSS-fixalpha", 
                   nsims=5000, this.M0=2, these.init.vals=cat2inits.ss, 
                   this.spat.area = intl.Nt.area, this.isAlphaFixed=T,
                   cat.Gevals.mat = pnw.Nt.cat.noCswarms.Gevals.mat, this.maxT=maxT.pnw.Nt, this.plot.loc="pnw/j-cats-pnw/north",
                   plot.title="PNW North Target (45-49N) JITTERED Events M2+, 1985-2018, Low Proposal Sds (=0.05), Self-Similar Alpha") 

# Crustal only
pnw.Nt.cat.noCswarms <- intl %>% filter(lat >= 45, lat < 49, lon >= -125, lon <= -116.5, 
                                        mag >= 2, date >= as.Date("1985-01-01"),
                                        is.na(remove), potential_duplicates_20210412 == F,
                                        !scheme2020 %in% c(1))
crustal.Nt.cat.noCswarms <- pnw.Nt.cat.noCswarms %>% filter(Depth2ModS > 10| is.na(Depth2ModS) )
which.crustal.Nt.cat.noCswarms <- which(intl.t$id %in% crustal.Nt.cat.noCswarms$id)
crustal.Nt.noCswarms.Gevals.mat <- pnw.Gevals.mat[, c(1:2, which.crustal.Nt.cat.noCswarms+2)]
jLMs1.crustal.Nt.noCswarms <- simCatSTBETASPNW(init.cat=intl.aux, region.label="north", swarms="noCswarms", target=T, 
                                               random.seed=1, crustal = T, file.label="crustal-Nt-noCswarms-Nt-cat2initsSS-fixalpha")
runSTBETASPlotsPNW(cat=jLMs1.crustal.Nt.noCswarms, 
                   cat.label = "jLM1crustal-Nt-noCswarms-Nt-cat2initsSS-fixalpha", 
                   nsims=5000, this.M0=2, these.init.vals=cat2inits.ss, 
                   this.spat.area = intl.Nt.area, this.isAlphaFixed=T,
                   cat.Gevals.mat = crustal.Nt.noCswarms.Gevals.mat, this.maxT=maxT.pnw.Nt, this.plot.loc="pnw/j-cats-pnw/north",
                   plot.title="PNW North Target (45-49N) JITTERED Events M2+, 1985-2018, Low Proposal Sds (=0.05), Self-Similar Alpha") 

jLMs2.crustal.Nt.noCswarms <- simCatSTBETASPNW(init.cat=intl.aux, region.label="north", swarms="noCswarms", target=T, 
                                               random.seed=2, crustal = T, file.label="crustal-Nt-noCswarms-Nt-cat2initsSS-fixalpha")
runSTBETASPlotsPNW(cat=jLMs2.crustal.Nt.noCswarms, 
                   cat.label = "jLM2crustal-Nt-noCswarms-Nt-cat2initsSS-fixalpha", 
                   nsims=5000, this.M0=2, these.init.vals=cat2inits.ss, 
                   this.spat.area = intl.Nt.area, this.isAlphaFixed=T,
                   cat.Gevals.mat = crustal.Nt.noCswarms.Gevals.mat, this.maxT=maxT.pnw.Nt, this.plot.loc="pnw/j-cats-pnw/north",
                   plot.title="PNW North Target (45-49N) JITTERED Events M2+, 1985-2018, Low Proposal Sds (=0.05), Self-Similar Alpha") 

jLMs3.crustal.Nt.noCswarms <- simCatSTBETASPNW(init.cat=intl.aux, region.label="north", swarms="noCswarms", target=T, 
                                               random.seed=3, crustal = T, file.label="crustal-Nt-noCswarms-Nt-cat2initsSS-fixalpha")
runSTBETASPlotsPNW(cat=jLMs3.crustal.Nt.noCswarms, 
                   cat.label = "jLM3crustal-Nt-noCswarms-Nt-cat2initsSS-fixalpha", 
                   nsims=5000, this.M0=2, these.init.vals=cat2inits.ss, 
                   this.spat.area = intl.Nt.area, this.isAlphaFixed=T,
                   cat.Gevals.mat = crustal.Nt.noCswarms.Gevals.mat, this.maxT=maxT.pnw.Nt, this.plot.loc="pnw/j-cats-pnw/north",
                   plot.title="PNW North Target (45-49N) JITTERED Events M2+, 1985-2018, Low Proposal Sds (=0.05), Self-Similar Alpha") 

jLMs4.crustal.Nt.noCswarms <- simCatSTBETASPNW(init.cat=intl.aux, region.label="north", swarms="noCswarms", target=T, 
                                               random.seed=4, crustal = T, file.label="crustal-Nt-noCswarms-Nt-cat2initsSS-fixalpha")
runSTBETASPlotsPNW(cat=jLMs4.crustal.Nt.noCswarms, 
                   cat.label = "jLM4crustal-Nt-noCswarms-Nt-cat2initsSS-fixalpha", 
                   nsims=5000, this.M0=2, these.init.vals=cat2inits.ss, 
                   this.spat.area = intl.Nt.area, this.isAlphaFixed=T,
                   cat.Gevals.mat = crustal.Nt.noCswarms.Gevals.mat, this.maxT=maxT.pnw.Nt, this.plot.loc="pnw/j-cats-pnw/north",
                   plot.title="PNW North Target (45-49N) JITTERED Events M2+, 1985-2018, Low Proposal Sds (=0.05), Self-Similar Alpha") 

jLMs5.crustal.Nt.noCswarms <- simCatSTBETASPNW(init.cat=intl.aux, region.label="north", swarms="noCswarms", target=T, 
                                               random.seed=5, crustal = T, file.label="crustal-Nt-noCswarms-Nt-cat2initsSS-fixalpha")
runSTBETASPlotsPNW(cat=jLMs5.crustal.Nt.noCswarms, 
                   cat.label = "jLM5crustal-Nt-noCswarms-Nt-cat2initsSS-fixalpha", 
                   nsims=5000, this.M0=2, these.init.vals=cat2inits.ss, 
                   this.spat.area = intl.Nt.area, this.isAlphaFixed=T,
                   cat.Gevals.mat = crustal.Nt.noCswarms.Gevals.mat, this.maxT=maxT.pnw.Nt, this.plot.loc="pnw/j-cats-pnw/north",
                   plot.title="PNW North Target (45-49N) JITTERED Events M2+, 1985-2018, Low Proposal Sds (=0.05), Self-Similar Alpha") 

jLMs6.crustal.Nt.noCswarms <- simCatSTBETASPNW(init.cat=intl.aux, region.label="north", swarms="noCswarms", target=T, 
                                               random.seed=6, crustal = T, file.label="crustal-Nt-noCswarms-Nt-cat2initsSS-fixalpha")
runSTBETASPlotsPNW(cat=jLMs6.crustal.Nt.noCswarms, 
                   cat.label = "jLM6crustal-Nt-noCswarms-Nt-cat2initsSS-fixalpha", 
                   nsims=5000, this.M0=2, these.init.vals=cat2inits.ss, 
                   this.spat.area = intl.Nt.area, this.isAlphaFixed=T,
                   cat.Gevals.mat = crustal.Nt.noCswarms.Gevals.mat, this.maxT=maxT.pnw.Nt, this.plot.loc="pnw/j-cats-pnw/north",
                   plot.title="PNW North Target (45-49N) JITTERED Events M2+, 1985-2018, Low Proposal Sds (=0.05), Self-Similar Alpha") 

jLMs7.crustal.Nt.noCswarms <- simCatSTBETASPNW(init.cat=intl.aux, region.label="north", swarms="noCswarms", target=T, 
                                               random.seed=7, crustal = T, file.label="crustal-Nt-noCswarms-Nt-cat2initsSS-fixalpha")
runSTBETASPlotsPNW(cat=jLMs7.crustal.Nt.noCswarms, 
                   cat.label = "jLM7crustal-Nt-noCswarms-Nt-cat2initsSS-fixalpha", 
                   nsims=5000, this.M0=2, these.init.vals=cat2inits.ss, 
                   this.spat.area = intl.Nt.area, this.isAlphaFixed=T,
                   cat.Gevals.mat = crustal.Nt.noCswarms.Gevals.mat, this.maxT=maxT.pnw.Nt, this.plot.loc="pnw/j-cats-pnw/north",
                   plot.title="PNW North Target (45-49N) JITTERED Events M2+, 1985-2018, Low Proposal Sds (=0.05), Self-Similar Alpha") 

jLMs8.crustal.Nt.noCswarms <- simCatSTBETASPNW(init.cat=intl.aux, region.label="north", swarms="noCswarms", target=T, 
                                               random.seed=8, crustal = T, file.label="crustal-Nt-noCswarms-Nt-cat2initsSS-fixalpha")
runSTBETASPlotsPNW(cat=jLMs8.crustal.Nt.noCswarms, 
                   cat.label = "jL85crustal-Nt-noCswarms-Nt-cat2initsSS-fixalpha", 
                   nsims=5000, this.M0=2, these.init.vals=cat2inits.ss, 
                   this.spat.area = intl.Nt.area, this.isAlphaFixed=T,
                   cat.Gevals.mat = crustal.Nt.noCswarms.Gevals.mat, this.maxT=maxT.pnw.Nt, this.plot.loc="pnw/j-cats-pnw/north",
                   plot.title="PNW North Target (45-49N) JITTERED Events M2+, 1985-2018, Low Proposal Sds (=0.05), Self-Similar Alpha") 

jLMs9.crustal.Nt.noCswarms <- simCatSTBETASPNW(init.cat=intl.aux, region.label="north", swarms="noCswarms", target=T, 
                                               random.seed=9, crustal = T, file.label="crustal-Nt-noCswarms-Nt-cat2initsSS-fixalpha")
runSTBETASPlotsPNW(cat=jLMs9.crustal.Nt.noCswarms, 
                   cat.label = "jLM9crustal-Nt-noCswarms-Nt-cat2initsSS-fixalpha", 
                   nsims=5000, this.M0=2, these.init.vals=cat2inits.ss, 
                   this.spat.area = intl.Nt.area, this.isAlphaFixed=T,
                   cat.Gevals.mat = crustal.Nt.noCswarms.Gevals.mat, this.maxT=maxT.pnw.Nt, this.plot.loc="pnw/j-cats-pnw/north",
                   plot.title="PNW North Target (45-49N) JITTERED Events M2+, 1985-2018, Low Proposal Sds (=0.05), Self-Similar Alpha") 
