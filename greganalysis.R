library(nlrx)
library(sensitivity)
library(ggplot2)

#######################################################
#making population density relative to resource abundance

netlogopath2 <- file.path("/Program Files/NetLogo 6.2.2")
modelpath2 <- file.path("C:\\Users\\Marcy\\Documents\\GitHub\\primate-social-groups-model\\basemodelprimatesocialgroups.nlogo")
outpath2 <- file.path("/Users/Marcy/Desktop/")



hypo0metrics <- c("mean-distance-traveled",
                  "mean [timemove] of primates / ticks",
                  "mean [timeeat] of primates / ticks",
                  "count primates")

hypo0constants <- list("group-recog-module?" = "false",
                       "go-tests-on?" = "false",
                       "resource-tests-on?" = "false",
                       "move-tests-on?" = "false") 

hypo0vars <- list(
  "abundance" = list(min = 200000, max = 1200000, qfun = "qunif"),
  "clump-size" = list(min = 1, max = 500, qfun = "qunif"), #new single param for clumps
  "energy-per-capita" = list(min = 4000, max = 9000, qfun = "qunif"),
  "qual-mean" = list(min = 25, max = 150, qfun = "qunif"),
  "qual-sd" = list(min = 1, max = 20, qfun = "qunif"),
  "extraction-rate-mean" = list(min = 2, max = 7, qfun = "qunif"),
  "extraction-rate-sd" = list(min = 0, max = 3, qfun = "qunif"),
  
  "tgt-dist" = list(min = 5, max = 40, qfun = "qunif"),
  "tgt-neighbor" = list(min = 0, max = 7, qfun = "qunif"),
  
  "other-primate-detection-radius" = list(min = 50, max = 100, qfun = "qunif"),
  "resource-detection-radius" = list(min = 50, max = 100, qfun = "qunif"),
  "regrowth-rate" = list(min = 0.5, max = 1.0, qfun = "qunif"),
  "patch-regrowth-interval" = list(min = 500, max = 3000, qfun = "qunif"),
  "movement-noise" = list(min = 10, max = 45, qfun = "qunif"),
  "max-move" = list(min = 10, max = 50, qfun = "qunif"))



nlhypo0 <- nl(nlversion = "6.2.2",
              nlpath = netlogopath2,
              modelpath = modelpath2,
              jvmmem = 12000)

nlhypo0@experiment <- experiment(expname = "gregchhypo0MEE",
                                 outpath = outpath2,
                                 repetition = 1,
                                 tickmetrics = "false",
                                 idsetup = "setup",
                                 idgo = "go",
                                 runtime = 4300,
                                 stopcond = "(ticks > 5000)",
                                 metrics = hypo0metrics,
                                 variables = hypo0vars,
                                 constants = hypo0constants)


nlhypo0@simdesign <- simdesign_morris(nl = nlhypo0,
                                      morristype = "oat",
                                      morrislevels = 8,
                                      morrisr = 32,
                                      morrisgridjump = 4, 
                                      nseeds = 1)


progressr::handlers("progress")
library(future)
plan(multisession)
resultsgregchhypo0 <- progressr::with_progress(run_nl_all(nlhypo0))

setsim(nlhypo0, "simoutput")<-resultsgregchhypo0
analysisgregchhypo0 <- analyze_nl(nlhypo0)


analysisgregchhypo0 <- split(analysisgregchhypo0, analysisgregchhypo0$metric)

library(reshape2)

musigma_processing <- function(data) {
  return(dcast(data, metric + parameter ~ index))
}


analysisgregchhypo0scatter <- data.frame()

for (df in analysisgregchhypo0) {
  analysisgregchhypo0scatter <- rbind(analysisgregchhypo0scatter, (musigma_processing(as.data.frame(df))))
}

 
library(ggrepel)
library(ggplot2)

metric.levels <- unique(analysisgregchhypo0scatter$metric)

for (metric in metric.levels) {
  print(ggplot(analysisgregchhypo0scatter[analysisgregchhypo0scatter$metric == metric,], aes(x = mustar, y = sigma, label = parameter)) +
          geom_point() + geom_label_repel(size = 3) +
          labs(title = paste0("Morris Elementary Effects, Greg ", metric)))
}

min(resultsgregchhypo0$`mean-distance-traveled`) * 10 / (4300 / 24)
max(resultsgregchhypo0$`mean-distance-traveled`) * 10 / (4300 / 24)

hist(resultsgregchhypo0$`mean [timemove] of primates / ticks`)
hist(resultsgregchhypo0$`mean-distance-traveled` * 10 / (4300 / 24))

a<-ggplot(resultsgregchhypo0, mapping = aes(x = `tgt-neighbor`, y = `mean-distance-traveled` * 10 / (4300 / 24))) + 
  geom_jitter() + 
  stat_summary(fun = "median", color = "red")

b<-ggplot(resultsgregchhypo0, mapping = aes(x = `abundance`, y = `mean-distance-traveled`)) + 
  geom_jitter(size = 1) + 
  stat_summary(fun = "median", color = "red")

c<-ggplot(resultsgregchhypo0, mapping = aes(x = `extraction-rate-mean`, y = `mean-distance-traveled`)) + 
  geom_jitter(size = 1) + 
  stat_summary(fun = "median", color = "red")

d<-ggplot(resultsgregchhypo0, mapping = aes(x = `patch-regrowth-freq`, y = `mean-distance-traveled`)) + 
  geom_jitter(size = 1) + 
  stat_summary(fun = "median", color = "red")

a<-ggplot(resultsgregchhypo0, mapping = aes(x = `max-move`, y = `mean-distance-traveled` * 10 / (4300 / 24))) + 
  geom_jitter() + 
  stat_summary(fun = "median", color = "red")
a


resultsgregchhypo0.edit <- resultsgregchhypo0[resultsgregchhypo0$`max-move` > 25,]
hist(resultsgregchhypo0.edit$`mean [timemove] of primates / ticks`)
hist(resultsgregchhypo0.edit$`mean-distance-traveled` * 10 / (4300 / 24))




multiplot(a, c, b, d, cols = 2)

#heatmaps for daily travel distance
a <- ggplot(resultsgregchhypo0, aes(x = `abundance`, y = `extraction-rate-mean`, fill = `mean-distance-traveled` *10 / (4300 / 24))) +
  geom_tile()+
  #scale_y_discrete(labels = y.scale) +
  scale_fill_gradient(low = "lightblue", high = "navy", guide = "colorbar")

b <- ggplot(resultsgregchhypo0, aes(x = `abundance`, y = `energy-per-capita`, fill = `mean-distance-traveled` *10 / (4300 / 24))) +
  geom_tile()+
  #scale_y_discrete(labels = y.scale) +
  scale_fill_gradient(low = "lightblue", high = "navy", guide = "colorbar")

c <- ggplot(resultsgregchhypo0, aes(x = `abundance`, y = `qual-mean`, fill = `mean-distance-traveled` *10 / (4300 / 24))) +
  geom_tile()+
  #scale_y_discrete(labels = y.scale) +
  scale_fill_gradient(low = "lightblue", high = "navy", guide = "colorbar")

d <- ggplot(resultsgregchhypo0, aes(x = `energy-per-capita`, y = `qual-mean`, fill = `mean-distance-traveled` *10 / (4300 / 24))) +
  geom_tile()+
  #scale_y_discrete(labels = y.scale) +
  scale_fill_gradient(low = "lightblue", high = "navy", guide = "colorbar")

d <- ggplot(resultsgregchhypo0, aes(x = `tgt-neighbor`, y = `tgt-dist`, fill = `mean-distance-traveled` *10 / (4300 / 24))) +
  geom_tile()+
  #scale_y_discrete(labels = y.scale) +
  scale_fill_gradient(low = "lightblue", high = "navy", guide = "colorbar")
d

a <- ggplot(resultsgregchhypo0, aes(x = `max-move`, y = `extraction-rate-mean`, fill = `mean-distance-traveled` *10 / (4300 / 24))) +
  geom_tile()+
  #scale_y_discrete(labels = y.scale) +
  scale_fill_gradient(low = "lightblue", high = "navy", guide = "colorbar")
a

a <- ggplot(resultsgregchhypo0, aes(x = `max-move`, y = `tgt-neighbor`, fill = `mean-distance-traveled` *10 / (4300 / 24))) +
  geom_tile()+
  #scale_y_discrete(labels = y.scale) +
  scale_fill_gradient(low = "lightblue", high = "navy", guide = "colorbar")
a

a <- ggplot(resultsgregchhypo0, aes(x = `tgt-dist`, y = `extraction-rate-mean`, fill = `mean-distance-traveled` *10 / (4300 / 24))) +
  geom_tile()+
  #scale_y_discrete(labels = y.scale) +
  scale_fill_gradient(low = "lightblue", high = "navy", guide = "colorbar")
a

a <- ggplot(resultsgregchhypo0, aes(x = `max-move`, y = `tgt-dist`, fill = `mean-distance-traveled` *10 / (4300 / 24))) +
  geom_tile()+
  #scale_y_discrete(labels = y.scale) +
  scale_fill_gradient(low = "lightblue", high = "navy", guide = "colorbar")
a


##heatmaps for activity budget
ggplot(resultsgregchhypo0, aes(x = `patch-regrowth-freq`, y = `extraction-rate-mean`, fill = `mean [timemove] of primates / ticks`)) +
  geom_tile()+
  #scale_y_discrete(labels = y.scale) +
  scale_fill_gradient(low = "lightblue", high = "navy", guide = "colorbar")

ggplot(resultsgregchhypo0, aes(x = `patch-regrowth-freq`, y = `qual-mean`, fill = `mean [timemove] of primates / ticks`)) +
  geom_tile()+
  #scale_y_discrete(labels = y.scale) +
  scale_fill_gradient(low = "lightblue", high = "navy", guide = "colorbar")

ggplot(resultsgregchhypo0, aes(x = `qual-mean`, y = `extraction-rate-mean`, fill = `mean [timemove] of primates / ticks`)) +
  geom_tile()+
  #scale_y_discrete(labels = y.scale) +
  scale_fill_gradient(low = "lightblue", high = "navy", guide = "colorbar")





####################################################

netlogopath2 <- file.path("/Program Files/NetLogo 6.2.2")
modelpath2 <- file.path("C:\\Users\\Marcy\\Documents\\GitHub\\primate-social-groups-model\\basemodelprimatesocialgroups.nlogo")
outpath2 <- file.path("/Users/Marcy/Desktop/")


hypo1metrics <- c("foraging-efficiency-time", 
                  "foraging-efficiency-dist", 
                  "mean-distance-traveled",
                  "var-energy-intake-biweek", #I'm throwing the H2 ones in here two
                  "var-energy-intake-monthly",
                  "mean [timemove] of primates / ticks",
                  "count primates")

hypo1constants <- list("group-recog-module?" = "false",
                       "go-tests-on?" = "false",
                       "resource-tests-on?" = "false",
                       "move-tests-on?" = "false",
                       "qual-sd" = 1,
                       "other-primate-detection-radius" = 50,
                       "resource-detection-radius" = 50,
                       "regrowth-rate" = 1.0,
                       "movement-noise" = 10) 

hypo1vars <- list(
  "abundance" = list(min = 200000, max = 1200000, qfun = "qunif"),
  "clump-size" = list(min = 1, max = 1000, qfun = "qunif"), #new single param for clumps
  "energy-per-capita" = list(min = 4000, max = 9000, qfun = "qunif"),
  "qual-mean" = list(min = 25, max = 150, qfun = "qunif"),
  "extraction-rate-mean" = list(min = 4, max = 11, qfun = "qunif"),
  "extraction-rate-sd" = list(min = 1, max = 5, qfun = "qunif"),
  
  "tgt-dist" = list(min = 1, max = 40, qfun = "qunif"),
  "tgt-neighbor" = list(min = 0, max = 11, qfun = "qunif"),
  
  "patch-regrowth-freq" = list(min = 1000, max = 3000, qfun = "qunif"),
  "max-move" = list(min = 25, max = 50, qfun = "qunif"))



nlhypo1 <- nl(nlversion = "6.2.2",
         nlpath = netlogopath2,
         modelpath = modelpath2,
         jvmmem = 12000)

nlhypo1@experiment <- experiment(expname = "gregchhypo1MEE",
                                 outpath = outpath2,
                                 repetition = 1,
                                 tickmetrics = "false",
                                 idsetup = "setup",
                                 idgo = "go",
                                 runtime = 4300,
                                 stopcond = "(ticks > 5000)",
                                 metrics = hypo1metrics,
                                 variables = hypo1vars,
                                 constants = hypo1constants)


nlhypo1@simdesign <- simdesign_morris(nl = nlhypo1,
                                 morristype = "oat",
                                 morrislevels = 12,
                                 morrisr = 72,
                                 morrisgridjump = 6, 
                                 nseeds = 1)


progressr::handlers("progress")
library(future)
plan(multisession)
resultsgregchhypo1 <- progressr::with_progress(run_nl_all(nlhypo1))

setsim(nlhypo1, "simoutput")<-resultsgregchhypo1
analysisgregchhypo1 <- analyze_nl(nlhypo1)


analysisgregchhypo1 <- split(analysisgregchhypo1, analysisgregchhypo1$metric)

library(reshape2)

musigma_processing <- function(data) {
  return(dcast(data, metric + parameter ~ index))
}


analysisgregchhypo1scatter <- data.frame()

for (df in analysisgregchhypo1) {
  analysisgregchhypo1scatter <- rbind(analysisgregchhypo1scatter, (musigma_processing(as.data.frame(df))))
}


library(ggrepel)
library(ggplot2)

metric.levels <- unique(analysisgregchhypo1scatter$metric)



pdf("gregChh1resultsnov.pdf")
dev.off()
hist(resultsgregchhypo1$`mean [timemove] of primates / ticks`)
hist(resultsgregchhypo1$`mean-distance-traveled` * 10 / (4300 / 24))
hist(resultsgregchhypo1$`foraging-efficiency-dist`)


ggplot(resultsgregchhypo1, aes(x = `patch-regrowth-freq`, y = `extraction-rate-mean`, fill = 10 *`mean-distance-traveled` / (4300 / 24))) +
  geom_tile()+
  #scale_y_discrete(labels = y.scale) +
  scale_fill_gradient(low = "lightblue", high = "navy", guide = "colorbar")

#################FIGURE 3 ###################################
figure1data <- analysisgregchhypo1scatter[analysisgregchhypo1scatter$metric == "foraging-efficiency-dist_mean",]
figure1data$number <- c(1, 4, 5, 8, 3, 2, 9, 10, 7, 6)
figure1data$readable.params <- c("abundance", #1
                                 "clump size", #4
                                 "energy per capita", #5
                                 "extraction rate (mean)", #8
                                 "extraction rate (SD)", #3
                                 "movement speed", #2
                                 "frequency of patch regrowth",#9
                                 "mean patch quality",#10
                                 "target distance", #7
                                 "target neighbors") #6

ggplot(figure1data, aes(x = mustar, y = sigma, label = parameter)) +
  geom_point() + geom_text_repel(aes(label=ifelse(mustar>1.1, as.character(readable.params),number)), size = 5) +
  #labs(title = "Sensitivity analysis, scenario A") +
  xlab(paste("\u03BC", "*"))+
  ylab("\u03C3 ") +
  scale_x_continuous(expand = c(0.1,0.1))+
  theme(
    axis.title = element_text(size = 20, color = "black"),
    axis.text = element_text(size = 14, color = "black"),
    
    panel.grid = element_line(color = "black"),
    panel.grid.major = element_line(color = "gray75"),
    panel.grid.minor = element_line(color = "gray90"), panel.background = element_rect(fill = "white", color = "gray50"),
  )




for (metric in metric.levels) {
  print(ggplot(analysisgregchhypo1scatter[analysisgregchhypo1scatter$metric == metric,], aes(x = mustar, y = sigma, label = parameter)) +
          geom_point() + geom_label_repel(size = 3) +
          labs(title = paste0("Morris Elementary Effects, Greg ", metric)))
}
dev.off()


x.scale <- mapply(function(x){as.character(signif(as.numeric(x), digits = 3) / 1000)}, levels(resultsgregchhypo1$abundance))

b<-ggplot(resultsgregchhypo1, mapping = aes(x = `abundance`, y = 100 * `mean-distance-traveled`/180)) + 
  geom_violin() +
  stat_summary(fun = "median", color = "red") +
  scale_x_discrete(labels = x.scale)+
  labs(title = "B.") +
  xlab("abundance x 1,000")+
  ylab("foraging efficiency") +
  theme(
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 10),
    
    panel.grid = element_line(color = "black"),
    panel.grid.major = element_line(color = "gray75"),
    panel.grid.minor = element_line(color = "gray90"), panel.background = element_rect(fill = "white", color = "gray50"),
  )
b

###########
a<-ggplot(resultsgregchhypo1, mapping = aes(x = `tgt-neighbor`, y = `foraging-efficiency-time`)) + 
  geom_jitter(size = 1) + 
  stat_summary(fun = "median", color = "red")

b<-ggplot(resultsgregchhypo1, mapping = aes(x = `abundance`, y = `foraging-efficiency-time`)) + 
  geom_jitter(size = 1) + 
  stat_summary(fun = "median", color = "red")

c<-ggplot(resultsgregchhypo1, mapping = aes(x = `tgt-dist`, y = `foraging-efficiency-time`)) + 
  geom_jitter(size = 1) + 
  stat_summary(fun = "median", color = "red")

d<-ggplot(resultsgregchhypo1, mapping = aes(x = `patch-regrowth-freq`, y = `foraging-efficiency-time`)) + 
  geom_jitter(size = 1) + 
  stat_summary(fun = "median", color = "red")

e<-ggplot(resultsgregchhypo1, mapping = aes(x = `extraction-rate-mean`, y = `foraging-efficiency-time`)) + 
  geom_jitter(size = 1) + 
  stat_summary(fun = "median", color = "red")

f<-ggplot(resultsgregchhypo1, mapping = aes(x = `extraction-rate-sd`, y = `foraging-efficiency-time`)) + 
  geom_jitter(size = 1) + 
  stat_summary(fun = "median", color = "red")

multiplot(a, c, b, d, e, cols = 1)






##############
a<-ggplot(resultsgregchhypo1, mapping = aes(x = `tgt-neighbor`, y = 100 * `mean-distance-traveled`/180)) + 
  geom_jitter(size = 1) + 
  stat_summary(fun = "median", color = "red")

b<-ggplot(resultsgregchhypo1, mapping = aes(x = `tgt-dist`, y = 100 * `mean-distance-traveled`/180)) + 
  geom_jitter(size = 1) + 
  stat_summary(fun = "median", color = "red")

c<-ggplot(resultsgregchhypo1, mapping = aes(x = `abundance`, y = 100 * `mean-distance-traveled`/180)) + 
  geom_jitter(size = 1) + 
  stat_summary(fun = "median", color = "red")
d<-ggplot(resultsgregchhypo1, mapping = aes(x = `patch-regrowth-freq`, y = `mean-distance-traveled`)) + 
  geom_jitter(size = 1) + 
  stat_summary(fun = "median", color = "red")
e<-ggplot(resultsgregchhypo1, mapping = aes(x = `extraction-rate-mean`, y = `mean-distance-traveled`)) + 
  geom_jitter(size = 1) + 
  stat_summary(fun = "median", color = "red")

f<-ggplot(resultsgregchhypo1, mapping = aes(x = `extraction-rate-sd`, y = `mean-distance-traveled`)) + 
  geom_jitter(size = 1) + 
  stat_summary(fun = "median", color = "red")

multiplot(a, b, c, d, e, f, cols = 1)


multiplot(a, b, c, cols = 1)


a<-ggplot(resultsgregchhypo1, mapping = aes(x = `tgt-neighbor`, y = `foraging-efficiency-dist`)) + 
  geom_jitter(size = 1) + 
  stat_summary(fun = "median", color = "red")

b<-ggplot(resultsgregchhypo1, mapping = aes(x = `abundance`, y = `foraging-efficiency-dist`)) + 
  geom_jitter(size = 1) + 
  stat_summary(fun = "median", color = "red")

c<-ggplot(resultsgregchhypo1, mapping = aes(x = `patch-regrowth-freq`, y = `foraging-efficiency-dist`)) + 
  geom_jitter(size = 1) + 
  stat_summary(fun = "median", color = "red")

d<-ggplot(resultsgregchhypo1, mapping = aes(x = `tgt-dist`, y = `foraging-efficiency-dist`)) + 
  geom_jitter(size = 1) + 
  stat_summary(fun = "median", color = "red")

e<-ggplot(resultsgregchhypo1, mapping = aes(x = `extraction-rate-mean`, y = `foraging-efficiency-dist`)) + 
  geom_jitter(size = 1) + 
  stat_summary(fun = "median", color = "red")

f<-ggplot(resultsgregchhypo1, mapping = aes(x = `clump-size`, y = `foraging-efficiency-dist`)) + 
  geom_jitter(size = 1) + 
  stat_summary(fun = "median", color = "red")

#######################################################

resultsgregchhypo1$`tgt-neighbor` <- as.factor(resultsgregchhypo1$`tgt-neighbor`)
resultsgregchhypo1$`abundance` <- as.factor(resultsgregchhypo1$`abundance`)
resultsgregchhypo1$`patch-regrowth-freq` <- as.factor(resultsgregchhypo1$`patch-regrowth-freq`)
resultsgregchhypo1$`tgt-dist` <- as.factor(resultsgregchhypo1$`tgt-dist`)
resultsgregchhypo1$`extraction-rate-mean` <- as.factor(resultsgregchhypo1$`extraction-rate-mean`)
resultsgregchhypo1$`clump-size` <- as.factor(resultsgregchhypo1$`clump-size`)
resultsgregchhypo1$`energy-per-capita` <- as.factor(resultsgregchhypo1$`energy-per-capita`)

resultsgregchhypo1$`qual-mean` <- as.factor(resultsgregchhypo1$`qual-mean`)

a<-ggplot(resultsgregchhypo1, aes(x = `tgt-neighbor`, y = log(`foraging-efficiency-dist`))) + 
  geom_violin() + 
  stat_summary(fun = "median", color = "red")+
  labs(title = "A.") +
  xlab("target neighbor")+
  ylab("foraging efficiency") +
  theme(
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 10),
    
    panel.grid = element_line(color = "black"),
    panel.grid.major = element_line(color = "gray75"),
    panel.grid.minor = element_line(color = "gray90"), panel.background = element_rect(fill = "white", color = "gray50"),
  )
  a

x.scale <- mapply(function(x){as.character(signif(as.numeric(x), digits = 3) / 1000)}, levels(resultsgregchhypo1$abundance))
  
b<-ggplot(resultsgregchhypo1, mapping = aes(x = `abundance`, y = `foraging-efficiency-dist`)) + 
  geom_violin() +
  stat_summary(fun = "median", color = "red") +
  scale_x_discrete(labels = x.scale)+
  labs(title = "B.") +
  xlab("abundance x 1,000")+
  ylab("foraging efficiency") +
  theme(
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 10),
    
    panel.grid = element_line(color = "black"),
    panel.grid.major = element_line(color = "gray75"),
    panel.grid.minor = element_line(color = "gray90"), panel.background = element_rect(fill = "white", color = "gray50"),
  )
b

x.scale <- mapply(function(x){as.character(round(as.numeric(x)))}, levels(resultsgregchhypo1$`patch-regrowth-freq`))

c<-ggplot(resultsgregchhypo1, mapping = aes(x = `patch-regrowth-freq`, y = log(`foraging-efficiency-dist`))) + 
  geom_violin() +
  stat_summary(fun = "median", color = "red") +
  scale_x_discrete(labels = x.scale)+
  labs(title = "C.") +
  xlab("frequency of patch regrowth")+
  ylab("foraging efficiency") +
  theme(
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 8),
    
    panel.grid = element_line(color = "black"),
    panel.grid.major = element_line(color = "gray75"),
    panel.grid.minor = element_line(color = "gray90"), panel.background = element_rect(fill = "white", color = "gray50"),
  )

x.scale <- mapply(function(x){as.character(round(as.numeric(x), digits = 2))}, levels(resultsgregchhypo1$`tgt-dist`))


d<-ggplot(resultsgregchhypo1, mapping = aes(x = `tgt-dist`, y = log(`foraging-efficiency-dist`))) + 
  geom_violin() +
  stat_summary(fun = "median", color = "red") +
  scale_x_discrete(labels = x.scale)+
  labs(title = "D.") +
  xlab("target distance")+
  ylab("foraging efficiency") +
  theme(
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 10),
    
    panel.grid = element_line(color = "black"),
    panel.grid.major = element_line(color = "gray75"),
    panel.grid.minor = element_line(color = "gray90"), panel.background = element_rect(fill = "white", color = "gray50"),
  )

x.scale <-  mapply(function(x){as.character(round(as.numeric(x), digits = 2))}, levels(resultsgregchhypo1$`energy-per-capita`))

d<-ggplot(resultsgregchhypo1, mapping = aes(x = `energy-per-capita`, y = log(`foraging-efficiency-dist`))) + 
  geom_violin() +
  stat_summary(fun = "median", color = "red") +
  scale_x_discrete(labels = x.scale)+
  labs(title = "D.") +
  xlab("energy per agent")+
  ylab("foraging efficiency") +
  theme(
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 10),
    
    panel.grid = element_line(color = "black"),
    panel.grid.major = element_line(color = "gray75"),
    panel.grid.minor = element_line(color = "gray90"), panel.background = element_rect(fill = "white", color = "gray50"),
  )
d

x.scale <- mapply(function(x){as.character(round(as.numeric(x), digits = 2))}, levels(resultsgregchhypo1$`extraction-rate-mean`))

e<-ggplot(resultsgregchhypo1, mapping = aes(x = `extraction-rate-mean`, y = `foraging-efficiency-dist`)) + 
  geom_violin() +
  stat_summary(fun = "median", color = "red") +
  scale_x_discrete(labels = x.scale)+
  labs(title = "E.") +
  xlab("mean extraction rate")+
  ylab("foraging efficiency") +
  theme(
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 10),
    
    panel.grid = element_line(color = "black"),
    panel.grid.major = element_line(color = "gray75"),
    panel.grid.minor = element_line(color = "gray90"), panel.background = element_rect(fill = "white", color = "gray50"),
  )

x.scale <- mapply(function(x){as.character(round(as.numeric(x)))}, levels(resultsgregchhypo1$`clump-size`))

f<-ggplot(resultsgregchhypo1, mapping = aes(x = `clump-size`, y = `foraging-efficiency-dist`)) + 
  geom_violin() + 
  stat_summary(fun = "median", color = "red") +
  scale_x_discrete(labels = x.scale)+
  labs(title = "F.") +
  xlab("clump size")+
  ylab("foraging efficiency") +
  theme(
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 10),
    
    panel.grid = element_line(color = "black"),
    panel.grid.major = element_line(color = "gray75"),
    panel.grid.minor = element_line(color = "gray90"), panel.background = element_rect(fill = "white", color = "gray50"),
  )


g<-ggplot(resultsgregchhypo1, mapping = aes(x = `qual-mean`, y = `foraging-efficiency-dist`)) + 
  geom_violin() + 
  stat_summary(fun = "median", color = "red") +
  scale_x_discrete(labels = x.scale)+
  labs(title = "g.") +
  xlab("mean patch size")+
  ylab("foraging efficiency") +
  theme(
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 10),
    
    panel.grid = element_line(color = "black"),
    panel.grid.major = element_line(color = "gray75"),
    panel.grid.minor = element_line(color = "gray90"), panel.background = element_rect(fill = "white", color = "gray50"),
  )
g

multiplot(a, b, c, d, e, f, cols = 2)

multiplot(a, b, cols = 1)

### scales for heatmap axes
tgtdist.scale <- mapply(function(x){as.character(signif(as.numeric(x), digits = 3))}, levels(as.factor(resultsgregchhypo1$`tgt-dist`)))
tgtneighbor.breaks <- c(0,NA,2,NA,4,NA,6,NA,8,NA,10,NA)
extractionrate.scale <- mapply(function(x){as.character(signif(as.numeric(x), digits = 3))}, levels(as.factor(resultsgregchhypo1$`extraction-rate-mean`)))
regrowthfreq.scale <- mapply(function(x){as.character(signif(as.numeric(x), digits = 3))}, levels(as.factor(resultsgregchhypo1$`patch-regrowth-freq`)))
qualmean.scale <- mapply(function(x){as.character(signif(as.numeric(x), digits = 3))}, levels(as.factor(resultsgregchhypo1$`qual-mean`)))

everyother <- function(x) x[seq_along(x) %% 2 == 0]


## heatmap attempt
library(ggplot2)
ggplot(resultsgregchhypo1, aes(x = `tgt-neighbor`, y = `patch-regrowth-freq`, fill = `foraging-efficiency-dist`)) +
  geom_tile()

y.scale <- mapply(function(x){as.character(signif(as.numeric(x), digits = 3) / 1000)}, levels(resultsgregchhypo1$abundance))

########################Figure 4###################################

aa <- ggplot(resultsgregchhypo1, aes(x = `tgt-neighbor`, y = `tgt-dist`, fill = log(`foraging-efficiency-dist`))) +
  geom_tile()+
  scale_y_discrete(labels = tgtdist.scale, breaks = everyother) +
  scale_x_discrete(breaks = tgtneighbor.breaks)+
  scale_fill_gradient(low = "lightblue", high = "navy", guide = "colorbar")+
  labs(title = "A.", fill = "log \nforaging\nefficiency") +
  xlab("target neighbors")+
  ylab("target distance") +  
  theme(
    axis.title = element_text(size = 10, color = "black"),
    axis.text = element_text(size = 8, color = "black"),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 10),
    
    panel.grid = element_line(color = "black"),
    panel.grid.major = element_line(color = "gray75"),
    panel.grid.minor = element_line(color = "gray90"), panel.background = element_rect(fill = "white", color = "gray50"),
  )
a <- ggplot(resultsgregchhypo1, aes(x = `extraction-rate-mean`, y = `patch-regrowth-freq`, fill = log(`foraging-efficiency-dist`))) +
  geom_tile()+
  scale_y_discrete(labels = regrowthfreq.scale, breaks = everyother) +
  scale_x_discrete(labels = extractionrate.scale, guide = guide_axis(n.dodge = 1, angle = 45), breaks = everyother) +
  scale_fill_gradient(low = "lightblue", high = "navy", guide = "colorbar")+
  labs(title = "B.", fill = "log \nforaging\nefficiency") +
  xlab("extraction rate (mean)")+
  ylab("patch regrowth interval") +
  theme(
    axis.title = element_text(size = 10, color = "black"),
    axis.text = element_text(size = 8, color = "black"),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 10),
    
    panel.grid = element_line(color = "black"),
    panel.grid.major = element_line(color = "gray75"),
    panel.grid.minor = element_line(color = "gray90"), panel.background = element_rect(fill = "white", color = "gray50"),
  )
b <- ggplot(resultsgregchhypo1, aes(x = `extraction-rate-mean`, y = `qual-mean`, fill = log(`foraging-efficiency-dist`))) +
  geom_tile()+
  scale_y_discrete(labels = qualmean.scale, breaks = everyother) +
  scale_x_discrete(labels = extractionrate.scale, guide = guide_axis(n.dodge = 1, angle = 45), breaks = everyother) +
  scale_fill_gradient(low = "lightblue", high = "navy", guide = "colorbar")+
  labs(title = "C.", fill = "log \nforaging\nefficiency") +
  xlab("extraction rate (mean)")+
  ylab("patch quality") +
  theme(
    axis.title = element_text(size = 10, color = "black"),
    axis.text = element_text(size = 8, color = "black"),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 10),
    
    panel.grid = element_line(color = "black"),
    panel.grid.major = element_line(color = "gray75"),
    panel.grid.minor = element_line(color = "gray90"), panel.background = element_rect(fill = "white", color = "gray50"),
  )
bb <- ggplot(resultsgregchhypo1, aes(x = `tgt-neighbor`, y = `extraction-rate-mean`, fill = log(`foraging-efficiency-dist`))) +
  geom_tile()+
  scale_y_discrete(labels = extractionrate.scale, breaks = everyother) +
  scale_x_discrete(breaks = tgtneighbor.breaks)+
  scale_fill_gradient(low = "lightblue", high = "navy", guide = "colorbar")+
  labs(title = "D.", fill = "log \nforaging\nefficiency") +
  xlab("target neighbors")+
  ylab("extraction rate (mean)") +  
  theme(
    axis.title = element_text(size = 10, color = "black"),
    axis.text = element_text(size = 8, color = "black"),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 10),
    
    panel.grid = element_line(color = "black"),
    panel.grid.major = element_line(color = "gray75"),
    panel.grid.minor = element_line(color = "gray90"), panel.background = element_rect(fill = "white", color = "gray50"),
  )
d <- ggplot(resultsgregchhypo1, aes(x = `tgt-neighbor`, y = `patch-regrowth-freq`, fill = log(`foraging-efficiency-dist`))) +
  geom_tile()+
  scale_y_discrete(labels = regrowthfreq.scale, breaks = everyother) +
  scale_x_discrete(breaks = tgtneighbor.breaks)+
  scale_fill_gradient(low = "lightblue", high = "navy", guide = "colorbar")+
  labs(title = "E.", fill = "log \nforaging\nefficiency") +
  xlab("target neighbors")+
  ylab("patch regrowth interval") +  
  theme(
    axis.title = element_text(size = 10, color = "black"),
    axis.text = element_text(size = 8, color = "black"),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 10),
    
    panel.grid = element_line(color = "black"),
    panel.grid.major = element_line(color = "gray75"),
    panel.grid.minor = element_line(color = "gray90"), panel.background = element_rect(fill = "white", color = "gray50"),
  )
e <- ggplot(resultsgregchhypo1, aes(x = `tgt-dist`, y = `patch-regrowth-freq`, fill = log(`foraging-efficiency-dist`))) +
  geom_tile()+
  scale_y_discrete(labels = regrowthfreq.scale, breaks = everyother) +
  scale_x_discrete(labels = tgtdist.scale, breaks = everyother, guide = guide_axis(n.dodge = 1, angle = 45)) +
  scale_fill_gradient(low = "lightblue", high = "navy", guide = "colorbar")+
  labs(title = "F.", fill = "log \nforaging\nefficiency") +
  xlab("target distance")+
  ylab("patch regrowth interval") +  
  theme(
    axis.title = element_text(size = 10, color = "black"),
    axis.text = element_text(size = 8, color = "black"),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 10),
    
    panel.grid = element_line(color = "black"),
    panel.grid.major = element_line(color = "gray75"),
    panel.grid.minor = element_line(color = "gray90"), panel.background = element_rect(fill = "white", color = "gray50"),
  )
multiplot(aa, a, b, bb, d, e, cols = 2)
# width 270, height 400


##################figure 5##############################



c <- ggplot(resultsgregchhypo1, aes(x = `tgt-dist`, y = `extraction-rate-mean`, fill = log(`foraging-efficiency-dist`))) +
  geom_tile()+
  scale_y_discrete(labels = extractionrate.scale, breaks = everyother) +
  scale_x_discrete(labels = tgtdist.scale, breaks = everyother, guide = guide_axis(n.dodge = 1, angle = 45)) +
  scale_fill_gradient(low = "lightblue", high = "navy", guide = "colorbar")+
  labs(title = "C.", fill = "log \nforaging\nefficiency") +
  xlab("target distance")+
  ylab("extraction rate (mean)") +  
  theme(
    axis.title = element_text(size = 10, color = "black"),
    axis.text = element_text(size = 8, color = "black"),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 10),
    
    panel.grid = element_line(color = "black"),
    panel.grid.major = element_line(color = "gray75"),
    panel.grid.minor = element_line(color = "gray90"), panel.background = element_rect(fill = "white", color = "gray50"),
  )

multiplot(a, b, c, d, e, cols = 2)
# width 540, height 600




ggplot(resultsgregchhypo1, aes(x = `tgt-neighbor`, y = `qual-mean`, fill = log(`foraging-efficiency-dist`))) +
  geom_tile()+
  #scale_y_discrete(labels = y.scale) +
  scale_fill_gradient(low = "lightblue", high = "navy", guide = "colorbar")




ggplot(resultsgregchhypo1, aes(x = `tgt-dist`, y = `qual-mean`, fill = log(`foraging-efficiency-dist`))) +
  geom_tile()+
  #scale_y_discrete(labels = y.scale) +
  scale_fill_gradient(low = "lightblue", high = "navy", guide = "colorbar")

ggplot(resultsgregchhypo1, aes(x = `tgt-neighbor`, y = `qual-mean`, fill = `foraging-efficiency-dist`)) +
  geom_tile() +
  scale_fill_gradient(low = "lightblue", high = "navy", guide = "colorbar")

ggplot(resultsgregchhypo1, aes(x = `extraction-rate-mean`, y = `qual-mean`, fill = log(`foraging-efficiency-dist`))) +
  geom_tile() +
  scale_fill_gradient(low = "lightblue", high = "navy", guide = "colorbar")


ggplot(resultsgregchhypo1, aes(x = `extraction-rate-mean`, y = `abundance`, fill = `foraging-efficiency-time`)) +
  geom_tile()+
  #scale_y_discrete(labels = y.scale) +
  scale_fill_gradient(low = "lightblue", high = "navy", guide = "colorbar")


ggplot(resultsgregchhypo1, aes(x = `extraction-rate-mean`, y = `patch-regrowth-freq`, fill = log(`foraging-efficiency-dist`))) +
  geom_tile()+
  #scale_y_discrete(labels = y.scale) +
  scale_fill_gradient(low = "lightblue", high = "navy", guide = "colorbar")

ggplot(resultsgregchhypo1, aes(x = `qual-mean`, y = `patch-regrowth-freq`, fill = log(`foraging-efficiency-dist`))) +
  geom_tile()+
  #scale_y_discrete(labels = y.scale) +
  scale_fill_gradient(low = "lightblue", high = "navy", guide = "colorbar")


ggplot(resultsgregchhypo1, aes(x = `tgt-neighbor`, y = `extraction-rate-mean`, fill = `var-energy-intake-monthly`)) +
  geom_tile()+
  scale_y_discrete(labels = y.scale) +
  scale_fill_gradient(low = "lightblue", high = "navy", guide = "colorbar")

ggplot(resultsgregchhypo1, aes(x = `tgt-neighbor`, y = `extraction-rate-mean`, fill = `var-energy-intake-biweekly`)) +
  geom_tile()+
  scale_y_discrete(labels = y.scale) +
  scale_fill_gradient(low = "lightblue", high = "navy", guide = "colorbar")


a<-ggplot(resultsgregchhypo1, mapping = aes(x = `patch-regrowth-freq`,
                                            y = `foraging-efficiency-time`, 
                                            color = `extraction-rate-mean`)) + 
  geom_jitter() + 
  stat_summary(fun = "median", color = "red")

b<-ggplot(resultsgregchhypo1, mapping = aes(x = `patch-regrowth-freq`, 
                                            y = `mean-distance-traveled`, 
                                            color = `abundance`)) + 
  geom_jitter(size = 1) + 
  stat_summary(fun = "median", color = "red")

c<-ggplot(resultsgregchhypo1, mapping = aes(x = `patch-regrowth-freq`, 
                                            y = `foraging-efficiency-dist`, 
                                            color = `abundance`)) + 
  geom_jitter(size = 1) + 
  stat_summary(fun = "median", color = "red")

multiplot(a, b, c, cols = 1)




#########################
#I'm going to just do the h1p2 and h1p3 using an MEE

hypo1metrics <-  c("foraging-efficiency-time", 
                  "foraging-efficiency-dist", 
                  "mean-distance-traveled",
                  "var-energy-intake-biweek", #I'm throwing the H2 ones in here two
                  "var-energy-intake-monthly",
                  "mean [timemove] of primates / ticks",
                  "count primates")

hypo1p2constants <- list("group-recog-module?" = "false",
                       "go-tests-on?" = "false",
                       "resource-tests-on?" = "false",
                       "move-tests-on?" = "false",
                       "qual-sd" = 1,
                       "other-primate-detection-radius" = 50,
                       "resource-detection-radius" = 50,
                       "regrowth-rate" = 1.0,
                       "movement-noise" = 10) 
                    

hypo1p2vars <- list(
  "abundance" = list(min = 200000, max = 1200000, qfun = "qunif"),
  "clump-size" = list(min = 1, max = 1000, qfun = "qunif"), #new single param for clumps
  "energy-per-capita" = list(min = 4000, max = 5000, qfun = "qunif"), #modified to exaggerate effect of gregariousness
  "qual-mean" = list(min = 25, max = 150, qfun = "qunif"),
  "extraction-rate-mean" = list(min = 8, max = 11, qfun = "qunif"), #modified to exaggerate effect of gregariousness
  "extraction-rate-sd" = list(min = 1, max = 5, qfun = "qunif"),
  
  "tgt-dist" = list(min = 1, max = 40, qfun = "qunif"),
  "tgt-neighbor" = list(min = 0, max = 11, qfun = "qunif"),
  
  "patch-regrowth-freq" = list(min = 2500, max = 3000, qfun = "qunif"), #modified to exaggerate effect of gregariousness
  "max-move" = list(min = 25, max = 50, qfun = "qunif"))



nlhypo1p2 <- nl(nlversion = "6.2.2",
              nlpath = netlogopath2,
              modelpath = modelpath2,
              jvmmem = 12000)

nlhypo1p2@experiment <- experiment(expname = "gregchhypo1p2MEE",
                                 outpath = outpath2,
                                 repetition = 1,
                                 tickmetrics = "false",
                                 idsetup = "setup",
                                 idgo = "go",
                                 runtime = 4300,
                                 stopcond = "(ticks > 8000)",
                                 metrics = hypo1metrics,
                                 variables = hypo1p2vars,
                                 constants = hypo1p2constants)


nlhypo1p2@simdesign <- simdesign_morris(nl = nlhypo1p2,
                                      morristype = "oat",
                                      morrislevels = 12,
                                      morrisr = 72,
                                      morrisgridjump = 6, 
                                      nseeds = 1)


progressr::handlers("progress")
library(future)
plan(multisession)
resultsgregchhypo1p2 <- progressr::with_progress(run_nl_all(nlhypo1p2))

setsim(nlhypo1p2, "simoutput")<-resultsgregchhypo1p2
analysisgregchhypo1p2 <- analyze_nl(nlhypo1p2)


analysisgregchhypo1p2 <- split(analysisgregchhypo1p2, analysisgregchhypo1p2$metric)



musigma_processing <- function(data) {
  return(dcast(data, metric + parameter ~ index))
}


analysisgregchhypo1p2scatter <- data.frame()
for (df in analysisgregchhypo1p2) {
  analysisgregchhypo1p2scatter <- rbind(analysisgregchhypo1p2scatter, (musigma_processing(as.data.frame(df))))
}


library(ggrepel)


metric.levels <- unique(analysisgregchhypo1p2scatter$metric)
 ### this is probably just for supplement
pdf("gregChH1p2MEEmusigmaplots.pdf")
for (metric in metric.levels) {
  print(ggplot(analysisgregchhypo1p2scatter[analysisgregchhypo1p2scatter$metric == metric,], aes(x = mustar, y = sigma, label = parameter)) +
          geom_point() + geom_label_repel(size = 3) +
          labs(title = paste0("Morris Elementary Effects H1P2, Greg ", metric)))
}
dev.off()

figure3data <- analysisgregchhypo1p2scatter[analysisgregchhypo1p2scatter$metric == "foraging-efficiency-dist_mean",]
figure3data$number <- 1:10
figure3data$readable.params <- c("abundance", 
                                 "clump size", 
                                 "energy per capita",
                                 "extraction rate (mean)", 
                                 "extraction rate (SD)", 
                                 "movement speed", 
                                 "frequency of patch regrowth",
                                 "mean patch quality",
                                 "target distance",
                                 "target neighbors")

ggplot(figure3data, aes(x = mustar, y = sigma, label = parameter)) +
  geom_point() + geom_text_repel(aes(label=ifelse(mustar>0.4, as.character(readable.params),number)), size = 5) +
  labs(title = "MEE for foraging efficiency, Scenario B") +
  xlab(paste("\u03BC", "*"))+
  ylab("\u03C3 ") +
  theme(
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 14),
    
    panel.grid = element_line(color = "black"),
    panel.grid.major = element_line(color = "gray75"),
    panel.grid.minor = element_line(color = "gray90"), panel.background = element_rect(fill = "white", color = "gray50"),
  )


a2<-ggplot(resultsgregchhypo1p2, mapping = aes(x = `tgt-neighbor`, y = log(`foraging-efficiency-dist`))) + 
  geom_jitter(size = 1) + 
  stat_summary(fun = "median", color = "red")

b2 <- ggplot(resultsgregchhypo1p2, mapping = aes(x = `tgt-dist`, y = log(`foraging-efficiency-dist`))) + 
  geom_jitter(size = 1) + 
  stat_summary(fun = "median", color = "red") 

c2 <- ggplot(resultsgregchhypo1p2, mapping = aes(x = `qual-mean`, y = log(`foraging-efficiency-dist`))) + 
  geom_jitter(size = 1) + 
  stat_summary(fun = "median", color = "red") 

d2 <- ggplot(resultsgregchhypo1p2, mapping = aes(x = `tgt-dist`, y = `var-energy-intake-monthly`)) + 
  geom_jitter(size = 1) + 
  stat_summary(fun = "median", color = "red") 


#for violin plot
a2<-ggplot(resultsgregchhypo1p2, mapping = aes(factor(x = `tgt-neighbor`), y = `foraging-efficiency-dist`)) + 
  geom_violin() + 
  stat_summary(fun = "median", color = "red")  +
  scale_x_discrete(labels = x.scale)+
  xlab("target neighbor")+
  ylab("foraging efficiency distance") +
  theme(
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 10),
    
    panel.grid = element_line(color = "black"),
    panel.grid.major = element_line(color = "gray75"),
    panel.grid.minor = element_line(color = "gray90"), panel.background = element_rect(fill = "white", color = "gray50"),
  )


b2 <- ggplot(resultsgregchhypo1p2, mapping = aes(factor(x = `tgt-dist`), y = `foraging-efficiency-dist`)) + 
  geom_violin() + 
  stat_summary(fun = "median", color = "red") +
  xlab("target distance")+
  ylab("foraging efficiency distance") +
  theme(
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 10),
    
    panel.grid = element_line(color = "black"),
    panel.grid.major = element_line(color = "gray75"),
    panel.grid.minor = element_line(color = "gray90"), panel.background = element_rect(fill = "white", color = "gray50"),
  )

c2 <- ggplot(resultsgregchhypo1p2, mapping = aes(factor(x = `patch-regrowth-freq`), y = `foraging-efficiency-dist`)) + 
  geom_violin() + 
  stat_summary(fun = "median", color = "red") +
  xlab("patch regrowth frequency")+
  ylab("foraging efficiency distance") +
  theme(
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 10),
    
    panel.grid = element_line(color = "black"),
    panel.grid.major = element_line(color = "gray75"),
    panel.grid.minor = element_line(color = "gray90"), panel.background = element_rect(fill = "white", color = "gray50"),
  )


multiplot(a, b, c, cols = 1)

ggplot(resultsgregchhypo1p2, aes(x = `tgt-neighbor`, y = `tgt-dist`, fill = log(`foraging-efficiency-dist`))) +
  geom_tile()+
 # scale_y_discrete(labels = tgtdist.scale) +
  scale_fill_gradient(low = "lightblue", high = "navy", guide = "colorbar")+
  labs(fill = "efficiency")




#####################################
#PREDICTION 3#
#####################################

hypo1p3constants <- hypo1p2constants

hypo1p3vars <- list(
  "abundance" = list(min = 200000, max = 1200000, qfun = "qunif"),
  "clump-size" = list(min = 1, max = 1000, qfun = "qunif"), #new single param for clumps
  "energy-per-capita" = list(min = 8000, max = 9000, qfun = "qunif"), #modified to minimize effect of gregariousness
  "qual-mean" = list(min = 25, max = 150, qfun = "qunif"),
  "extraction-rate-mean" = list(min = 4, max = 7, qfun = "qunif"), #modified to minimize effect of gregariousness
  "extraction-rate-sd" = list(min = 1, max = 5, qfun = "qunif"),
  
  "tgt-dist" = list(min = 1, max = 40, qfun = "qunif"),
  "tgt-neighbor" = list(min = 0, max = 11, qfun = "qunif"),
  
  "patch-regrowth-freq" = list(min = 1000, max = 2000, qfun = "qunif"), #modified to minimize effect of gregariousness
  "max-move" = list(min = 25, max = 50, qfun = "qunif"))


nlhypo1p3 <- nl(nlversion = "6.2.2",
                nlpath = netlogopath2,
                modelpath = modelpath2,
                jvmmem = 12000)

nlhypo1p3@experiment <- experiment(expname = "gregchhypo1p3MEE",
                                   outpath = outpath2,
                                   repetition = 1,
                                   tickmetrics = "false",
                                   idsetup = "setup",
                                   idgo = "go",
                                   runtime = 4300,
                                   stopcond = "(ticks > 8000)",
                                   metrics = hypo1metrics,
                                   variables = hypo1p3vars,
                                   constants = hypo1p3constants)


nlhypo1p3@simdesign <- simdesign_morris(nl = nlhypo1p3,
                                        morristype = "oat",
                                        morrislevels = 12,
                                        morrisr = 72,
                                        morrisgridjump = 6, 
                                        nseeds = 1)


progressr::handlers("progress")
library(future)
plan(multisession)
resultsgregchhypo1p3 <- progressr::with_progress(run_nl_all(nlhypo1p3))

setsim(nlhypo1p3, "simoutput")<-resultsgregchhypo1p3
analysisgregchhypo1p3 <- analyze_nl(nlhypo1p3)

analysisgregchhypo1p3 <- split(analysisgregchhypo1p3, analysisgregchhypo1p3$metric)



musigma_processing <- function(data) {
  return(dcast(data, metric + parameter ~ index))
}
library(reshape2)

analysisgregchhypo1p3scatter <- data.frame()
for (df in analysisgregchhypo1p3) {
  analysisgregchhypo1p3scatter <- rbind(analysisgregchhypo1p3scatter, (musigma_processing(as.data.frame(df))))
}


library(ggrepel)


metric.levels <- unique(analysisgregchhypo1p3scatter$metric)
 
### only for supplement
pdf("gregChH1p3MEEmusigmaplots.pdf")
for (metric in metric.levels) {
  print(ggplot(analysisgregchhypo1p3scatter[analysisgregchhypo1p3scatter$metric == metric,], aes(x = mustar, y = sigma, label = parameter)) +
          geom_point() + geom_label_repel(size = 3) +
          labs(title = paste0("Morris Elementary Effects H1P3, Greg ", metric)))
}
dev.off()

figure5data <- analysisgregchhypo1p3scatter[analysisgregchhypo1p3scatter$metric == "foraging-efficiency-dist_mean",]
figure5data$number <- 1:10
figure5data$readable.params <- c("abundance", 
                                 "clump size", 
                                 "energy per capita",
                                 "extraction rate (mean)", 
                                 "extraction rate (SD)", 
                                 "movement speed", 
                                 "frequency of patch regrowth",
                                 "mean patch quality",
                                 "target distance",
                                 "target neighbors")

ggplot(figure5data, aes(x = mustar, y = sigma, label = parameter)) +
  geom_point() + geom_text_repel(aes(label=ifelse(mustar>1, as.character(readable.params),number)), size = 5) +
  labs(title = "MEE for foraging efficiency, Scenario C") +
  xlab(paste("\u03BC", "*"))+
  ylab("\u03C3 ") +
  theme(
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 14),
    
    panel.grid = element_line(color = "black"),
    panel.grid.major = element_line(color = "gray75"),
    panel.grid.minor = element_line(color = "gray90"), panel.background = element_rect(fill = "white", color = "gray50"),
  ) 

a3<-ggplot(resultsgregchhypo1p3, mapping = aes(x = `tgt-neighbor`, y = log(`foraging-efficiency-dist`))) + 
  geom_jitter(size = 1) + 
  stat_summary(fun = "median", color = "red")

b3<-ggplot(resultsgregchhypo1p3, mapping = aes(x = `tgt-dist`, y = log(`foraging-efficiency-dist`))) + 
  geom_jitter(size = 1) + 
  stat_summary(fun = "median", color = "red")

c<-ggplot(resultsgregchhypo1p3, mapping = aes(x = `clump-size`, y = `foraging-efficiency-dist`)) + 
  geom_jitter(size = 1) + 
  stat_summary(fun = "median", color = "red")


d3<-ggplot(resultsgregchhypo1p3, mapping = aes(x = `qual-mean`, y = `foraging-efficiency-dist`)) + 
  geom_jitter(size = 1) + 
  stat_summary(fun = "median", color = "red")

e3 <- ggplot(resultsgregchhypo1p3, mapping = aes(x = `tgt-dist`, y = `var-energy-intake-monthly`)) + 
  geom_jitter(size = 1) + 
  stat_summary(fun = "median", color = "red") 


#for violin plot
a3<-ggplot(resultsgregchhypo1p3, mapping = aes(factor(x = `tgt-neighbor`), y = `foraging-efficiency-dist`)) + 
  geom_violin() + 
  stat_summary(fun = "median", color = "red") +
  xlab("target neighbor")+
  ylab("foraging efficiency distance") +
  theme(
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 10),
    
    panel.grid = element_line(color = "black"),
    panel.grid.major = element_line(color = "gray75"),
    panel.grid.minor = element_line(color = "gray90"), panel.background = element_rect(fill = "white", color = "gray50"),
  )

b3<-ggplot(resultsgregchhypo1p3, mapping = aes(factor(x = `tgt-dist`), y = `foraging-efficiency-dist`)) + 
  geom_violin() + 
  stat_summary(fun = "median", color = "red") +
  xlab("target distance")+
  ylab("foraging efficiency distance") +
  theme(
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 10),
    
    panel.grid = element_line(color = "black"),
    panel.grid.major = element_line(color = "gray75"),
    panel.grid.minor = element_line(color = "gray90"), panel.background = element_rect(fill = "white", color = "gray50"),
  )


d3<-ggplot(resultsgregchhypo1p3, mapping = aes(factor(x = `qual-mean`), y = `foraging-efficiency-dist`)) + 
  geom_violin() + 
  stat_summary(fun = "median", color = "red") +
  xlab("qual mean")+
  ylab("foraging efficiency distance") +
  theme(
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 10),
    
    panel.grid = element_line(color = "black"),
    panel.grid.major = element_line(color = "gray75"),
    panel.grid.minor = element_line(color = "gray90"), panel.background = element_rect(fill = "white", color = "gray50"),
  )


multiplot(a, b, c, d, cols = 1)

a <- ggplot(resultsgregchhypo1p3, aes(x = `tgt-neighbor`, y = `tgt-dist`, fill = log(`foraging-efficiency-dist`))) +
  geom_tile()+
  scale_y_discrete(labels = tgtdist.scale) +
  scale_fill_gradient(low = "lightblue", high = "navy", guide = "colorbar") +
  labs(fill = "efficiency")
b <- ggplot(resultsgregchhypo1p3, aes(x = `clump-size`, y = `qual-mean`, fill = log(`foraging-efficiency-dist`))) +
  geom_tile()+
  #scale_y_discrete(labels = y.scale) +
  scale_fill_gradient(low = "lightblue", high = "navy", guide = "colorbar") +
  labs(fill = "efficiency")
multiplot(a, b, cols = 1)
dev.off()

#### Figure , both H1P2 and H1 P3###################
multiplot(a2, b2, c2, a3, b3, d3, cols = 2)




###########Figure 6 ###################
resultsgregchhypo1p2$`tgt-neighbor` <- as.factor(resultsgregchhypo1p2$`tgt-neighbor`)
resultsgregchhypo1p2$`abundance` <- as.factor(resultsgregchhypo1p2$`abundance`)
resultsgregchhypo1p2$`patch-regrowth-freq` <- as.factor(resultsgregchhypo1p2$`patch-regrowth-freq`)
resultsgregchhypo1p2$`tgt-dist` <- as.factor(resultsgregchhypo1p2$`tgt-dist`)
resultsgregchhypo1p2$`extraction-rate-mean` <- as.factor(resultsgregchhypo1p2$`extraction-rate-mean`)
resultsgregchhypo1p2$`clump-size` <- as.factor(resultsgregchhypo1p2$`clump-size`)
resultsgregchhypo1p2$`energy-per-capita` <- as.factor(resultsgregchhypo1p2$`energy-per-capita`)
resultsgregchhypo1p2$`qual-mean` <- as.factor(resultsgregchhypo1p2$`qual-mean`)

resultsgregchhypo1p3$`tgt-neighbor` <- as.factor(resultsgregchhypo1p3$`tgt-neighbor`)
resultsgregchhypo1p3$`abundance` <- as.factor(resultsgregchhypo1p3$`abundance`)
resultsgregchhypo1p3$`patch-regrowth-freq` <- as.factor(resultsgregchhypo1p3$`patch-regrowth-freq`)
resultsgregchhypo1p3$`tgt-dist` <- as.factor(resultsgregchhypo1p3$`tgt-dist`)
resultsgregchhypo1p3$`extraction-rate-mean` <- as.factor(resultsgregchhypo1p3$`extraction-rate-mean`)
resultsgregchhypo1p3$`clump-size` <- as.factor(resultsgregchhypo1p3$`clump-size`)
resultsgregchhypo1p3$`energy-per-capita` <- as.factor(resultsgregchhypo1p3$`energy-per-capita`)
resultsgregchhypo1p3$`qual-mean` <- as.factor(resultsgregchhypo1p3$`qual-mean`)

clump.scale <- mapply(function(x){as.character(signif(as.numeric(x), digits = 2))}, levels(as.factor(resultsgregchhypo1p3$`clump-size`)))


a <- ggplot(resultsgregchhypo1p2, aes(x = `tgt-neighbor`, y = `tgt-dist`, fill = log(`foraging-efficiency-dist`))) +
  geom_tile()+
  scale_y_discrete(labels = tgtdist.scale, breaks = everyother) +
  scale_x_discrete(breaks = tgtneighbor.breaks)+
  scale_fill_gradient(low = "lightblue", high = "navy", guide = "colorbar")+
  labs(title = "A.", fill = "log \nforaging\nefficiency") +
  xlab("target neighbors")+
  ylab("target distance") +  
  theme(
    axis.title = element_text(size = 10, color = "black"),
    axis.text = element_text(size = 8, color = "black"),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 10),
    
    panel.grid = element_line(color = "black"),
    panel.grid.major = element_line(color = "gray75"),
    panel.grid.minor = element_line(color = "gray90"), panel.background = element_rect(fill = "white", color = "gray50"),
  )

b <- ggplot(resultsgregchhypo1p3, aes(x = `tgt-neighbor`, y = `tgt-dist`, fill = log(`foraging-efficiency-dist`))) +
  geom_tile()+
  scale_y_discrete(labels = tgtdist.scale, breaks = everyother) +
  scale_x_discrete(breaks = tgtneighbor.breaks)+
  scale_fill_gradient(low = "lightblue", high = "navy", guide = "colorbar")+
  labs(title = "B.", fill = "log \nforaging\nefficiency") +
  xlab("target neighbors")+
  ylab("target distance") +  
  theme(
    axis.title = element_text(size = 10, color = "black"),
    axis.text = element_text(size = 8, color = "black"),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 10),
    
    panel.grid = element_line(color = "black"),
    panel.grid.major = element_line(color = "gray75"),
    panel.grid.minor = element_line(color = "gray90"), panel.background = element_rect(fill = "white", color = "gray50"),
  )
c <- ggplot(resultsgregchhypo1p3, aes(x = `clump-size`, y = `qual-mean`, fill = log(`foraging-efficiency-dist`))) +
  geom_tile()+
  scale_y_discrete(labels = qualmean.scale, breaks = everyother) +
  scale_x_discrete(labels = clump.scale, breaks = everyother, guide = guide_axis(n.dodge = 1, angle = 45)) +
  xlab("clump size") +
  ylab("patch quality")+
  scale_fill_gradient(low = "lightblue", high = "navy", guide = "colorbar") +
  labs(title = "C.", fill = "log\nforaging\nefficiency")+  
  theme(
    axis.title = element_text(size = 10, color = "black"),
    axis.text = element_text(size = 8, color = "black"),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 10),
    
    panel.grid = element_line(color = "black"),
    panel.grid.major = element_line(color = "gray75"),
    panel.grid.minor = element_line(color = "gray90"), panel.background = element_rect(fill = "white", color = "gray50"),
  )
multiplot(a, b, c, layout = matrix(c(1, 2, NA, 3), nrow = 2, byrow = TRUE))











### hypothesis 2 - variation over time
#### this does not have its own MEE anymore - the data was collected at the same time as H1
hist(resultsgregchhypo1$`var-energy-intake-monthly`)

##this is actually figure 7
figure5data <- analysisgregchhypo1scatter[analysisgregchhypo1scatter$metric == "var-energy-intake-monthly_mean",]
figure5data$number <- c(2, 4, 6, 8, 7, 1, 9, 5, 10, 3)
figure5data$readable.params <- c("abundance", #2
                                 "clump size", #4
                                 "energy per capita",#6
                                 "extraction rate (mean)", #8
                                 "extraction rate (SD)", #7
                                 "movement speed",#1 
                                 "frequency of patch regrowth",#9
                                 "mean patch quality", #5
                                 "target distance",#10
                                 "target neighbors") #3

ggplot(figure5data, aes(x = mustar, y = sigma, label = parameter)) +
  geom_point() + geom_text_repel(aes(label=ifelse(mustar>0.3, as.character(readable.params),number)), size = 5) +
  #labs(title = "Sensitivity analysis for variance in intake") +
  xlab(paste("\u03BC", "*"))+
  ylab("\u03C3 ") +
  theme(
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 14, color = "black"),
    
    panel.grid = element_line(color = "black"),
    panel.grid.major = element_line(color = "gray75"),
    panel.grid.minor = element_line(color = "gray90"), panel.background = element_rect(fill = "white", color = "gray50"),
  ) 


##################figure 7########################

a <- ggplot(resultsgregchhypo1, aes(x = `extraction-rate-mean`, y = `patch-regrowth-freq`, fill = (`var-energy-intake-monthly`))) +
  geom_tile()+
  scale_y_discrete(labels = regrowthfreq.scale, breaks = everyother) +
  scale_x_discrete(labels = extractionrate.scale, breaks = everyother, guide = guide_axis(n.dodge = 1, angle = 45)) +
  scale_fill_gradient(low = "lightblue", high = "navy", guide = "colorbar") +
  labs(title = "A.", fill = "variance\nin intake") +
  xlab("extraction rate (mean)")+
  ylab("patch regrowth interval") +
  theme(
    axis.title = element_text(size = 10, color = "black"),
    axis.text = element_text(size = 8, color = "black"),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 10),
    
    panel.grid = element_line(color = "black"),
    panel.grid.major = element_line(color = "gray75"),
    panel.grid.minor = element_line(color = "gray90"), panel.background = element_rect(fill = "white", color = "gray50"),
  )

b <- ggplot(resultsgregchhypo1, aes(x = `extraction-rate-mean`, y = `tgt-dist`, fill = (`var-energy-intake-monthly`))) +
  geom_tile()+
  scale_y_discrete(labels = tgtdist.scale, breaks = everyother) +
  scale_x_discrete(labels = extractionrate.scale, breaks = everyother, guide = guide_axis(n.dodge = 1, angle = 45)) +
  scale_fill_gradient(low = "lightblue", high = "navy", guide = "colorbar") +
  labs(title = "B.", fill = "variance\nin intake") +
  xlab("extraction rate (mean)")+
  ylab("target distance") +
  theme(
    axis.title = element_text(size = 10, color = "black"),
    axis.text = element_text(size = 8, color = "black"),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 10),
    
    panel.grid = element_line(color = "black"),
    panel.grid.major = element_line(color = "gray75"),
    panel.grid.minor = element_line(color = "gray90"), panel.background = element_rect(fill = "white", color = "gray50"),
  )
c <- ggplot(resultsgregchhypo1, aes(x = `patch-regrowth-freq`, y = `tgt-dist`, fill = (`var-energy-intake-monthly`))) +
  geom_tile()+
  scale_y_discrete(labels = tgtdist.scale, breaks = everyother) +
  scale_x_discrete(labels = regrowthfreq.scale, breaks = everyother, guide = guide_axis(n.dodge = 1, angle = 45)) +
  scale_fill_gradient(low = "lightblue", high = "navy", guide = "colorbar") +
  labs(title = "C.", fill = "variance\nin intake") +
  xlab("patch regrowth interval")+
  ylab("target distance") +
  theme(
    axis.title = element_text(size = 10, color = "black"),
    axis.text = element_text(size = 8, color = "black"),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 10),
    
    panel.grid = element_line(color = "black"),
    panel.grid.major = element_line(color = "gray75"),
    panel.grid.minor = element_line(color = "gray90"), panel.background = element_rect(fill = "white", color = "gray50"),
  )

multiplot(a, b, c, cols = 2)





figure9data <- analysisgregchhypo1p2scatter[analysisgregchhypo1p2scatter$metric == "var-energy-intake-monthly_mean",]
figure9data$number <- 1:10
figure9data$readable.params <- c("abundance", 
                                 "clump size", 
                                 "energy per capita",
                                 "extraction rate (mean)", 
                                 "extraction rate (SD)", 
                                 "movement speed", 
                                 "frequency of patch regrowth",
                                 "mean patch quality",
                                 "target distance",
                                 "target neighbors")

ggplot(figure9data, aes(x = mustar, y = sigma, label = parameter)) +
  geom_point() + geom_text_repel(aes(label=ifelse(mustar>0.2, as.character(readable.params),number)), size = 5) +
  labs(title = "MEE for variance in intake, Scenario B") +
  xlab(paste("\u03BC", "*"))+
  ylab("\u03C3 ") +
  theme(
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 14),
    
    panel.grid = element_line(color = "black"),
    panel.grid.major = element_line(color = "gray75"),
    panel.grid.minor = element_line(color = "gray90"), panel.background = element_rect(fill = "white", color = "gray50"),
  ) 





##################figure 8############################
harsh.regrowthfreq.scale <- mapply(function(x){as.character(signif(as.numeric(x), digits = 3))}, levels(as.factor(resultsgregchhypo1p2$`patch-regrowth-freq`)))


ggplot(resultsgregchhypo1p2, aes(x = `patch-regrowth-freq`, y = `tgt-dist`, fill = (`var-energy-intake-monthly`))) +
  geom_tile()+
  scale_y_discrete(labels = tgtdist.scale, breaks = everyother) +
  scale_x_discrete(labels = harsh.regrowthfreq.scale, breaks = everyother, guide = guide_axis(n.dodge = 1, angle = 45)) +
  scale_fill_gradient(low = "lightblue", high = "navy", guide = "colorbar") +
  labs(fill = "variance\nin intake")+
  xlab("patch regrowth interval")+
  ylab("target distance") +
  theme(
    axis.title = element_text(size = 10, color = "black"),
    axis.text = element_text(size = 8, color = "black"),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 10),
    
    panel.grid = element_line(color = "black"),
    panel.grid.major = element_line(color = "gray75"),
    panel.grid.minor = element_line(color = "gray90"), panel.background = element_rect(fill = "white", color = "gray50"),
  )
















figure11data <- analysisgregchhypo1p3scatter[analysisgregchhypo1p3scatter$metric == "var-energy-intake-monthly_mean",]
figure11data$number <- 1:10
figure11data$readable.params <- c("abundance", 
                                 "clump size", 
                                 "energy per capita",
                                 "extraction rate (mean)", 
                                 "extraction rate (SD)", 
                                 "movement speed", 
                                 "frequency of patch regrowth",
                                 "mean patch quality",
                                 "target distance",
                                 "target neighbors")

ggplot(figure11data, aes(x = mustar, y = sigma, label = parameter)) +
  geom_point() + geom_text_repel(aes(label=ifelse(mustar>0.02, as.character(readable.params),number)), size = 5) +
  labs(title = "MEE for variance in intake, Scenario C") +
  xlab(paste("\u03BC", "*"))+
  ylab("\u03C3 ") +
  theme(
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 14),
    
    panel.grid = element_line(color = "black"),
    panel.grid.major = element_line(color = "gray75"),
    panel.grid.minor = element_line(color = "gray90"), panel.background = element_rect(fill = "white", color = "gray50"),
  ) 




##### directly compare scenario B and C, in terms of both efficiency and variance?


############figure 9######################

resultsgregchhypo1$scenario <- rep("general", times = 792)
resultsgregchhypo1p2$scenario <- rep("fast-extracting, slow-renewing", times = 792)
resultsgregchhypo1p3$scenario <- rep("slow-extracting, fast-renewing", times = 792)

compare.scenarios.data <- data.frame()

compare.scenarios.data <- rbind(resultsgregchhypo1, resultsgregchhypo1p2, resultsgregchhypo1p3)




a <- ggplot(compare.scenarios.data, aes(x = factor(scenario, level =c("general", "fast-extracting, slow-renewing", "slow-extracting, fast-renewing")), 
                                        y = log(`foraging-efficiency-dist`))) +
  geom_violin(fill = "gray50") +
  labs(title = "A.")+
  xlab("resource scenario")+
  ylab("log foraging efficiency") +
  scale_x_discrete(breaks = c("general", 
                              "fast-extracting, slow-renewing", 
                              "slow-extracting, fast-renewing"), 
                   labels=(c("general", "fast-extract\nslow-renew", "slow-extract\nfast-renew")))+
  theme(
    axis.title = element_text(size = 10, color = "black"),
    axis.text = element_text(size = 7, color = "black"),
    
    panel.grid = element_line(color = "black"),
    panel.grid.major = element_line(color = "gray75"),
    panel.grid.minor = element_line(color = "gray90"), panel.background = element_rect(fill = "white", color = "gray50"),
  ) 




b <- ggplot(compare.scenarios.data, aes(x = factor(scenario, level = c("general", "fast-extracting, slow-renewing", "slow-extracting, fast-renewing")), 
                                        y = log(`var-energy-intake-monthly`))) +
  geom_violin(fill = "gray50")+
  labs(title = "B.")+
  xlab("resource scenario")+
  ylab("log variance in intake") +
  scale_x_discrete(breaks = c("general", 
                              "fast-extracting, slow-renewing", 
                              "slow-extracting, fast-renewing"), 
                   labels=(c("general", "fast-extract\nslow-renew", "slow-extract\nfast-renew")))+
  theme(
    axis.title = element_text(size = 10, color = "black"),
    axis.text = element_text(size = 7, color = "black"),
    
    panel.grid = element_line(color = "black"),
    panel.grid.major = element_line(color = "gray75"),
    panel.grid.minor = element_line(color = "gray90"), panel.background = element_rect(fill = "white", color = "gray50"),
  ) 

multiplot(a, b, cols = 2)















hypo2metrics <- c("var-energy-intake-biweek",
                  "var-energy-intake-quarterly") #something here is not working

hypo2constants <- list("starting-pop-primates" = 200,
                       "group-recog-module?" = "false",
                       "go-tests-on?" = "false",
                       "resource-tests-on?" = "false",
                       "move-tests-on?" = "false" ) 

hypo2vars <- list(
  "abundance" = list(min = 300000, max = 800000, qfun = "qunif"),
       "clump-size" = list(min = 1, max = 1000, qfun = "qunif"),
       "qual-mean" = list(min = 5, max = 20, qfun = "qunif"),
       "qual-sd" = list(min = 1, max = 20, qfun = "qunif"),
       "extraction-rate-mean" = list(min = 1, max = 3, qfun = "qunif"),
       "extraction-rate-sd" = list(min = 0, max = 3, qfun = "qunif"),
       "tgt-dist" = list(min = 0.5, max = 5.0, qfun = "qunif"),
       "tgt-neighbor" = list(min = 0, max = 9, qfun = "qunif"),
       "other-primate-detection-radius" = list(min = 5, max = 10, qfun = "qunif"),
       "resource-detection-radius" = list(min = 5, max = 10, qfun = "qunif"),
       "regrowth-rate" = list(min = 0.05, max = 1.0, qfun = "qunif"),
       "patch-regrowth-freq" = list(min = 1000, max = 3000, qfun = "qunif"),
       "movement-noise" = list(min = 10, max = 45, qfun = "qunif"))



nlhypo2 <- nl(nlversion = "6.2.2",
              nlpath = netlogopath2,
              modelpath = modelpath2,
              jvmmem = 8000)

nlhypo2@experiment <- experiment(expname = "gregchhypo2MEE",
                                 outpath = outpath2,
                                 repetition = 1,
                                 tickmetrics = "false",
                                 idsetup = "setup",
                                 idgo = "go",
                                 runtime = 10000,
                                 stopcond = "(ticks > 8000)",
                                 metrics = hypo2metrics,
                                 variables = hypo2vars,
                                 constants = hypo2constants)


nlhypo2@simdesign <- simdesign_morris(nl = nlhypo2,
                                      morristype = "oat",
                                      morrislevels = 8,
                                      morrisr = 32,
                                      morrisgridjump = 4, 
                                      nseeds = 1)

library(future)
plan(multisession)

progressr::handlers("progress")
resultsgregchhypo2 <- progressr::with_progress(run_nl_all(nlhypo2))

setsim(nlhypo2, "simoutput")<-resultsgregchhypo2

analysisgregchhypo2 <- analyze_nl(nlhypo2)

analysisgregchhypo2 <- split(analysisgregchhypo2, analysisgregchhypo2$metric)

musigma_processing <- function(data) {
  return(dcast(data, metric + parameter ~ index))
}

analysisgregchhypo2scatter <- data.frame()

for (df in analysisgregchhypo2) {
  analysisgregchhypo2scatter <- rbind(analysisgregchhypo2scatter, (musigma_processing(as.data.frame(df))))
}


library(ggrepel)


metric.levels <- unique(analysisgregchhypo2scatter$metric)

pdf("gregChH2MEEmusigmaplots-fixed.pdf")
for (metric in metric.levels) {
  print(ggplot(analysisgregchhypo2scatter[analysisgregchhypo2scatter$metric == metric,], aes(x = mustar, y = sigma, label = parameter)) +
          geom_point() + geom_label_repel(size = 3) +
          labs(title = paste0("Morris Elementary Effects, Greg ", metric)))
}
dev.off()
 

a<-ggplot(resultsgregchhypo2, mapping = aes(x = `extraction-rate-mean`, y = `var-energy-intake-biweek`)) + 
  geom_jitter(size = 1) + 
  stat_summary(fun = "median", color = "red")

b<-ggplot(resultsgregchhypo2, mapping = aes(x = `abundance`, y = `var-energy-intake-biweek`)) + 
  geom_jitter(size = 1) + 
  stat_summary(fun = "median", color = "red")

c<-ggplot(resultsgregchhypo2, mapping = aes(x = `tgt-neighbor`, y = `var-energy-intake-biweek`)) + 
  geom_jitter(size = 1) + 
  stat_summary(fun = "median", color = "red")

d<-ggplot(resultsgregchhypo2, mapping = aes(x = `tgt-dist`, y = `var-energy-intake-biweek`)) + 
  geom_jitter(size = 1) + 
  stat_summary(fun = "median", color = "red")

e<-ggplot(resultsgregchhypo2, mapping = aes(x = `qual-mean`, y = `var-energy-intake-biweek`)) + 
  geom_jitter(size = 1) + 
  stat_summary(fun = "median", color = "red")

f<-ggplot(resultsgregchhypo2, mapping = aes(x = `qual-sd`, y = `var-energy-intake-biweek`)) + 
  geom_jitter(size = 1) + 
  stat_summary(fun = "median", color = "red")

multiplot(a, b, c, d, e, f, cols = 1)



a<-ggplot(resultsgregchhypo2, mapping = aes(x = `extraction-rate-mean`, y = `var-energy-intake-quarterly`)) + 
  geom_jitter(size = 1) + 
  stat_summary(fun = "median", color = "red")

b<-ggplot(resultsgregchhypo2, mapping = aes(x = `abundance`, y = `var-energy-intake-quarterly`)) + 
  geom_jitter(size = 1) + 
  stat_summary(fun = "median", color = "red")

c<-ggplot(resultsgregchhypo2, mapping = aes(x = `tgt-neighbor`, y = `var-energy-intake-quarterly`)) + 
  geom_jitter(size = 1) + 
  stat_summary(fun = "median", color = "red")

d<-ggplot(resultsgregchhypo2, mapping = aes(x = `tgt-dist`, y = `var-energy-intake-quarterly`)) + 
  geom_jitter(size = 1) + 
  stat_summary(fun = "median", color = "red")

#for violin plot

a<-ggplot(resultsgregchhypo2, mapping = aes(factor(x = `extraction-rate-mean`), y = `var-energy-intake-quarterly`)) + 
  geom_violin() + 
  stat_summary(fun = "median", color = "red")

b<-ggplot(resultsgregchhypo2, mapping = aes(factor(x = `abundance`), y = `var-energy-intake-quarterly`)) + 
  geom_violin() + 
  stat_summary(fun = "median", color = "red")

c<-ggplot(resultsgregchhypo2, mapping = aes(factor(x = `tgt-neighbor`), y = `var-energy-intake-quarterly`)) + 
  geom_violin() + 
  stat_summary(fun = "median", color = "red")

d<-ggplot(resultsgregchhypo2, mapping = aes(factor(x = `tgt-dist`), y = `var-energy-intake-quarterly`)) + 
  geom_violin() + 
  stat_summary(fun = "median", color = "red")
multiplot(a, b, c, d, cols = 1)

###old figure 5!####
multiplot(c, d, a, b, cols = 1)

#####











#new heatmaps
############################new heatmaps

#Figure 4, foraging efficiency general scenario
#####
newheatmap4Adata <- output4A #tgt-neighbors and tgt-dist
newheatmap4Bdata <- output4B7A #extraction rate x patch regrowth int
newheatmap4Cdata <- output4C #extraction rate x patch quality
newheatmap4Ddata <- output4D #tgt-neighbors x extraction rate
newheatmap4Edata <- output4E #tgt-neighbors x patch regrwoth int
newheatmap4Fdata <- output4F7C #tgt-dist x patch regrowth int

fig4A <- ggplot(newheatmap4Adata, aes(x = tgt.neighbor, y = tgt.dist, 
                           fill = (foraging.efficiency.dist))) +
  geom_tile()+
  #scale_y_discrete(labels = tgtdist.scale, breaks = everyother) +
  #scale_x_discrete(labels = tgtneighbor.breaks, breaks = everyother, guide = guide_axis(n.dodge = 1, angle = 45)) +
  scale_fill_distiller(palette = "YlGn", direction = 1) +
  labs(fill = "foraging\nefficiency", title = "A.")+
  xlab("target neighbors")+
  ylab("target distance") +
  theme(
    axis.title = element_text(size = 10, color = "black"),
    axis.text = element_text(size = 8, color = "black"),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 10),
    panel.grid = element_line(color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), 
    panel.background = element_rect(fill = "white", color = "gray50")
  )



fig4B <- ggplot(newheatmap4Bdata, aes(x = extraction.rate.mean, y = patch.regrowth.interval, 
                             fill = (foraging.efficiency.dist))) +
  geom_tile()+
  #scale_y_discrete(labels = tgtdist.scale, breaks = everyother) +
  #scale_x_discrete(labels = tgtneighbor.breaks, breaks = everyother, guide = guide_axis(n.dodge = 1, angle = 45)) +
  scale_fill_distiller(palette = "YlGn", direction = 1) +
  labs(fill = "foraging\nefficiency", title = "B.")+
  xlab("extraction rate (mean)")+
  ylab("patch regrowth interval") +
  theme(
    axis.title = element_text(size = 10, color = "black"),
    axis.text = element_text(size = 8, color = "black"),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 10),
    panel.grid = element_line(color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), 
    panel.background = element_rect(fill = "white", color = "gray50")
  )

#extraction rate x patch quality
fig4C <- ggplot(newheatmap4Cdata, aes(x = extraction.rate.mean, y = qual.mean, 
                                      fill = (foraging.efficiency.dist))) +
  geom_tile()+
  #scale_y_discrete(labels = tgtdist.scale, breaks = everyother) +
  #scale_x_discrete(labels = tgtneighbor.breaks, breaks = everyother, guide = guide_axis(n.dodge = 1, angle = 45)) +
  scale_fill_distiller(palette = "YlGn", direction = 1) +
  labs(fill = "foraging\nefficiency", title = "C.")+
  xlab("extraction rate (mean)")+
  ylab("patch quality") +
  theme(
    axis.title = element_text(size = 10, color = "black"),
    axis.text = element_text(size = 8, color = "black"),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 10),
    panel.grid = element_line(color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), 
    panel.background = element_rect(fill = "white", color = "gray50")
  )

#tgt-neighbors x extraction rate
fig4D <- ggplot(newheatmap4Ddata, aes(x = tgt.neighbor, y = extraction.rate.mean, 
                                      fill = (foraging.efficiency.dist))) +
  geom_tile()+
  #scale_y_discrete(labels = tgtdist.scale, breaks = everyother) +
  #scale_x_discrete(labels = tgtneighbor.breaks, breaks = everyother, guide = guide_axis(n.dodge = 1, angle = 45)) +
  scale_fill_distiller(palette = "YlGn", direction = 1) +
  labs(fill = "foraging\nefficiency", title = "D.")+
  xlab("target neighbors")+
  ylab("extraction rate (mean)") +
  theme(
    axis.title = element_text(size = 10, color = "black"),
    axis.text = element_text(size = 8, color = "black"),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 10),
    panel.grid = element_line(color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), 
    panel.background = element_rect(fill = "white", color = "gray50")
  )

#tgt-neighbors x patch regrwoth int
fig4E <- ggplot(newheatmap4Edata, aes(x = tgt.neighbor, y = patch.regrowth.interval, 
                                      fill = (foraging.efficiency.dist))) +
  geom_tile()+
  #scale_y_discrete(labels = tgtdist.scale, breaks = everyother) +
  #scale_x_discrete(labels = tgtneighbor.breaks, breaks = everyother, guide = guide_axis(n.dodge = 1, angle = 45)) +
  scale_fill_distiller(palette = "YlGn", direction = 1) +
  labs(fill = "foraging\nefficiency", title = "E.")+
  xlab("target neighbors")+
  ylab("patch regrowth interval") +
  theme(
    axis.title = element_text(size = 10, color = "black"),
    axis.text = element_text(size = 8, color = "black"),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 10),
    panel.grid = element_line(color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), 
    panel.background = element_rect(fill = "white", color = "gray50")
  )

fig4F <- ggplot(newheatmap4Fdata, aes(x = tgt.dist, y = patch.regrowth.interval, 
                             fill = (foraging.efficiency.dist))) +
  geom_tile()+
  #scale_y_discrete(labels = tgtdist.scale, breaks = everyother) +
  #scale_x_discrete(labels = tgtneighbor.breaks, breaks = everyother, guide = guide_axis(n.dodge = 1, angle = 45)) +
  scale_fill_distiller(palette = "YlGn", direction = 1) +
  labs(fill = "foraging\nefficiency", title = "F.")+
  xlab("target distance")+
  ylab("patch regrowth interval") +
  theme(
    axis.title = element_text(size = 10, color = "black"),
    axis.text = element_text(size = 8, color = "black"),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 10),
    panel.grid = element_line(color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), 
    panel.background = element_rect(fill = "white", color = "gray50")
  )

multiplot(fig4A, fig4B, fig4C, fig4D, fig4E, fig4F, cols = 2)
#####

#Figure 5, foraging efficiency in alt scenarios
#####
newheatmap5Adata <- output5A #tgt-neighbors x tgt-dist
newheatmap5Bdata <- output5B #tgt-neighbors x tgt-dist
newheatmap5Cdata <- output5C #clump size x patch quality

fig5A <- ggplot(newheatmap5Adata, aes(x = tgt.neighbor, y = tgt.dist, 
                                      fill = foraging.efficiency.dist)) +
  geom_tile()+
  #scale_y_discrete(labels = tgtdist.scale, breaks = everyother) +
  #scale_x_discrete(labels = tgtneighbor.breaks, breaks = everyother, guide = guide_axis(n.dodge = 1, angle = 45)) +
  scale_fill_distiller(palette = "YlGn", direction = 1) +
  labs(fill = "foraging\nefficiency", title = "A.") + 
  xlab("target neighbors")+
  ylab("target distance") +
  theme(
    axis.title = element_text(size = 10, color = "black"),
    axis.text = element_text(size = 8, color = "black"),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 10),
    panel.grid = element_line(color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), 
    panel.background = element_rect(fill = "white", color = "gray50")
  )

fig5B <- ggplot(newheatmap5Bdata, aes(x = tgt.neighbor, y = tgt.dist, 
                                      fill = foraging.efficiency.dist)) +
  geom_tile()+
  #scale_y_discrete(labels = tgtdist.scale, breaks = everyother) +
  #scale_x_discrete(labels = tgtneighbor.breaks, breaks = everyother, guide = guide_axis(n.dodge = 1, angle = 45)) +
  scale_fill_distiller(palette = "YlGn", direction = 1) +
  labs(fill = "foraging\nefficiency", title = "B.") + 
  xlab("target neighbors")+
  ylab("target distance") +
  theme(
    axis.title = element_text(size = 10, color = "black"),
    axis.text = element_text(size = 8, color = "black"),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 10),
    panel.grid = element_line(color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), 
    panel.background = element_rect(fill = "white", color = "gray50")
  )

fig5C <- ggplot(newheatmap5Cdata, aes(x = clump.size, y = qual.mean, 
                                      fill = foraging.efficiency.dist)) +
  geom_tile()+
  #scale_y_discrete(labels = tgtdist.scale, breaks = everyother) +
  #scale_x_discrete(labels = tgtneighbor.breaks, breaks = everyother, guide = guide_axis(n.dodge = 1, angle = 45)) +
  scale_fill_distiller(palette = "YlGn", direction = 1) +
  labs(fill = "foraging\nefficiency", title = "C.") + 
  xlab("clump size")+
  ylab("patch quality") +
  theme(
    axis.title = element_text(size = 10, color = "black"),
    axis.text = element_text(size = 8, color = "black"),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 10),
    panel.grid = element_line(color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), 
    panel.background = element_rect(fill = "white", color = "gray50")
  )
multiplot(fig5A, fig5B, fig5C, cols = 2)

#####

#Figure 7: variance in intake, general scenario
#####
newheatmap7Adata <- output4B7A
newheatmap7Bdata <- output7B
newheatmap7Cdata <- output4F7C


fig7A <- ggplot(newheatmap7Adata, aes(x = extraction.rate.mean, y = patch.regrowth.interval, 
                           fill = var.energy.intake.monthly)) +
  geom_tile()+
  #scale_y_discrete(labels = tgtdist.scale, breaks = everyother) +
  #scale_x_discrete(labels = tgtneighbor.breaks, breaks = everyother, guide = guide_axis(n.dodge = 1, angle = 45)) +
  scale_fill_distiller(palette = "YlGn", direction = 1) +
  labs(fill = "var in\nintake", title = "A.") + 
  xlab("extraction rate (mean)")+
  ylab("patch regrowth interval") +
  theme(
    axis.title = element_text(size = 10, color = "black"),
    axis.text = element_text(size = 8, color = "black"),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 10),
    panel.grid = element_line(color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), 
    panel.background = element_rect(fill = "white", color = "gray50")
  )

fig7B <- ggplot(newheatmap7Bdata, aes(x = extraction.rate.mean, y = tgt.dist, 
                                      fill = var.energy.intake.monthly)) +
  geom_tile()+
  #scale_y_discrete(labels = tgtdist.scale, breaks = everyother) +
  #scale_x_discrete(labels = tgtneighbor.breaks, breaks = everyother, guide = guide_axis(n.dodge = 1, angle = 45)) +
  scale_fill_distiller(palette = "YlGn", direction = 1) +
  labs(fill = "var in\nintake", title = "B.") + 
  xlab("extraction rate (mean)")+
  ylab("target distance") +
  theme(
    axis.title = element_text(size = 10, color = "black"),
    axis.text = element_text(size = 8, color = "black"),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 10),
    panel.grid = element_line(color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), 
    panel.background = element_rect(fill = "white", color = "gray50")
  )

fig7C <- ggplot(newheatmap7Cdata, aes(x = tgt.dist, y = patch.regrowth.interval, 
                                      fill = var.energy.intake.monthly)) +
  geom_tile()+
  #scale_y_discrete(labels = tgtdist.scale, breaks = everyother) +
  #scale_x_discrete(labels = tgtneighbor.breaks, breaks = everyother, guide = guide_axis(n.dodge = 1, angle = 45)) +
  scale_fill_distiller(palette = "YlGn", direction = 1) +
  labs(fill = "var in\nintake", title = "C.") + 
  xlab("target distance")+
  ylab("patch regrowth interval") +
  theme(
    axis.title = element_text(size = 10, color = "black"),
    axis.text = element_text(size = 8, color = "black"),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 10),
    
    panel.grid = element_line(color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), 
    panel.background = element_rect(fill = "white", color = "gray50")
  )

multiplot(fig7A, fig7B, fig7C, cols = 2)

#####

#Figure 8: variance in intake, alt scenarios
######
newheatmap8data <- output8

ggplot(newheatmap8data, aes(x = tgt.dist, y = patch.regrowth.interval, 
                             fill = var.energy.intake.monthly)) +
  geom_tile()+
  #scale_y_discrete(labels = tgtdist.scale, breaks = everyother) +
  #scale_x_discrete(labels = tgtneighbor.breaks, breaks = everyother, guide = guide_axis(n.dodge = 1, angle = 45)) +
  scale_fill_distiller(palette = "YlGn", direction = 1) +
  labs(fill = "var in\nintake", title = "") + 
  xlab("target distance")+
  ylab("patch regrowth interval") +
  theme(
    axis.title = element_text(size = 10, color = "black"),
    axis.text = element_text(size = 8, color = "black"),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 10),
    panel.grid = element_line(color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), 
    panel.background = element_rect(fill = "white", color = "gray50")
  )
#####


#pattern-matching figures (supplements)
######


hist(resultsgregchhypo1$`mean [timemove] of primates / ticks`)
hist(resultsgregchhypo1$`mean-distance-traveled`* 10 / (4300 / 24))


days <- 4300 / 24
resultsgregchhypo1$distinkms <- (resultsgregchhypo1$`mean-distance-traveled`* 10) / days 
resultsgregchhypo1$distinkms <- resultsgregchhypo1$distinkms / 1000


     
#daily path length, comparison Vidal-Cordasco data and model output
a<- ggplot(VidalCardaso, aes(x = `DMD (km/day)`)) +
  geom_histogram(aes(y=..count../sum(..count..))) +
  xlim(0, 10) +
  ylim(0, 0.2) +
  labs(title = "A.") + 
  xlab("mean daily movment distance (km/day)\nVidal-Cordasco et al. 2020") +
  ylab("proportion of species")+
  theme(
    axis.title = element_text(size = 10, color = "black"),
    axis.text = element_text(size = 8, color = "black"),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 10),
    
    panel.grid = element_line(color = "black"),
    panel.grid.major = element_line(color = "gray75"),
    panel.grid.minor = element_line(color = "gray90"), 
    panel.background = element_rect(fill = "white", color = "gray50"),
  )



b<- ggplot(resultsgregchhypo1, aes(x = distinkms)) +
  geom_histogram(aes(y=..count../sum(..count..))) +
  xlim(0, 10) +
  ylim(0, 0.2) +
  labs(title = "B.") +
  xlab("mean daily movement distance\n(converted)") +
  ylab("proportion of simulations")+
  theme(
    axis.title = element_text(size = 10, color = "black"),
    axis.text = element_text(size = 8, color = "black"),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 10),
    
    panel.grid = element_line(color = "black"),
    panel.grid.major = element_line(color = "gray75"),
    panel.grid.minor = element_line(color = "gray90"), 
    panel.background = element_rect(fill = "white", color = "gray50"),
  )


multiplot(a, b, cols = 2)




#activity budget, comparison Kamilar & Cooper data and model output
trimmedKamilarCooper <- KamilarCooperactivitybudgetdata[complete.cases(KamilarCooperactivitybudgetdata[,37:38]),]
trimmedKamilarCooper$`move/move + feed` <- as.numeric(trimmedKamilarCooper$`move/move + feed`)

a<- ggplot(trimmedKamilarCooper, aes(x = `move/move + feed`)) +
  geom_histogram(aes(y=..count../sum(..count..))) +
  xlim(0, 1.0) +
  ylim(0, 0.2) +
  labs(title = "A.") + 
  xlab("time moving out of time moving or feeding\nKamilar & Cooper 2013") +
  ylab("proportion of species")+
  theme(
    axis.title = element_text(size = 10, color = "black"),
    axis.text = element_text(size = 8, color = "black"),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 10),
    
    panel.grid = element_line(color = "black"),
    panel.grid.major = element_line(color = "gray75"),
    panel.grid.minor = element_line(color = "gray90"), 
    panel.background = element_rect(fill = "white", color = "gray50"),
  )



b<- ggplot(resultsgregchhypo1, aes(x = `mean [timemove] of primates / ticks`)) +
  geom_histogram(aes(y=..count../sum(..count..))) +
  xlim(0, 1.0) +
  ylim(0, 0.2) +
  labs(title = "B.") +
  xlab("time moving of all timesteps") +
  ylab("proportion of simulations")+
  theme(
    axis.title = element_text(size = 10, color = "black"),
    axis.text = element_text(size = 8, color = "black"),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 10),
    
    panel.grid = element_line(color = "black"),
    panel.grid.major = element_line(color = "gray75"),
    panel.grid.minor = element_line(color = "gray90"), 
    panel.background = element_rect(fill = "white", color = "gray50"),
  )


multiplot(a, b, cols = 2)
#####