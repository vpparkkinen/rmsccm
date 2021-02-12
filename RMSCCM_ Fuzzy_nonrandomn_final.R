##  Title  : Robustness and model selection in configurational causal modeling

# Reproduces the benchmark results in Appendix section A2.2

#Preamble

# Comments after library() calls specify package versions used to run the script
library(cna) # V 3.0.1
library(QCA) # V 3.6
library(plyr) # V 1.8.5
library(dplyr) # V 0.8.5
library(gdata) # V 2.18.0
library(ggplot2) # V 3.3.0
library(reshape2) # V 0.8.8
library(gridExtra) # V 2.3
library(cowplot) # V 1.0.0
library(cnaOpt) # V 0.2.0 

# setwd("yourWD") #replace yourWD with correct path to funcs_revised.R
source("funcs_revised.R")
 

suppressWarnings(RNGversion("3.6.1"))






n.sets <- 1000 # number of datasets
tcomp <- c(1,2) # target complexity, 1 or 2 outcomes
samplesize <- c(1,3) # sample size multipliers

# noise ratios (= percentage of cases incompatible with the ground truth)
# noise <- c(0.05, 0.15, 0.25) 
noise <- 0.1 
fuzzvals <- seq(0, 0.4, 0.01)

# fit thresholds:
# first element is the lower bound for con/cov used in the re-analysis series from which  
# FRscore and MaxFit models are selected 
# second and third element are con and cov thresholds for Conv0.75 and Conv0.8, respectively 
fitfloor <- c(0.7, 0.75, 0.8) 

# prepare test settings 
setlist <- list(tcomp, samplesize, noise)
settings <- expand.grid(setlist) 
colnames(settings) <- c("tcomp", "ssize", "noise") 

# lists for storing results
results <- vector("list", nrow(settings)) # aggregate results stored here
inspect <- vector("list", nrow(settings)) # selected models by each approach stored here

set.seed(21) # seed for selecting seeds, for replicability
seeds <- sample(.Machine$integer.max, nrow(settings))

nfacs <- 6 # no. of factors in the data sets
comple <- 2:3 # max complexity of targets

tarops <- list(c("E"), c("D", "E")) # target outcome options

tquant <- 0.98 # percentile above which scored models are selected

# loop over test parameter combinations 

for(s in 1:nrow(settings)){
  
cat("test", s, "of", nrow(settings),"--", "target outcomes:", settings[s,1], " ",
    "sample size multiplier:", settings[s,2] , " ",
    "noise ratio:", settings[s,3], "\n")

set.seed(seeds[s])


# select outcome for the target DGSs
outcome <- tarops[[settings[s,1]]]

# create the ground truth DGSs that are the search targets
target <- replicate(n.sets, randomCsf(nfacs, outcome = outcome, compl = comple))

# create ideal data sets that match the targets
alldata <- allCombs(c(rep(2,nfacs)))-1 #all logically possible configurations
testdata <- lapply(lapply(target, function(x) selectCases(x, alldata)), ct2df) # select rows that match the targets

# noise_rows <- setdiff(alldata, testdata[[1]])
# 
# 
# weights <- rnorm(nrow(noise_rows), 10, 3)
# bias_data <- 

#rowno_alldata <- seq(1, nrow(alldata))

# adjust sample size 
testdatadup <- testdata 
t <- 1
while(t<settings[s,2]){
  testdata <- mapply(rbind, testdata, testdatadup, SIMPLIFY = F)
  t<-t+1
}

# set noise ratio
nratio <- settings[s,3]


# bias the sampling of noise rows
#noiseprobs <- rep(c(1, 40), length.out = round(nrow(testdata[[1]])*nratio))

noiserepeats <- round(nrow(testdata[[1]]) * nratio)

Noisedata <- vector("list", length(testdata))
for (i in 1:length(testdata)){
 Noisedata[[i]] <- bring_theNoiseV2(testdata[[i]], alldata,
                                     no.replace = round(nrow(testdata[[i]])*nratio),
                                     cleanbias = NULL,
                                     noisebias = NULL,
                                     rep.noise = noiserepeats,
                                     rep.clean = NULL,
                                     add = F)
 }

Noisedata <- lapply(Noisedata, function(x) makeFuzzy(x, fuzzvalues = fuzzvals))

# ----analyses----
 
# quote calls for performing the analyses (in order to save them with results later)

# complete reanalysis series
call_reanalysis <- quote(frscore(x = Noisedata[[i]], type = c(1, fitfloor[1], -0.1), inus.only = T, rm.dup.factors = F, rm.const.factors=F,
                         strict = T, ordering = as.list(outcome),
                         allconcov = T, normalize =T))

# conventional fit threshold option 1 (con=cov=0.75)
call_convention1 <- quote(csf1(cna(Noisedata[[i]], inus.only = T,  rm.dup.factors = F,
                      rm.const.factors = F, con = fitfloor[2], cov = fitfloor[2], ordering = as.list(outcome),
                      strict = T)))


# convention, fit option 2 (con=cov=0.8)
call_convention2 <- quote(csf1(cna(Noisedata[[i]], inus.only = T, rm.dup.factors = F,
                      rm.const.factors = F, con = fitfloor[3], cov = fitfloor[3], ordering = as.list(outcome),
                      strict = T)))


# perform a reanalysis series and calculate robustness scores
llrobu <- vector("list", length(Noisedata))
for (i in 1:length(Noisedata)){
  llrobu[[i]] <- eval(call_reanalysis)  
}


# analyze data sets with the two conventional fit threshold options

llconv <- vector("list", length(Noisedata))
for (i in 1:length(Noisedata)){
  llconv[[i]] <- eval(call_convention1)
}
  
llconv2 <- vector("list", length(Noisedata))
for (i in 1:length(Noisedata)){
  llconv2[[i]] <- eval(call_convention2)
}


ress <- list(llrobu, llconv, llconv2)

# replace comma with star in non-cohering models to allow is.submodel
# to be called on those -- does not affect causal interpretation of models 
# (note, frscore does this for the results of the re-analysis series)
ress[2:3] <- lapply(ress[2:3], function(y) lapply(y, function(x) {x$condition <- gsub("\\),\\(", "\\)*\\(", 
                                                         as.character(x$condition));return(x)}))


# calculate total number of models returned for each approach across all data sets
total_no_sols <- lapply(ress, function(x) sum(unlist(lapply(x, nrow))))
avgsols_per_data <- lapply(total_no_sols, function(x) x / n.sets) #average per data set


# add con*cov column to results
ress <- lapply(ress, function(x) lapply(x, function(y){y$concov <- y[,3]*y[,4];return(y)}))

# check correctness
ress <- lapply(ress, function(x){
  for (i in 1:length(x)){
    if (nrow(x[[i]])==0 || is.na(x[[i]]$condition)) {x[[i]] <- cbindX(x[[i]], data.frame(correct = NA))  
    }else 
      {x[[i]]$correct <- is.submodel(x[[i]]$condition, target[[i]])}
  };return(x)})

# for all data sets, check if a correct model is found in _all_ the results returned for each data set
correct_in_raw <- lapply(ress, function(y) 
  lapply(y, function(x)if (is.na(x$correct)){FALSE}else{any(x$correct)}))

# for all data sets, check if correct found or no solutions at all returned; i.e are the results fallacy free 
fallacy_free_raw <- lapply(ress, function(y) 
  lapply(y, function(x) any(x$correct) || is.na(x$correct) || nrow(x)==0))

fallacy_free_ratio_raw <- lapply(fallacy_free_raw, function(x) sum(unlist(x), na.rm = T) / length(Noisedata))

# prepare results from re-analysis series for maxfit selection
pre_MaxFit <- ress[[1]]

# select maxfit top con*cov models from the results of the re-analysis series
# maxfit selects models with con*cov above 98 percentile, from all models produced by the highest
# fit threshold that produced non-null results

# first select only models that were produced with the highest fit threshold
# setting that produced a result
pre_MaxFit <- lapply(pre_MaxFit, function(x) x %>% mutate(thresh=cnacon*cnacov))
maxFit_selected <- lapply(pre_MaxFit, function(y) if (!is.null(y)){
            y[y$thresh == max(y$thresh),]})

# from the above, select models with con*cov above the 98 percentile 
maxFit_selected <- lapply(maxFit_selected, function(x) if (!is.null(x)){
  x[x$concov >= quantile(x$concov, tquant, na.rm = T),]})

# drop columns not needed anymore and select unique models
maxFit_selected <- lapply(maxFit_selected, function(x)
  {if (is.na(x$condition)) {x<-x} else{
  x <- x[,c('condition','consistency','coverage','complexity', 'concov',
            'correct')] %>% distinct(condition, .keep_all = T)};return(x)})

# select models with robustness score in the top 98 percentile from results of the re-analysis series
FRscore_selected <- lapply(ress[[1]], function(y) if (!is.null(y)){
            y[y$score >= quantile(y$score, tquant, na.rm = T),]})

# drop columns not needed anymore and select unique models
FRscore_selected <- lapply(FRscore_selected, function(x)
  {if (is.na(x$condition)) {x<-x} else{
  x <- x[,c('condition','consistency','coverage','complexity', 'concov',
            'correct', 'score')] %>% distinct(condition, .keep_all = T)};return(x)})

# correctness / fallacy free statistics for MaxFit
topthrbb <- rbind.fill(maxFit_selected)
avg_sols_MaxFit <- nrow(topthrbb[!is.na(topthrbb$condition),]) / n.sets
cor_in_MaxFit <- lapply(maxFit_selected, function(x) any(x$correct))
avg_cor_MaxFit <- sum(unlist(cor_in_MaxFit), na.rm = T) / n.sets
ffree_in_MaxFit <- lapply(maxFit_selected, function(x) is.na(x$condition) || any(x$correct))
avg_ffree_MaxFit <- sum(unlist(ffree_in_MaxFit), na.rm = T) / n.sets


# correctness / fallacy free statistics for the two conventional settings

convs_selected <- lapply(ress[2:length(ress)], function(x) lapply(x, function(y)
  {if(is.na(y$condition)){y <- y}else{y <- y[y$concov >= quantile(y$concov, tquant, na.rm = T),]};return(y)}))

convs_selected_noNA <- lapply(convs_selected, function(y) lapply(y, function(x) if(!is.na(x$condition)){x}))
avg_sols_convs <- lapply(convs_selected_noNA, function(x) nrow(rbind.fill(x)) / n.sets)
avg_sols_convs <- lapply(avg_sols_convs, function(x){
       if(identical(x, numeric(0))){x<-0};return(x)}) 

cor_in_convs <- lapply(convs_selected, function(x) lapply(x, function(y) any(y$correct)))
sum_cor_in_convs <- lapply(cor_in_convs, function(x) sum(unlist(x), na.rm = T))
avg_cor_convs <- lapply(sum_cor_in_convs, function(x) x/length(Noisedata))
ffree_in_convs <- lapply(convs_selected, function(x) lapply(x, function(y) 
  if(is.na(y$condition)){TRUE}else{any(y$correct)}))
sum_ffree_in_convs <- lapply(ffree_in_convs, function(x) sum(unlist(x), na.rm = T))
avg_ffree_convs <- lapply(sum_ffree_in_convs, function(x) x/length(Noisedata))


# correctness / fallacy free statistics for robustness scoring
cor_in_FRscore <- lapply(FRscore_selected, function(y) if (is.null(y)){FALSE}else{any(y$correct)})
avg_cor_FRscore <- sum(unlist(cor_in_FRscore), na.rm = TRUE) / length(cor_in_FRscore)
ffree_in_FRscore <- lapply(FRscore_selected, function(y) is.null(y) || any(y$correct) || is.na(y$condition))
avg_ffree_FRscore <- sum(unlist(ffree_in_FRscore), na.rm = TRUE) / length(ffree_in_FRscore)

avg_sols_FRscore <- nrow(rbind.fill(FRscore_selected)[!is.na(rbind.fill(FRscore_selected)$condition),]) / n.sets

# store selected solutions and test results for all approaches
allscores <- mapply(list, ress[[1]], ress[[2]], ress[[3]], FRscore_selected,
                           convs_selected[[1]], convs_selected[[2]], maxFit_selected,
                            target, 
                    c(call_reanalysis),  c(nratio), SIMPLIFY = F)

allscores <- lapply(allscores, function(x){ names(x) <- c("solutions whole series",
                                                          "solutions conventional setting1","solutions conventional setting2", 
                                                          "top sols robustness", "top sols conv1",
                                                          "top sols conv2",
                                                          'top sols MaxFit', "target", 
                                                               "reanalysis type", "noise ratio"); return(x)})


# add completeness ratios to selected solutions
for (n in 4:7){
allscores <- lapply(allscores, function(x) {x[[n]] <- add_completeness(x,n);return(x)})  
}

# save selected results plus additional info for each approach
inspect[[s]] <- allscores 

# average completeness
avg_completeness <- vector("numeric", 4)
j <- 1
for (i in 4:7){
  avg_completeness[j] <- avg_comp(allscores, i)
  j<-j+1
}

# ratio of correct model found per data set in _all_ models returned in the re-analysis series'  
ratiocor_raw_reanalysis <- sum(unlist(correct_in_raw[[1]]), na.rm = TRUE) / length(correct_in_raw[[1]])

# same for conventional settings
ratiocor_raw_conv <- lapply(correct_in_raw[-1], function(x)
                            sum(unlist(x), na.rm = TRUE) / length(x))

# store average benchmark results for each test type
# order : robustness, convention1, convention2,

# total number of returned models
soltotals <- c(total_no_sols[[1]], total_no_sols[[2]], total_no_sols[[3]], 
               total_no_sols[[1]])

# average number of returned models
sols_avgs <- c(avgsols_per_data[[1]], avgsols_per_data[[2]], avgsols_per_data[[3]],
               avgsols_per_data[[1]])

# average number of models in the set of selected models for each approach
solstop_avgs <- c(avg_sols_FRscore, 
                  avg_sols_convs[[1]], avg_sols_convs[[2]],
                    avg_sols_MaxFit)

# raw correctness ratios 
raw_cor_ratios <- c(ratiocor_raw_reanalysis, 
                    ratiocor_raw_conv[[1]], ratiocor_raw_conv[[2]],
                    ratiocor_raw_reanalysis)
# raw fallacy free ratios
raw_ffree_ratios <- c(fallacy_free_ratio_raw[[1]], fallacy_free_ratio_raw[[2]], 
                      fallacy_free_ratio_raw[[3]], fallacy_free_ratio_raw[[1]])

# correctness ratios for selected models
sel_cor_ratios <- c(avg_cor_FRscore, 
                    avg_cor_convs[[1]], avg_cor_convs[[2]],
                    avg_cor_MaxFit)

# fallacy free ratios for selected models
sel_ffree_ratios <- c(avg_ffree_FRscore, avg_ffree_convs[[1]], 
                      avg_ffree_convs[[2]], 
                       avg_ffree_MaxFit)


# store fit threshold floor values
fitfloors <- c(fitfloor, fitfloor[1])

#store results
results[[s]] <- data.frame("selection"=c("robustness", "convention 0.75",
                                         'convention 0.8', 
                                        'max fit threshold'),
           "dsets" = n.sets, "noi ratio" = nratio, 
           "N" = nrow(Noisedata[[1]]), 
           "t.asfs" = settings[s,1], "fit floor"=fitfloors,
           "total sols"=soltotals, "avg sols" = sols_avgs, 
           "corr.all sols" = raw_cor_ratios, "ffree all sols" = raw_ffree_ratios, 
           "avg sel sols" = solstop_avgs, 
           "corr.selected" = sel_cor_ratios, 
           "ffree selected"=sel_ffree_ratios,
           "completeness" = avg_completeness)

}


  # summarize results across all tests

results_all <- rbind.fill(results)




results_all_means <- results_all %>% group_by(selection) %>% 
  summarise_at(c("total.sols", "avg.sols", "corr.all.sols", "ffree.all.sols",
                 "avg.sel.sols", "corr.selected", "ffree.selected", "completeness"), .funs = mean)

results_all_means

# plots

toplot <- results_all

# plot1 -- noise levels
# tp_noise_cor <- toplot[,c(1,3,12,13,14)]
# 
# tp_noise_cor_mean <- tp_noise_cor %>% group_by(selection, noi.ratio) %>% 
#  summarise(corr.selected=mean(corr.selected), ffree = mean(ffree.selected), 
#            completeness = mean(completeness))
# 
# k1 <- melt(tp_noise_cor_mean, id.vars = c("noi.ratio",  "selection")) 
# k1$pos <- rep(tp_noise_cor_mean$ffree, 3)
# k1$lab <- rep(c(rep("conv 0.75", 3), rep("conv 0.8", 3),  rep("maxfit", 3), rep("robust", 3)), 3)
# k1$noi.ratio <- as.character(k1$noi.ratio)
# 
# plot1 <- ggplot(k1, aes(x = noi.ratio, y = value, group = selection, alpha = variable)) +
#   geom_bar(stat="identity", position = position_dodge(width = .8), width=.7) +
#   geom_text(aes(x = noi.ratio, y = pos, label = lab, vjust = -3), 
#             position = position_dodge(0.7), size = 3.1, show.legend = F) +
#   scale_fill_grey(start = 0, end = 0.6) +
#   theme_bw()+
#   theme(plot.title = element_text(size = 6))+
#   scale_x_discrete(name ="noise ratio")+ylim(0, 1.1) + 
#   ylab("correctness / fallacy freeness") +
#   scale_alpha_manual(values=c(0.6, .06, .9)) +
#   ggtitle("results with different noise levels") +
#   theme(plot.title = element_text(size = 20, face="bold"))  
# 
# plot1
# 

# plot2 -- target asfs

tp_asfs <- toplot[,c(1,5,12,13,14)]

tp_asfs_mean <- tp_asfs %>% group_by(selection, t.asfs) %>% 
 summarise(corr.selected=mean(corr.selected), ffree = mean(ffree.selected), 
           completeness = mean(completeness))

k2 <- melt(tp_asfs_mean, id.vars = c("t.asfs",  "selection"))
k2$pos <- rep(tp_asfs_mean$ffree, 3)
k2$lab <- rep(c(rep("conv 0.75", 2), rep("conv 0.8", 2),  rep("maxfit", 2), rep("robust", 2)), 3)
k2$t.asfs <- as.character(k2$t.asfs)
plot2 <- ggplot(k2, aes(x = t.asfs, y = value, group = selection, alpha = variable)) +
  geom_bar(stat="identity", position = position_dodge(width = .8), width=.7) +
  geom_text(aes(x = t.asfs, y = pos, label = lab, vjust = -3), position = position_dodge(0.7), show.legend = F) +
  scale_fill_grey(start = 0, end = 0.6) +
  theme_bw()+
  theme(plot.title = element_text(size = 6))+
  scale_x_discrete(name ="target asfs")+ylim(0, 1.1) + 
  ylab("correctness / fallacy freeness") +
  scale_alpha_manual(values=c(0.6, .06, .9)) +
  ggtitle("results with different target complexities") +
  theme(plot.title = element_text(size = 20, face="bold"))

plot2

# plot3 -- sample size multiplier  
# NOTE: actual N per target asf is displayed in a small table in the plot --
# where this is rendered will depend on the size of the plot window, change
# draw_grob arguments accordingly

tp_N <- toplot[,c(1,4,12,13,14)]

tp_N_mean <- tp_N %>% group_by(selection, N) %>% 
 summarise(corr.selected=mean(corr.selected), ffree = mean(ffree.selected), 
           completeness = mean(completeness))

tp_N_mean$sp <- rep(c(rep(samplesize[1], 2), rep(samplesize[2], 2)), nrow(results_all_means))

low <- unique(tp_N_mean[tp_N_mean$sp==samplesize[1],]$N) 
high <- unique(tp_N_mean[tp_N_mean$sp==samplesize[2],]$N)

ns <- data.frame("target" = c("1 asf", "2 asf"), "small N"=c(low[2], low[1]), "large N"=c(high[2], high[1]))

ngrob <- draw_grob(tableGrob(ns, rows = NULL), x=0.224, y=0.8, width=0.01, height=0.01)


tp_N_mean <- tp_N_mean %>% group_by(selection, sp) %>% 
   summarise(corr.selected=mean(corr.selected), ffree = mean(ffree), 
           completeness = mean(completeness))


k3 <- melt(tp_N_mean, id.vars = c("sp",  "selection"))
k3$pos <- rep(tp_N_mean$ffree, 3)
k3$lab <- rep(c(rep("conv 0.75", 2), rep("conv 0.8", 2),  rep("maxfit", 2), rep("robust", 2)),3)
k3$sp <- as.character(k3$sp)
nplot <- ggplot(k3, aes(x = sp, y = value, group = selection, alpha = variable)) +
  geom_bar(stat="identity", position = position_dodge(width = .8), 
           width=.7) +
  geom_text(aes(x = sp, y = pos, label = lab, vjust = -3.), position = position_dodge(0.8),
            show.legend = F) +

  scale_fill_grey(start = 0, end = 0.6) +
  theme_bw()+
  theme(plot.title = element_text(size = 6))+
  scale_x_discrete(name ="sample size multiplier")+ylim(0, 1.2) + 
  ylab("correctness / fallacy freeness") +
  scale_alpha_manual(values=c(0.6, .06, .9)) +
  ggtitle("results with different sample sizes") +
  theme(plot.title = element_text(size = 20, face="bold"))

plot3 <- ggdraw(nplot) + ngrob 
plot3



# plot 4 -- average results

selsols <- results_all_means[,c(1,6)]
toplot_all <- results_all_means


selsols$selection <- c("conv 0.75", "conv 0.8", "maxfit", "robust")
colnames(selsols)[colnames(selsols)=="avg.sel.sols"] <- "avg sols"
selsols$`avg sols` <- signif(selsols$`avg sols`, 3)
tbl <- draw_grob(tableGrob(selsols, rows=NULL, theme = ttheme_default(base_size = 7.7)), 
                 x=0.75, y=0.73 , width=.01, height=.01)

tp_agg_res <- toplot_all[,c(1,7,8,9)]

k4 <- melt(tp_agg_res, id.vars = "selection")
k4$lab <- rep(c("conv 0.75", "conv 0.8",  "maxfit", "robust"), 3)

k4 <- k4[c(which(k4$variable=="ffree.selected"),which(k4$variable=="corr.selected"),which(k4$variable=="completeness")) ,]
k4$variable <- as.character(k4$variable)
k4$variable  <- factor(k4$variable , levels=unique(k4$variable ))

plotsel <- ggplot(k4, aes(x = variable, y = value, fill = selection)) +
  geom_bar(stat="identity", position = position_dodge(width = .8), width=.7) +
  scale_fill_manual(values = c("grey87", "grey68", "grey49", "grey20")) +
  theme_bw()+
  theme(plot.title = element_text(size = 6))+
  scale_x_discrete("performance measure")+ylim(0, 1) + 
  ggtitle("Average results with each selection method") +
  theme(plot.title = element_text(size = 20, face="bold"))
  
plot4 <- ggdraw(plotsel) + tbl
plot4



