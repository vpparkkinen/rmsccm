##################################################################################
##                                                                              ##
##  Title  : Robustness and model selection in configurational causal modeling  ##
##    Authors: Veli-Pekka Parkkinen & Michael Baumgartner                       ##
##                                                                              ##
##                                                                              ##
##################################################################################


# Code for replicating the examples discussed in sections Overfitting and Robustness,
# and the benchmarking results discussed in section Benchmarking. 


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

setwd("yourWD") #replace yourWD with correct path to 'funcs_revised.R'
source("funcs_revised.R")


suppressWarnings(RNGversion("3.6.1"))


#------Section Overfitting-----------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------

#-Analysis of Table 1(a)----------------------------------------------------------------------------
tab1a <- structure(list(A = c(0, 0, 1, 1, 1, 0, 0, 1, 0, 0, 1, 0, 1, 0, 
1, 1), B = c(0, 0, 1, 1, 1, 1, 0, 0, 1, 1, 1, 0, 0, 1, 0, 0), 
    C = c(0, 1, 1, 1, 0, 0, 1, 0, 1, 1, 0, 0, 0, 0, 1, 1), D = c(1, 
    1, 1, 0, 1, 1, 0, 1, 0, 1, 0, 0, 0, 0, 1, 0), E = c(0, 0, 
    1, 1, 1, 0, 0, 1, 1, 1, 1, 0, 1, 0, 1, 0)), row.names = c("c1", 
"c2", "c3", "c4", "c5", "c6", "c7", "c8", "c9", "c10", "c11", "c12", 
"c13", "c14", "c15", "c16"), class = "data.frame")
tab1a
# CNA analysis at maximal fit thresholds
cna(tab1a)
# QCA analysis at maximal fit thresholds
tt <- truthTable(tab1a, outcome="E")
minimize(tt, details=T, include="?")

#-Table 1(b)----------------------------------------------------------------------------------------

# reanalyze data at systematically varied con/cov thresholds
tab1bmods <- cnalow4(tab1a, ordering = list("E"), attempt = seq(1, 0.75, -0.05), allconcov = TRUE)

# bind the results into single df and select columns matching table 1(b)
tab1bmods <- do.call(rbind, tab1bmods) %>% 
  select(condition, consistency, coverage, cnacon, cnacov) 

# filter out duplicates of solutions, select only those shown in Table 1(b)
tab1b <- tab1bmods %>% group_by(condition) %>%  filter(cnacon*cnacov==max(cnacon*cnacov))

# add correctness column
tab1b$correct <- is.submodel(tab1b$condition, "A+B*C<->E")

tab1b


#-Reproduction of the results in Figure 1. ----------------------------------------------------------

set.seed(3) # for replicability
# generate the configurations compatible with the target (3) (p. 9)
fullData <- ct2df(full.ct(5))
target <- "A + B*C <-> E"
dat <- ct2df(selectCases(target, fullData))
# generate the configurations incompatible with the target (3)
dat.neg <- setdiff(fullData,dat)

size <- 1:5 # increase sample size, from factor 1 to factor 5
remove <- 1:3 # incompatibility shares
n <- 200 # number of iterations of each trial

# Increase sample size loop
size.score <- vector("list", length(size))
for(j in 1:length(size)){
    cat(j, "\n")
    dat.test <- do.call(rbind, rep(list(dat),size[j]))
    
    # Increase number of incompatible cases loop 
    remove.score <- vector("list", length(remove))
    for(k in 1:length(remove)){
          z <-remove[k]*size[j]
          
          # Test loop
          score <- vector("list", n)
          for(i in 1:n){
               dat.lim <- dat.test[-sample(1:nrow(dat.test), z), ] 
               dat.real <- rbind(dat.neg[sample(1:nrow(dat.neg), z), ], dat.lim)
               x <- cnaOpt(dat.real, "E",rm.dup.factors = F,rm.const.factors = F)$condition
               y <- any(is.submodel(x, target))
               score[[i]] <- list(x,y)
          } # end test loop
           
     remove.score[[k]] <- score
    
    } # end remove loop
    
    size.score[[j]] <- remove.score
} # end size loop


result <- lapply(size.score, function(x) lapply(x, function(y) (n-sum(sapply(y, "[[", 2)))/n))

result.size <- lapply(result, function(x) do.call(rbind, x))
overall <- as.data.frame(do.call(cbind, result.size))
colnames(overall) <- c("size1","size2","size3","size4","size5")
overall$id <- row.names(overall)

overall



# PLOT
ofK1 <- melt(overall, id.vars = "id")
colnames(ofK1) <- c("factor", "size", "value")
of_plot <- ggplot(ofK1, aes(x = size, y = value,group = 1, fill = factor)) + 
  geom_bar(stat = "identity", position = position_dodge2(), width=.7) + scale_fill_grey(start = 0, end = .8)+
  theme_bw()+
  theme(plot.title = element_text(size = 6))+
  scale_x_discrete(name ="")+ylim(0, 1)

of_plot


#------Section Robustness-----------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------

# Table 2
# The following code reproduces contents of Table 2. Creates an object 'tb2' - a list of two.
# The first element of the list tb2 displays the models, their counts, and raw- and normalized
# FR-scores, as in Table 2.
# The second element is a list, where each named element contains counts of sub- and supermodels
# a model that appears as a name of that list element

tb2 <- frscore(x = tab1a, type = c(1, 0.75, -.05), allconcov = T, verbose = T, ordering = list("E"),
               strict = T, normalize = F)

tb2$score_raw <- tb2$score
temp_tb2 <- frscore(x = tab1a, type = c(1, 0.75, -.05), allconcov = T, verbose = F, ordering = list("E"),
               strict = T, normalize = T)

tb2[[1]]$score_norm <- temp_tb2$score
tb2[[1]]$condition <- as.character(tb2[[1]]$condition)
tb2[[1]] <- tb2[[1]] %>% group_by(condition) %>% mutate(q = n())
tb2[[1]] <- tb2[[1]][, c(2,11,9,10)] %>% distinct(condition, .keep_all = T) 
colnames(tb2[[1]])[colnames(tb2[[1]])=="score"] <- "score_raw"

tb2 #contents of Table 2

# footnote 11: a complete listing of model tokens 
# and the fit threshold settings (columns 'cnacon', 'cnacov', and 'conmsc') 
# at which they were returned

footnote11 <- rbind.fill(cnalow4(tab1a, attempt = seq(1, 0.75, -.05), allconcov = T, ordering = list("E"),
                 strict=T))

footnote11 # a complete listing of model tokens 
           # and the fit threshold settings (columns 'cnacon') at which they were returned



#------Section Benchmarking---------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------

# Reproduces the benchmark results in section Benchmarking and the plots in figures 2. and 3.
# Benchmark tests for three approaches to model selection:
# 1./2. Selection by best model fit from results obtained with conventional fit threshold settings 
# con=cov=0.75 and con=cov=0.8 (Conv0.75, Conv0.8),
# 3. scanning model fit thresholds for largest values that produces non-null results and selecting by best model fit (MaxFit), and 
# 4. selecting by FR-score (FRscore).

# preparations for test series

# NOTE: running the full test battery with 1000 data sets will take some time. To get a 
# result faster, reduce the number assigned to n.sets.


n.sets <- 1000 # number of datasets
tcomp <- c(1,2) # target complexity, 1 or 2 outcomes
samplesize <- c(1,3) # sample size multipliers

# noise ratios (= percentage of cases incompatible with the ground truth)
noise <- c(0.05, 0.15, 0.25) 

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


# adjust sample size 
testdatadup <- testdata 
t <- 1
while(t<settings[s,2]){
  testdata <- mapply(rbind, testdata, testdatadup, SIMPLIFY = F)
  t<-t+1
}

# add noise
nratio <- settings[s,3]

Noisedata <- vector("list", length(testdata))
for (i in 1:length(testdata)){
 Noisedata[[i]] <- bring_theNoise(testdata[[i]], alldata,
                                     no.replace = round(nrow(testdata[[i]])*nratio), 
                                     add = F)
 }

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
tp_noise_cor <- toplot[,c(1,3,12,13,14)]

tp_noise_cor_mean <- tp_noise_cor %>% group_by(selection, noi.ratio) %>% 
 summarise(corr.selected=mean(corr.selected), ffree = mean(ffree.selected), 
           completeness = mean(completeness))

k1 <- melt(tp_noise_cor_mean, id.vars = c("noi.ratio",  "selection")) 
k1$pos <- rep(tp_noise_cor_mean$ffree, 3)
k1$lab <- rep(c(rep("conv 0.75", 3), rep("conv 0.8", 3),  rep("maxfit", 3), rep("robust", 3)), 3)
k1$noi.ratio <- as.character(k1$noi.ratio)

plot1 <- ggplot(k1, aes(x = noi.ratio, y = value, group = selection, alpha = variable)) +
  geom_bar(stat="identity", position = position_dodge(width = .8), width=.7) +
  geom_text(aes(x = noi.ratio, y = pos, label = lab, vjust = -3), 
            position = position_dodge(0.7), size = 3.1, show.legend = F) +
  scale_fill_grey(start = 0, end = 0.6) +
  theme_bw()+
  theme(plot.title = element_text(size = 6), legend.title = element_blank())+
  scale_x_discrete(name ="noise ratio")+ylim(0, 1.1) + 
  ylab(element_blank()) +

  scale_alpha_manual(values=c(0.6, .06, .9), label = c("correctness", "fallacy-freeness", "completeness")) +
  ggtitle("results with different noise levels") +
  theme(plot.title = element_text(size = 20, face="bold"))  

plot1


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
  theme(plot.title = element_text(size = 6), legend.title = element_blank())+
  scale_x_discrete(name ="target asfs")+ylim(0, 1.1) + ylab(element_blank()) +
  scale_alpha_manual(values=c(0.6, .06, .9), label = c("correctness", "fallacy-freeness", "completeness")) +
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

ns <- data.frame("target" = c("1 out", "2 out"), "small"=c(low[2], low[1]), "large"=c(high[2], high[1]))


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
  theme(plot.title = element_text(size = 6), legend.title = element_blank())+
  scale_x_discrete(name ="sample size multiplier")+ylim(0, 1.2) + ylab(element_blank()) +
  scale_alpha_manual(values=c(0.6, .06, .9), label = c("correctness", "fallacy-freeness", "completeness")) +
  ggtitle("results with different sample sizes") +
  theme(plot.title = element_text(size = 20, face="bold"))

plot3 <- ggdraw(nplot) + ngrob 
plot3

# plot 4 -- average results

selsols <- results_all_means[,c(1,6)]
toplot_all <- results_all_means

selsols <- rename(selsols, approach=selection)
selsols$approach <- c("conv 0.75", "conv 0.8", "maxfit", "robust")
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
  theme(plot.title = element_text(size = 6), legend.title = element_blank())+
  scale_x_discrete(labels = c("fallacy-freeness", "correctness", "completeness"))+
  xlab(element_blank()) +
  ylab(element_blank()) +
  ylim(0, 1) + 
  ggtitle("Average results with each approach") +
  theme(plot.title = element_text(size = 20, face="bold"))

plot4 <- ggdraw(plotsel) + tbl
plot4



