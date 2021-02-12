
##################################################################################
##                                                                              ##
##  Title  : Robustness and model selection in configurational causal modeling  ##
##           Auxiliary functions (to be sourced in the replication script)      ##
##                                                                              ##
##################################################################################
library(plyr)
library(dplyr)
library(gdata) 


# helper function for calculating completeness scores -- false solutions get 0
add_completeness <- function(a,b){

  cor <- a[[b]]
  
  if (is.na(cor$condition)){
    cor$completeness <- 0
    out <- cor
  } else {
  
  if (all(cor$correct)==TRUE){
    cor$completeness <- cor$complexity / getcomp(a$target)
    out <- cor
  } else{
    tr <- cor[cor$correct == TRUE,]
    fa <- cor[cor$correct == FALSE,]
    tr$completeness <- tr$complexity / getcomp(a$target)
    fa$completeness <- 0
    out <- rbind(tr, fa)
  }
  }  
    return(out)
}

# helper function for calculating average completeness scores
avg_comp <- function(scores, num){
  se <- lapply(scores, '[[',num)
  tops <- lapply(se, function(x) max(x$completeness))
  topsum <- sum(unlist(tops))
  out <- topsum / length(se)
  return(out)
}


# Introduce random noise to a data set. Use this to replicate 
# the results in 'RMSCCM replication.R'
# 'cleandata': data frame of clean data
#'allcombos': data frame of all logically possible configurations for factors in cleandata
#'no.replace' = how many rows to replace / add
#'add': if TRUE, add the noise rows to clean data, if FALSE, replace rows
#'bias': bias the noise toward particular configurations; passed to 'some' as 'prob' argument
bring_theNoise <- function(cleandata, allcombos, no.replace, add=FALSE, bias=NULL){
  if (add){
    
    a <- cleandata
    b <- some(setdiff(allcombos,cleandata), no.replace, replace = FALSE, prob=bias)
    Noisedata <- rbind(a, b)
    
  } else {
    
    a <- cleandata[sample(nrow(cleandata), nrow(cleandata)-no.replace, replace = FALSE),]
    b <- some(setdiff(allcombos,cleandata), no.replace, replace = FALSE, prob=bias)
    Noisedata <- rbind(a, b)
  }
  return(Noisedata)
}

# add or replace rows with noise -- new function for creating
# noisy data sets that allows biasing the samplinf of both clean and noise rows,
# or setting
# a fixed value of duplicates of randomly selected clean or noisy rows
bring_theNoiseV2 <- function(cleandata, 
                           allcombos, 
                           no.replace, 
                           add = FALSE, 
                           noisebias = NULL,
                           cleanbias = NULL,
                           rep.noise = NULL,
                           rep.clean = NULL){
  
  pnoise <- setdiff(allcombos, cleandata)
  
  if (add){
    a <- makedat(cleandata, bias = cleanbias, rep.rows = rep.clean)
  } else {
    a <- makedat(some(cleandata, nrow(cleandata)-no.replace), 
                 bias = cleanbias, rep.rows = rep.clean)
  }
  if (nrow(pnoise) < no.replace) {
     pnoise <- rbind(pnoise, some(pnoise, no.replace - nrow(pnoise)))
   }
  b <- makedat(some(pnoise, no.replace, replace = FALSE), bias = noisebias, rep.rows = rep.noise)
  

  return(rbind(a, b))
}

#bring_theNoiseV2 helper function
makedat <- function(x, bias = NULL, rep.rows = NULL){
  if (!is.null(bias)){
    return(x[sample(nrow(x), nrow(x), prob = bias, replace = TRUE),])
  }
  if (!is.null(rep.rows)){
    
    duprow <- rep(sample(1:nrow(x), 1), rep.rows)
    if (length(duprow) >= nrow(x)){
      out <- x[duprow,]
    } else {
      tempout <- setdiff(x, x[duprow,])
      out <- rbind(tempout[1:(nrow(tempout)-rep.rows+1),], x[duprow,])
    }
  } else {
    out <- x
  } 
  return(out)
}






# Perform a reanalysis series on data set x and and calculate FR-scores for the returned models.
# Reanalysis type is defined in argument 'type', a numeric vector of three values from the unit interval.
# Other arguments are passed to cnalow4 and subScore.
# Takes cna arguments as additional arguments.
# frscore calls cnalow4 with allconcov = TRUE such that all con/cov 
# combinations defined by 'type' are scanned.
# Calls subScore with normalize = TRUE to calculate normalized FR-scores.
# Returns a data frame of cna solutions and their FR-scores.
frscore <- function(x, type, allconcov=TRUE, output = "csf", normalize = TRUE, verbose = F,...){
  cl <- match.call()
  cl$type <- cl$normalize <- cl$verbose <- cl$scoretype <- NULL
  cl[[1]] <- as.name("cnalow4")
  attempt <- seq(type[1], type[2], type[3])
  cl$attempt <- attempt
  cl$output <- output
  cl$allconcov <- allconcov
  clres <- eval.parent(cl)
  rescomb <- rbind.fill(clres)
  rescomb <- subset(rescomb, is.na(rescomb[,1])==F)
  rescomb$condition <- gsub("\\),\\(", "\\)*\\(", as.character(rescomb$condition))
  out <- subScore(rescomb, normalize = normalize, verbose = verbose)
  return(out)
  }

# cnalow4 is a helper function used in frscore() to perform a re-analysis series on a data set.
# Can also be used separately to automate multiple analyses of a data set with varying con/cov thresholds.
# attempt determines values to be used as con/cov thresholds in the analyses.
# confirst = TRUE will start from con=cov=max(attempt) and lowers con by a number of steps 
# determined by lowstep, along the sequence determined by attempt,
# then lowers cov equal amount of steps, and repeats until con=cov=min(attempt).
# covfirst = TRUE does the same, but lowers cov first.  
# conmsc can be used to set a con.msc in the analyses to differ from con -- each cna run 
# is performed with con.msc = con - conmsc. 
# allconcov = TRUE will run cna with all unique combinations of con/cov values given in attempt -- i.e. performs
# a full reanalysis series that can be used for FR-robustness scoring.
# Returns a list of either asfs (atomic models) or csfs (complex models), 
# depending on argument 'output', which defaults to "csf", 
# and the con/cov/conmsc threshold values that produced each model. 
# attempt must be a declining sequence. 
# output must take either value "asf" or "csf". 
# Accepts cna arguments as additional arguments.

cnalow4 <- function(..., what = "c", 
                   attempt = seq(0.95, 0.75, -0.05), conmsc = 0, lowstep = 0, confirst=FALSE, 
                   covfirst=FALSE, allconcov=FALSE, ncsf = 20, output = "csf"){
  if (lowstep==0 & allconcov==FALSE){stop("lowstep must be greater than zero")}
  if (lowstep>length(attempt)){stop("lowstep must be smaller than length of attempt")}
  cl <- match.call()
  cl$attempt <- cl$conmsc <- cl$asf <- cl$ncsf <- cl$csf <- cl$max <- cl$lowstep <- cl$confirst <- cl$covfirst <- cl$allconcov <- cl$output <- NULL
  cl[[1]] <- as.name("cna")
  
  #create a data frame of values to be used as con/cov/con.msc
  if (allconcov==TRUE){ccargs <- as.data.frame(expand.grid(attempt, attempt))
    colnames(ccargs)<-c("lowfirst", "lowsec")
    ccargs$conmsc <- ccargs$lowfirst-conmsc} else {
  ccargs <- ccvalues2(attempt = attempt, lstep = lowstep)
  if (confirst==TRUE){ccargs$conmsc <- ccargs$lowfirst-conmsc}
  if (covfirst==TRUE){ccargs$conmsc <- ccargs$lowsec-conmsc}}
  
  #run cna with con/cov values from ccargs
  sols <- vector("list", length = nrow(ccargs))
  for (i in 1:length(sols)){
    if (confirst | allconcov){
      cl$con <- ccargs[i,"lowfirst"]
      cl$cov <- ccargs[i, "lowsec"]}
    if (covfirst){
      cl$con <- ccargs[i, "lowsec"]
      cl$cov <- ccargs[i, "lowfirst"]}
    cl$con.msc <- ccargs[i, "conmsc"]
    #if (output=="csf"){sols[[i]] <- csf1(eval.parent(cl), n = ncsf)} #for cna versions < 2
    if (output=="csf"){sols[[i]] <- csf(eval.parent(cl), n = ncsf)} #using csf() from cna package
    if (output=="asf"){sols[[i]] <- asf(eval.parent(cl))}
    dt <- data.frame(cnacon = rep(cl$con, nrow(sols[[i]])), 
                                  cnacov = rep(cl$cov, nrow(sols[[i]])), 
                                  conmsc = rep(cl$con.msc, nrow(sols[[i]])))
    sols[[i]] <- cbind(sols[[i]], dt)
  }
  
  return(sols)
}

#helper function for cnalow4 
ccvalues2 <- function(attempt, lstep){
  if (lstep>=length(attempt)){stop("lstep must be an integer smaller than length of attempt")}
  
  #create a sequence that starts descending first
  by <- attempt[1]-attempt[2]
  colo <- split(attempt, attempt)
  repind <- seq(length(colo), 1, -lstep)
  colo[repind[2:length(repind)]] <- lapply(colo[repind[2:length(repind)]], function(x) rep(x, lstep+1))
  mod <- as.integer(length(colo)%%lstep)
  if (mod==0){mod<-lstep}
  if (length(colo[[1]])==1){colo[[1]] <- rep(colo[[1]], mod)}
  colo <- unname(unlist(colo))
  colo <- rev(colo)
  cohi <- lohelp(colo = colo, by = by) #create the other sequence
  return(data.frame(lowfirst=colo, lowsec=cohi))
  }

#helper function for cnalow4
lohelp <- function(colo, by){
  out <- vector("numeric", length(colo))
  out[1] <- colo[1]
  for (i in 2:length(colo)){
    if (colo[i]<colo[i-1]){
      out[i]<-out[i-1]}
    else{out[i]<-out[i-1]-by}
    }
  out
  }




#helper function for subScore that tests for submodel-hood, for those cases when needed in subScore
subCounter <- function(s,p){
  if (is.submodel(s[,1],p[,1])){
    return(data.frame(mod=s[,1], subsc=p[,2], supmod=p[,1], supsc=s[,2], stringsAsFactors=FALSE))
  } else
    return(data.frame(mod=s[,1], subsc=0, supmod=p[,1], supsc=0, stringsAsFactors=FALSE))
  }

# Function that calculates FRscores for models extracted from cna results objects.
# Input must be a data frame of csfs or asfs outputted by cna().
# Returns the same data frame with added 'score' column.
# verbose = TRUE returns in addition a list displaying for each model the unique sub- supermodel types 
# of the model, and their contribution (i.e. number of model tokens) to the model's score
# scoretype -argument specifies how the score is calculated:
# "full" counts sub- and supermodels, "submodel" counts submodels only, and
# "supermodel" counts supermodels only.

subScore <- function(sols, weigh = FALSE, normalize = F, verbose = F, scoretype = "full"){
  if (NA %in% sols$condition){sols <- sols[!is.na(sols$condition),]}
  if (nrow(sols)==0){warning("nothing to test") 
    sols <- cbindX(sols, data.frame("score" = NA)); return(sols)}else if(nrow(sols)==1){
      sols$score <- 0
      if(verbose){
        scsums <- list(NULL)
        names(scsums) <- sols$condition
        return(list(sols, scsums))  
      }else{return(sols)}
    }else{
      
      out <- sols
      out$condition <- as.character(out$condition)
      out <- out[order(out$condition),]
      mf <- as.data.frame(table(sols$condition), stringsAsFactors = FALSE)
      
      if (nrow(mf)==1){
        if(scoretype %in% c("submodel", "supermodel")){
          out$score <- rep((mf$Freq-1)*2/2, mf$Freq)
        }else{
          out$score <- rep((mf$Freq-1)*2, mf$Freq)
        }
        if(weigh){out$score <- out$score*out$concov}
        if(normalize){if (max(out$score>=1)){out$score <- out$score / max(out$score)}}
        if(verbose){
          elems <- (mf$Freq-1)*2
          if (scoretype %in% c("supermodel", "submodel")){elems <- elems / 2}
          names(elems) <- out$condition[1]
          
          scsums <- list(elems)
          names(scsums) <- out$condition[1]
          return(list(out, scsums))
        }else{
          return(out)
        }
      }else{
        
        mf <- mf[order(mf[,1]),]
        sscore <- vector("list", nrow(mf))
        
        for (m in 1:nrow(mf)){
          subres <- vector("list", nrow(mf[-m,]))
          for(mo in 1:nrow(mf[-m,])){
            subres[[mo]] <- if (nchar(mf[,1][m]) > nchar(mf[-m,][,1][mo])){
              data.frame(mod=mf[,1][m], subsc=0, supmod=mf[-m,][,1][mo], supsc=0, stringsAsFactors=FALSE)
            }else{
              subCounter(mf[m,], mf[-m,][mo,])
            }  
          } 
          sscore[[m]] <- subres
        }
        
        sc <- rbind.fill(lapply(sscore, function(y) rbind.fill(y)))
        
        if(verbose){
          bs <- sc[, c(1,3,4)]
          colnames(bs)[colnames(bs) == "supsc"] <- "sub.frequency"
          colnames(bs)[colnames(bs) == "mod"] <- "model"
          bysup <- bs %>% group_split(supmod)
          supnames <- unlist(lapply(bysup, function(x) unique(x$supmod)))
          names(bysup) <- supnames
          subspermod <- lapply(bysup, function(x) x[,c(1,3)])
          subspermod <- lapply(subspermod, function(x) as.data.frame(x, stringsAsFactors = FALSE))
          
          sps <- sc[, c(1,2,3)]
          sps <- data.frame(supermodel = sps[,3], sup.frequency = sps[,2], mod = sps[,1])
          bysub <- sps %>% group_split(mod)
          subnames <- supnames <- unlist(lapply(bysub, function(x) unique(x$mod)))
          names(bysub) <- subnames
          superpermod <- lapply(bysub, function(x) x[,2])
          superpermod <- lapply(superpermod, function(x) as.data.frame(x, stringsAsFactors = FALSE))
          
          robbasis <- mapply(cbind, subspermod, superpermod, SIMPLIFY = F)
          mfs <- mf
          colnames(mfs)[colnames(mfs) == "Var1"] <- "model"
          dups <- lapply(names(robbasis), function(x) mfs[mfs[,1]==x,])
          dupscores <- lapply(dups, function(x) x %>% 
                                mutate(sub.frequency=Freq-1, sup.frequency = Freq-1, Freq = NULL))
          dupscores <- lapply(dupscores, function(x) if(x[,2] == 0){x[-1,]}else{x})
          robbasis <- mapply(rbind, robbasis, dupscores, SIMPLIFY = F)
          robred <- lapply(robbasis, function(x) x[x[,2] + x[,3] > 0,])
          
          if (scoretype == "full") {
            scsums <- lapply(robred, function(x) 
              if(nrow(x) == 0){x<-NULL}else{apply(x[,c(2,3)], 1, sum)})
          }
          
          if (scoretype == "supermodel") {
            scsums <- lapply(robred, function(x) 
              if(nrow(x) == 0){x<-NULL}else{x[,3]})
          }
          
          if (scoretype == "submodel") {
            scsums <- lapply(robred, function(x) 
              if(nrow(x) == 0){x<-NULL}else{x[,2]})
          }
          
          
          for (i in 1:length(scsums)){
            if(!is.null(scsums[[i]])){names(scsums[[i]]) <- robred[[i]][,1]}
          }
          
          scsums <- lapply(scsums, function(x) x[x>0])
          scsums <- lapply(scsums, function(x) if (length(x)<1){NULL}else{x})
        }
        
        
        
        pre.ssc <- sc[,c(1,2)] %>% group_by(mod) %>% mutate_at(vars("subsc"), sum) %>% distinct  
        pre.susc <- sc[,c(3,4)] %>% group_by(supmod) %>% mutate_at(vars("supsc"), sum) %>% distinct  
        pre.ssc <- pre.ssc[order(pre.ssc$mod),]
        pre.susc <- pre.susc[order(pre.susc$supmod),]
        
        if (scoretype == "full") {out$score <- rep(pre.ssc$subsc, mf$Freq) + rep(pre.susc$supsc, mf$Freq) +
          (rep(mf$Freq, mf$Freq)-1)*2}
        
        if (scoretype == "supermodel") {out$score <- rep(pre.ssc$subsc, mf$Freq)  +
          (rep(mf$Freq, mf$Freq)-1)*2/2}
        
        if (scoretype == "submodel") {out$score <- rep(pre.susc$supsc, mf$Freq) +
          (rep(mf$Freq, mf$Freq)-1)*2/2}
        
        if(weigh){out$score<-out$score*out$concov}
        if(normalize){if (max(out$score>=1)){out$score <- out$score / max(out$score)}}
        if(verbose==TRUE){return(list(out, scsums))}else{
          return(out)}}
    }
}





# calculates complexity of a solution -- slightly modified version of original by Mathias Ambuehl
getcomp <- function (cond){
    lhsides <- lapply(cna:::extract_asf(cond), cna:::lhs)
    ll <- lengths(strsplit(unlist(lhsides), "[\\+\\*]"))
    vapply(cna:::C_relist_Int(ll, lengths(lhsides)), sum, integer(1))
}

# minimizes and extracts csfs from cna solution objects -- slightly modified version of original by Mathias Ambuehl
# NOTE: only needed if running the scripts with cna versions < 2.0
csf1 <- function(x, n = 20, 
                 details = intersect(x$details,
                                      c("exhaustiveness", "faithfulness", "coherence"))){
 if(nrow(asf(x))==0){a <- csf(x,n=n)} else {
  a <- minimalizeCsf(x, n = n)
  a$n.asf <- a$redundantParts <- NULL
  nms <- names(a)
  names(a)[match("con", nms)] <- "consistency"
  names(a)[match("cov", nms)] <- "coverage"
  class(a) <- c("condTbl", "data.frame")
  class(a$outcome) <- "outcomeString"
  attributes(a$condition) <- list(class = "stdComplex")
  
  getComplexity <- function(cond){
    cond <- a$condition
    lhsides <- lapply(cna:::extract_asf(cond), cna:::lhs)
    ll <- lengths(strsplit(unlist(lhsides), "[\\+\\*]"))
    vapply(cna:::C_relist_Int(ll, lengths(lhsides)), sum, integer(1))
  }
  a$complexity <- getComplexity(a$condition)
  a <- a[with(a, order(complexity, -consistency * coverage)), , drop = FALSE]
  
  dd <- c("exhaustiveness", "faithfulness", "coherence")
  details <- dd[pmatch(details, dd, 0L)]

  rownames(a) <- NULL
                                     }
  a
}


### Modified versions of randomConds functions from the cna package,
### for creating random csf's with a common cause structure 

cc_randomCsf <- function(x, outcome = NULL, n.asf = NULL, compl = NULL, ccause = FALSE) 
{
  x <- full.ct(x)
  stopifnot(length(x) >= 4)
  type <- attr(x, "type")
  if (is.null(outcome)) {
    if (is.null(n.asf)) 
      n.asf <- 2:pmin(ncol(x) - 2, 4)
    if (length(n.asf) > 1) 
      n.asf <- sample(n.asf, 1)
    outcome <- sample(names(x), n.asf)
  }
  else {
    if (length(outcome) > ncol(x) - 2) 
      stop("The number of outcomes is limited to number of factors minus 2.")
    n.asf <- length(outcome)
  }
  if (is.null(compl)) 
    compl <- 2:pmin(ncol(x) - 1, 4)
  lhs_factors <- setdiff(names(x), outcome)
  outCsf <- ""
  for (i in seq_along(outcome)) {
    comcause <- NULL
    if (i > 1) {
      # if (ccause){
      #   comcause <- sample(lhs_factors, 1)
      #   }
      lhs_factors <- c(lhs_factors, outc)
    }
    outc <- outcome[[i]]
    xx <- x[c(lhs_factors, outc)]
    if (nzchar(outCsf)) {
      xx <- selectCases(outCsf, xx)
    }
    # comcause <- NULL
    repeat {
      if (i>1){
        if (ccause){
          precc <- gsub("\\)\\*\\(|\\(|\\)|<->.|[^A-Za-z]", "", outCsf)
          comcause <- sample(unlist(strsplit(precc, "")), 1)
          #comcause <- sample(unlist(strsplit(gsub("[^A-Za-z]", "", outCsf), "")), 1)
        }
      }
      rasf <- cc_randomAsf(xx, outcome = outc, compl = compl, 
                           how = "minimal", cc = comcause)
      if (!(cna:::lhs(rasf) %in% c("0", "1"))) 
        break
    }
    outCsf <- paste0(outCsf, if (nzchar(outCsf)) 
      "*", "(", rasf, ")")
  }
  if (ccause) {list(structure(outCsf, class = c("stdComplex", "character")), comcause)} else {
    structure(outCsf, class = c("stdComplex", "character")) 
  }
  
}


cc_randomAsf <- function (x, outcome = NULL, compl = NULL, how = c("inus", 
                                                                   "minimal"), cc = NULL) 
{
  how <- match.arg(how)
  x <- configTable(x, rm.dup.factors = FALSE, rm.const.factors = FALSE, 
                verbose = FALSE)
  if (how == "inus") 
    x <- full.ct(x)
  stopifnot(length(x) >= 3)
  if (!is.null(outcome)) 
    stopifnot(length(outcome) == 1, outcome %in% names(x))
  type <- attr(x, "type")
  if (is.null(compl)) {
    compl <- pmin(ncol(x) - 1, 4)
    if (compl > 2) 
      compl <- 2:compl
  }
  else {
    compl <- compl[compl <= length(x) - 1]
    compl <- compl[compl >= 1]
  }
  if (length(compl) < 1) 
    stop("Invalid specification of compl")
  if (length(compl) > 1) {
    n.msc <- sample(compl, 1)
    len.msc <- sample(compl, n.msc, replace = TRUE)
  }
  else {
    n.msc <- compl
    len.msc <- rep(compl, n.msc)
  }
  p <- ncol(x)
  tti <- cna:::ctInfo(x)
  if (type %in% c("cs", "fs")) {
    vals <- mapply(c, tti$resp_nms, tolower(tti$resp_nms), 
                   USE.NAMES = TRUE, SIMPLIFY = FALSE)
  }
  else {
    nvals <- lengths(tti$uniqueValues)
    vals <- split(tti$resp_nms, rep(seq_len(p), nvals))
    names(vals) <- names(nvals)
  }
  rhsFactor <- if (is.null(outcome)) {
    sample(names(vals), 1)
  }
  else {
    outcome
  }
  rhs <- if (type %in% c("cs", "fs")) {
    rhsFactor
  }
  else {
    sample(vals[[rhsFactor]], 1)
  }
  #vals1 <- vals[-match(rhsFactor, names(vals))]
  #vals1 <- vals[-match(c(rhsFactor, toupper(cc)), names(vals))]
  vals1 <- vals[-match(rhsFactor, names(vals))]
  repeat {
    if (!is.null(cc)){
      #prelhs <- lapply(len.msc-length(cc), cna:::msamp, vals1)
      prelhs <- lapply(len.msc, cna:::msamp, vals1)
      if(sample(c(TRUE, FALSE), 1))prelhs[[length(prelhs)+1]] <- cc
      if (!cc %in% unlist(prelhs)){prelhs[[length(prelhs)]][[length(prelhs[[length(prelhs)]])]] <- cc}
      # prelhs <- lapply(prelhs, function(x) {if (!cc %in% x){
      #   x[length(x)] <- cc
      # } else {x <- x};return(x)})
      
      
      lhs <- cna:::hconcat(list(prelhs), c("+", "*"))  
    } else {
      lhs <- cna:::hconcat(list(lapply(len.msc, cna:::msamp, vals1)), c("+", "*"))
    }
    
    
    lhs <- rreduce(lhs, x, full = how == "inus")
    #if (!lhs %in% c("0", "1") & (cc %in% unlist(strsplit(lhs, "")) | is.null(cc))) 
    if (!lhs %in% c("0", "1") & if (!is.null(cc)){cc %in% unlist(strsplit(lhs, ""))}else{TRUE}) 
      break
  }
  structure(paste0(lhs, "<->", rhs), class = c("stdAtomic", 
                                               "character"))
}


  