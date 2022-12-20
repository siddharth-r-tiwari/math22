library(ggplot2)
library(Metrics)
library(markovchain)
library(tidyverse)
library(patchwork)


get_pl2 <- function(){
  ratings_pl2 <- read.csv("placeboII_ratings.csv", strip.white = TRUE, stringsAsFactors = FALSE)
  
  pain_pl2<-c()
  pain_pl2 <- (sapply(unique(ratings_pl2$subject), function(sub) {
    tmp = subset(ratings_pl2, subject == sub)
    tmp <- na.omit(tmp)
    
    if(sum(as.numeric(tmp$phase == 1)) > 0){
      cbind(subject = sub,
            time = tmp$day,
            pain = tmp$pain) 
    }
  }))
  pain_pl2 <- do.call(rbind, pain_pl2)
  pain_pl2 <- na.omit(pain_pl2)
  pain_pl2 <- as.data.frame(pain_pl2)
  colnames(pain_pl2) <- c("subject","time", "pain")
  
  return(pain_pl2)
}

get_ldopa <- function(){
  ratings_ldopa <- read.csv("ldopa_ratings.csv", strip.white = TRUE, stringsAsFactors = FALSE)
  ratings_ldopa <- na.omit(ratings_ldopa)
  
  pain_ldopa <- (lapply(unique(ratings_ldopa$id), function(sub) {
    tmp = subset(ratings_ldopa, sub == id)
    tmp <- na.omit(tmp)
    
    start <- as.numeric(as.POSIXct(tmp[1,4], format="%d-%b-%Y"))
    
    cbind(
      subject = sub,
      time = (as.numeric(as.POSIXct(tmp[,4], format="%d-%b-%Y")) - start)/86400,
      pain = tmp$pain
    ) 
    
    
  }))
  
  pain_ldopa <- do.call(rbind, pain_ldopa)
  pain_ldopa <- na.omit(pain_ldopa)
  pain_ldopa <- as.data.frame(pain_ldopa)
  colnames(pain_ldopa) <- c("subject", "time", "pain") 
  
  return(pain_ldopa)
}

softmax <- function(x) {x / sum(x)}

pain_pl2 <- get_pl2()
pain_ldopa <- get_ldopa()

sub = "SAT009"

#extract sequence of pain ratings from subject
sequence <- round(as.numeric(subset(pain_ldopa, subject == sub)$pain), digits = 0)

#fit markov model
mcFitMLE <- markovchainFit(data = sequence)

transition_matrix <- mcFitMLE[["estimate"]]@transitionMatrix

#round to two digits for demonstration
round(transition_matrix, digits = 2)

#obtain steady state vector
#left eigenspace decomposition (requires transpose)
#extract real value and scale to components to sum to 1 due to floating errors
eigenvector <- softmax(abs(Re(eigen(t(transition_matrix))[["vectors"]][,1])))

round(eigenvector, digits=2)

end_behavior <- (as.numeric(rownames(transition_matrix)) %*% eigenvector)

end_behavior

metric <- sequence[1] - end_behavior

metric