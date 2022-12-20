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

change_in_pain <- function(pain, n){
  
  #iterate through subjects to return matrix of scores
  scores <- sapply(unique(pain$subject), function(sub) {
    
    #extract sequence of pain ratings from subject
    full_sequence <- round(as.numeric(subset(pain, subject == sub)$pain), digits = 0)
    
    #calculate len of n% of sequence
    len <- round(n*length(full_sequence), digits = 0)
    
    #n% of sequence
    sequence <- full_sequence[1:len]
    
    #fit markov model
    transition_matrix <- markovchainFit(data = sequence)[["estimate"]]@transitionMatrix
    
    #obtain steady state vector
    #left eigenspace decomposition (requires transpose)
    #extract real value and scale to components to sum to 1 due to floating errors
    eigenvector <- softmax(abs(Re(eigen(t(transition_matrix))[["vectors"]][,1])))
    
    #calculate metrics
    data.frame(subject = sub,
               #post - pre
               metric_tot = full_sequence[length(full_sequence)] - full_sequence[1],
               #end = weighted average of all possible pain values
               metric_exp = (as.numeric(colnames(transition_matrix)) %*% eigenvector) - sequence[1],
               metric_std = sequence[length(sequence)] - sequence[1],
               
               #absolute difference between the metrics
               diff_exp = abs(full_sequence[length(full_sequence)] - full_sequence[1] - (as.numeric(colnames(transition_matrix)) %*% eigenvector) - sequence[1]),
               diff_std = abs(full_sequence[length(full_sequence)] - full_sequence[1] - (sequence[length(sequence)] - sequence[1])))
    
  }
  )
  
  #reformat scores
  return(data.frame(t(scores)))
}

n_s <- c(0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9)

#iterate through n's
rmses <- sapply(n_s, function(n) {
  
  #obtain change-in-pain metrics for both studies
  ldopa <- change_in_pain(pain_ldopa, n)
  pl2 <- change_in_pain(pain_pl2, n)  
  
  
  #calculate rmse for both metrics in both studies
  data.frame(
    n_s = n,
    ldopa_rmse_exp = as.numeric(rmse(as.numeric(c(ldopa$metric_tot)), as.numeric(c(ldopa$metric_exp)))),
    ldopa_rmse_std = as.numeric(rmse(as.numeric(c(ldopa$metric_tot)), as.numeric(c(ldopa$metric_std)))),
    pl2_rmse_exp = as.numeric(rmse(as.numeric(c(pl2$metric_tot)), as.numeric(c(pl2$metric_exp)))),
    pl2_rmse_std = as.numeric(rmse(as.numeric(c(pl2$metric_tot)), as.numeric(c(pl2$metric_std))))
  )
  
})

rmses <- as.data.frame(t(rmses))

ldopa_df <- rmses %>%
  select(n_s, ldopa_rmse_exp, ldopa_rmse_std) %>%
  gather(key = "metrictype", value = "rmse", -n_s)

ldopa_df <- as.data.frame(lapply(ldopa_df, unlist))

ldopa_plot <- ggplot(ldopa_df, aes(x = n_s, y = rmse)) +
  geom_line(aes(color = metrictype, linetype = metrictype)) + 
  scale_color_manual(values = c("darkred", "steelblue"))

pl2_df <- rmses %>%
  select(n_s, pl2_rmse_exp, pl2_rmse_std) %>%
  gather(key = "metrictype", value = "rmse", -n_s)

pl2_df <- as.data.frame(lapply(pl2_df, unlist))

pl2_plot <- ggplot(pl2_df, aes(x = n_s, y = rmse)) +
  geom_line(aes(color = metrictype, linetype = metrictype)) + 
  scale_color_manual(values = c("darkred", "steelblue"))


ldopa_plot + ggtitle("RMSE Between Metrics by Percen. of Seq. for LDopa") + xlab("Percentage of Sequence") + ylab("RMSE")

ggsave("ldopa.png", width=8, height=4)

pl2_plot + ggtitle("RMSE Between Metrics by Percen. of Seq. for PL2") + xlab("Percentage of Sequence") + ylab("RMSE")

ggsave("pl2.png", width=8, height=4)