
#Libraries
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(pander)
library(lubridate)
library(rbi)
library(rbi.helpers)

#SIR model

model_str <- "
model SIR_deterministic {
  const N = 41;
  obs Infected;
  state S,R,I1,I2,I3,I4,I5,I;
  param i_beta;
  param i_gamma;

  sub parameter {
    i_beta ~ gaussian(2.4, .5)
    i_gamma ~ gaussian(0.66, .1)
  }
  sub initial {
    inline I = I1 + I2 + I3 + I4 + I5 
    S <- N-1
    I <- 1
    I1 <- 1
    I2 <- 0
    I3 <- 0
    I4 <- 0
    I5 <- 0
  }
  sub transition {
    ode{
      dS/dt  = -(i_beta * S * I) / N
      dI1/dt =  (i_beta * S * I) / N - 5 * i_gamma * I1
      dI2/dt =  5 * i_gamma * (I1 - I2)
      dI3/dt =  5 * i_gamma * (I2 - I3)
      dI4/dt =  5 * i_gamma * (I3 - I4)
      dI5/dt =  5 * i_gamma * (I4 - I5)
      dI/dt = (i_beta * S * I) / N - 5 * i_gamma * I5
      dR/dt  =  5 * i_gamma * I5 
    }
   }
  sub observation {
  Infected ~ poisson(I)
  }

}"


# Create a libbi object
model <- bi_model(lines = stringi::stri_split_lines(model_str)[[1]])
bi_model <- libbi(model)



# Load the data
y <- data.frame(time = c(0,1,2,3,4,5,6), value = c(1,10,24,22,5,1,0))

#Define number of cores and minimum particles
ncores <- 2
minParticles <- max(ncores, 4)


input_lst <- list(N = 52196381)
bi <- sample(bi_model, end_time = 6,input = input_lst, obs = list(Incidence=y), nsamples = 10000, nparticles = minParticles, nthreads = ncores, proposal = 'prior') %>% 
  adapt_particles(min = minParticles, max = minParticles*200) %>%
 adapt_proposal(min = 0.05, max = 0.4) %>%
  sample(nsamples = 10000, thin = 5) %>% # burn in 
  sample(nsamples = 10000, thin = 5)
bi_lst <- bi_read(bi %>% sample_obs)


#Retrieve results and plot

  #Parameters  
       library('coda')
       traces <- mcmc(get_traces(bi))
       plot(traces)


   #Infected
   obs_bi <- sample_obs(bi)
   os <- summary(obs_bi, type="obs")

   ggplot(os, aes(x=time))+
       geom_line(aes(y=Median)) +
       geom_ribbon(aes(ymin=`1st Qu.`, ymax=`3rd Qu.`), alpha=0.5) +
       geom_point(aes(y=value), y, color="darkred") +
       ylab("Infected")
