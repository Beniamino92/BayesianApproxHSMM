# BayesianApproxHSMM

BayesianApproxHSMM contains the stan software (as well as R utilities) to model time series and sequential data using a  Bayesian  HMM that is a reformulation of any given HSMM as detailed in "Bayesian Approximations to Hidden Semi-Markov Models for Telemetric Monitoring of Physical Activity" (2022) by B.Hadj-Amar, J.Jewson, M.Fiecas. The article and supplemental material can be found at https://www.researchgate.net/publication/361145205_Bayesian_Approximations_to_Hidden_Semi-Markov_Models_for_Telemetric_Monitoring_of_Physical_Activity

* Main Inference Scripts:
```
tutorial.Rmd, case_study.Rmd, case_study_frequency.jl
```



## Overview

This software is illustrated in two studies.  First, in `tutorial.Rmd` we generate data from a three state HSMM with Gaussian emission and Poisson dwell durations. Estimation procedures are carried out using both frequentist (EM) and Bayesian approaches, where we validate each model using several diagnostic tools. We also illustrate the use of Bayes factors for model selection between our proposed approach and the HMM. Second, in `case_study.Rmd` we analyze the physical activity time series that was investigated using a frequentist HMM by Huang et al. (2018). We seek to conduct a similar study but within a Bayesian framework and considering the 
extra flexibility afforded by our proposed methodology to investigate departures from the HMM. 

<!-- Finally, `case_study_frequency.jl` illustrates the use of our Metropolis-within-Gibbs sampler to estimate the frequency that drives the overall variation in the physical activity time series. 
 -->

We implemented code for the following dwell durations: geometric (i.e. HMM), poisson, negative-binomial and the unstructured start geometric tail distribution (Sansom and Thomson, 2001). The probabilistic programming framework associated with stan makes it easy for practitioners to consider further dwell distributions to the ones considered here. Users need only to change the corresponding function in our stan files. Moreover, the stan modelling language implements compressed row storage sparse matrix representation and multiplication, which provides considerable speed up when the sparsity is greater than 90%.
In our applied scenario we consider dwell-approximation thresholds as big as a = (150,10,10) with sparsity of greater than 98% allowing us to take considerable advantage of this formulation.

## Example - BayesApproxHSMM

Here is an example for using our (stan) software in R from `tutorial.Rmd`.

```r
K <- 3 # n states 
m <- rep(5, K) # dwell threshold

data.stan <- list(N = length(obs), K = K, y = obs,
                  m = m,  mu_0 = rep(mean(obs), K), 
                  sigma_0 = 2, a_0 = rep(0.01, K), b_0 = rep(0.01, K),
                  alpha_0 = matrix(1, nrow = K, ncol = K-1))

if ((K / sum(m) < 0.1)) {
  stan_path <- "stan/bayesHSMMapprox_GaussEmis_PoissDur.stan"
} else {
  stan_path <- "stan/bayesHSMMapprox_GaussEmis_PoissDur_sparse.stan"
}

HSMM.stan <- stan(file = stan_path, data = data.stan, 
                 init = function(){HSMM.init.stan(K, obs, rep(10, K))}, 
                 warmup = 1000, chains = 1, iter = (1+5)*1000, cores = 1, 
                 control = list(adapt_delta=0.99, stepsize=0.01, max_treedepth = 20))
```                 
   
   
        
```r
sims <- extract(HSMM.stan)
HSMM.predictive.plot(sims, obs, m, ndraw = 50)
```


<p align="center">
<img src="https://github.com/Beniamino92/BayesianApproxHSMM/blob/master/figures/example_postpred.png" width="500" heigth="250"/> 
</p>


## Identifying the Periodicity 

We have also developed a (high-performance) software written  in  Julia  v1.6  to identify the frequency that drives the overall variation in the  time series. We ran  the  sampler  for  5000 iterations  which took around 3 seconds on an Intel®CoreTMi5 2 GHz Processor 16 GB RAM. An example of our implementation `case_study_frequency.jl` is given below, where we also highlight the visual output of our software: (left) the trace plot (after burn-in) of the posterior sample of the frequency, noting that the acceptance rate was 28%; (right) 20 draws from the posterior predictive distribution of the stationary periodic model, as well as the posterior mean of the oscillatory signal.



```julia
M = 1 # number of frequencies
GibbsSamplerOscillatory(obs, M, N_MCMC, hyperparms; plt = true)
```

<p align="center">
<img src="https://github.com/Beniamino92/BayesianApproxHSMM/blob/master/figures/posterior_freq.png" width="400" heigth="400"/> <img src="https://github.com/Beniamino92/BayesianApproxHSMM/blob/master/figures/posterior_pred_freq.png" width="400" heigth="400"/>
</p>


### Contact Information and Acknowledgments

BayesianApproxHSMM was developed by Beniamino Hadj-Amar (Beniamino.Hadj-Amar@rice.edu) and Jack Jewson (jack.jewson@barcelonagse.eu). 
Parts of this software were adapted from "Hidden Markov Models for Time Series" (Zucchini et al., 2017), "Hidden Markov Models" (Stan User's Guide, Sec.2.6), "Bayesian Model Search for Nonstationary Periodic Time Series" (Hadj-Amar et al., 2019).



