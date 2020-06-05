# BayesianApproxHSMM

BayesianApproxHSMM contains the stan software (as well as R utilities) to model time series and sequential data via a Bayesian  HMM that is a reformulation of any given HSMM as detailed in "Bayesian Approximation of Hidden Semi-Markov Models" (2020) by B.Hadj-Amar, J.Jewson, M.Fiecas. 

## Contents

This software is illustrated in two studies.  First, we generate data from a three state HSMM with Gaussian emission and Poisson dwell durations. Estimation procedures are carried out using both frequentist (EM) and Bayesian approaches, where we validate each model using several diagnostic tools. We also illustrate the use of Bayes factors for model selection between our proposed approach and the HMM. Second, we analyze the physical activity time series that was investigated using a frequentist HMM by Huang et al. (2018). We seek to conduct a similar study but within a Bayesian framework and considering the 
extra flexibility afforded by our proposed methodology to investigate departures from the HMM. 


The probabilistic programming framework associated with stan makes it easy for practitioners to consider further dwell distributions to the ones considered here. Users need only to change the corresponding function in our stan files.

* Main Inference Scripts:
```
tutorial.Rmd, case_study.Rmd
```


