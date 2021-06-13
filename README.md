# BayesianApproxHSMM

BayesianApproxHSMM contains the stan software (as well as R utilities) to model time series and sequential data using a  Bayesian  HMM that is a reformulation of any given HSMM as detailed in "Bayesian Approximations to Hidden Semi-Markov Models" (2020) by B.Hadj-Amar, J.Jewson, M.Fiecas. The article and supplemental material can be found at https://www.researchgate.net/publication/342200521_Bayesian_Approximations_to_Hidden_Semi-Markov_Models

* Main Inference Scripts:
```
tutorial.Rmd, case_study.Rmd
```



## Contents

This software is illustrated in two studies.  First, we generate data from a three state HSMM with Gaussian emission and Poisson dwell durations. Estimation procedures are carried out using both frequentist (EM) and Bayesian approaches, where we validate each model using several diagnostic tools. We also illustrate the use of Bayes factors for model selection between our proposed approach and the HMM. Second, we analyze the physical activity time series that was investigated using a frequentist HMM by Huang et al. (2018). We seek to conduct a similar study but within a Bayesian framework and considering the 
extra flexibility afforded by our proposed methodology to investigate departures from the HMM. 


We implemented code for the following dwell durations: geometric (i.e. HMM), poisson, negative-binomial and the unstructured start geometric tail distribution (Sansom and Thomson, 2001). The probabilistic programming framework associated with stan makes it easy for practitioners to consider further dwell distributions to the ones considered here. Users need only to change the corresponding function in our stan files. Moreover, the stan modelling language implements compressed row storage sparse matrix representation and multiplication, which provides considerable speed up when the sparsity is greater than 90%.
In our applied scenario we consider dwell-approximation thresholds as big as a = (150,10,10) with sparsity of greater than 98% allowing us to take considerable advantage of this formulation.


## boh 1
![GitHub Logo](/figures/posterior_pred_freq.png){:height="50%" width="50%"}

## boh 2
![alt text](/figures/posterior_pred_freq.png){:height="50%" width="50%"}

## boh 3
![plot](/figures/posterior_pred_freq.png){:height="50%" width="50%"}

### Contact Information and Acknowledgments

BayesianApproxHSMM was developed by Beniamino Hadj-Amar (Beniamino.Hadj-Amar@rice.edu) and Jack Jewson (jack.jewson@barcelonagse.eu). 
Parts of this software were adapted from "Hidden Markov Models for Time Series" (Zucchini et al., 2017) and  "Hidden Markov Models" (Stan User's Guide, Sec.2.6).




