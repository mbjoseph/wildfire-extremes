# Extreme value theory sandbox

Earth Lab Data Meetup, July 29, 2016

## What's here? 

- Toy data in `data/d.csv` that contains timesteps and observations
- Script to explore the data: `R/eda.R`
- Script to estimate parameters for the generalized extreme value distribution using maximum likelihood `R/mle.R` or Bayesian approaches `R/bayes.R`
- Script to simulate block maxima as a function of within block sample size: `R/block_maxima_varying_n.R`
- Some helper functions `R/helpers.R` 
- [Stan](http://mc-stan.org/) model statement for the Bayesian implementation `stan/gev.stan`

To get a local copy of this repository:

```
git clone https://github.com/mbjoseph/demo_evt.git
```
