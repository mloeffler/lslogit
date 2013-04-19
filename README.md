lslogit
=======

*lslogit* is a [Stata](http://www.stata.com/) command that fits mixed logit models with particular focus on the estimation of structural labor supply models. It uses maximum simulated likelihood methods ([Train, 2009](#train_2009)). It makes use of Stata's maximum likelihood framework and is implemented as d2-evaluator ([Gould et al., 2010](#gould_etal_2010)). *lslogit* is written in Mata, Stata's matrix programming language.

## Features
- Flexible utility specifications (translog, quadratic, Box-Cox).
- Observed heterogeneity in preferences as well as alternative specific taste shifters.
- Unobserved preference heterogeneity via multivariate normally distributed random coefficients.
- Integrating out the wage prediction error in order to use of the full wage distribution.
- Automatic calculation of the marginal utility with respect to consumption.
- Allows to impose contraints on the marginal utility of consumption.
- Simultaneous estimation of wages and preferences (work in progress).

## Thanks
- [Arne Risa Hole](http://www.shef.ac.uk/economics/people/hole) for his excellent Stata command [mixlogit](http://www.shef.ac.uk/economics/people/hole/stata).

## Reference
- <a name="train_2009"></a>Train, K. E. (2009). [Discrete Choice Methods with Simulation](http://elsa.berkeley.edu/books/choice2.html), second edn, Cambridge University Press.
- <a name="gould_etal_2010"></a>Gould, W., Pitblado, J., and Poi, B. (2010). [Maximum Likelihood Estimation with Stata](http://www.stata.com/bookstore/maximum-likelihood-estimation-stata/), fourth edn, Stata Press.
