lslogit
=======

*lslogit* is a [Stata](http://www.stata.com/) command that fits mixed logit models with particular focus on the estimation of structural labor supply models. It uses maximum simulated likelihood methods ([Train, 2009](#train_2009)). It makes use of Stata's maximum likelihood framework and is implemented as d2-evaluator. *lslogit* is written in Mata, Stata's matrix programming language.

## Features
- Flexible utility specifications (translog, quadratic, Box-Cox).
- Observed heterogeneity in preferences as well as alternative specific taste shifters.
- Unobserved preference heterogeneity via multivariate normally distributed random coefficients.
- Integrating out the wage prediction error in order to use of the full wage distribution.
- Automatic calculation of the marginal utility with respect to consumption.
- Simultaneous estimation of wages and preferences (work in progress).

## Thanks
- [Arne Risa Hole](http://www.shef.ac.uk/economics/people/hole) for his excellent Stata command [mixlogit](http://www.shef.ac.uk/economics/people/hole/stata).

## Reference
- <a id="train_2009">Train, K. E. (2009)</a>. *Discrete Choice Methods with Simulation*, second edn, Cambridge University Press.
