lslogit
=======

*lslogit* is a [Stata](http://www.stata.com/) command that fits mixed logit models with particular focus on the estimation of structural labor supply models. It uses maximum simulated likelihood methods (Train, 2009). *lslogit* is written in Mata, Stata's matrix programming language.

## Features
- Flexible utility specifications (translog, quadratic, Box-Cox).
- Observed and unobserved heterogeneity in preferences via multivariate normally distributed random coefficients.
- Integrating out the wage prediction error.
- Simultaneous estimation of wages and preferences.
- Automatic calculation of the marginal utility with respect to consumption.

## Thanks
- [Arne Risa Hole](http://www.shef.ac.uk/economics/people/hole) for his excellent Stata command [mixlogit](http://www.shef.ac.uk/economics/people/hole/stata).

## Reference
- Train, K. E. (2009). *Discrete Choice Methods with Simulation*, second edn, Cambridge University Press.