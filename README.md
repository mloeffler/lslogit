Stata package lslogit
=====================

*lslogit* is a [Stata](http://www.stata.com/) command that fits mixed logit models with particular focus on the estimation of structural labor supply models. It uses maximum simulated likelihood methods ([Train, 2009](#train_2009)). It makes use of Stata's maximum likelihood framework and is implemented as d2-evaluator ([Gould et al., 2010](#gould_etal_2010)). *lslogit* is written in Mata, Stata's matrix programming language.

I presentend a preview of *lslogit* at the [Stata Conference in New Orleans](http://www.stata.com/meeting/new-orleans13/) in July 2013. You can find the slides of this presentation [here](https://ideas.repec.org/p/boc/norl13/8.html).


## Features
- Flexible utility specifications (translog, quadratic, Box-Cox).
- Observed heterogeneity in preferences as well as alternative specific taste shifters.
- Unobserved preference heterogeneity via multivariate normally distributed random coefficients.
- Integrating out the wage prediction error in order to use of the full wage distribution.
- Automatic calculation of the marginal utility with respect to consumption.
- Allows to impose contraints on the marginal utility of consumption.
- Simultaneous estimation of wages and preferences (work in progress).


## Installation

You can install the latest version of *lslogit* it by typing:

	. net from https://mloeffler.github.io/stata/
	. net install lslogit

Done.


## Thanks
- [Arne Risa Hole](http://www.shef.ac.uk/economics/people/hole) for his excellent Stata command [mixlogit](http://www.shef.ac.uk/economics/people/hole/stata).


## Reference
- <a name="train_2009"></a>Train, K. E. (2009). [Discrete Choice Methods with Simulation](http://elsa.berkeley.edu/books/choice2.html), second edn, Cambridge University Press.
- <a name="gould_etal_2010"></a>Gould, W., Pitblado, J., and Poi, B. (2010). [Maximum Likelihood Estimation with Stata](http://www.stata.com/bookstore/maximum-likelihood-estimation-stata/), fourth edn, Stata Press.


## Info

Copyright (C) 2012-2015, [Max Löffler](http://www.zew.de/en/staff/mlo)

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License along
with this program; if not, write to the Free Software Foundation, Inc.,
51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

