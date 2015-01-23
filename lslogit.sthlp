{smcl}
{* *! version 0.4, 23jan2015}{...}
{vieweralsosee "[R] clogit" "help clogit"}{...}
{vieweralsosee "[R] mixlogit" "help mixlogit"}{...}
{vieweralsosee "lslogit on GitHub" "browse https://github.com/mloeffler/lslogit/"}{...}
{viewerjumpto "Syntax" "lslogit##syntax"}{...}
{viewerjumpto "Description" "lslogit##description"}{...}
{viewerjumpto "Options" "lslogit##options"}{...}
{viewerjumpto "References" "lslogit##references"}{...}
{viewerjumpto "Author" "lslogit##author"}{...}
{title:Title}

{phang}
{bf:lslogit} {hline 2} Estimate mixed logit labor supply models


{marker syntax}{...}
{title:Syntax}

{p 8 17 2}
{cmdab:lslogit}
    {depvar}
    {ifin}
    {weight}
    {cmd:,}
    {cmdab:gr:oup:(}{varname}{cmd:)}
    {cmdab:c:onsumption:(}{varname}{cmd:)}
    {cmdab:l:eisure:(}{varlist}{cmd:)}
    [{it:options}]

{p 8 17 2}
{cmdab:lslpred}
    {it:{help newvarlist}}
    {ifin}
    [{cmd:,} {cmd:pc1} {cmd:xb}]


{synoptset 21 tabbed}{...}
{synopthdr:lslogit options}
{synoptline}
{syntab:Model}
{p2coldent:* {opth gr:oup(varname)}}matched group variable{p_end}
{synopt:{opt quad:ratic}}quadratic utility function (default){p_end}
{synopt:{opt tran:slog}}translog utility function{p_end}
{synopt:{opt boxc:ox}}Box-Cox utility function{p_end}
{synopt:{opt boxcc(#)}}normalize consumption for Box-Cox transformation{p_end}
{synopt:{opt boxcl(#)}}normalize leisure for Box-Cox transformation{p_end}
{synopt:{opt dr:aws(#)}}number of Halton draws used to approximate{p_end}
{synopt:{opt burn(#)}}number of initial Halton draws to burn{p_end}
{synopt:{opt [no]round}}enable/disable rounding of hourly and monthly wages{p_end}

{syntab:Right hand side}
{p2coldent:* {opth c:onsumption(varname)}}consumption variable{p_end}
{p2coldent:* {opth l:eisure(varlist)}}leisure variable(s){p_end}
{synopt:{opth cx(varlist)}}interactions with consumption{p_end}
{synopt:{opth lx1(varlist)}}interactions with leisure (person 1){p_end}
{synopt:{opth lx2(varlist)}}interactions with leisure (person 2){p_end}
{synopt:{opth c2x(varlist)}}interactions with consumption squared{p_end}
{synopt:{opth l2x1(varlist)}}interactions with leisure squared (person 1){p_end}
{synopt:{opth l2x2(varlist)}}interactions with leisure squared (person 2){p_end}
{synopt:{opth ind:eps(varlist)}}taste shifters{p_end}

{syntab:Random preferences}
{synopt:{opt rand:vars(coeflist)}}indices of random coefficients{p_end}
{synopt:{opt corr}}correlation between random coefficients{p_end}

{syntab:Wage estimation}
{synopt:{opth hw:age(varlist)}}hourly wage rate(s){p_end}
{synopt:{opth heckm:an(varlist)}}wage equation{p_end}
{synopt:{opth select(varlist)}}selection equation{p_end}

{syntab:Wage prediction}
{synopt:{opth wagep:red(varlist)}}integrate out wage prediction error(s){p_end}
{synopt:{opth day:s(varname)}}days per year{p_end}
{synopt:{opt totalt:ime(#)}}time endowment per week{p_end}
{synopt:{opt hecksig:ma(sigma)}}standard error(s) of the wage regression(s){p_end}
{synopt:{opt taxr:eg(name)}}estimates of the tax regrssion{p_end}
{synopt:{opth tria1(varlist)}}interactions with earnings (person 1){p_end}
{synopt:{opth tria2(varlist)}}interactions with earnings (person 2){p_end}
{synoptline}
{p2colreset}{...}
{p 4 6 2}
* {opt group(varname)}, {opt consumption(varname)} and {opt leisure(varlist)} are required.{p_end}

{p 4 6 2}
{cmd:fweight}s are allowed (see {help weight}), but they are interpreted to
apply to the group, not to individual observations. See
{mansection R clogitRemarksUseofweights:{it:Use of weights}} in {bf:[R] clogit}.{p_end}


{synoptset 21 tabbed}{...}
{synopthdr:lslpred options}
{synoptline}
{synopt:{opt pc1}}predict choice probabilities (default){p_end}
{synopt:{opt xb}}predict systematic utility{p_end}
{synoptline}
{p2colreset}{...}
{p 4 6 2}
If two variables are specified in {it:{help newvarlist}}, the first will contain
the choice probabilities, the latter will contain the systematic utilities.


{marker description}{...}
{title:Description}

{pstd}
{cmd:lslogit} fits mixed logit models with particular focus on the estimation
of structural labor supply models in the discrete choice context. It allows for
several different functional forms, very flexible utility specifications,
both observed and unobserved heterogeneity in preferences and can be easily
augmented with flexible dummy refinement. Moreover, it allows to integrate out
the wage prediction error during the estimation.

{pstd}
Technically, it is an extension to the Stata command {help clogit} and the
user-written routine {help mixlogit} ({help lslogit##hole_2007:Hole, 2007}). {cmd:lslogit}
makes use of maximum simulated likelihood methods
({help lslogit##train_2009:Train, 2009}).

{pstd}
{cmd:lslpred} predicts the systematic utility and/or the choice probabilities
according to the estimated coefficients.


{marker options}{...}
{title:Options}

{dlgtab:Model specification}

{phang}
{opth group(varname)} specifies an identifier variable for the matched groups.

{phang}
{opt quadratic} specifies the functional form of the systematic utility to be
quadratic ({help lslogit##keane_moffitt_1998:Keane and Moffitt, 1998}). If no
functional form is specified, {cmd:lslogit} will assume a quadratic utility
function.

{phang}
{opt translog} specifies a translog utility function ({help lslogit##vansoest_1995:van Soest, 1995}).

{phang}
{opt boxcox} specifies a Box-Cox transformed utility function
({help lslogit##aaberge_etal_1995:Aaberge et al., 1995}). The likelihood
function converges more easily when consumption and leisure variables
have a smaller range. Therefore, consumption is divided by {opt boxcc(#)} and
leisure is divided by {opt boxcl(#)} (defaults are {bf:boxcl(1000)} and
{bf:boxcl(80)}).

{phang}
{opt draws(#)} is the number of Halton draws that are used to approximate
the likelihood contribution of each group (default is {bf:draws(50)}).

{phang}
{opt burn(#)} denotes the number of highly correlated initial Halton draws
to be burned (default is {bf:burn(15)}).

{phang}
{opt [no]round} enables/disables rounding of hourly and monthly wages. Rounding
may cause problems when estimating wages and working hours jointly because of 
non-concave parts in the likelihood function.

{dlgtab:Right hand side variables}

{phang}
{opth consumption(varname)} specifies the variable which contains the consumption
or disposable income, respectively, for every choice for every household. Do not
transform the consumption variable according to the chosen functional form
of the utility function as this will be done by {cmd:lslogit} automatically.

{pmore}
Observed heterogeneity in the preferences for consumption can be introduced by
adding interaction variables using the {opth cx(varlist)} option. Interaction with
consumption squared can be specified by {opth c2x(varlist)}.

{phang}
{opth leisure(varlist)} specifies one or two variables containing the amount of
leisure time for the choice alternatives. Again, do not transform the variables
according to the functional form of the utility function. The first leisure
variable is referred to as {it:L1}. If {cmd:lslogit} finds a second leisure
variable, it assumes a two-decision-maker household and denotes the second leisure
term as {it:L2}.

{pmore}
By specifying {opth lx1(varlist)} and {opth l2x1(varlist)} one can introduce
preference heterogeneity in leisure and leisure squared for the first decision
maker. {opth lx2(varlist)} and {opth l2x2(varlist)} denote the interation variables
for the second decision maker in the household.

{phang}
{opth indeps(varlist)} allows to specify taste shifters that differ over choice
alternatives, e.g., fixed costs, disutility from welfare participation or hours
restriction.

{dlgtab:Random preferences}

{phang}
{opth randvars(coeflist)} specifies the indices of the coefficients that are
assumed to be random. Thus, in a model without observed heterogeneity in
preferences for consumption, specifying {bf:randvars(1)} says that the coefficient
on consumption follows a normal distribution. The mean is given by the corresponding
estimate, the estimated standard deviation will be denoted as {it:sd_1}. Indices
are separated by spaces, e.g., {bf:randvars(1 4 6)}. The standard deviations will
be denoted as {it:sd_1}, {it:sd_4} and {it:sd_6}.

{phang}
{opt corr} tells the estimation routine to allow for correlation between all random
coefficients. In the example above, this will lead to additional covariance estimates
{it:sd_1_4}, {it:sd_1_6} and {it:sd_4_6}.

{pmore}
Please note that the sign of the estimates is irrelevant, just interpret them as
being positive ({help lslogit##hole_2007:Hole, 2007}).

{dlgtab:Wage prediction handling}

{phang}
{opth wagepred(varlist)} specifies one or two dummy variables (depending on the
number of decision makers in the household) denoting whether or not the wage
prediction error has to be integrated out during the estimation.

{phang}
{opth days(varname)} allows to specify a variable containing the number of days
within the tax year. This may be helpful when using pooled cross-sectional data
of different years. (By default, {cmd:lslogit} assumes that a year has 365 days.)

{phang}
{opt totaltime(#)} denotes the time endowment per week (default is
{bf:totaltime(80)}). Working hours are calculated as difference between the time
endowment and the amout of leisure for every choice alternative.

{phang}
{opt hecksigma(sigma)} specifies the estimated standard deviation(s) of the wage
equation (separated by spaces for two decision makers in the household). This is
needed in order to integrate the wage prediction error out.

{phang}
{opt taxreg(name)} specifies the name of the stored "tax regression"
estimates. {cmd:lslogit} makes use of these estimates to predict the consumption
for every possible draw from the wage distribution. The "tax regression" has
to follow a specific setup:

{pmore}
{it:dpi = a*mw + b*mw^2 + c*mw*ia + d*mw*ia^2 + e*f(x,ia) + f + g}

{pmore}
with {it:dpi} as disposable income, {it:mw} denotes monthly earnings, {it:ia}
refers to a set of interaction variables, {it:x} is a set of characteristics
unrelated to earnings, e.g., non-labor income, children or age. {it:a, b, c, d, e}
and {it:f} denote coefficients to be estimated, {it:g} is a normally distributed
error term.

{pmore}
Interactions for decision maker one can be specified using {opth tria1(varlist)}. If
a second decision maker lives in the household, interactions can be specified
using {opth tria2(varlist)}. His/her monthly earnings, squared earnings and interacted
earnings enter the regression equation between {it:d*mw*ia^2} and
{it:e*f(x,ia)} in the example above, the sorting order remains the same as for
the first decision maker.


{marker references}{...}
{title:References}

{marker aaberge_etal_1995}{...}
{phang}
Aaberge, R., Dagsvik, J. K. and Strøm, S. (1995). Labor supply responses and welfare effects of tax reforms,
{it:The Scandinavian Journal of Economics} 97(4): 635-659.

{marker hole_2007}{...}
{phang}
Hole, A. R. (2007). Fitting mixed logit models by using maximum simulated likelihood,
{it:Stata Journal} 7(3), 388-401.

{marker keane_moffitt_1998}{...}
{phang}
Keane, M. P. and Moffitt, R. (1998). A structural model of multiple welfare program participation and labor supply,
{it:International Economic Review} 39(3): 553-589.

{marker train_2009}{...}
{phang}
Train, K. E. (2009). {it:Discrete Choice Methods with Simulation}. 2nd ed. Cambridge University Press.

{marker vansoest_1995}{...}
{phang}
van Soest, A. (1995). Structural models of family labor supply -- a discrete choice approach,
{it:The Journal of Human Resources} 30(1): 63-88.
{p_end}


{marker author}{...}
{title:Author}

{pstd}
{cmd:lslogit} was written by Max Löffler ({browse "mailto:loeffler@zew.de":loeffler@zew.de}),
Centre for European Economic Research (ZEW), Mannheim, Germany. Comments and
suggestions are welcome.


{marker license}{...}
{title:License}

{pstd}
Copyright (C) 2012-2015, {browse "mailto:loeffler@zew.de":Max Löffler}

{pstd}
This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

{pstd}
This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

{pstd}
You should have received a copy of the GNU General Public License along
with this program; if not, write to the Free Software Foundation, Inc.,
51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

