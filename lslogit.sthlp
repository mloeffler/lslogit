{smcl}
{* *! version 1.2.0  28dec2012}{...}
{findalias asfradohelp}{...}
{vieweralsosee "" "--"}{...}
{vieweralsosee "[R] help" "help help"}{...}
{viewerjumpto "Syntax" "lslogit##syntax"}{...}
{viewerjumpto "Description" "lslogit##description"}{...}
{viewerjumpto "Options" "lslogit##options"}{...}
{viewerjumpto "Remarks" "lslogit##remarks"}{...}
{title:Title}

{phang}
{bf:lslogit} {hline 2} Mixed logit labor supply model


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
{newvar}
{ifin}
[{cmd:,} {cmd:xb} {cmd:pc1}]


{synoptset 21 tabbed}{...}
{synopthdr:lslogit options}
{synoptline}
{syntab:Model}
{p2coldent:* {opth gr:oup(varname)}}matched group variable{p_end}
{synopt:{opt quad:ratic}}quadratic utility function{p_end}
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
{cmd:fweight}s, {cmd:iweight}s, and {cmd:pweight}s are allowed (see {help weight}),
but they are interpreted to apply to the group, not to individual observations. See
{mansection R clogitRemarksUseofweights:{it:Use of weights}} in {bf:[R] clogit}.{p_end}


{synoptset 21 tabbed}{...}
{synopthdr:lslpred options}
{synoptline}
{synopt:{opt xb}}predict systematic utility{p_end}
{synopt:{opt pc1}}predict choice probabilities{p_end}
{synoptline}
{p2colreset}{...}


{marker description}{...}
{title:Description}

{pstd}
{cmd:lslogit} fits mixed logit models with particular focus on the estimation
of structural labor supply models. It allows for several different utility
specifications as well as observed and unobserved heterogeneity in
preferences. Maximum simulated likelihood methods are used to incorporate
normally distributed random coefficients (Train, 2009).

{pstd}
{cmd:lslpred} calculates the whatever statistic for the variables in
{varlist} when the data are not stratified.


{marker options}{...}
{title:Options}

{dlgtab:Blab}

{phang}
{opt detail} displays detailed output of the calculation.

{phang}
{opt meanonly} restricts the calculation to be based on only the means.  The
default is to use a trimmed mean.

{phang}
{opt format} requests that the summary statistics be displayed using
the display formats associated with the variables, rather than the default
{cmd:g} display format; see
{findalias frformats}.

{phang}
{opt separator(#)} specifies how often to insert separation lines
into the output.  The default is {cmd:separator(5)}, meaning that a
line is drawn after every 5 variables.  {cmd:separator(10)} would draw a line
after every 10 variables.  {cmd:separator(0)} suppresses the separation line.

{phang}
{opth generate(newvar)} creates {it:newvar} containing the whatever
values.


{marker remarks}{...}
{title:Remarks}

{pstd}
For detailed information on the whatever statistic, see {bf:[R] intro}.


{marker author}{...}
{title:Author}

{phang}
{cmd:lslogit} was written by Max Löffler (loeffler@iza.org), Institute for the
Study of Labor (IZA), Bonn, Germany.
    

