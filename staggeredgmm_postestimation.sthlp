{smcl}
{* *! version 0.1.0  30oct2026}{...}
{vieweralsosee "staggeredgmm" "help staggeredgmm"}{...}
{viewerjumpto "Postestimation commands" "staggeredgmm_postestimation##description"}{...}
{viewerjumpto "estat jtest" "staggeredgmm_postestimation##jtest"}{...}
{viewerjumpto "Examples" "staggeredgmm_postestimation##examples"}{...}
{viewerjumpto "Stored results" "staggeredgmm_postestimation##results"}{...}
{title:Title}

{phang}
{bf:staggeredgmm postestimation} {hline 2} Postestimation tools for staggeredgmm


{marker description}{...}
{title:Postestimation commands}

{pstd}
The following postestimation command is of special interest after {cmd:staggeredgmm}:

{synoptset 20}{...}
{p2coldent:Command}Description{p_end}
{synoptline}
{synopt:{helpb staggeredgmm_postestimation##jtest:estat jtest}}specification test of parallel trends and no anticipation{p_end}
{synoptline}

{pstd}
The following standard postestimation commands are also available:

{synoptset 20}{...}
{p2coldent:Command}Description{p_end}
{synoptline}
{synopt:{helpb estat vce}}variance-covariance matrix of the estimators{p_end}
{synopt:{helpb estimates}}cataloging estimation results{p_end}
{synopt:{helpb lincom}}point estimates, standard errors, testing, and inference for linear combinations of coefficients{p_end}
{synopt:{helpb nlcom}}point estimates, standard errors, testing, and inference for nonlinear combinations of coefficients{p_end}
{synopt:{helpb test}}Wald tests of simple and composite linear hypotheses{p_end}
{synoptline}


{marker jtest}{...}
{title:Syntax for estat jtest}

{p 8 17 2}
{cmd:estat} {cmd:jtest} [{cmd:,} {opt full} {opt w:indow(#)}]

{synoptset 20}{...}
{synopthdr}
{synoptline}
{synopt:{opt full}}use every pre-treatment placebo moment{p_end}
{synopt:{opt w:indow(#)}}width of the local pre-window, in periods; default {cmd:3}; ignored with {opt full}{p_end}
{synoptline}


{title:Description of estat jtest}

{pstd}
{cmd:estat jtest} is a serial-correlation robust joint test of the identifying
assumptions, following Section 4.4 of Arora and Bijani (2026). It is built on
the never-treated event-study basis: fixing a base pre-period t0, each treated
cohort g and period t != t0 defines

{p 8 8 2}
D(g,t) = (Ybar(g,t) - Ybar(g,t0)) - (Ybar(inf,t) - Ybar(inf,t0)),

{pstd}
where inf denotes the never-treated group. Those with t >= g identify the
cohort-by-time effects. Those with t < g are placebo restrictions whose
expectation is zero only under parallel trends and no anticipation. Collecting
J placebo moments into p with covariance Omega_p implied by the fitted
serial-correlation model, the statistic p' inv(Omega_p) p is asymptotically
chi-squared with J degrees of freedom under the null. Because the placebo
moments do not involve the treatment effects, no degrees of freedom are
subtracted for estimated parameters. If Omega_p is rank deficient, a truncated
inverse is used and the degrees of freedom are reduced accordingly, with a
warning.

{pstd}
The base period is t0 = g_min - 1, the last period before any cohort adopts, so
that every moment has an untreated reference.

{pstd}
A large number of restrictions degrades the chi-squared approximation in small
panels. By default the test therefore uses only the placebo moments within
{opt window()} periods before each cohort's adoption, which holds its nominal
size across sample sizes at essentially no cost in power, since a differential
trend or an anticipation effect is already visible just before treatment.
{opt full} uses every placebo moment and is well sized only when the number of
units is comfortably large.

{pstd}
The test requires a never-treated control group.


{marker examples}{...}
{title:Examples}

{phang2}{cmd:. sysuse sgmm_sim_panel, clear}{p_end}
{phang2}{cmd:. staggeredgmm y, ivar(unit_id) time(year) gvar(cohort)}{p_end}
{phang2}{cmd:. estat jtest}{p_end}
{phang2}{cmd:. estat jtest, window(5)}{p_end}
{phang2}{cmd:. estat jtest, full}{p_end}
{phang2}{cmd:. matrix list r(moments)}{p_end}


{marker results}{...}
{title:Stored results}

{pstd}
{cmd:estat jtest} stores the following in {cmd:r()}:

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Scalars}{p_end}
{synopt:{cmd:r(J)}}test statistic{p_end}
{synopt:{cmd:r(df)}}degrees of freedom{p_end}
{synopt:{cmd:r(p)}}p-value{p_end}
{synopt:{cmd:r(base_period)}}base period t0{p_end}
{synopt:{cmd:r(n_moments)}}number of placebo moments{p_end}
{synopt:{cmd:r(n_dropped)}}placebo moments dropped for referencing an unobserved group-period{p_end}
{synopt:{cmd:r(window)}}width of the local pre-window (local test only){p_end}

{p2col 5 20 24 2: Macros}{p_end}
{synopt:{cmd:r(type)}}{cmd:local} or {cmd:full}{p_end}
{synopt:{cmd:r(weighting)}}weighting scheme of the fitted model{p_end}

{p2col 5 20 24 2: Matrices}{p_end}
{synopt:{cmd:r(moments)}}placebo moments: {cmd:g}, {cmd:t}, {cmd:event_time}, {cmd:moment}{p_end}


{title:Reference}

{phang}
Arora, P. and R. Bijani. 2026. Estimating treatment effects under staggered
timing and non-spherical errors. Available at SSRN,
{browse "https://doi.org/10.2139/ssrn.6558759":doi:10.2139/ssrn.6558759}.
{p_end}
