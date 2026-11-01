{smcl}
{* *! version 0.1.0  30oct2026}{...}
{vieweralsosee "staggeredgmm postestimation" "help staggeredgmm_postestimation"}{...}
{viewerjumpto "Syntax" "staggeredgmm##syntax"}{...}
{viewerjumpto "Description" "staggeredgmm##description"}{...}
{viewerjumpto "Options" "staggeredgmm##options"}{...}
{viewerjumpto "Data requirements" "staggeredgmm##data"}{...}
{viewerjumpto "Remarks" "staggeredgmm##remarks"}{...}
{viewerjumpto "Examples" "staggeredgmm##examples"}{...}
{viewerjumpto "Stored results" "staggeredgmm##results"}{...}
{viewerjumpto "References" "staggeredgmm##references"}{...}
{viewerjumpto "Authors" "staggeredgmm##authors"}{...}
{title:Title}

{phang}
{bf:staggeredgmm}: Generalised Method of Moments Estimation of Treatment Effects under Staggered Adoption


{marker syntax}{...}
{title:Syntax}

{p 8 17 2}
{cmd:staggeredgmm}
{depvar}
[{it:covariates}]
{ifin}{cmd:,}
{opth iv:ar(varname)}
{opth ti:me(varname)}
{opth gv:ar(varname)}
[{it:options}]

{synoptset 28 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:Required}
{synopt:{opth iv:ar(varname)}}unit identifier; numeric or string{p_end}
{synopt:{opth ti:me(varname)}}calendar time; integer-valued{p_end}
{synopt:{opth gv:ar(varname)}}first treated period of each unit; {cmd:0} for never treated{p_end}

{syntab:Model}
{synopt:{opt weight:ing(scheme)}}covariance model for the optimal weight: {cmd:pooled_toeplitz} (default),
{cmd:cohort_toeplitz} or {cmd:unrestricted}{p_end}
{synopt:{opt never:treated}}require a never-treated control group{p_end}
{synopt:{opt nonever:treated}}do not use never-treated units as controls{p_end}

{syntab:Iteration}
{synopt:{opt max:iter(#)}}maximum number of iterated-GMM reweighting steps; default {cmd:100}{p_end}
{synopt:{opt tol:erance(#)}}convergence tolerance on the largest change in any CATT; default {cmd:1e-6}{p_end}

{syntab:Reporting}
{synopt:{opt l:evel(#)}}set confidence level; default is {cmd:level(95)}{p_end}
{synoptline}
{p 4 6 2}
{it:covariates} are baseline (pre-treatment, time-invariant) covariates and may contain factor variables; see {help fvvarlist}.{p_end}
{p 4 6 2}
See {help staggeredgmm_postestimation} for features available after estimation, including {cmd:estat jtest} and {cmd:estat aggregate}.{p_end}


{marker description}{...}
{title:Description}

{pstd}
{cmd:staggeredgmm} estimates cohort-by-time average treatment effects (CATTs) under
staggered treatment adoption by the generalized method of moments, following
Arora and Bijani (2026). Each CATT is identified by every clean 2x2
difference-in-differences comparison available for it: against never-treated
units and against cohorts that are not yet treated. The clean comparisons are
stacked into an over-identified moment system and combined by iterated optimal
GMM, with the weight built from a model of the serial correlation of the
errors.

{pstd}
Three covariance models are available through {opt weighting()}. They
correspond to GMM-T, GMM-HT and GMM-U in the paper.

{pstd}
Two aggregate ATTs are reported. The {it:treated-observation} ATT weights every
treated unit-period equally, so that cohorts observed for longer after
treatment carry more weight; this is the aggregation used in the paper and
coincides with the "simple" aggregation of Callaway and Sant'Anna (2021). The
{it:cohort-equal} ATT gives every treated cohort the same total weight, spread
equally over its own post-treatment periods.


{marker options}{...}
{title:Options}

{dlgtab:Required}

{phang}
{opth ivar(varname)} identifies units. It may be numeric or string.

{phang}
{opth time(varname)} gives the calendar period. It must be integer-valued, and
the set of periods present in the estimation sample, pooled across units, must
be consecutive.

{phang}
{opth gvar(varname)} gives each unit's first treated period, with {cmd:0} for
units that are never treated. It must be constant within unit. Units whose
{cmd:gvar} falls after the last observed period are never treated within the
window and are pooled with the never-treated group as clean controls.

{dlgtab:Model}

{phang}
{opt weighting(scheme)} chooses the covariance model used to form the optimal
GMM weight. {cmd:pooled_toeplitz} (the default) uses a single stationary
autocovariance sequence shared by all groups (GMM-T). {cmd:cohort_toeplitz}
uses a separate stationary sequence for each cohort (GMM-HT). {cmd:unrestricted}
uses a full, unrestricted within-cohort covariance (GMM-U). Unique
abbreviations are accepted, as are the paper's labels {cmd:t}, {cmd:ht} and
{cmd:u}.

{phang}
{opt nevertreated} requires that at least one never-treated unit exists and
exits with an error otherwise. {opt nonevertreated} excludes never-treated units
from the controls, relying on not-yet-treated cohorts alone, and warns if any are
present. By default, never-treated units are detected automatically and used
when present. The two options may not be combined.

{dlgtab:Iteration}

{phang}
{opt maxiter(#)} sets the maximum number of iterated-GMM reweighting steps. The
iteration can converge linearly at a slow rate; the default of 100 is generous.

{phang}
{opt tolerance(#)} sets the convergence tolerance on the largest absolute change
in any CATT between successive iterations.

{dlgtab:Reporting}

{phang}
{opt level(#)}; see {helpb estimation options##level():[R] Estimation options}.


{marker data}{...}
{title:Data requirements}

{pstd}
The data must have one observation per unit-period. The following conditions
are checked before any estimation, and each exits with an error naming the
offending observations, units, periods or cohorts:

{p 8 10 2}- the outcome is numeric; {cmd:time()} and {cmd:gvar()} are integer-valued;{p_end}
{p 8 10 2}- there are no missing values in {cmd:ivar()}, {cmd:time()} or {cmd:gvar()} (missing outcomes are permitted);{p_end}
{p 8 10 2}- each unit-period combination appears exactly once;{p_end}
{p 8 10 2}- {cmd:gvar()} is constant within unit;{p_end}
{p 8 10 2}- no period is absent from the whole panel (individual units may miss individual periods);{p_end}
{p 8 10 2}- no unit is treated at or before the first observed period;{p_end}
{p 8 10 2}- if never-treated units are coded 0, period 0 does not lie after the first observed period (otherwise 0 would also mean "first treated in period 0");{p_end}
{p 8 10 2}- at least one cohort adopts strictly inside the observed window.{p_end}

{pstd}
Units treated at or before the first observed period have no pre-treatment
period, so none of their effects is identified, and including them silently
would contaminate the aggregate. Restrict the sample before estimating; see the
examples below. Never-treated units may also be coded with any value after the
last observed period; such units are pooled with those coded 0.

{pstd}
Stored results include a weight matrix with one row per group-period and
matrices with one row per cohort-by-time effect. Stata/BE limits matrices to
800 rows, Stata/SE to 11,000 and Stata/MP to 65,534 (see {help limits}); a
design that exceeds the limit is rejected before estimation with an error that
reports the size required.


{marker remarks}{...}
{title:Remarks}

{pstd}
{bf:Unbalanced panels.} Missing outcomes are supported. Group-period means are
taken over the units observed at that period; autocovariances are divided by
the number of residual pairs actually observed at each lag; and the variance of
a group-period mean accounts for the overlap between the unit sets contributing
at each pair of periods. When the panel is balanced all three reduce exactly to
the complete-data formulas. Validity requires that missingness be independent
of the outcome given group and period.

{pstd}
{bf:Partial identification.} A cohort-by-time effect with no clean comparison
is not estimable. It is excluded from {cmd:e(b)} and reported with a missing
estimate and {cmd:identified} = 0 in {cmd:e(catt)}. Where an unidentified cell
carries positive aggregation weight, the corresponding aggregate is missing,
and an identified-subset aggregate renormalised over the estimable cells is
reported alongside it, together with the share of weight it covers.

{pstd}
{bf:Convergence.} Three states are reported separately. {cmd:e(converged)} is 1
only when the largest change in any CATT fell below {opt tolerance()}.
{cmd:e(solve_ok)} records the weaker fact that at least one weighted solve
succeeded. {cmd:e(termination)} gives the reason the iteration stopped:
{cmd:converged}, {cmd:max_iter_reached}, {cmd:singular_weight_matrix},
{cmd:singular_normal_equations} or {cmd:non_finite_update}.

{pstd}
When {cmd:e(solve_ok)} is 0, no weighted step was ever completed and the
reported estimates are the identity-weighted seed (the per-cell average of the
clean comparisons), with sandwich standard errors, rather than a GMM estimate
under the requested weighting. That always raises a warning. A cohort-specific
or unrestricted covariance carries more parameters than a short panel with few
units per cohort can support. When {cmd:e(solve_ok)} is 1 but
{cmd:e(converged)} is 0, the estimates are the last valid iterate; this warns
for the two Toeplitz weightings but not for {cmd:unrestricted}, where slow
convergence is routine.

{pstd}
{bf:Covariates.} When {it:covariates} are supplied, the outcome-regression
adjustment of Section 4.5 of the paper is applied, so that parallel trends need
hold only conditional on the covariates. For each clean comparison, the change
in the control group's outcome is replaced by the treated cohort's mean
covariate vector times the coefficients of a regression, within the control
group, of the outcome change on the covariates and an intercept. Where a
covariate varies within a unit, its value at the unit's earliest observed
period is used, with a warning. Where the regression is not identified (too few
complete-case control units, or collinear covariates), that comparison reverts
to the unconditional form, with a warning reporting how many did.

{pstd}
{bf:Covariates and the specification test.} {cmd:estat jtest} does {bf:not}
adjust for covariates. Its placebo moments are unadjusted group means, so it
tests {bf:unconditional} parallel trends and no anticipation, even after a
covariate-adjusted fit. It is therefore not a test of the conditional parallel
trends assumption under which a covariate-adjusted model is estimated: when
covariates are distributed differently across cohorts, the test can reject
although conditional parallel trends hold, and a non-rejection does not
validate the conditional assumption. Both the estimation output and
{cmd:estat jtest} print a reminder whenever covariates are used.

{pstd}
{bf:Exporting results.} {cmd:e(b)} and {cmd:e(V)} hold the cohort-by-time
effects, so {cmd:esttab}, {cmd:outreg2}, {cmd:lincom} and {cmd:test} work on
them directly. To export the aggregate ATTs as coefficients with standard
errors, use {cmd:estat aggregate, post}, which replaces {cmd:e(b)} and
{cmd:e(V)} with the two aggregates and their covariance and keeps all other
results; see {help staggeredgmm_postestimation##aggregate:estat aggregate}.


{marker examples}{...}
{title:Examples}

{pstd}Simulated panel: 60 units, 33 periods, five cohorts and a never-treated group{p_end}
{phang2}{cmd:. sysuse sgmm_sim_panel, clear}{p_end}
{phang2}{cmd:. staggeredgmm y, ivar(unit_id) time(year) gvar(cohort)}{p_end}

{pstd}Specification test of parallel trends and no anticipation{p_end}
{phang2}{cmd:. estat jtest}{p_end}
{phang2}{cmd:. estat jtest, full}{p_end}

{pstd}Aggregate ATTs, and posting them for {cmd:esttab} or {cmd:outreg2}{p_end}
{phang2}{cmd:. estat aggregate}{p_end}
{phang2}{cmd:. estimates store catt}{p_end}
{phang2}{cmd:. estat aggregate, post}{p_end}
{phang2}{cmd:. estimates store att}{p_end}

{pstd}Cohort-specific Toeplitz weighting{p_end}
{phang2}{cmd:. staggeredgmm y, ivar(unit_id) time(year) gvar(cohort) weighting(cohort_toeplitz)}{p_end}

{pstd}Adjusting for baseline covariates (the J-test remains unconditional){p_end}
{phang2}{cmd:. staggeredgmm y x1 i.x2, ivar(unit_id) time(year) gvar(cohort)}{p_end}
{phang2}{cmd:. estat jtest}{p_end}

{pstd}No never-treated group: not-yet-treated controls only{p_end}
{phang2}{cmd:. staggeredgmm y if cohort != 0, ivar(unit_id) time(year) gvar(cohort)}{p_end}

{pstd}
Without never-treated units, the effects of the last cohort, and of the other
cohorts once every later cohort is treated, have no clean comparison. Only 30
of the 90 effects are identified; the full aggregates are missing, and the
identified-subset aggregates are reported with the share of weight they cover.

{pstd}The full table of cohort-by-time effects, including any that are not identified{p_end}
{phang2}{cmd:. matrix list e(catt)}{p_end}

{pstd}Always-treated units are rejected{p_end}
{phang2}{cmd:. replace cohort = 1 if cohort == 10}{p_end}
{phang2}{cmd:. staggeredgmm y, ivar(unit_id) time(year) gvar(cohort)}{p_end}

{pstd}
Cohort 1 is treated in the first observed period, so it has no pre-treatment
period and none of its effects is identified. The command exits with an error
naming the cohort; drop such units, for example with {cmd:if cohort != 1}.


{marker results}{...}
{title:Stored results}

{pstd}
{cmd:staggeredgmm} stores the following in {cmd:e()}:

{synoptset 24 tabbed}{...}
{p2col 5 24 28 2: Scalars}{p_end}
{synopt:{cmd:e(N)}}number of observations with a non-missing outcome{p_end}
{synopt:{cmd:e(N_rows)}}number of unit-period rows in the estimation sample, including rows with a missing outcome{p_end}
{synopt:{cmd:e(N_units)}}number of units{p_end}
{synopt:{cmd:e(T_periods)}}number of periods{p_end}
{synopt:{cmd:e(T_min)}, {cmd:e(T_max)}}first and last period{p_end}
{synopt:{cmd:e(n_cohorts)}}number of treated cohorts{p_end}
{synopt:{cmd:e(N_beta)}}number of cohort-by-time effects{p_end}
{synopt:{cmd:e(N_identified)}}number of identified cohort-by-time effects{p_end}
{synopt:{cmd:e(N_2x2)}}number of clean comparisons in the moment system{p_end}
{synopt:{cmd:e(n_dropped_2x2)}}comparisons dropped for referencing an unobserved group-period{p_end}
{synopt:{cmd:e(moment_rank)}}rank of the moment space{p_end}
{synopt:{cmd:e(n_kept)}}eigen-directions of the weight retained at the last weighted solve{p_end}
{synopt:{cmd:e(n_zero_lag)}}autocovariance lags with no observed residual pair{p_end}
{synopt:{cmd:e(converged)}}1 if the iteration converged{p_end}
{synopt:{cmd:e(solve_ok)}}1 if at least one weighted solve succeeded{p_end}
{synopt:{cmd:e(n_iter)}}number of iterations{p_end}
{synopt:{cmd:e(max_change)}}largest change in any CATT at the last iteration{p_end}
{synopt:{cmd:e(tol)}, {cmd:e(max_iter)}}convergence settings{p_end}
{synopt:{cmd:e(has_never)}}1 if never-treated units are used as controls{p_end}
{synopt:{cmd:e(n_never_zero)}}units with {cmd:gvar} = 0{p_end}
{synopt:{cmd:e(n_never_late)}}units with {cmd:gvar} after the last period{p_end}
{synopt:{cmd:e(cov_n_fallback)}}comparisons whose covariate adjustment fell back to the unconditional form{p_end}
{synopt:{cmd:e(balanced)}}1 if no outcome is missing{p_end}
{synopt:{cmd:e(n_missing_cells)}}number of missing unit-period outcomes{p_end}
{synopt:{cmd:e(att_cw)}, {cmd:e(se_cw)}}treated-observation ATT and its standard error{p_end}
{synopt:{cmd:e(att_cw_id)}, {cmd:e(se_cw_id)}}identified-subset version{p_end}
{synopt:{cmd:e(share_cw)}}share of treated-observation weight on identified cells{p_end}
{synopt:{cmd:e(att_ew)}, {cmd:e(se_ew)}}cohort-equal ATT and its standard error{p_end}
{synopt:{cmd:e(att_ew_id)}, {cmd:e(se_ew_id)}}identified-subset version{p_end}
{synopt:{cmd:e(share_ew)}}share of cohort-equal weight on identified cells{p_end}
{synopt:{cmd:e(level)}}confidence level{p_end}

{p2col 5 24 28 2: Macros}{p_end}
{synopt:{cmd:e(cmd)}}{cmd:staggeredgmm}{p_end}
{synopt:{cmd:e(subcmd)}}{cmd:aggregate} after {cmd:estat aggregate, post}{p_end}
{synopt:{cmd:e(cmdline)}}command as typed{p_end}
{synopt:{cmd:e(depvar)}}name of dependent variable{p_end}
{synopt:{cmd:e(ivar)}, {cmd:e(timevar)}, {cmd:e(gvar)}}unit, time and cohort variables{p_end}
{synopt:{cmd:e(weighting)}}weighting scheme{p_end}
{synopt:{cmd:e(termination)}}reason the iteration stopped{p_end}
{synopt:{cmd:e(covariates)}}covariates as specified{p_end}
{synopt:{cmd:e(xnames)}}covariate terms after factor-variable expansion{p_end}
{synopt:{cmd:e(title)}}title in estimation output{p_end}
{synopt:{cmd:e(properties)}}{cmd:b V}{p_end}
{synopt:{cmd:e(estat_cmd)}}program used to implement {cmd:estat}{p_end}

{p2col 5 24 28 2: Matrices}{p_end}
{synopt:{cmd:e(b)}}identified cohort-by-time effects; {cmd:_b[g}{it:G}{cmd::t}{it:T}{cmd:]} is the effect for cohort {it:G} at period {it:T}{p_end}
{synopt:{cmd:e(V)}}variance-covariance matrix of {cmd:e(b)}{p_end}
{synopt:{cmd:e(catt)}}all cohort-by-time cells: {cmd:g}, {cmd:t}, {cmd:event_time}, {cmd:estimate}, {cmd:std_error}, {cmd:identified}, {cmd:n_comparisons}{p_end}
{synopt:{cmd:e(w_cw)}, {cmd:e(w_ew)}}aggregation weights over all cells{p_end}
{synopt:{cmd:e(N_g)}}units per group (never-treated pool first, if used){p_end}
{synopt:{cmd:e(cohorts)}}treated cohorts{p_end}
{synopt:{cmd:e(sgmm_W)}, {cmd:e(sgmm_ymeans)}, {cmd:e(sgmm_Ngt)}, {cmd:e(sgmm_groups)}}internal objects used by {cmd:estat jtest}{p_end}

{p2col 5 24 28 2: Functions}{p_end}
{synopt:{cmd:e(sample)}}marks estimation sample; rows with a missing outcome are included, because they define the panel{p_end}


{marker references}{...}
{title:References}

{phang}
Arora, P. and R. Bijani. 2026. Estimating treatment effects under staggered
timing and non-spherical errors. Available at SSRN,
{browse "https://doi.org/10.2139/ssrn.6558759":doi:10.2139/ssrn.6558759}.

{phang}
Callaway, B. and P. H. C. Sant'Anna. 2021. Difference-in-differences with
multiple time periods. {it:Journal of Econometrics} 225(2): 200-230.


{marker authors}{...}
{title:Authors}

{pstd}
Parush Arora, Ashoka University, parush.arora@ashoka.edu.in{break}
Rishabh Bijani, Ashoka University, rishabhbijani@gmail.com

{pstd}
Please cite the paper above and the software:
Arora, P. and R. Bijani. 2026. staggeredgmm: Generalised Method of Moments Estimation of Treatment Effects under Staggered Adoption.
Stata package version 0.1.0.
{browse "https://github.com/RishabhBijani/staggeredGMM_STATA"}

{pstd}
Released under the MIT licence.
{p_end}
