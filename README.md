# staggeredgmm (Stata)

Generalised Method of Moments Estimation of Treatment Effects under Staggered Adoption.

`staggeredgmm` estimates cohort-by-time average treatment effects (CATTs) under
staggered treatment adoption by the generalized method of moments, following
Arora and Bijani (2026). Each CATT is identified by every clean 2x2
difference-in-differences comparison available for it, against never-treated
units and against cohorts that are not yet treated. The clean comparisons are
stacked into an over-identified moment system. They are combined by iterated
optimal GMM, with the weight built from a model of the serial correlation of
the errors.

- [Installation](#installation)
- [Quick start](#quick-start)
- [Syntax](#syntax)
- [Options](#options)
- [What is estimated](#what-is-estimated)
- [Specification test: estat jtest](#specification-test-estat-jtest)
- [Exporting results: esttab and outreg2](#exporting-results-esttab-and-outreg2)
- [Stored results](#stored-results)
- [Examples](#examples)
- [Data requirements](#data-requirements)
- [Unbalanced panels, partial identification and convergence](#unbalanced-panels-partial-identification-and-convergence)
- [Files](#files)
- [Citation](#citation)
- [Authors and licence](#authors-and-licence)

## Installation

From GitHub:

```stata
net install staggeredgmm, from("https://raw.githubusercontent.com/RishabhBijani/staggeredGMM_STATA/main/") replace
```

Once the package is on SSC:

```stata
ssc install staggeredgmm
```

Requires Stata 16 or later. Installation also places the example dataset
`sgmm_sim_panel` where `sysuse` can find it.

## Quick start

```stata
sysuse sgmm_sim_panel, clear
staggeredgmm y, ivar(unit_id) time(year) gvar(cohort)
estat jtest
```

## Syntax

```
staggeredgmm depvar [covariates] [if] [in], ivar(varname) time(varname) gvar(varname) [options]

estat jtest [, full window(#)]

estat aggregate [, post level(#)]
```

`covariates` are baseline (pre-treatment, time-invariant) covariates. They may
contain factor variables, as in `csdid`.

## Options

| Option | Description |
|---|---|
| `ivar(varname)` | unit identifier, numeric or string (required) |
| `time(varname)` | calendar period, integer-valued (required) |
| `gvar(varname)` | first treated period of each unit, `0` if never treated (required) |
| `weighting(scheme)` | `pooled_toeplitz` (default), `cohort_toeplitz` or `unrestricted`; unique abbreviations and the paper's labels `t`, `ht`, `u` are accepted |
| `nevertreated` | require a never-treated control group, and exit with an error if there is none |
| `nonevertreated` | do not use never-treated units as controls; rely on not-yet-treated cohorts alone |
| `maxiter(#)` | maximum number of iterated-GMM reweighting steps; default `100` |
| `tolerance(#)` | convergence tolerance on the largest change in any CATT; default `1e-6` |
| `level(#)` | confidence level; default `95` |

By default, never-treated units are detected automatically and used as
controls when present. Units whose `gvar` falls after the last observed period
are never treated within the window, so they are pooled with the never-treated
group.

### Weighting schemes

| Scheme | Paper label | Covariance model for the optimal weight |
|---|---|---|
| `pooled_toeplitz` | GMM-T | one stationary autocovariance sequence shared by all groups |
| `cohort_toeplitz` | GMM-HT | a separate stationary sequence for each cohort |
| `unrestricted` | GMM-U | a full, unrestricted within-cohort covariance |

`pooled_toeplitz` has the fewest parameters and is the safest choice in short
panels or with few units per cohort. `unrestricted` needs many units per cohort
relative to the number of periods.

## What is estimated

**Cohort-by-time effects.** For each treated cohort $g$ and period $t \ge g$,
the CATT is the average effect at $t$ on units first treated at $g$. Estimates
are posted in `e(b)` with names `g<G>:t<T>`, so `_b[g10:t12]` is the effect for
cohort 10 in period 12. `lincom`, `test` and `nlcom` work as usual.

**Aggregate ATTs.** Two aggregates are reported, with standard errors:

- The *treated-observation* ATT (`e(att_cw)`) weights every treated unit-period
  equally, so cohorts observed for longer after treatment carry more weight.
  This is the aggregation used in the paper and coincides with the "simple"
  aggregation of Callaway and Sant'Anna (2021).
- The *cohort-equal* ATT (`e(att_ew)`) gives every treated cohort the same total
  weight, spread equally over its own post-treatment periods.

**Covariates.** When covariates are supplied, the outcome-regression adjustment
of Section 4.5 of the paper is applied, so parallel trends need only hold
conditional on the covariates. For each clean comparison, the change in the
control group's outcome is replaced by the treated cohort's mean covariate
vector times the coefficients of a within-control regression of the outcome
change on the covariates and an intercept. Where that regression is not
identified (too few complete-case control units, or collinear covariates), the
comparison reverts to the unconditional form, with a warning reporting how many
did (`e(cov_n_fallback)`). Where a covariate varies within a unit, its value in
the unit's earliest observed period is used, with a warning.

> **Covariates and the specification test.** `estat jtest` does **not** adjust
> for covariates. Even after a covariate-adjusted fit, it tests
> **unconditional** parallel trends and no anticipation, not the
> covariate-conditional parallel trends assumption under which the effects are
> estimated. If covariates are distributed differently across cohorts, the test
> can reject although conditional parallel trends hold, and a non-rejection
> does not validate the conditional assumption. The estimation output and
> `estat jtest` print a reminder whenever covariates are used.

## Specification test: estat jtest

`estat jtest` is a serial-correlation robust joint test of **unconditional**
parallel trends and no anticipation (Section 4.4 of the paper). Its placebo
moments are never covariate-adjusted; see the box under
[Covariates](#what-is-estimated). It fixes a base pre-period
$t_0 = g_{\min} - 1$, the last period before any cohort adopts. For each
treated cohort $g$ and pre-period $t < g$ it forms the placebo moment

$$D(g,t) = \bigl(\bar Y_{g,t} - \bar Y_{g,t_0}\bigr) - \bigl(\bar Y_{\infty,t} - \bar Y_{\infty,t_0}\bigr),$$

where $\infty$ denotes the never-treated group. The statistic
$J = p'\,\Omega_p^{-1}\,p$ uses the covariance $\Omega_p$ of the placebo
moments implied by the fitted serial-correlation model. It is asymptotically
$\chi^2$ with as many degrees of freedom as placebo moments. No degrees of
freedom are subtracted for estimated parameters, because the placebo moments do
not involve the treatment effects.

| Option | Description |
|---|---|
| `window(#)` | use only the placebo moments within `#` periods before each cohort's adoption; default `3` |
| `full` | use every pre-treatment placebo moment |

The default local test holds its nominal size across sample sizes. A
differential trend or an anticipation effect is already visible just before
treatment, so restricting to the window costs little power. `full` is well
sized only when the number of units is comfortably large. The test requires a
never-treated control group. Results are returned in `r(J)`, `r(df)`, `r(p)`,
`r(base_period)`, `r(n_moments)`, `r(n_dropped)`, `r(window)`, `r(type)`,
`r(scope)` (always `unconditional`), `r(weighting)` and the matrix `r(moments)`.

## Exporting results: esttab and outreg2

`e(b)` and `e(V)` hold the cohort-by-time effects, so `esttab`, `outreg2`,
`lincom` and `test` work on them directly. The aggregate ATTs are exported with
`estat aggregate, post`. It replaces `e(b)` and `e(V)` with the two aggregates,
`att_cw` (treated-observation) and `att_ew` (cohort-equal), and their exact
joint covariance. It keeps every other result, so `estat jtest` still works.

```stata
staggeredgmm y, ivar(unit_id) time(year) gvar(cohort)
estimates store catt               // cohort-by-time effects
estat aggregate, post
estimates store att                // aggregate ATTs
esttab att, se                     // or: outreg2 using results, replace
esttab catt, se
```

Without `post`, `estat aggregate` displays the same table and returns it in
`r(b)` and `r(V)` without changing `e()`. When some cohort-by-time effects are
not identified, both aggregates are taken over the identified effects only, and
the output reports the share of weight they cover (`r(share_cw)`,
`r(share_ew)`).

## Stored results

`staggeredgmm` is an e-class command. The main results are listed below; `help
staggeredgmm` has the full list.

| Result | Content |
|---|---|
| `e(b)`, `e(V)` | identified CATTs and their variance-covariance matrix |
| `e(catt)` | every cohort-by-time cell: `g`, `t`, `event_time`, `estimate`, `std_error`, `identified`, `n_comparisons` |
| `e(att_cw)`, `e(se_cw)` | treated-observation ATT and standard error |
| `e(att_ew)`, `e(se_ew)` | cohort-equal ATT and standard error |
| `e(att_cw_id)`, `e(att_ew_id)` (and `se_`, `share_`) | identified-subset aggregates and the share of weight they cover |
| `e(w_cw)`, `e(w_ew)` | aggregation weights over all cells |
| `e(N)` | observations with a non-missing outcome |
| `e(N_rows)`, `e(N_units)`, `e(T_periods)`, `e(n_cohorts)` | unit-period rows (including missing outcomes), units, periods, treated cohorts |
| `e(N_beta)`, `e(N_identified)`, `e(N_2x2)`, `e(moment_rank)` | design dimensions |
| `e(converged)`, `e(solve_ok)`, `e(n_iter)`, `e(max_change)`, `e(termination)` | convergence diagnostics |
| `e(has_never)`, `e(n_never_zero)`, `e(n_never_late)` | control-group information |
| `e(balanced)`, `e(n_missing_cells)`, `e(n_dropped_2x2)` | missing-data information |

## Examples

```stata
* Simulated panel: 60 units, 33 periods, five cohorts and a never-treated group
sysuse sgmm_sim_panel, clear
staggeredgmm y, ivar(unit_id) time(year) gvar(cohort)

* Specification tests
estat jtest
estat jtest, window(5)
estat jtest, full

* Aggregate ATTs, posted for export
estat aggregate
estat aggregate, post

* Cohort-specific Toeplitz weighting
staggeredgmm y, ivar(unit_id) time(year) gvar(cohort) weighting(cohort_toeplitz)

* Baseline covariates, including a factor variable (the J-test stays unconditional)
staggeredgmm y x1 i.x2, ivar(unit_id) time(year) gvar(cohort)
estat jtest

* Not-yet-treated controls only
staggeredgmm y, ivar(unit_id) time(year) gvar(cohort) nonevertreated

* No never-treated group at all: 30 of the 90 effects are identified
staggeredgmm y if cohort != 0, ivar(unit_id) time(year) gvar(cohort)

* Every cohort-by-time cell, and inference on a combination of effects
staggeredgmm y, ivar(unit_id) time(year) gvar(cohort)
matrix list e(catt)
lincom _b[g10:t12] - _b[g10:t10]

* Always-treated units are rejected with an error naming the cohort
replace cohort = 1 if cohort == 10
staggeredgmm y, ivar(unit_id) time(year) gvar(cohort)
```

`unrestricted` weighting needs many units per cohort relative to the number of
periods. With 10 units per cohort and 33 periods, `sgmm_sim_panel` is too small
for it: the command warns that no weighted step could be completed and reports
the identity-weighted starting estimate.

## Data requirements

The data must contain one observation per unit-period. Each of the following is
checked before estimation, and a failure exits with an error naming the
offending observations, units, periods or cohorts:

- the outcome is numeric, and `time()` and `gvar()` are integer-valued;
- `ivar()`, `time()` and `gvar()` have no missing values (missing outcomes are allowed);
- each unit-period combination appears exactly once;
- `gvar()` is constant within unit;
- no period is absent from the whole panel (individual units may miss individual periods);
- no unit is treated at or before the first observed period;
- if never-treated units are coded 0, period 0 does not lie after the first observed period (otherwise 0 would also mean "first treated in period 0"; code never-treated units with a value after the last period instead);
- at least one cohort adopts strictly inside the observed window.

Units treated at or before the first observed period have no pre-treatment
period, so none of their effects is identified. Including them silently would
contaminate the aggregates, so they must be removed before estimation.

Stored results include a weight matrix with one row per group-period. Stata/BE
limits matrices to 800 rows, Stata/SE to 11,000 and Stata/MP to 65,534. A
design that exceeds the limit is rejected before estimation, with an error
reporting the size required.

## Unbalanced panels, partial identification and convergence

**Unbalanced panels.** Missing outcomes are supported.
- Group-period means are taken over the units observed in that period.
- Autocovariances are divided by the number of residual pairs actually observed at each lag.
- The variance of a group-period mean accounts for the overlap between the unit sets observed in each pair of periods.

On a balanced panel all three reduce exactly to the complete-data formulas.
Validity requires missingness to be independent of the outcome given group and
period. A clean comparison that references a group-period with no observed unit
is dropped from the moment system, with a warning (`e(n_dropped_2x2)`).

**Partial identification.** A CATT with no clean comparison cannot be
estimated. This happens, for example, to the last cohort's effects when there is
no never-treated group. Such a CATT is excluded from `e(b)` and appears in
`e(catt)` with a missing estimate and `identified = 0`. If an unidentified cell
carries positive aggregation weight, the corresponding aggregate is missing.
An identified-subset aggregate, renormalised over the estimable cells, is then
reported alongside it with the share of weight it covers.

**Convergence.** Three states are reported separately:
- `e(converged)` is 1 only when the largest change in any CATT fell below `tolerance()`.
- `e(solve_ok)` records the weaker fact that at least one weighted solve succeeded.
- `e(termination)` gives the reason the iteration stopped: `converged`, `max_iter_reached`, `singular_weight_matrix`, `singular_normal_equations` or `non_finite_update`.

When `e(solve_ok)` is 0, the reported estimates are the identity-weighted seed
(the per-cell average of the clean comparisons) with sandwich standard errors.
They are not a GMM estimate under the requested weighting, and a warning is
always issued.

## Files

| File | Purpose |
|---|---|
| `staggeredgmm.ado` | estimation command (ado and Mata) |
| `staggeredgmm_estat.ado` | `estat jtest` and `estat aggregate` |
| `staggeredgmm.sthlp`, `staggeredgmm_postestimation.sthlp` | help files |
| `staggeredgmm.pkg`, `stata.toc` | `net install` metadata |
| `sgmm_sim_panel.dta` | simulated panel used in the examples |
| `test/staggeredgmm_test.do` | test script: exact recovery, invariances, error handling, postestimation |

## Citation

Please cite the paper and the software:

> Arora, P. and Bijani, R. (2026). "Estimating Treatment Effects under
> Staggered Timing and Non-Spherical Errors." Available at SSRN.
> [doi:10.2139/ssrn.6558759](https://doi.org/10.2139/ssrn.6558759)

> Arora, P. and Bijani, R. (2026). staggeredgmm: Generalised Method of Moments Estimation of Treatment Effects under Staggered Adoption.
> Stata package version 0.1.0.
> https://github.com/RishabhBijani/staggeredGMM_STATA

```bibtex
@unpublished{arora_bijani_2026,
  title  = {Estimating Treatment Effects under Staggered Timing and Non-Spherical Errors},
  author = {Arora, Parush and Bijani, Rishabh},
  year   = {2026},
  note   = {Available at SSRN},
  doi    = {10.2139/ssrn.6558759}
}

@manual{staggeredgmm_stata,
  title  = {staggeredgmm: Generalised Method of Moments Estimation of Treatment Effects under Staggered Adoption},
  author = {Arora, Parush and Bijani, Rishabh},
  year   = {2026},
  note   = {Stata package version 0.1.0},
  url    = {https://github.com/RishabhBijani/staggeredGMM_STATA}
}
```

## Authors and licence

Parush Arora, Ashoka University, parush.arora@ashoka.edu.in
Rishabh Bijani, Ashoka University, rishabhbijani@gmail.com

Released under the MIT licence.
