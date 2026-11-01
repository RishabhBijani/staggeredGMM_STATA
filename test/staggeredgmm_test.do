* Test script for staggeredgmm.
*
* Run from this test/ directory:   do staggeredgmm_test.do
*
* Uses the ado files in the parent directory, not any installed copy, and the
* example dataset sgmm_sim_panel.dta from the parent directory. Every check
* is self-contained: exact recovery on noiseless data, the J-test restriction
* counts of the simulated design, invariances, error handling, stored results
* and postestimation. Runs with varabbrev off, as required for SSC.

version 16
clear all
set more off
set varabbrev off
set linesize 120

local here "`c(pwd)'"
adopath ++ "`here'/.."
discard
local sim "`here'/../sgmm_sim_panel.dta"

global sgmm_fail 0
global sgmm_pass 0
global sgmm_failed ""

program sgmm_record
	args case what ok detail
	if `ok' {
		global sgmm_pass = $sgmm_pass + 1
	}
	else {
		global sgmm_fail = $sgmm_fail + 1
		global sgmm_failed "$sgmm_failed `case':`what'"
		di as err "  FAIL `case' `what' `detail'"
	}
end

program sgmm_expect_error
	gettoken name 0 : 0
	gettoken rc cmd : 0
	capture `cmd'
	local got = _rc
	sgmm_record E `name' `=`got' == `rc'' "expected rc `rc', got `got'"
end

* A noiseless panel: cohorts 3 and 4 (two units each), two never-treated units,
* four periods, effects (3,3) = 1, (3,4) = 2, (4,4) = 3.
program sgmm_noiseless
	clear
	quietly set obs 24
	quietly gen long unit_id = ceil(_n / 4)
	quietly gen int year = _n - (unit_id - 1) * 4
	quietly gen int cohort = cond(unit_id <= 2, 3, cond(unit_id <= 4, 4, 0))
	quietly gen double tau = 0
	quietly replace tau = 1 if cohort == 3 & year == 3
	quietly replace tau = 2 if cohort == 3 & year == 4
	quietly replace tau = 3 if cohort == 4 & year == 4
	quietly gen double y = unit_id / 10 + year^2 / 7 + tau
	drop tau
end

di _n as txt "{hline 70}" _n "Estimation" _n "{hline 70}"

* exact recovery on a noiseless panel, under every weighting
foreach w in pooled_toeplitz cohort_toeplitz unrestricted {
	sgmm_noiseless
	quietly staggeredgmm y, ivar(unit_id) time(year) gvar(cohort) weighting(`w')
	tempname C
	matrix `C' = e(catt)
	local ok = abs(`C'[1,4] - 1) < 1e-8 & abs(`C'[2,4] - 2) < 1e-8 & abs(`C'[3,4] - 3) < 1e-8
	local ok = `ok' & `C'[1,1] == 3 & `C'[2,1] == 3 & `C'[3,1] == 4 & `C'[3,2] == 4 & `C'[2,3] == 1
	sgmm_record A exact_recovery_`w' `ok' ""
}

* J-test restriction counts on the simulated panel (70 full, 14 local, base 9)
quietly use "`sim'", clear
quietly staggeredgmm y, ivar(unit_id) time(year) gvar(cohort)
quietly estat jtest, full
sgmm_record A jtest_full_70 `=r(n_moments) == 70 & r(base_period) == 9' ""
quietly estat jtest
sgmm_record A jtest_local_14 `=r(n_moments) == 14' ""
tempname m
matrix `m' = r(moments)
local ok 1
forvalues i = 1/`=rowsof(`m')' {
	if `m'[`i', 2] >= `m'[`i', 1] | `m'[`i', 2] == 9 local ok 0
}
sgmm_record A jtest_placebo_pretreatment `ok' ""
sgmm_record A jtest_p_formula `=reldif(r(p), chi2tail(r(df), r(J))) < 1e-12' ""
sgmm_record A jtest_scope `="`r(scope)'" == "unconditional"' ""

* _b[] naming and agreement with e(catt)
tempname Cc agg
matrix `Cc' = e(catt)
matrix `agg' = e(w_cw) * `Cc'[1..., 4]
sgmm_record A b_naming `=abs(_b[g10:t10] - `Cc'[1, 4]) < 1e-12' ""
sgmm_record A aggregate_identity `=reldif(e(att_cw), `agg'[1, 1]) < 1e-12' ""

* e(N) counts non-missing outcomes, e(N_rows) every row
preserve
quietly replace y = . in 5/9
quietly staggeredgmm y, ivar(unit_id) time(year) gvar(cohort)
sgmm_record A N_nonmissing `=e(N) == 1975 & e(N_rows) == 1980' "N=`=e(N)' N_rows=`=e(N_rows)'"
restore

di _n as txt "{hline 70}" _n "Invariances" _n "{hline 70}"

quietly use "`sim'", clear
quietly staggeredgmm y, ivar(unit_id) time(year) gvar(cohort)
tempname base_b
matrix `base_b' = e(b)

preserve
set seed 101
tempvar u
quietly gen double `u' = runiform()
sort `u'
quietly staggeredgmm y, ivar(unit_id) time(year) gvar(cohort)
sgmm_record B row_order `=mreldif(e(b), `base_b') < 1e-10' ""
restore

preserve
quietly gen str7 sid = "unit" + string(unit_id, "%03.0f")
quietly staggeredgmm y, ivar(sid) time(year) gvar(cohort)
sgmm_record B string_id `=mreldif(e(b), `base_b') < 1e-10' ""
restore

preserve
quietly replace year = year + 1000
quietly replace cohort = cohort + 1000 if cohort != 0
quietly staggeredgmm y, ivar(unit_id) time(year) gvar(cohort)
tempname bs
matrix `bs' = e(b)
sgmm_record B time_shift `=colsof(`bs') == colsof(`base_b') & mreldif(`bs', `base_b') < 1e-10' ""
restore

preserve
quietly gen junk = "x"
quietly gen junk2 = .
quietly staggeredgmm y, ivar(unit_id) time(year) gvar(cohort)
sgmm_record B extra_vars `=mreldif(e(b), `base_b') < 1e-10' ""
restore

preserve
quietly replace cohort = 99 if cohort == 0
quietly staggeredgmm y, ivar(unit_id) time(year) gvar(cohort)
sgmm_record B late_cohort_as_never `=mreldif(e(b), `base_b') < 1e-10 & e(n_never_late) == 10' ""
restore

* time running through 0: never-treated coded after the last period is accepted
preserve
quietly replace year = year - 5
quietly replace cohort = cohort - 5 if cohort != 0
quietly replace cohort = 99 if cohort == 0
quietly staggeredgmm y, ivar(unit_id) time(year) gvar(cohort)
tempname bz
matrix `bz' = e(b)
sgmm_record B time_through_zero `=mreldif(`bz', `base_b') < 1e-10' ""
restore

preserve
quietly staggeredgmm y if unit_id != 7, ivar(unit_id) time(year) gvar(cohort)
tempname bif
matrix `bif' = e(b)
quietly drop if unit_id == 7
quietly staggeredgmm y, ivar(unit_id) time(year) gvar(cohort)
sgmm_record B if_equals_drop `=mreldif(e(b), `bif') < 1e-12' ""
restore

di _n as txt "{hline 70}" _n "Options" _n "{hline 70}"

quietly use "`sim'", clear
quietly staggeredgmm y, ivar(unit_id) time(year) gvar(cohort) weighting(cohort_toeplitz)
tempname bht
matrix `bht' = e(b)
quietly staggeredgmm y, ivar(unit_id) time(year) gvar(cohort) weighting(ht)
local ok = mreldif(e(b), `bht') < 1e-12 & "`e(weighting)'" == "cohort_toeplitz"
quietly staggeredgmm y, ivar(unit_id) time(year) gvar(cohort) weighting(coh)
local ok = `ok' & mreldif(e(b), `bht') < 1e-12
sgmm_record C weighting_abbrev `ok' ""
quietly staggeredgmm y, ivar(unit_id) time(year) gvar(cohort) nevertreated
sgmm_record C nevertreated_explicit `=mreldif(e(b), `base_b') < 1e-12' ""
quietly staggeredgmm y, ivar(unit_id) time(year) gvar(cohort) nonever
sgmm_record C nonevertreated_abbrev `=e(has_never) == 0 & e(n_never_zero) > 0' ""
quietly staggeredgmm y, ivar(unit_id) time(year) gvar(cohort) level(90)
local ok = e(level) == 90
capture staggeredgmm, level(80)
sgmm_record C level_option `=`ok' & _rc == 0' ""

preserve
quietly gen byte xc = (x1 > 0) + (x2 > 0)
quietly tab xc, gen(xcd)
quietly staggeredgmm y i.xc, ivar(unit_id) time(year) gvar(cohort)
tempname bfv
matrix `bfv' = e(b)
quietly staggeredgmm y xcd2 xcd3, ivar(unit_id) time(year) gvar(cohort)
sgmm_record C factor_covariates `=mreldif(e(b), `bfv') < 1e-10' ""
quietly gen byte one = 1
capture staggeredgmm y i.one, ivar(unit_id) time(year) gvar(cohort)
sgmm_record C covariates_all_base `=_rc == 459' ""
restore

di _n as txt "{hline 70}" _n "estat aggregate" _n "{hline 70}"

quietly use "`sim'", clear
quietly staggeredgmm y, ivar(unit_id) time(year) gvar(cohort)
local att_cw = e(att_cw)
local att_ew = e(att_ew)
local se_cw  = e(se_cw)
local se_ew  = e(se_ew)
local nobs   = e(N)
quietly estat aggregate
tempname rb rV
matrix `rb' = r(b)
matrix `rV' = r(V)
local ok = reldif(`rb'[1,1], `att_cw') < 1e-12 & reldif(`rb'[1,2], `att_ew') < 1e-12
local ok = `ok' & reldif(sqrt(`rV'[1,1]), `se_cw') < 1e-10 & reldif(sqrt(`rV'[2,2]), `se_ew') < 1e-10
local ok = `ok' & abs(r(share_cw) - 1) < 1e-12 & abs(r(share_ew) - 1) < 1e-12 & "`e(subcmd)'" == ""
sgmm_record D aggregate_values `ok' ""
estimates store catt
quietly estat aggregate, post
local names : colnames e(b)
local ok = "`names'" == "att_cw att_ew" & "`e(subcmd)'" == "aggregate" & "`e(cmd)'" == "staggeredgmm"
local ok = `ok' & reldif(_b[att_cw], `att_cw') < 1e-12 & reldif(_se[att_ew], `se_ew') < 1e-10
local ok = `ok' & e(N) == `nobs' & e(N_units) == 60
sgmm_record D aggregate_post `ok' ""
capture estat jtest
sgmm_record D jtest_after_post `=_rc == 0' ""
capture staggeredgmm
sgmm_record D replay_after_post `=_rc == 0' ""
capture estat aggregate, post
sgmm_record D post_twice_refused `=_rc == 498' ""
quietly estimates restore catt
sgmm_record D restore_catt `=colsof(e(b)) == 90' ""

* partial identification: aggregates over the identified effects only
quietly staggeredgmm y if cohort != 0, ivar(unit_id) time(year) gvar(cohort)
local ok = e(att_cw) >= . & e(N_identified) == 30
quietly estat aggregate
tempname rb2
matrix `rb2' = r(b)
local ok = `ok' & reldif(`rb2'[1,1], e(att_cw_id)) < 1e-12 & reldif(r(share_cw), e(share_cw)) < 1e-12
sgmm_record D aggregate_identified_subset `ok' ""

di _n as txt "{hline 70}" _n "Error handling" _n "{hline 70}"

local s staggeredgmm y, ivar(unit_id) time(year) gvar(cohort)
sgmm_noiseless

preserve
quietly replace year = year + 0.5 in 3
sgmm_expect_error noninteger_time 459 `s'
restore
preserve
quietly replace year = . in 3
sgmm_expect_error missing_time 459 `s'
restore
preserve
quietly replace cohort = . in 5
sgmm_expect_error missing_cohort 459 `s'
restore
preserve
quietly replace y = . in 4
capture `s'
sgmm_record E missing_outcome_ok `=_rc == 0' ""
restore
preserve
quietly expand 2 in 1
sgmm_expect_error duplicate_rows 459 `s'
restore
preserve
quietly replace cohort = 4 if unit_id == 1 & year == 2
sgmm_expect_error cohort_varies 459 `s'
restore
preserve
quietly drop if year == 2
sgmm_expect_error period_gap 459 `s'
restore
preserve
quietly replace cohort = 1 if cohort == 3
sgmm_expect_error always_treated 459 `s'
restore
preserve
quietly drop if cohort == 0
sgmm_expect_error nevertreated_none 459 `s' nevertreated
restore
preserve
quietly keep if year == 1
sgmm_expect_error single_period 459 `s'
restore
preserve
quietly replace year = year - 2
quietly replace cohort = cohort - 2 if cohort != 0
sgmm_expect_error period_zero_in_window 459 `s'
restore
sgmm_expect_error maxiter_zero 198 `s' maxiter(0)
sgmm_expect_error maxiter_fraction 198 `s' maxiter(2.5)
sgmm_expect_error tol_negative 198 `s' tolerance(-1)
sgmm_expect_error bad_weighting 198 `s' weighting(nonsense)
sgmm_expect_error both_never_options 198 `s' nevertreated nonevertreated
sgmm_expect_error unknown_option 198 `s' bogus
sgmm_expect_error covariate_is_key 198 staggeredgmm y year, ivar(unit_id) time(year) gvar(cohort)
sgmm_expect_error same_var_twice 198 staggeredgmm y, ivar(unit_id) time(unit_id) gvar(cohort)
preserve
quietly tostring y, replace force
sgmm_expect_error string_outcome 109 `s'
restore

quietly use "`sim'", clear
preserve
quietly replace x1 = . if cohort == 13
sgmm_expect_error covariate_missing_in_cohort 459 staggeredgmm y x1 x2, ivar(unit_id) time(year) gvar(cohort)
restore

* estat errors
preserve
quietly drop if cohort == 0
quietly staggeredgmm y, ivar(unit_id) time(year) gvar(cohort)
capture estat jtest
sgmm_record E jtest_needs_never `=_rc == 459' ""
restore
quietly staggeredgmm y, ivar(unit_id) time(year) gvar(cohort)
capture estat jtest, window(0)
sgmm_record E jtest_window_zero `=_rc == 198' ""
capture estat jtest, window(2.5)
sgmm_record E jtest_window_fraction `=_rc == 198' ""

* replay
capture staggeredgmm
sgmm_record E replay `=_rc == 0' ""

di _n as txt "{hline 70}"
di as txt "Checks passed: " as res $sgmm_pass as txt "   failed: " as res $sgmm_fail
if $sgmm_fail > 0 {
	di as err "FAILED:$sgmm_failed"
	exit 9
}
di as res "All checks passed."
