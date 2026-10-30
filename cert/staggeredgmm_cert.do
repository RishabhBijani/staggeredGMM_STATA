* Certification script for staggeredgmm.
*
* Run from this cert/ directory:   do staggeredgmm_cert.do
*
* Part A compares staggeredgmm against reference values produced by the R
* package staggeredGMM 0.1.0 (CRAN) on 23 designs (sgmm_reference.dta, built
* by make_reference.R from the same input files in data/). Part B mirrors the
* deterministic tests of the R package's test suite: exact recovery, comparison
* counts, J-test restriction counts, input validation and invariances.
*
* The script uses the ado files in the parent directory, not any installed copy.

version 16
clear all
set more off
set linesize 120

local here "`c(pwd)'"
adopath ++ "`here'/.."
discard
mata: mata clear

global sgmm_fail 0
global sgmm_pass 0
global sgmm_failed ""
foreach k in coef se vcov {
	global sgmm_maxd_`k' 0
}

frame create ref
frame ref: use "`here'/sgmm_reference.dta", clear

mata:
real scalar sgmm_reldiff(real colvector a, real colvector b)
{
	real colvector ok

	if (rows(a) != rows(b)) return(.)
	if (any((a :>= .) :!= (b :>= .))) return(.)
	ok = selectindex(b :< .)
	if (rows(ok) == 0) return(0)
	return(max(abs(a[ok] - b[ok]) :/ rowmax((J(rows(ok), 1, 1), abs(b[ok])))))
}
end

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

* Compare a Stata matrix (any shape, compared column-major) with reference rows.
program sgmm_vec
	args case key mat tol
	tempname R
	frame ref: quietly count if rcase == "`case'" & key == "`key'"
	local nr = r(N)
	if `nr' == 0 exit
	frame ref: mkmat value if rcase == "`case'" & key == "`key'", matrix(`R')
	tempname dd
	mata: st_numscalar("`dd'", sgmm_reldiff(vec(st_matrix("`mat'")), st_matrix("`R'")))
	local d = `dd'
	local ok = (`d' < . & `d' <= `tol')
	sgmm_record `case' `key' `ok' "reldiff=`d'"
	if `tol' < 1 & inlist("`key'", "coef", "se", "vcov") {
		global sgmm_maxd_`key' = max(cond(`d' < ., `d', 1e300), ${sgmm_maxd_`key'})
	}
end

* Compare a scalar with a reference value; tol 0 means exact.
program sgmm_sca
	args case key valexp tol
	local val = `valexp'
	frame ref: quietly count if rcase == "`case'" & key == "`key'"
	local nrow = r(N)
	frame ref: quietly summarize value if rcase == "`case'" & key == "`key'", meanonly
	local n  = r(N)
	local rv = r(mean)
	if `nrow' == 0 exit
	if `n' == 0 local rv = .
	if `rv' >= . {
		local ok = (`val' >= .)
	}
	else if `tol' == 0 {
		local ok = (`val' == `rv')
	}
	else {
		local ok = (abs(`val' - `rv') / max(1, abs(`rv')) <= `tol')
	}
	sgmm_record `case' `key' `ok' "stata=`val' R=`rv'"
end

program sgmm_str
	args case key val
	frame ref: quietly levelsof svalue if rcase == "`case'" & key == "`key'", local(rv) clean
	local ok = ("`val'" == "`rv'")
	sgmm_record `case' `key' `ok' "stata=`val' R=`rv'"
end

program sgmm_case
	syntax, case(string) data(string) cmd(string asis) [JTests JFULLonly NOASSERT NOKEPT]
	local here "`c(pwd)'"
	quietly use "`here'/data/`data'.dta", clear
	capture noisily `cmd'
	if _rc {
		sgmm_record `case' run 0 "rc=`=_rc'"
		exit
	}
	local tol = cond("`noassert'" != "", 1e300, 1e-7)

	tempname C est se idf nc V W
	matrix `C'   = e(catt)
	matrix `est' = `C'[1..., 4]
	matrix `se'  = `C'[1..., 5]
	matrix `idf' = `C'[1..., 6]
	matrix `nc'  = `C'[1..., 7]
	matrix `V'   = e(V)
	sgmm_vec `case' coef `est' `tol'
	sgmm_vec `case' se `se' `tol'
	sgmm_vec `case' identified `idf' 0
	sgmm_vec `case' n_comparisons `nc' 0
	sgmm_vec `case' vcov `V' `tol'
	foreach a in CW EW {
		local l = strlower("`a'")
		sgmm_sca `case' `a'_estimate e(att_`l') `tol'
		sgmm_sca `case' `a'_std_error e(se_`l') `tol'
		sgmm_sca `case' `a'_estimate_identified e(att_`l'_id) `tol'
		sgmm_sca `case' `a'_std_error_identified e(se_`l'_id) `tol'
		sgmm_sca `case' `a'_share e(share_`l') 1e-12
		tempname w
		matrix `w' = e(w_`l')'
		sgmm_vec `case' w_`a' `w' 1e-12
	}
	matrix `W' = e(sgmm_W)
	sgmm_sca `case' W_trace `=trace(`W')' `tol'
	tempname wsum
	mata: st_numscalar("`wsum'", sum(st_matrix("`W'")))
	local ws = `wsum'
	sgmm_sca `case' W_sum `ws' `tol'
	if "`noassert'" == "" {
		foreach k in converged solve_ok n_iter {
			sgmm_sca `case' `k' e(`k') 0
		}
		sgmm_str `case' termination "`e(termination)'"
		sgmm_sca `case' max_change e(max_change) `tol'
		if "`nokept'" == "" sgmm_sca `case' n_kept e(n_kept) 0
	}
	foreach k in N_units T_periods n_cohorts N_beta N_identified N_2x2 n_dropped_2x2 ///
	             moment_rank has_never n_never_zero n_never_late cov_n_fallback      ///
	             balanced n_missing_cells {
		sgmm_sca `case' `k' e(`k') 0
	}
	sgmm_sca `case' n_obs e(N) 0
	if "`noassert'" == "" sgmm_sca `case' n_zero_lag e(n_zero_lag) 0
	tempname Ng
	matrix `Ng' = e(N_g)'
	sgmm_vec `case' N_g `Ng' 0

	local specs `""local 3" "full 3" "local 5""'
	if "`jfullonly'" != "" {
		local jtests jtests
		local specs `""full 3""'
	}
	if "`jtests'" != "" {
		foreach spec of local specs {
			gettoken typ win : spec
			local win = strtrim("`win'")
			local p J_`typ'`win'_
			local opt = cond("`typ'" == "full", "full", "window(`win')")
			capture noisily estat jtest, `opt'
			if _rc {
				sgmm_record `case' `p'run 0 "rc=`=_rc'"
				continue
			}
			tempname m mm
			matrix `m' = r(moments)
			matrix `mm' = `m'[1..., 4]
			local rJ  = r(J)
			local rp  = r(p)
			local rdf = r(df)
			local rb  = r(base_period)
			local rn  = r(n_moments)
			local rd  = r(n_dropped)
			sgmm_sca `case' `p'stat `rJ' `tol'
			sgmm_sca `case' `p'p `rp' `tol'
			sgmm_sca `case' `p'df `rdf' 0
			sgmm_sca `case' `p'base `rb' 0
			sgmm_sca `case' `p'n_moments `rn' 0
			sgmm_sca `case' `p'n_dropped `rd' 0
			sgmm_vec `case' `p'moments `mm' 1e-10
		}
	}
	di as txt "  `case' done"
end

di _n as txt "{hline 70}" _n "Part A: agreement with the R package (23 designs)" _n "{hline 70}"

local base staggeredgmm y, ivar(unit_id) time(year) gvar(cohort)
quietly {
noisily sgmm_case, case(S1)  data(sim_panel) cmd(`base') jtests
noisily sgmm_case, case(S2)  data(sim_panel) cmd(`base' weighting(cohort_toeplitz)) jtests
noisily sgmm_case, case(S3)  data(sim_panel) cmd(`base' weighting(unrestricted)) jtests
noisily sgmm_case, case(S4)  data(sim_panel) cmd(staggeredgmm y x1 x2, ivar(unit_id) time(year) gvar(cohort)) jtests
noisily sgmm_case, case(S5)  data(s5) cmd(`base') jtests
noisily sgmm_case, case(S5b) data(s5) cmd(`base' weighting(cohort_toeplitz)) jtests
noisily sgmm_case, case(S5c) data(s5) cmd(staggeredgmm y x1 x2, ivar(unit_id) time(year) gvar(cohort))
noisily sgmm_case, case(T1_pooled_toeplitz) data(t1) cmd(`base') nokept
noisily sgmm_case, case(T1_cohort_toeplitz) data(t1) cmd(`base' weighting(cohort_toeplitz)) nokept
noisily sgmm_case, case(T1_unrestricted)    data(t1) cmd(`base' weighting(unrestricted)) nokept
noisily sgmm_case, case(T2)  data(t2) cmd(`base')
noisily sgmm_case, case(T3)  data(t3) cmd(staggeredgmm y xa xb, ivar(unit_id) time(year) gvar(cohort))
noisily sgmm_case, case(T4_pooled_toeplitz) data(t4) cmd(`base') jfullonly
noisily sgmm_case, case(T4_cohort_toeplitz) data(t4) cmd(`base' weighting(cohort_toeplitz)) jfullonly
noisily sgmm_case, case(T4_unrestricted)    data(t4) cmd(`base' weighting(unrestricted)) jfullonly
noisily sgmm_case, case(T5)  data(t5) cmd(`base' nonevertreated)
noisily sgmm_case, case(T6)  data(t6) cmd(`base')
noisily sgmm_case, case(T7)  data(sim_panel) cmd(`base' maxiter(1))
noisily sgmm_case, case(T8)  data(t8) cmd(`base')
noisily sgmm_case, case(T9)  data(t9) cmd(staggeredgmm y x, ivar(unit_id) time(year) gvar(cohort))
noisily sgmm_case, case(T10) data(t10) cmd(`base' weighting(unrestricted) maxiter(1))
}

di _n as txt "Informational only (not asserted): beck_banks, not-yet-treated controls."
di as txt "The iterated weighting diverges on this design; R itself changes by up to"
di as txt "5e-6 (pooled) and 3e+1 (cohort) under a 1e-13 perturbation of the outcome."
quietly {
noisily sgmm_case, case(S6)  data(beck_banks_1977) cmd(staggeredgmm ln_gini, ivar(state) time(wrkyr) gvar(branch_reform)) noassert
noisily sgmm_case, case(S6b) data(beck_banks_1977) cmd(staggeredgmm ln_gini, ivar(state) time(wrkyr) gvar(branch_reform) weighting(cohort_toeplitz)) noassert
}

di _n as txt "{hline 70}" _n "Part B: deterministic tests" _n "{hline 70}"

* exact recovery on a noiseless panel, under every weighting
foreach w in pooled_toeplitz cohort_toeplitz unrestricted {
	quietly use "`here'/data/t1.dta", clear
	quietly staggeredgmm y, ivar(unit_id) time(year) gvar(cohort) weighting(`w')
	tempname C
	matrix `C' = e(catt)
	local ok = abs(`C'[1,4] - 1) < 1e-8 & abs(`C'[2,4] - 2) < 1e-8 & abs(`C'[3,4] - 3) < 1e-8
	local ok = `ok' & `C'[1,1] == 3 & `C'[2,1] == 3 & `C'[3,1] == 4 & `C'[3,2] == 4 & `C'[2,3] == 1
	sgmm_record B exact_recovery_`w' `ok' ""
}

* J-test restriction counts on sim_panel (Appendix E: 70 full, 14 local)
quietly use "`here'/data/sim_panel.dta", clear
quietly staggeredgmm y, ivar(unit_id) time(year) gvar(cohort)
quietly estat jtest, full
sgmm_record B jtest_full_70 `=r(n_moments) == 70 & r(base_period) == 9' ""
quietly estat jtest
sgmm_record B jtest_local_14 `=r(n_moments) == 14' ""
tempname m
matrix `m' = r(moments)
local ok 1
forvalues i = 1/`=rowsof(`m')' {
	if `m'[`i', 2] >= `m'[`i', 1] | `m'[`i', 2] == 9 local ok 0
}
sgmm_record B jtest_placebo_pretreatment `ok' ""
sgmm_record B jtest_p_formula `=reldif(r(p), chi2tail(r(df), r(J))) < 1e-12' ""

* _b[] naming and agreement with e(catt)
tempname Cc agg
matrix `Cc' = e(catt)
matrix `agg' = e(w_cw) * `Cc'[1..., 4]
sgmm_record B b_naming `=abs(_b[g10:t10] - `Cc'[1, 4]) < 1e-12' ""
sgmm_record B aggregate_identity `=reldif(e(att_cw), `agg'[1, 1]) < 1e-12' ""

* invariances
quietly staggeredgmm y, ivar(unit_id) time(year) gvar(cohort)
tempname base_b
matrix `base_b' = e(b)

preserve
set seed 101
tempvar u
quietly gen double `u' = runiform()
sort `u'
quietly staggeredgmm y, ivar(unit_id) time(year) gvar(cohort)
sgmm_record B invariance_row_order `=mreldif(e(b), `base_b') < 1e-10' ""
restore

preserve
quietly gen str7 sid = "unit" + string(unit_id, "%03.0f")
quietly staggeredgmm y, ivar(sid) time(year) gvar(cohort)
sgmm_record B invariance_string_id `=mreldif(e(b), `base_b') < 1e-10' ""
restore

preserve
quietly replace year = year + 1000
quietly replace cohort = cohort + 1000 if cohort != 0
quietly staggeredgmm y, ivar(unit_id) time(year) gvar(cohort)
tempname bs
matrix `bs' = e(b)
sgmm_record B invariance_time_shift `=colsof(`bs') == colsof(`base_b') & mreldif(`bs', `base_b') < 1e-10' ""
restore

preserve
quietly gen junk = "x"
quietly gen junk2 = .
quietly staggeredgmm y, ivar(unit_id) time(year) gvar(cohort)
sgmm_record B invariance_extra_vars `=mreldif(e(b), `base_b') < 1e-10' ""
restore

preserve
quietly replace cohort = 99 if cohort == 0
quietly staggeredgmm y, ivar(unit_id) time(year) gvar(cohort)
sgmm_record B late_cohort_as_never `=mreldif(e(b), `base_b') < 1e-10 & e(n_never_late) == 10' ""
restore

* if/in restriction equals dropping the observations
preserve
quietly staggeredgmm y if unit_id != 7, ivar(unit_id) time(year) gvar(cohort)
tempname bif
matrix `bif' = e(b)
quietly drop if unit_id == 7
quietly staggeredgmm y, ivar(unit_id) time(year) gvar(cohort)
sgmm_record B if_equals_drop `=mreldif(e(b), `bif') < 1e-12' ""
restore

* weighting() abbreviations and the explicit nevertreated option
quietly staggeredgmm y, ivar(unit_id) time(year) gvar(cohort) weighting(cohort_toeplitz)
tempname bht
matrix `bht' = e(b)
quietly staggeredgmm y, ivar(unit_id) time(year) gvar(cohort) weighting(ht)
local ok = mreldif(e(b), `bht') < 1e-12 & "`e(weighting)'" == "cohort_toeplitz"
quietly staggeredgmm y, ivar(unit_id) time(year) gvar(cohort) weighting(coh)
local ok = `ok' & mreldif(e(b), `bht') < 1e-12
sgmm_record B weighting_abbrev `ok' ""
quietly staggeredgmm y, ivar(unit_id) time(year) gvar(cohort) nevertreated
sgmm_record B nevertreated_explicit `=mreldif(e(b), `base_b') < 1e-12' ""
quietly staggeredgmm y, ivar(unit_id) time(year) gvar(cohort) nonever
sgmm_record B nonevertreated_abbrev `=e(has_never) == 0 & e(n_never_zero) > 0' ""
capture staggeredgmm y, ivar(unit_id) time(year) gvar(cohort) bogus
sgmm_record B unknown_option `=_rc == 198' ""

* level() on estimation and replay
quietly staggeredgmm y, ivar(unit_id) time(year) gvar(cohort) level(90)
local ok = e(level) == 90
capture staggeredgmm, level(80)
sgmm_record B level_option `=`ok' & _rc == 0' ""

* factor-variable covariates equal hand-made indicators
preserve
quietly gen byte xc = (x1 > 0) + (x2 > 0)
quietly tab xc, gen(xcd)
quietly staggeredgmm y i.xc, ivar(unit_id) time(year) gvar(cohort)
tempname bfv
matrix `bfv' = e(b)
quietly staggeredgmm y xcd2 xcd3, ivar(unit_id) time(year) gvar(cohort)
sgmm_record B factor_covariates `=mreldif(e(b), `bfv') < 1e-10' ""
quietly gen byte one = 1
capture staggeredgmm y i.one, ivar(unit_id) time(year) gvar(cohort)
sgmm_record B covariates_all_base `=_rc == 459' ""
restore

* input validation: each must fail with the stated return code
program sgmm_expect_error
	gettoken name 0 : 0
	gettoken rc cmd : 0
	capture `cmd'
	local got = _rc
	sgmm_record B err_`name' `=`got' == `rc'' "expected rc `rc', got `got'"
end

local s staggeredgmm y, ivar(unit_id) time(year) gvar(cohort)
quietly use "`here'/data/t1.dta", clear

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
sgmm_record B missing_outcome_ok `=_rc == 0' ""
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
sgmm_expect_error maxiter_zero 198 `s' maxiter(0)
sgmm_expect_error maxiter_fraction 198 `s' maxiter(2.5)
sgmm_expect_error tol_negative 198 `s' tolerance(-1)
sgmm_expect_error bad_weighting 198 `s' weighting(nonsense)
sgmm_expect_error both_never_options 198 `s' nevertreated nonevertreated
sgmm_expect_error covariate_is_key 198 staggeredgmm y year, ivar(unit_id) time(year) gvar(cohort)
sgmm_expect_error same_var_twice 198 staggeredgmm y, ivar(unit_id) time(unit_id) gvar(cohort)
preserve
quietly tostring y, replace force
sgmm_expect_error string_outcome 109 `s'
restore

* estat jtest errors
quietly use "`here'/data/t8.dta", clear
quietly staggeredgmm y, ivar(unit_id) time(year) gvar(cohort)
capture estat jtest
sgmm_record B jtest_needs_never `=_rc == 459' ""
quietly use "`here'/data/sim_panel.dta", clear
quietly staggeredgmm y, ivar(unit_id) time(year) gvar(cohort)
capture estat jtest, window(0)
sgmm_record B jtest_window_zero `=_rc == 198' ""
capture estat jtest, window(2.5)
sgmm_record B jtest_window_fraction `=_rc == 198' ""

* replay
capture staggeredgmm
sgmm_record B replay `=_rc == 0' ""

di _n as txt "{hline 70}"
di as txt "Largest relative difference from R across asserted cases:"
foreach k in coef se vcov {
	di as txt "  `k'" _col(12) as res %9.2e ${sgmm_maxd_`k'}
}
di as txt "{hline 70}"
di as txt "Checks passed: " as res $sgmm_pass as txt "   failed: " as res $sgmm_fail
if $sgmm_fail > 0 {
	di as err "FAILED:$sgmm_failed"
	exit 9
}
di as res "All checks passed."
