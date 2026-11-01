*! version 0.1.0  30oct2026  Parush Arora, Rishabh Bijani
*! Generalised Method of Moments Estimation of Treatment Effects under Staggered Adoption

program staggeredgmm, eclass
	version 16

	if replay() {
		if "`e(cmd)'" != "staggeredgmm" error 301
		syntax [, Level(cilevel)]
		_sgmm_display, level(`level')
		exit
	}

	local cmdline `"staggeredgmm `0'"'

	syntax varlist(min=1 fv) [if] [in], IVar(varname) TIme(varname) GVar(varname) ///
		[ WEIGHTing(string) MAXiter(real 100) TOLerance(real 1e-6)          ///
		  Level(cilevel) * ]

	gettoken depvar indepvars : varlist
	local indepvars `indepvars'
	_fv_check_depvar `depvar'

	if `maxiter' >= . | `maxiter' < 1 | `maxiter' != floor(`maxiter') {
		di as err "maxiter() must be a single whole number of at least 1."
		exit 198
	}
	if `tolerance' >= . | `tolerance' <= 0 {
		di as err "tolerance() must be a single positive finite number."
		exit 198
	}

	local wt = strlower(strtrim(`"`weighting'"'))
	if `"`wt'"' == "" local wt pooled_toeplitz
	else if inlist(`"`wt'"', "t", "gmm-t")  local wt pooled_toeplitz
	else if inlist(`"`wt'"', "ht", "gmm-ht") local wt cohort_toeplitz
	else if inlist(`"`wt'"', "u", "gmm-u")  local wt unrestricted
	else if !inlist(`"`wt'"', "pooled_toeplitz", "cohort_toeplitz", "unrestricted") {
		local hit
		local nhit 0
		foreach c in pooled_toeplitz cohort_toeplitz unrestricted {
			if strpos("`c'", `"`wt'"') == 1 {
				local hit `c'
				local ++nhit
			}
		}
		if `nhit' != 1 {
			di as err `"weighting() must be one of pooled_toeplitz, cohort_toeplitz or unrestricted (or t, ht, u); got "`weighting'"."'
			exit 198
		}
		local wt `hit'
	}

	* nevertreated and nonevertreated are parsed by hand: syntax would read a
	* typed "nonevertreated" as switching off a declared nevertreated option.
	local never_on  0
	local never_off 0
	foreach opt of local options {
		local L = strlen(`"`opt'"')
		if `L' >= 5 & `"`opt'"' == substr("nevertreated", 1, `L') {
			local never_on 1
		}
		else if `L' >= 7 & `"`opt'"' == substr("nonevertreated", 1, `L') {
			local never_off 1
		}
		else {
			di as err `"option `opt' not allowed"'
			exit 198
		}
	}
	if `never_on' & `never_off' {
		di as err "nevertreated and nonevertreated may not be combined."
		exit 198
	}
	local never = -1
	if `never_on'  local never = 1
	if `never_off' local never = 0

	if "`depvar'" == "`ivar'" | "`depvar'" == "`time'" | "`depvar'" == "`gvar'" | ///
	   "`ivar'" == "`time'" | "`ivar'" == "`gvar'" | "`time'" == "`gvar'" {
		di as err "depvar, ivar(), time() and gvar() must name four distinct variables."
		exit 198
	}

	capture confirm numeric variable `depvar'
	if _rc {
		di as err "`depvar' (the outcome) must be numeric; it is a string variable."
		di as err "Convert it explicitly rather than relying on automatic coercion."
		exit 109
	}
	foreach v in time gvar {
		capture confirm numeric variable ``v''
		if _rc {
			local role = cond("`v'" == "time", "time", "cohort")
			di as err "``v'' (`role') must be an integer-valued numeric variable; it is a string variable."
			exit 109
		}
	}

	marksample touse, novarlist
	quietly count if `touse'
	if r(N) == 0 error 2000

	local cvars
	local cnames
	if `"`indepvars'"' != "" {
		fvrevar `indepvars', list
		local cbase `r(varlist)'
		foreach v of local cbase {
			if inlist("`v'", "`depvar'", "`ivar'", "`time'", "`gvar'") {
				di as err "Covariates must not name the outcome, ivar(), time() or gvar() variable."
				exit 198
			}
			capture confirm numeric variable `v'
			if _rc {
				di as err "Covariate `v' must be numeric; string covariates are not supported."
				exit 109
			}
		}
		fvexpand `indepvars' if `touse'
		local terms `r(varlist)'
		local fvars
		foreach term of local terms {
			local parts : subinstr local term "#" " ", all
			foreach part of local parts {
				if ustrregexm("`part'", "^[0-9]+[a-z]*\.(.+)$") {
					local fvars `fvars' `=ustrregexs(1)'
				}
			}
		}
		local fvars : list uniq fvars
		foreach v of local fvars {
			quietly summarize `v' if `touse', meanonly
			if r(N) == 0 | r(min) == r(max) {
				di as err "Factor covariate `v' takes fewer than two distinct values in the estimation sample."
				exit 459
			}
		}
		foreach term of local terms {
			if ustrregexm("`term'", "(^|#)([0-9]*(b|o|bo)|co)\.") continue
			fvrevar `term'
			local cvars `cvars' `r(varlist)'
			local cnames `cnames' `term'
		}
		if "`cvars'" == "" {
			di as err "Every covariate term is a base or omitted level; no covariate is left to adjust for."
			di as err "Drop the covariates or specify terms that vary."
			exit 459
		}
	}

	tempname b V catt wcw wew Ng coh W ym Ngt grp sc
	mata: _sgmm_main("`depvar'", "`ivar'", "`time'", "`gvar'", "`cvars'",  ///
		"`cnames'", "`touse'", "`wt'", `never', `maxiter', `tolerance',     ///
		"`b'", "`V'", "`catt'", "`wcw'", "`wew'", "`Ng'", "`coh'", "`W'",    ///
		"`ym'", "`Ngt'", "`grp'", "`sc'")

	quietly count if `touse' & !missing(`depvar')
	local nobs = r(N)
	ereturn post `b' `V', esample(`touse') depname(`depvar') obs(`nobs')

	local scn : colnames `sc'
	local j 0
	foreach nm of local scn {
		local ++j
		ereturn scalar `nm' = `sc'[1, `j']
	}
	ereturn scalar level = `level'

	ereturn matrix catt         = `catt'
	ereturn matrix w_cw         = `wcw'
	ereturn matrix w_ew         = `wew'
	ereturn matrix N_g          = `Ng'
	ereturn matrix cohorts      = `coh'
	ereturn matrix sgmm_W       = `W'
	ereturn matrix sgmm_ymeans  = `ym'
	ereturn matrix sgmm_Ngt     = `Ngt'
	ereturn matrix sgmm_groups  = `grp'

	ereturn local termination   "`termination'"
	ereturn local covariates    `"`indepvars'"'
	ereturn local xnames        `"`cnames'"'
	ereturn local weighting     "`wt'"
	ereturn local gvar          "`gvar'"
	ereturn local timevar       "`time'"
	ereturn local ivar          "`ivar'"
	ereturn local depvar        "`depvar'"
	ereturn local title         "Staggered-adoption GMM estimator"
	ereturn local estat_cmd     "staggeredgmm_estat"
	ereturn local properties    "b V"
	ereturn local cmdline       `"`cmdline'"'
	ereturn local cmd           "staggeredgmm"

	_sgmm_display, level(`level')
end


program _sgmm_display
	syntax [, Level(cilevel)]

	local w "`e(weighting)'"
	if "`w'" == "pooled_toeplitz"      local wl "pooled Toeplitz (GMM-T)"
	else if "`w'" == "cohort_toeplitz" local wl "cohort-specific Toeplitz (GMM-HT)"
	else if "`w'" == "unrestricted"    local wl "unrestricted within-cohort (GMM-U)"
	else local wl "`w'"

	di _n as txt "Staggered-adoption GMM estimator"
	di as txt "Weighting: " as res "`wl'" _n
	di as txt "  Observations          " as res %8.0f e(N)
	di as txt "  Units                 " as res %8.0f e(N_units)
	di as txt "  Periods               " as res %8.0f e(T_periods)
	di as txt "  Treated cohorts       " as res %8.0f e(n_cohorts)
	local cid = string(e(N_beta)) + " (" + string(e(N_identified)) + ")"
	di as txt "  CATTs (identified)    " as res %8s "`cid'"
	di as txt "  Clean comparisons     " as res %8.0f e(N_2x2)
	di as txt "  Moment-space rank     " as res %8.0f e(moment_rank)
	di as txt "  Iterations            " as res %8.0f e(n_iter)
	di as txt "  Converged             " as res %8s cond(e(converged) == 1, "yes", "no")
	if e(converged) != 1 {
		di as txt "  Termination           " as res "`e(termination)'"
	}
	if e(balanced) != 1 {
		local mc = string(e(n_missing_cells)) + " missing cell(s)"
		di as txt "  Unbalanced            " as res "`mc'"
	}
	di ""
	_sgmm_aggline "Treated-observation ATT" cw
	_sgmm_aggline "Cohort-equal ATT" ew
	di ""

	ereturn display, level(`level')

	if e(N_identified) < e(N_beta) {
		tempname C
		matrix `C' = e(catt)
		local list
		local nlist 0
		local nun = e(N_beta) - e(N_identified)
		forvalues i = 1/`=rowsof(`C')' {
			if `C'[`i', 6] == 0 & `nlist' < 8 {
				local sep = cond(`nlist' > 0, ",", "")
				local list `list'`sep' (g=`=`C'[`i',1]', t=`=`C'[`i',2]')
				local ++nlist
			}
		}
		if `nun' > 8 local list `list', ... (`nun' total)
		di as txt "{p 0 6 2}Note: `nun' cohort-by-time effect(s) have no clean comparison and are not estimated: `list'. Their rows in e(catt) are missing.{p_end}"
	}
	if `"`e(covariates)'"' != "" {
		di as txt "{p 0 6 2}Baseline covariates: `e(xnames)'{p_end}"
		if e(cov_n_fallback) > 0 & e(cov_n_fallback) < . {
			di as txt "  `=e(cov_n_fallback)' comparison(s) fell back to the unconditional form."
		}
		di as txt "{p 0 6 2}Note: {cmd:estat jtest} does not adjust for these covariates. It tests " ///
			"{bf:unconditional} parallel trends, not the covariate-conditional parallel trends " ///
			"assumed by this fit; see {help staggeredgmm_postestimation##jtest:estat jtest}.{p_end}"
	}
end


program _sgmm_aggline
	args label which
	if e(att_`which') >= . {
		di as txt %-28s "`label'" as res "      NA" as txt "  (identified subset: " ///
			as res %6.4f e(att_`which'_id) as txt ", se " as res %6.4f e(se_`which'_id) ///
			as txt ", " as res %4.1f 100 * e(share_`which') as txt "% of weight)"
	}
	else {
		di as txt %-28s "`label'" " " as res %9.4f e(att_`which') ///
			as txt "  (se " as res %6.4f e(se_`which') as txt ")"
	}
end


version 16
mata:

struct _sgmm_panel {
	real scalar      N_units, T_min, T_max, TT, has_never
	real scalar      n_never_zero, n_never_late, n_missing_cells, n_obs, has_cov
	real colvector   unit_row, time_idx, y, cohort_of, grp_of, in_design
	real colvector   treated, all_groups, N_g
	real matrix      Y, obs, N_gt, X_unit
	string colvector unit_lab
}


string scalar _sgmm_num(real scalar x)
{
	return(strtrim(strofreal(x, "%20.0g")))
}


string colvector _sgmm_numlab(real colvector x)
{
	string colvector s
	real scalar      i

	s = J(rows(x), 1, "")
	for (i = 1; i <= rows(x); i++) s[i] = _sgmm_num(x[i])
	return(s)
}


string colvector _sgmm_slab(real colvector x)
{
	string colvector out
	real scalar      i

	out = J(rows(x), 1, "")
	for (i = 1; i <= rows(x); i++) {
		out[i] = (x[i] < 0 ? "m" + _sgmm_num(-x[i]) : _sgmm_num(x[i]))
	}
	return(out)
}


string scalar _sgmm_preview(string colvector x, real scalar n)
{
	string colvector u
	string scalar    out
	real scalar      i, j, seen, total

	if (rows(x) == 0) return("")
	total = rows(uniqrows(x))
	u = J(0, 1, "")
	for (i = 1; i <= rows(x); i++) {
		if (rows(u) >= n) break
		seen = 0
		for (j = 1; j <= rows(u); j++) {
			if (u[j] == x[i]) {
				seen = 1
				break
			}
		}
		if (!seen) u = u \ x[i]
	}
	out = invtokens(u', ", ")
	if (total > n) out = out + ", ... (" + _sgmm_num(total) + " total)"
	return(out)
}


void _sgmm_abort(string scalar msg, real scalar rc)
{
	errprintf("%s\n", msg)
	exit(rc)
}


void _sgmm_warn(string scalar msg)
{
	printf("{p 0 9 2}{err}Warning:{txt} %s{p_end}\n", msg)
}


real scalar _sgmm_integerish(real colvector x)
{
	real colvector f

	f = select(x, x :< .)
	if (rows(f) == 0) return(0)
	return(all(f :== trunc(f)))
}


void _sgmm_check_missing(real colvector miss, real colvector obsno,
                         string scalar nm, string scalar role)
{
	real colvector w

	if (any(miss)) {
		w = select(obsno, miss)
		_sgmm_abort(nm + " (" + role + ") contains " + _sgmm_num(rows(w)) +
			" missing value(s); observations: " +
			_sgmm_preview(_sgmm_numlab(w), 5) + ".", 459)
	}
}


real colvector _sgmm_group(transmorphic colvector x, transmorphic colvector uniq)
{
	real colvector p, grp, newg
	real scalar    i, n

	n = rows(x)
	p = order(x, 1)
	newg = J(n, 1, 0)
	newg[1] = 1
	for (i = 2; i <= n; i++) {
		if (x[p[i]] != x[p[i - 1]]) newg[i] = 1
	}
	grp = J(n, 1, .)
	grp[p] = runningsum(newg)
	uniq = x[select(p, newg)]
	return(grp)
}


real colvector _sgmm_match(real colvector x, real colvector table)
{
	real colvector out
	real scalar    i, j

	out = J(rows(x), 1, .)
	for (i = 1; i <= rows(x); i++) {
		for (j = 1; j <= rows(table); j++) {
			if (x[i] == table[j]) {
				out[i] = j
				break
			}
		}
	}
	return(out)
}


struct _sgmm_panel scalar _sgmm_validate(string scalar yname, string scalar idname,
	string scalar tname, string scalar gname, string rowvector xvars,
	string scalar touse, real scalar never)
{
	struct _sgmm_panel scalar P
	real colvector obsno, y, t, g, idn, miss, unit_row, time_idx, ord, dup
	real colvector gs, bad, cohort_of, present, isp, gaps, ucoh, neg, always
	real colvector treated, is_never, all_groups, in_design, grp_of, N_g, ix
	real colvector first, badu, w
	real matrix    Y, obs, N_gt, Xf, firstval
	string colvector ids, ulab
	transmorphic colvector uniq
	real scalar    n, isstr, i, j, u, N, T_min, T_max, TT, n0, nl, nn
	real scalar    has_never, G, gi, nc, v, n_always

	obsno = selectindex(st_data(., touse) :!= 0)
	n = rows(obsno)
	if (n == 0) _sgmm_abort("no observations", 2000)

	y = st_data(obsno, yname)
	t = st_data(obsno, tname)
	g = st_data(obsno, gname)
	isstr = st_isstrvar(idname)

	if (!_sgmm_integerish(t)) {
		_sgmm_abort(tname + " (time) must be an integer-valued numeric variable. " +
			"Non-integer values would be silently truncated.", 459)
	}
	if (!_sgmm_integerish(g)) {
		_sgmm_abort(gname + " (cohort) must be an integer-valued numeric variable. " +
			"Code never-treated units as 0.", 459)
	}

	if (isstr) {
		ids  = st_sdata(obsno, idname)
		miss = (ids :== "")
	}
	else {
		idn  = st_data(obsno, idname)
		miss = (idn :>= .)
	}
	_sgmm_check_missing(miss, obsno, idname, "unit id")
	_sgmm_check_missing((t :>= .), obsno, tname, "time")
	_sgmm_check_missing((g :>= .), obsno, gname, "cohort")

	uniq = J(0, 1, .)
	if (isstr) {
		unit_row = _sgmm_group(ids, uniq)
		ulab = uniq
	}
	else {
		unit_row = _sgmm_group(idn, uniq)
		ulab = _sgmm_numlab(uniq)
	}
	N = rows(ulab)

	ord = order((unit_row, t, (1::n)), (1, 2, 3))
	dup = J(n, 1, 0)
	for (i = 2; i <= n; i++) {
		if (unit_row[ord[i]] == unit_row[ord[i - 1]] & t[ord[i]] == t[ord[i - 1]]) {
			dup[ord[i]] = 1
		}
	}
	if (any(dup)) {
		w = select(obsno, dup)
		_sgmm_abort("Each unit-period combination must appear exactly once; found " +
			_sgmm_num(rows(w)) + " duplicate observation(s). Offending observations: " +
			_sgmm_preview(_sgmm_numlab(w), 5) + ".", 459)
	}

	T_min = min(t)
	T_max = max(t)
	TT    = T_max - T_min + 1
	if (TT < 2) {
		_sgmm_abort("At least two distinct time periods are required; " + tname +
			" takes only one value.", 459)
	}
	time_idx = t :- T_min :+ 1

	ord      = order((unit_row, time_idx), (1, 2))
	unit_row = unit_row[ord]
	time_idx = time_idx[ord]
	y        = y[ord]
	gs       = g[ord]

	cohort_of = J(N, 1, .)
	bad       = J(N, 1, 0)
	for (i = 1; i <= n; i++) {
		u = unit_row[i]
		if (cohort_of[u] >= .) cohort_of[u] = gs[i]
		else if (cohort_of[u] != gs[i]) bad[u] = 1
	}
	if (any(bad)) {
		ix = selectindex(bad)
		_sgmm_abort(gname + " (cohort) must be constant within unit; " +
			_sgmm_num(rows(ix)) + " unit(s) have more than one cohort value. Affected: " +
			_sgmm_preview(ulab[ix], 5) + ".", 459)
	}

	present = uniqrows(t)
	if (rows(present) != TT) {
		isp = J(TT, 1, 0)
		isp[present :- T_min :+ 1] = J(rows(present), 1, 1)
		gaps = select((T_min::T_max), (isp :== 0))
		_sgmm_abort("The set of periods present in the data must be consecutive: " +
			"no period may be absent for every unit. Period(s) absent from the whole panel: " +
			_sgmm_preview(_sgmm_numlab(gaps), 5) + ". The stationary weighting indexes " +
			"autocovariances by calendar distance, so an empty period would create lag " +
			"distances that correspond to nothing in the data. Individual units may be " +
			"missing individual periods; that is the unbalanced case and is supported.", 459)
	}

	ucoh = uniqrows(cohort_of)
	neg  = select(ucoh, (ucoh :< 0))
	if (rows(neg) > 0) {
		_sgmm_abort(gname + " (cohort) contains negative value(s): " +
			_sgmm_preview(_sgmm_numlab(neg), 5) +
			". Cohorts must be 0 (never treated) or a treated period.", 459)
	}
	if (T_min < 0 & T_max >= 0 & any(cohort_of :== 0)) {
		_sgmm_abort(gname + " (cohort) codes never-treated units as 0, but period 0 lies " +
			"inside the observed window (" + _sgmm_num(T_min) + " to " + _sgmm_num(T_max) +
			"), so 0 could also mean first treated in period 0. Recode never-treated units " +
			"to a value after the last period (for example " + _sgmm_num(T_max + 1) +
			"), which is read as never treated within the window, or shift " + tname +
			" and the treated cohort values by the same amount so that every period " +
			"is positive.", 459)
	}
	always = select(ucoh, ((ucoh :> 0) :& (ucoh :<= T_min)))
	if (rows(always) > 0) {
		n_always = sum((cohort_of :> 0) :& (cohort_of :<= T_min))
		_sgmm_abort(_sgmm_num(n_always) + " unit(s) in cohort(s) " +
			_sgmm_preview(_sgmm_numlab(always), 5) +
			" are treated at or before the first observed period (" + _sgmm_num(T_min) +
			"). Such units have no pre-treatment period, so none of their " +
			"cohort-by-time effects is identified, and including them would " +
			"silently contaminate the aggregate.", 459)
	}
	treated = select(ucoh, ((ucoh :> T_min) :& (ucoh :<= T_max)))
	if (rows(treated) == 0) {
		_sgmm_abort("No treated cohort is identified: " + gname +
			" contains no value strictly inside the observed period range (" +
			_sgmm_num(T_min) + ", " + _sgmm_num(T_max) + "].", 459)
	}

	is_never = ((cohort_of :== 0) :| (cohort_of :> T_max))
	n0 = sum(cohort_of :== 0)
	nl = sum(cohort_of :> T_max)
	nn = n0 + nl

	if (never == 1 & nn == 0) {
		_sgmm_abort("nevertreated was requested, but no unit is never treated within " +
			"the observed window: no unit has " + gname + " == 0 or " + gname + " > " +
			_sgmm_num(T_max) + ". Omit the option to detect automatically, or specify " +
			"nonevertreated to rely on not-yet-treated controls alone.", 459)
	}
	has_never = (never == -1 ? (nn > 0) : never)
	if (never == 0 & nn > 0) {
		_sgmm_warn("nonevertreated: " + _sgmm_num(nn) + " never-treated unit(s) are " +
			"present but will not be used as controls. Precision will be lower and some " +
			"cohort-by-time effects may become unidentified.")
	}

	grp_of     = cohort_of :* (is_never :== 0)
	all_groups = (has_never ? (0 \ treated) : treated)
	in_design  = (has_never ? J(N, 1, 1) : (is_never :== 0))

	G   = rows(all_groups)
	N_g = J(G, 1, 0)
	for (gi = 1; gi <= G; gi++) {
		N_g[gi] = sum((grp_of :== all_groups[gi]) :& in_design)
	}
	if (any(N_g :== 0)) {
		_sgmm_abort("Internal inconsistency: a group contains no units. " +
			"Please report this with a reproducible example.", 9)
	}

	Y = J(N, TT, .)
	for (i = 1; i <= n; i++) Y[unit_row[i], time_idx[i]] = y[i]
	obs = (Y :< .)

	N_gt = J(G, TT, 0)
	for (gi = 1; gi <= G; gi++) {
		ix = selectindex((grp_of :== all_groups[gi]) :& in_design)
		N_gt[gi, .] = colsum(obs[ix, .])
	}

	P.has_cov = (cols(xvars) > 0)
	if (P.has_cov) {
		Xf = (J(n, 1, 1), st_data(obsno, xvars))
		Xf = Xf[ord, .]
		first = J(N, 1, .)
		for (i = n; i >= 1; i--) first[unit_row[i]] = i
		P.X_unit = Xf[first, .]

		nc       = cols(Xf)
		firstval = J(N, nc, .)
		badu     = J(N, 1, 0)
		for (i = 1; i <= n; i++) {
			u = unit_row[i]
			for (j = 1; j <= nc; j++) {
				v = Xf[i, j]
				if (v < .) {
					if (firstval[u, j] >= .) firstval[u, j] = v
					else if (firstval[u, j] != v) badu[u] = 1
				}
			}
		}
		if (any(badu)) {
			ix = selectindex(badu)
			_sgmm_warn(_sgmm_num(rows(ix)) + " unit(s) have a covariate value that is " +
				"not constant over time. Covariates are documented as baseline " +
				"(pre-treatment, time-invariant); each affected unit's value at its " +
				"earliest observed period has been used. Affected units: " +
				_sgmm_preview(ulab[ix], 5) + ".")
		}
	}

	P.N_units         = N
	P.T_min           = T_min
	P.T_max           = T_max
	P.TT              = TT
	P.has_never       = has_never
	P.n_never_zero    = n0
	P.n_never_late    = nl
	P.n_obs           = n
	P.unit_row        = unit_row
	P.time_idx        = time_idx
	P.y               = y
	P.cohort_of       = cohort_of
	P.grp_of          = grp_of
	P.in_design       = in_design
	P.treated         = treated
	P.all_groups      = all_groups
	P.N_g             = N_g
	P.Y               = Y
	P.obs             = obs
	P.N_gt            = N_gt
	P.n_missing_cells = sum(obs :== 0)
	P.unit_lab        = ulab
	return(P)
}


real matrix _sgmm_cells(real colvector treated, real scalar T_max)
{
	real matrix C
	real scalar k, s, pos

	C   = J(sum(T_max :- treated :+ 1), 2, .)
	pos = 0
	for (k = 1; k <= rows(treated); k++) {
		for (s = treated[k]; s <= T_max; s++) {
			pos++
			C[pos, 1] = treated[k]
			C[pos, 2] = s
		}
	}
	return(C)
}


real matrix _sgmm_contrasts(real matrix C, real colvector treated,
                            real scalar T_min, real scalar has_never)
{
	real matrix    CD
	real colvector ctrl
	real scalar    nb, total, i, a, s, pos, g, tt

	nb    = rows(C)
	total = 0
	for (i = 1; i <= nb; i++) {
		total = total + (C[i, 1] - T_min) * (sum(treated :> C[i, 2]) + has_never)
	}
	if (total == 0) return(J(0, 5, .))

	CD  = J(total, 5, .)
	pos = 0
	for (i = 1; i <= nb; i++) {
		g  = C[i, 1]
		tt = C[i, 2]
		ctrl = select(treated, (treated :> tt))
		if (has_never) ctrl = (0 \ ctrl)
		for (a = 1; a <= rows(ctrl); a++) {
			for (s = T_min; s <= g - 1; s++) {
				pos++
				CD[pos, 1] = g
				CD[pos, 2] = tt
				CD[pos, 3] = ctrl[a]
				CD[pos, 4] = s
				CD[pos, 5] = i
			}
		}
	}
	return(CD)
}


real matrix _sgmm_gmeans(real matrix M, real colvector grp_of,
                         real colvector in_design, real colvector groups)
{
	real matrix    out, Mg
	real colvector ix
	real scalar    gi

	out = J(rows(groups), cols(M), .)
	for (gi = 1; gi <= rows(groups); gi++) {
		ix = selectindex((grp_of :== groups[gi]) :& in_design)
		if (rows(ix) > 0) {
			Mg = M[ix, .]
			out[gi, .] = colsum(Mg) :/ colsum(Mg :< .)
		}
	}
	return(out)
}


real scalar _sgmm_qr_rank(real matrix X0, real scalar tol)
{
	real matrix    x
	real rowvector qraux, work1, work2, perm
	real scalar    n, p, lup, k, l, j, jj, m, nrmxl, t, tt

	x = X0
	n = rows(x)
	p = cols(x)
	qraux = sqrt(colsum(x :^ 2))
	work1 = qraux
	work2 = qraux
	for (j = 1; j <= p; j++) {
		if (work2[j] == 0) work2[j] = 1
	}
	lup = min((n, p))
	k   = p + 1
	for (l = 1; l <= lup; l++) {
		while (!(l >= k | qraux[l] >= work2[l] * tol)) {
			perm = J(1, p, .)
			m = 0
			for (jj = 1; jj <= p; jj++) {
				if (jj != l) {
					m++
					perm[m] = jj
				}
			}
			perm[p] = l
			x     = x[., perm]
			qraux = qraux[perm]
			work1 = work1[perm]
			work2 = work2[perm]
			k--
		}
		if (l == n) continue
		nrmxl = sqrt(sum(x[|l, l \ n, l|] :^ 2))
		if (nrmxl == 0) continue
		if (x[l, l] != 0) nrmxl = (x[l, l] < 0 ? -abs(nrmxl) : abs(nrmxl))
		x[|l, l \ n, l|] = x[|l, l \ n, l|] :/ nrmxl
		x[l, l] = 1 + x[l, l]
		for (j = l + 1; j <= p; j++) {
			t = -sum(x[|l, l \ n, l|] :* x[|l, j \ n, j|]) / x[l, l]
			x[|l, j \ n, j|] = x[|l, j \ n, j|] + t :* x[|l, l \ n, l|]
			if (qraux[j] != 0) {
				tt = 1 - (abs(x[l, j]) / qraux[j])^2
				if (tt < 0) tt = 0
				if (abs(tt) < 1e-6) {
					qraux[j] = sqrt(sum(x[|l + 1, j \ n, j|] :^ 2))
					work1[j] = qraux[j]
				}
				else {
					qraux[j] = qraux[j] * sqrt(tt)
				}
			}
		}
		qraux[l] = x[l, l]
		x[l, l]  = -nrmxl
	}
	return(min((k - 1, n)))
}


real colvector _sgmm_delta_cov(real matrix CD, struct _sgmm_panel scalar P,
	real matrix ym, real matrix cm, real colvector g_pos, real colvector c_pos,
	real colvector tp, real colvector tr, real scalar n_fb)
{
	real colvector ctrl, o, grp, units, dY, keep, gam
	real matrix    Xc
	real scalar    n2, ncov, i, j, k, r1, ok

	n2   = rows(CD)
	ncov = cols(P.X_unit)
	ctrl = J(n2, 1, .)
	n_fb = 0
	o = order((CD[., 3], CD[., 4], CD[., 2], (1::n2)), (1, 2, 3, 4))

	i = 1
	while (i <= n2) {
		j = i
		while (j < n2) {
			if (CD[o[j + 1], 3] == CD[o[i], 3] & CD[o[j + 1], 4] == CD[o[i], 4] &
			    CD[o[j + 1], 2] == CD[o[i], 2]) j++
			else break
		}
		grp = o[|i \ j|]
		r1  = grp[1]

		units = selectindex((P.grp_of :== CD[r1, 3]) :& P.in_design)
		dY    = P.Y[units, tp[r1]] - P.Y[units, tr[r1]]
		Xc    = P.X_unit[units, .]
		keep  = ((rowmissing(Xc) :== 0) :& (dY :< .))

		ok = 0
		if (sum(keep) > ncov) {
			Xc = select(Xc, keep)
			dY = select(dY, keep)
			if (_sgmm_qr_rank(Xc, 1e-7) == ncov) {
				gam = qrsolve(Xc, dY)
				if (!hasmissing(gam)) ok = 1
			}
		}

		if (ok) {
			ctrl[grp] = cm[g_pos[grp], .] * gam
		}
		else {
			n_fb = n_fb + rows(grp)
			for (k = 1; k <= rows(grp); k++) {
				ctrl[grp[k]] = ym[c_pos[grp[k]], tp[grp[k]]] - ym[c_pos[grp[k]], tr[grp[k]]]
			}
		}
		i = j + 1
	}
	return(ctrl)
}


real matrix _sgmm_fe_resid(real colvector yadj, real colvector unit_row,
                           real colvector time_idx, real scalar N, real scalar TT)
{
	real matrix    R, D, Dk, info
	real colvector ok, u, tt, yy, delta, res
	real rowvector keepc
	real scalar    n, i, p, a, b

	R  = J(N, TT, .)
	ok = selectindex(yadj :< .)
	if (rows(ok) == 0) return(R)
	u  = unit_row[ok]
	tt = time_idx[ok]
	yy = yadj[ok]
	n  = rows(yy)

	D = J(n, TT, 0)
	for (i = 1; i <= n; i++) D[i, tt[i]] = 1
	keepc = selectindex(colsum(D) :> 0)

	info = panelsetup(u, 1)
	for (p = 1; p <= rows(info); p++) {
		a = info[p, 1]
		b = info[p, 2]
		yy[|a \ b|] = yy[|a \ b|] :- mean(yy[|a \ b|])
		D[|a, 1 \ b, TT|] = D[|a, 1 \ b, TT|] :- mean(D[|a, 1 \ b, TT|])
	}

	if (cols(keepc) <= 1) {
		res = yy
	}
	else {
		keepc = keepc[2..cols(keepc)]
		Dk    = D[., keepc]
		delta = qrsolve(Dk, yy)
		res   = yy - Dk * delta
	}
	for (i = 1; i <= n; i++) R[u[i], tt[i]] = res[i]
	return(R)
}


real colvector _sgmm_autocov(real matrix R, real colvector npair)
{
	real matrix    O, Rz
	real colvector sigma
	real scalar    TT, nr, d, np

	TT = cols(R)
	nr = rows(R)
	O  = (R :< .)
	Rz = editmissing(R, 0)
	sigma = J(TT, 1, 0)
	npair = J(TT, 1, 0)
	for (d = 0; d <= TT - 1; d++) {
		np = sum(O[|1, 1 \ nr, TT - d|] :* O[|1, 1 + d \ nr, TT|])
		npair[d + 1] = np
		if (np > 0) {
			sigma[d + 1] = sum(Rz[|1, 1 \ nr, TT - d|] :* Rz[|1, 1 + d \ nr, TT|]) / np
		}
	}
	return(sigma)
}


real matrix _sgmm_unrcov(real matrix R)
{
	real matrix O, Rz

	O  = (R :< .)
	Rz = editmissing(R, 0)
	return(editmissing(cross(Rz, Rz) :/ cross(O, O), 0))
}


real matrix _sgmm_overlap(real matrix obs_g)
{
	real rowvector na

	na = colsum(obs_g)
	return(editmissing(cross(obs_g, obs_g) :/ (na' * na), 0))
}


real matrix _sgmm_toeplitz(real colvector s)
{
	real matrix T
	real scalar n, i, j

	n = rows(s)
	T = J(n, n, .)
	for (i = 1; i <= n; i++) {
		for (j = 1; j <= n; j++) T[i, j] = s[abs(i - j) + 1]
	}
	return(T)
}


real matrix _sgmm_weight(real matrix R, real matrix obs, real colvector grp_of,
	real colvector in_design, real colvector groups, string scalar weighting,
	real scalar nzl)
{
	real matrix    W, S
	real colvector sp, sg, np, ix
	real scalar    TT, G, gi, a, b

	TT  = cols(R)
	G   = rows(groups)
	W   = J(G * TT, G * TT, 0)
	nzl = 0
	np  = J(0, 1, .)
	if (weighting == "pooled_toeplitz") {
		sp  = _sgmm_autocov(R[selectindex(in_design), .], np)
		nzl = sum(np :== 0)
	}
	for (gi = 1; gi <= G; gi++) {
		ix = selectindex((grp_of :== groups[gi]) :& in_design)
		if (rows(ix) == 0) continue
		if (weighting == "pooled_toeplitz") {
			S = _sgmm_toeplitz(sp)
		}
		else if (weighting == "cohort_toeplitz") {
			sg  = _sgmm_autocov(R[ix, .], np)
			nzl = nzl + sum(np :== 0)
			S   = _sgmm_toeplitz(sg)
		}
		else {
			S = _sgmm_unrcov(R[ix, .])
		}
		a = (gi - 1) * TT + 1
		b = gi * TT
		W[|a, a \ b, b|] = S :* _sgmm_overlap(obs[ix, .])
	}
	return(W)
}


real matrix _sgmm_reduced_M(real colvector d, real matrix v, real matrix W)
{
	real matrix M

	M = (d * d') :* (v' * W * v)
	return((M + M') :/ 2)
}


real scalar _sgmm_reduced_solve(real colvector d, real matrix v, real matrix uQ,
	real colvector uDe, real matrix W, real scalar tol_rel, real matrix QtAQ,
	real colvector QtAD, real scalar n_kept)
{
	real matrix    M, X, Vk
	real vector    L
	real colvector lam, keep, dk
	real scalar    lmax

	M = _sgmm_reduced_M(d, v, W)
	symeigensystem(M, X, L)
	lam  = colshape(L, 1)
	lmax = max(lam)
	if (lmax >= . | lmax <= 0) return(0)
	keep = selectindex(lam :> tol_rel * lmax)
	if (rows(keep) == 0) return(0)
	Vk   = X[., keep]
	dk   = 1 :/ lam[keep]
	QtAQ = uQ' * (Vk * (dk :* (Vk' * uQ)))
	QtAD = uQ' * (Vk * (dk :* (Vk' * uDe)))
	n_kept = rows(keep)
	return(1)
}


real scalar _sgmm_rsolve(real matrix A, real matrix B, real matrix X)
{
	real matrix Ainv
	real scalar rc

	Ainv = lusolve(A, I(rows(A)), 1e-300)
	if (hasmissing(Ainv)) return(0)
	rc = 1 / (max(colsum(abs(A))) * max(colsum(abs(Ainv))))
	if (rc >= . | rc < epsilon(1)) return(0)
	if (rows(B) == 0) X = Ainv
	else X = lusolve(A, B, 1e-300)
	return(1)
}


real matrix _sgmm_ginv(real matrix A)
{
	real matrix    U, Vt
	real colvector s, pos
	real vector    s0

	svd(A, U, s0, Vt)
	s   = colshape(s0, 1)
	pos = selectindex(s :> max((sqrt(epsilon(1)) * max(s), 0)))
	if (rows(pos) == 0) return(J(cols(A), rows(A), 0))
	return(Vt[pos, .]' * diag(1 :/ s[pos]) * U[., pos]')
}


void _sgmm_agg(real colvector w, real colvector beta, real matrix V,
               real colvector id_idx, real scalar nb, real rowvector out)
{
	real colvector ident, wid, ws
	real scalar    full, share, es, ss

	ident = J(nb, 1, 0)
	ident[id_idx] = J(rows(id_idx), 1, 1)
	full  = (sum(w :* (ident :== 0)) <= 0)
	wid   = w[id_idx]
	share = sum(wid)
	if (share > 0) {
		ws = wid :/ share
		es = sum(ws :* beta)
		ss = sqrt(max((0, ws' * V * ws)))
	}
	else {
		es = .
		ss = .
	}
	out = ((full ? es : .), (full ? ss : .), es, ss, share)
}


void _sgmm_main(string scalar yname, string scalar idname, string scalar tname,
	string scalar gname, string scalar xvars_s, string scalar xlab_s,
	string scalar touse, string scalar weighting, real scalar never,
	real scalar max_iter, real scalar tol,
	string scalar bname, string scalar Vname, string scalar cattname,
	string scalar wcwname, string scalar wewname, string scalar Ngname,
	string scalar cohname, string scalar Wname, string scalar ymname,
	string scalar Ngtname, string scalar grpname, string scalar scname)
{
	struct _sgmm_panel scalar P
	real matrix    C, CD, ym, cm, Q, QS, uQ, U, Uu, Vt, v, W_last, QtAQ, QtAQ_out
	real matrix    V_id, M, R, catt, beta_new
	real colvector g_pos, c_pos, tp, tr, ok, keep, tchg, ctrl, Delta, ncomp
	real colvector id_idx, ident, beta, d, uDe, QtAD, beta_old
	real colvector beta_full, tau, yadj, cell_id, offset, se_full, Tg_of
	real colvector Ng_of, wcw, wew, cpos, Ng_treat, Tg_treat, ix, qtq
	real vector    s0
	real rowvector aggcw, aggew, scv
	string matrix  bstripe, cstripe
	string rowvector scn
	real scalar    T_min, T_max, TT, G, nb, n2, n_dropped, cov_fb, nid, i, k, r
	real scalar    rank, converged, solve_ok, n_iter, max_change, n_kept, nzl
	real scalar    it, kk, c, K, denom, solved, nzl_it, nk_it
	real scalar    jj, maxdim
	real colvector bg
	string rowvector xl
	string colvector bl
	string scalar  termination

	P = _sgmm_validate(yname, idname, tname, gname, tokens(xvars_s), touse, never)

	T_min = P.T_min
	T_max = P.T_max
	TT    = P.TT
	G     = rows(P.all_groups)

	C  = _sgmm_cells(P.treated, T_max)
	nb = rows(C)

	// Stored results include a (G*TT) by (G*TT) weight matrix and nb by nb
	// matrices; check them against the matrix limit of this Stata edition first.
	maxdim = (c("MP") ? 65534 : (c("SE") ? 11000 : 800))
	if (max((G * TT, nb)) > maxdim) {
		_sgmm_abort("This design needs Stata matrices of dimension " +
			_sgmm_num(max((G * TT, nb))) + " (" + _sgmm_num(G) + " groups x " +
			_sgmm_num(TT) + " periods for the weight matrix; " + _sgmm_num(nb) +
			" cohort-by-time effects), but this edition of Stata allows at most " +
			_sgmm_num(maxdim) + ". Use Stata/SE or Stata/MP, or shorten the panel.", 908)
	}

	CD = _sgmm_contrasts(C, P.treated, T_min, P.has_never)
	if (rows(CD) == 0) {
		_sgmm_abort("No clean comparison exists for this design: every treated cohort " +
			"lacks both a never-treated control and a later not-yet-treated cohort. " +
			"No treatment effect is identified.", 459)
	}

	ym = _sgmm_gmeans(P.Y, P.grp_of, P.in_design, P.all_groups)

	g_pos = _sgmm_match(CD[., 1], P.all_groups)
	c_pos = _sgmm_match(CD[., 3], P.all_groups)
	tp    = CD[., 2] :- T_min :+ 1
	tr    = CD[., 4] :- T_min :+ 1

	n2 = rows(CD)
	ok = J(n2, 1, 0)
	for (i = 1; i <= n2; i++) {
		ok[i] = (ym[g_pos[i], tp[i]] < .) & (ym[g_pos[i], tr[i]] < .) &
		        (ym[c_pos[i], tp[i]] < .) & (ym[c_pos[i], tr[i]] < .)
	}
	n_dropped = sum(ok :== 0)
	if (n_dropped > 0) {
		keep  = selectindex(ok)
		if (rows(keep) == 0) {
			_sgmm_warn(_sgmm_num(n_dropped) + " of " + _sgmm_num(n_dropped) +
				" clean comparison(s) reference a group-period with no observed unit and " +
				"have been dropped from the moment system.")
			_sgmm_abort("Every clean comparison references a group-period with no " +
				"observed unit; no treatment effect is identified.", 459)
		}
		CD    = CD[keep, .]
		g_pos = g_pos[keep]
		c_pos = c_pos[keep]
		tp    = tp[keep]
		tr    = tr[keep]
		_sgmm_warn(_sgmm_num(n_dropped) + " of " + _sgmm_num(n_dropped + rows(CD)) +
			" clean comparison(s) reference a group-period with no observed unit and " +
			"have been dropped from the moment system.")
	}
	n2 = rows(CD)

	tchg = J(n2, 1, .)
	for (i = 1; i <= n2; i++) tchg[i] = ym[g_pos[i], tp[i]] - ym[g_pos[i], tr[i]]

	cov_fb = 0
	if (!P.has_cov) {
		ctrl = J(n2, 1, .)
		for (i = 1; i <= n2; i++) ctrl[i] = ym[c_pos[i], tp[i]] - ym[c_pos[i], tr[i]]
	}
	else {
		cm   = _sgmm_gmeans(P.X_unit, P.grp_of, P.in_design, P.all_groups)
		ctrl = _sgmm_delta_cov(CD, P, ym, cm, g_pos, c_pos, tp, tr, cov_fb)
		if (cov_fb > 0) {
			_sgmm_warn(_sgmm_num(cov_fb) + " of " + _sgmm_num(n2) + " comparison(s) " +
				"could not identify the covariate adjustment (too few complete-case " +
				"control units, or a rank-deficient fit) and reverted to the " +
				"unconditional comparison for those contrasts only.")
		}
	}
	Delta = tchg - ctrl
	if (hasmissing(Delta)) {
		bg = uniqrows(select(CD[., 1], (Delta :>= .)))
		xl = tokens(xlab_s)
		bl = J(0, 1, "")
		for (i = 1; i <= rows(bg); i++) {
			k = _sgmm_match(bg[i], P.all_groups)
			for (jj = 2; jj <= cols(cm); jj++) {
				if (cm[k, jj] >= .) bl = bl \ (xl[jj - 1] + " in cohort " + _sgmm_num(bg[i]))
			}
		}
		_sgmm_abort("A covariate has no observed value for any unit of a treated cohort, " +
			"so that cohort's covariate mean is undefined and its covariate-adjusted " +
			"comparisons (" + _sgmm_num(sum(Delta :>= .)) + " of " + _sgmm_num(n2) +
			") cannot be formed. Affected: " + _sgmm_preview(bl, 8) + ". Drop the " +
			"covariate, or supply it for at least one unit of each treated cohort.", 459)
	}

	Q = J(n2, nb, 0)
	for (i = 1; i <= n2; i++) Q[i, CD[i, 5]] = 1
	ncomp  = colsum(Q)'
	id_idx = selectindex(ncomp :> 0)
	if (rows(id_idx) == 0) {
		_sgmm_abort("No cohort-by-time effect is anchored by a clean comparison; " +
			"nothing is identified.", 459)
	}
	nid   = rows(id_idx)
	ident = (ncomp :> 0)

	if (nid < nb) {
		ix = selectindex(ident :== 0)
		cstripe = J(rows(ix), 1, "")
		for (i = 1; i <= rows(ix); i++) {
			cstripe[i] = "(g=" + _sgmm_num(C[ix[i], 1]) + ", t=" + _sgmm_num(C[ix[i], 2]) + ")"
		}
		_sgmm_warn(_sgmm_num(nb - nid) + " of " + _sgmm_num(nb) + " cohort-by-time " +
			"effect(s) have no clean comparison and are reported as missing, not as an " +
			"estimated zero. Affected: " + _sgmm_preview(cstripe, 8) + "." +
			(P.has_never ? "" : " A never-treated control group, if available, would identify more."))
	}

	QS   = Q[., id_idx]
	qtq  = 1 :/ colsum(QS)'
	beta = qtq :* (QS' * Delta)

	U = J(n2, G * TT, 0)
	for (i = 1; i <= n2; i++) {
		U[i, (g_pos[i] - 1) * TT + tp[i]] = U[i, (g_pos[i] - 1) * TT + tp[i]] + 1
		U[i, (g_pos[i] - 1) * TT + tr[i]] = U[i, (g_pos[i] - 1) * TT + tr[i]] - 1
		U[i, (c_pos[i] - 1) * TT + tp[i]] = U[i, (c_pos[i] - 1) * TT + tp[i]] - 1
		U[i, (c_pos[i] - 1) * TT + tr[i]] = U[i, (c_pos[i] - 1) * TT + tr[i]] + 1
	}
	if (rows(U) >= cols(U)) {
		svd(U, Uu, s0, Vt)
		v  = Vt'
	}
	else {
		svd(U', v, s0, Vt)
		Uu = Vt'
	}
	d    = colshape(s0, 1)
	uQ   = Uu' * QS
	uDe  = Uu' * Delta
	rank = sum(d :> max(d) * epsilon(1) * max((rows(U), cols(U))))
	U    = J(0, 0, .)
	Uu   = J(0, 0, .)

	K      = rows(P.treated)
	offset = J(K, 1, 0)
	for (k = 2; k <= K; k++) offset[k] = offset[k - 1] + (T_max - P.treated[k - 1] + 1)
	cell_id = J(P.n_obs, 1, .)
	for (r = 1; r <= P.n_obs; r++) {
		c = P.cohort_of[P.unit_row[r]]
		if (c > T_min & c <= T_max) {
			kk = _sgmm_match(c, P.treated)
			if (P.time_idx[r] + T_min - 1 >= c) {
				cell_id[r] = offset[kk] + (P.time_idx[r] + T_min - 1 - c) + 1
			}
		}
	}

	converged   = 0
	solve_ok    = 0
	termination = "max_iter_reached"
	n_iter      = 0
	max_change  = .
	n_kept      = .
	nzl         = 0
	QtAQ_out    = J(0, 0, .)
	W_last      = J(0, 0, .)
	QtAQ        = J(0, 0, .)
	QtAD        = J(0, 1, .)
	beta_new    = J(0, 1, .)
	nzl_it      = 0
	nk_it       = .
	V_id        = J(0, 0, .)

	for (it = 1; it <= max_iter; it++) {
		n_iter   = it
		beta_old = beta

		beta_full = J(nb, 1, .)
		beta_full[id_idx] = beta
		tau = J(P.n_obs, 1, 0)
		for (r = 1; r <= P.n_obs; r++) {
			if (cell_id[r] < .) tau[r] = beta_full[cell_id[r]]
		}
		yadj = P.y - tau

		R      = _sgmm_fe_resid(yadj, P.unit_row, P.time_idx, P.N_units, TT)
		W_last = _sgmm_weight(R, P.obs, P.grp_of, P.in_design, P.all_groups,
		                      weighting, nzl_it)
		nzl    = nzl_it

		solved = _sgmm_reduced_solve(d, v, uQ, uDe, W_last, 1e-8, QtAQ, QtAD, nk_it)
		if (!solved) {
			termination = "singular_weight_matrix"
			break
		}
		n_kept = nk_it

		if (!_sgmm_rsolve(QtAQ, QtAD, beta_new)) {
			termination = "singular_normal_equations"
			break
		}
		if (hasmissing(beta_new)) {
			termination = "non_finite_update"
			break
		}

		solve_ok   = 1
		QtAQ_out   = QtAQ
		max_change = max(abs(beta_new - beta_old))
		beta       = beta_new

		if (max_change < . & max_change < tol) {
			converged   = 1
			termination = "converged"
			break
		}
	}

	if (solve_ok) {
		if (!_sgmm_rsolve(QtAQ_out, J(0, 0, .), V_id)) V_id = _sgmm_ginv(QtAQ_out)
	}
	else {
		M    = _sgmm_reduced_M(d, v, W_last)
		V_id = (qtq :* (uQ' * M * uQ)) :* qtq'
	}
	V_id = (V_id + V_id') :/ 2

	if (!solve_ok) {
		_sgmm_warn("No weighted GMM step could be completed (termination: " +
			termination + "). The reported estimates are the identity-weighted seed, " +
			"not a " + weighting + "-weighted estimator, and standard errors use the " +
			"sandwich fallback. " +
			(n_kept < . ? "The estimated weighting matrix supports " + _sgmm_num(n_kept) +
			  " independent moment direction(s), but " + _sgmm_num(nid) +
			  " effect(s) must be identified. " : "") +
			"A cohort-specific covariance carries more parameters than a short panel " +
			"with few units per cohort can support; a more parsimonious weighting(), " +
			"or more units, is needed for this design.")
	}
	else if (!converged & weighting != "unrestricted") {
		_sgmm_warn("The efficient GMM reweighting did not converge (termination: " +
			termination + ") after " + _sgmm_num(n_iter) + " of at most " +
			_sgmm_num(max_iter) + " iteration(s), with a final change of " +
			strtrim(sprintf("%9.2e", max_change)) + ". Estimates are the last valid " +
			"iterate. Consider a larger maxiter(), a looser tolerance(), or a more " +
			"parsimonious weighting().")
	}

	beta_full = J(nb, 1, .)
	beta_full[id_idx] = beta
	se_full = J(nb, 1, .)
	se_full[id_idx] = sqrt(rowmax((J(nid, 1, 0), diagonal(V_id))))

	cpos     = _sgmm_match(C[., 1], P.all_groups)
	Ng_of    = P.N_g[cpos]
	Tg_of    = T_max :- C[., 1] :+ 1
	Ng_treat = P.N_g[_sgmm_match(P.treated, P.all_groups)]
	Tg_treat = T_max :- P.treated :+ 1
	denom    = sum(Ng_treat :* Tg_treat)
	wcw = Ng_of :/ denom
	wcw = wcw :/ sum(wcw)
	wew = (1 / K) :/ Tg_of
	wew = wew :/ sum(wew)

	aggcw = J(1, 0, .)
	aggew = J(1, 0, .)
	_sgmm_agg(wcw, beta, V_id, id_idx, nb, aggcw)
	_sgmm_agg(wew, beta, V_id, id_idx, nb, aggew)

	bstripe = J(nid, 2, "")
	for (i = 1; i <= nid; i++) {
		bstripe[i, 1] = "g" + _sgmm_slab(C[id_idx[i], 1])
		bstripe[i, 2] = "t" + _sgmm_slab(C[id_idx[i], 2])
	}
	st_matrix(bname, beta')
	st_matrixcolstripe(bname, bstripe)
	st_matrix(Vname, V_id)
	st_matrixrowstripe(Vname, bstripe)
	st_matrixcolstripe(Vname, bstripe)

	catt = (C, C[., 2] - C[., 1], beta_full, se_full, ident, ncomp)
	st_matrix(cattname, catt)
	cstripe = J(nb, 2, "")
	for (i = 1; i <= nb; i++) {
		cstripe[i, 1] = "g" + _sgmm_slab(C[i, 1])
		cstripe[i, 2] = "t" + _sgmm_slab(C[i, 2])
	}
	st_matrixrowstripe(cattname, cstripe)
	st_matrixcolstripe(cattname, (J(7, 1, ""), ("g" \ "t" \ "event_time" \ "estimate" \
		"std_error" \ "identified" \ "n_comparisons")))

	st_matrix(wcwname, wcw')
	st_matrixcolstripe(wcwname, cstripe)
	st_matrix(wewname, wew')
	st_matrixcolstripe(wewname, cstripe)

	st_matrix(Ngname, P.N_g')
	st_matrixcolstripe(Ngname, (J(G, 1, ""), ("g" :+ _sgmm_slab(P.all_groups))))
	st_matrix(cohname, P.treated')
	st_matrixcolstripe(cohname, (J(K, 1, ""), ("g" :+ _sgmm_slab(P.treated))))
	st_matrix(Wname, W_last)
	st_matrix(ymname, ym)
	st_matrix(Ngtname, P.N_gt)
	st_matrix(grpname, P.all_groups')
	st_matrixcolstripe(grpname, (J(G, 1, ""), ("g" :+ _sgmm_slab(P.all_groups))))

	scn = ("N_rows", "N_units", "T_periods", "T_min", "T_max", "n_cohorts", "N_beta",
	       "N_identified", "N_2x2", "n_dropped_2x2", "moment_rank", "n_kept",
	       "n_zero_lag", "converged", "solve_ok", "n_iter", "max_change", "tol",
	       "max_iter", "has_never", "n_never_zero", "n_never_late", "cov_n_fallback",
	       "balanced", "n_missing_cells",
	       "att_cw", "se_cw", "att_cw_id", "se_cw_id", "share_cw",
	       "att_ew", "se_ew", "att_ew_id", "se_ew_id", "share_ew")
	scv = (P.n_obs, P.N_units, TT, T_min, T_max, K, nb,
	       nid, n2, n_dropped, rank, n_kept,
	       nzl, converged, solve_ok, n_iter, max_change, tol,
	       max_iter, P.has_never, P.n_never_zero, P.n_never_late,
	       cov_fb,
	       (P.n_missing_cells == 0), P.n_missing_cells,
	       aggcw, aggew)
	st_matrix(scname, scv)
	st_matrixcolstripe(scname, (J(cols(scn), 1, ""), scn'))

	st_local("termination", termination)
}

end
