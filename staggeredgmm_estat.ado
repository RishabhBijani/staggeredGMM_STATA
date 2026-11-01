*! version 0.1.0  30oct2026  Rishabh Bijani, Parush Arora
*! Postestimation for staggeredgmm: estat jtest, estat aggregate

program staggeredgmm_estat, rclass
	version 16
	if "`e(cmd)'" != "staggeredgmm" error 301

	gettoken sub rest : 0, parse(" ,")
	local lsub = strlen(`"`sub'"')
	if `lsub' >= 2 & `"`sub'"' == substr("jtest", 1, `lsub') {
		_sgmm_jtest `rest'
		return add
	}
	else if `lsub' >= 3 & `"`sub'"' == substr("aggregate", 1, `lsub') {
		_sgmm_aggregate `rest'
		return add
	}
	else {
		estat_default `0'
		return add
	}
end


program _sgmm_jtest, rclass
	syntax [, FULL Window(real 3)]

	if `window' >= . | `window' < 1 | `window' != floor(`window') {
		di as err "window() must be a single whole number of at least 1."
		exit 198
	}
	if e(has_never) != 1 {
		di as err "The specification test requires a never-treated control group: the "  ///
			"event-study basis it is built on differences every cohort against the same " ///
			"never-treated reference. This design has none."
		exit 459
	}
	capture confirm matrix e(sgmm_W)
	if _rc {
		di as err "The estimation results carry no weighting matrix, so the covariance " ///
			"of the placebo moments cannot be formed."
		exit 459
	}

	local type = cond("`full'" != "", "full", "local")
	tempname J df p base nmom ndrop mom
	mata: _sgmmj_run("`type'", `window', "`J'", "`df'", "`p'", "`base'", "`nmom'", "`ndrop'", "`mom'")

	local w "`e(weighting)'"
	if "`w'" == "pooled_toeplitz"      local wl "pooled Toeplitz (GMM-T)"
	else if "`w'" == "cohort_toeplitz" local wl "cohort-specific Toeplitz (GMM-HT)"
	else if "`w'" == "unrestricted"    local wl "unrestricted within-cohort (GMM-U)"
	else local wl "`w'"

	di _n as txt "Specification test for parallel trends and no anticipation"
	di as txt "Hansen J on the pre-treatment placebo restrictions" _n
	if "`type'" == "local" {
		di as txt "  Restriction set   " as res "local pre-window (`window' period(s) before adoption)"
	}
	else {
		di as txt "  Restriction set   " as res "all pre-treatment placebo moments"
	}
	di as txt "  Base period       " as res %1.0f scalar(`base')
	di as txt "  Weighting         " as res "`wl'"
	di as txt "  Moments           " as res %1.0f scalar(`nmom')
	di ""
	di as txt "  J = " as res %9.4f scalar(`J') as txt ", df = " as res %1.0f scalar(`df') ///
		as txt ", p = " as res %6.4f scalar(`p')

	local scope "unconditional"
	if `"`e(covariates)'"' != "" {
		di ""
		di as txt "{hline 78}"
		di as err "Important: " as txt "this model adjusts for covariates, but this test does not."
		di as txt "{p 0 0 2}The placebo moments are unadjusted group means, so the test checks " ///
			"{bf:unconditional} parallel trends and no anticipation. It does {bf:not} test the " ///
			"covariate-conditional parallel trends assumption under which the effects were " ///
			"estimated. A rejection can reflect differences in covariates across cohorts rather " ///
			"than a failure of that assumption, and a non-rejection does not validate it.{p_end}"
		di as txt "{hline 78}"
	}

	return scalar J           = scalar(`J')
	return scalar df          = scalar(`df')
	return scalar p           = scalar(`p')
	return scalar base_period = scalar(`base')
	return scalar n_moments   = scalar(`nmom')
	return scalar n_dropped   = scalar(`ndrop')
	if "`type'" == "local" return scalar window = `window'
	return local  type        "`type'"
	return local  scope       "`scope'"
	return local  weighting   "`e(weighting)'"
	return matrix moments     = `mom'
end


program _sgmm_aggregate, rclass
	syntax [, POST Level(cilevel)]

	if "`e(subcmd)'" == "aggregate" {
		if "`post'" != "" {
			di as err "The aggregate ATTs are already posted in e(b); rerun staggeredgmm first."
			exit 498
		}
		ereturn display, level(`level')
		tempname bP VP
		matrix `bP' = e(b)
		matrix `VP' = e(V)
		return matrix b = `bP'
		return matrix V = `VP'
		return scalar share_cw = e(share_cw)
		return scalar share_ew = e(share_ew)
		exit
	}
	foreach m in b V catt w_cw w_ew {
		capture confirm matrix e(`m')
		if _rc {
			di as err "e(`m') not found; estat aggregate needs the results of staggeredgmm."
			exit 301
		}
	}

	tempname bA VA sh
	mata: _sgmma_run("`bA'", "`VA'", "`sh'")
	local share_cw = `sh'[1, 1]
	local share_ew = `sh'[1, 2]

	di _n as txt "Aggregate average treatment effects on the treated"
	di as txt "  att_cw  treated-observation ATT" ///
		cond(`share_cw' < 1 - 1e-12, " (identified cells only: " + string(100 * `share_cw', "%4.1f") + "% of its weight)", "")
	di as txt "  att_ew  cohort-equal ATT" ///
		cond(`share_ew' < 1 - 1e-12, " (identified cells only: " + string(100 * `share_ew', "%4.1f") + "% of its weight)", "")

	if "`post'" == "" {
		_sgmm_aggshow `bA' `VA', level(`level')
	}
	else {
		_sgmm_aggpost `bA' `VA', level(`level')
	}

	return matrix b = `bA'
	return matrix V = `VA'
	return scalar share_cw = `share_cw'
	return scalar share_ew = `share_ew'
end


* Display the aggregates as a coefficient table, leaving e() unchanged.
program _sgmm_aggshow, eclass
	syntax anything(name=mats), Level(cilevel)
	gettoken bA VA : mats
	tempname hold b V
	matrix `b' = `bA'
	matrix `V' = `VA'
	_estimates hold `hold', copy restore
	ereturn post `b' `V'
	ereturn display, level(`level')
end


* Replace e(b) and e(V) by the aggregates, keeping every other result.
program _sgmm_aggpost, eclass
	syntax anything(name=mats), Level(cilevel)
	gettoken bA VA : mats
	tempname b V
	matrix `b' = `bA'
	matrix `V' = `VA'

	local depvar `"`e(depvar)'"'
	local nobs = e(N)
	tempvar es
	quietly gen byte `es' = e(sample)
	local scal : e(scalars)
	local macs : e(macros)
	local mats : e(matrices)
	foreach x of local scal {
		tempname s_`x'
		scalar `s_`x'' = e(`x')
	}
	local keepm
	foreach x of local mats {
		if !inlist("`x'", "b", "V") {
			tempname m_`x'
			matrix `m_`x'' = e(`x')
			local keepm `keepm' `x'
		}
	}
	foreach x of local macs {
		local c_`x' `"`e(`x')'"'
	}

	ereturn post `b' `V', esample(`es') depname(`depvar') obs(`nobs')
	foreach x of local scal {
		if !inlist("`x'", "N", "rank") ereturn scalar `x' = scalar(`s_`x'')
	}
	foreach x of local keepm {
		ereturn matrix `x' = `m_`x''
	}
	foreach x of local macs {
		if "`x'" != "cmd" ereturn local `x' `"`c_`x''"'
	}
	ereturn local subcmd "aggregate"
	ereturn local cmd `"`c_cmd'"'
	ereturn display, level(`level')
	di as txt "{p 0 6 2}e(b) and e(V) now hold the aggregate ATTs; all other results of " ///
		"staggeredgmm are kept. Rerun staggeredgmm to restore the cohort-by-time effects.{p_end}"
end

version 16
mata:

void _sgmmj_warn(string scalar msg)
{
	printf("{p 0 9 2}{err}Warning:{txt} %s{p_end}\n", msg)
}


void _sgmmj_abort(string scalar msg, real scalar rc)
{
	errprintf("%s\n", msg)
	exit(rc)
}


void _sgmma_run(string scalar bn, string scalar Vn, string scalar shn)
{
	real rowvector b
	real matrix    V, C, A, VA
	real colvector id, wc, we

	b  = st_matrix("e(b)")
	V  = st_matrix("e(V)")
	C  = st_matrix("e(catt)")
	wc = st_matrix("e(w_cw)")'
	we = st_matrix("e(w_ew)")'
	id = selectindex(C[., 6] :== 1)
	if (rows(id) != cols(b)) {
		_sgmmj_abort("e(b) does not match the identified cells of e(catt).", 498)
	}
	A  = ((wc[id] :/ sum(wc[id]))' \ (we[id] :/ sum(we[id]))')
	VA = A * V * A'
	VA = (VA + VA') :/ 2
	st_matrix(bn, (A * b')')
	st_matrixcolstripe(bn, (J(2, 1, ""), ("att_cw" \ "att_ew")))
	st_matrix(Vn, VA)
	st_matrixrowstripe(Vn, (J(2, 1, ""), ("att_cw" \ "att_ew")))
	st_matrixcolstripe(Vn, (J(2, 1, ""), ("att_cw" \ "att_ew")))
	st_matrix(shn, (sum(wc[id]), sum(we[id])))
}


string scalar _sgmmj_num(real scalar x)
{
	return(strtrim(strofreal(x, "%20.0g")))
}


void _sgmmj_run(string scalar type, real scalar window, string scalar Jn,
	string scalar dfn, string scalar pn, string scalar basen, string scalar nmomn,
	string scalar ndropn, string scalar momn)
{
	real matrix    ym, Ngt, W, Up, Om, X, Vk, mom
	real colvector groups, cohorts, mom_g, mom_t, tt, g_pos, t_idx, ok, p_hat
	real colvector lam, keep, proj
	real vector    L
	real scalar    T_min, TT, nt_pos, t0, t0_idx, k, g, s, n_mom, n_dropped
	real scalar    i, lmax, df, stat, pval, G, j

	ym      = st_matrix("e(sgmm_ymeans)")
	Ngt     = st_matrix("e(sgmm_Ngt)")
	W       = st_matrix("e(sgmm_W)")
	groups  = st_matrix("e(sgmm_groups)")'
	cohorts = st_matrix("e(cohorts)")'
	T_min   = st_numscalar("e(T_min)")
	TT      = st_numscalar("e(T_periods)")
	G       = rows(groups)

	nt_pos = 0
	for (j = 1; j <= G; j++) {
		if (groups[j] == 0) nt_pos = j
	}
	if (nt_pos == 0) {
		_sgmmj_abort("The specification test requires a never-treated control group.", 459)
	}

	t0 = min(cohorts) - 1
	if (t0 < T_min) {
		_sgmmj_abort("The earliest treated cohort adopts in the first observed period, " +
			"so no common untreated base period exists and no placebo restriction can " +
			"be formed.", 459)
	}
	t0_idx = t0 - T_min + 1

	mom_g = J(0, 1, .)
	mom_t = J(0, 1, .)
	for (k = 1; k <= rows(cohorts); k++) {
		g  = cohorts[k]
		tt = J(0, 1, .)
		for (s = T_min; s <= g - 1; s++) {
			if (type == "local" & s < g - window) continue
			if (s == t0) continue
			tt = tt \ s
		}
		if (rows(tt) == 0) continue
		mom_g = mom_g \ J(rows(tt), 1, g)
		mom_t = mom_t \ tt
	}
	if (rows(mom_g) == 0) {
		_sgmmj_abort("No placebo restriction is available" +
			(type == "local" ? " with window(" + _sgmmj_num(window) +
			  "). Try a wider window or the full option." : "."), 459)
	}

	g_pos = J(rows(mom_g), 1, .)
	for (i = 1; i <= rows(mom_g); i++) {
		for (j = 1; j <= G; j++) {
			if (groups[j] == mom_g[i]) {
				g_pos[i] = j
				break
			}
		}
	}
	t_idx = mom_t :- T_min :+ 1

	ok = J(rows(mom_g), 1, 0)
	for (i = 1; i <= rows(mom_g); i++) {
		ok[i] = (Ngt[g_pos[i], t_idx[i]] > 0) & (Ngt[g_pos[i], t0_idx] > 0) &
		        (Ngt[nt_pos, t_idx[i]] > 0) & (Ngt[nt_pos, t0_idx] > 0)
	}
	n_dropped = sum(ok :== 0)
	if (n_dropped > 0) {
		keep = selectindex(ok)
		if (rows(keep) == 0) {
			_sgmmj_warn(_sgmmj_num(n_dropped) + " placebo moment(s) reference a group-period " +
				"with no observed unit and have been dropped from the test.")
			_sgmmj_abort("Every placebo moment references a group-period with no " +
				"observed unit; the test cannot be formed.", 459)
		}
		mom_g = mom_g[keep]
		mom_t = mom_t[keep]
		g_pos = g_pos[keep]
		t_idx = t_idx[keep]
		_sgmmj_warn(_sgmmj_num(n_dropped) + " placebo moment(s) reference a group-period " +
			"with no observed unit and have been dropped from the test.")
	}
	n_mom = rows(mom_g)

	p_hat = J(n_mom, 1, .)
	Up    = J(n_mom, G * TT, 0)
	for (i = 1; i <= n_mom; i++) {
		p_hat[i] = (ym[g_pos[i], t_idx[i]] - ym[g_pos[i], t0_idx]) -
		           (ym[nt_pos, t_idx[i]] - ym[nt_pos, t0_idx])
		Up[i, (g_pos[i] - 1) * TT + t_idx[i]]  = 1
		Up[i, (g_pos[i] - 1) * TT + t0_idx]    = -1
		Up[i, (nt_pos - 1) * TT + t_idx[i]]    = Up[i, (nt_pos - 1) * TT + t_idx[i]] - 1
		Up[i, (nt_pos - 1) * TT + t0_idx]      = Up[i, (nt_pos - 1) * TT + t0_idx] + 1
	}

	Om = Up * W * Up'
	Om = (Om + Om') :/ 2
	symeigensystem(Om, X, L)
	lam  = colshape(L, 1)
	lmax = max(lam)
	if (lmax >= . | lmax <= 0) {
		_sgmmj_abort("The estimated covariance of the placebo moments is not positive " +
			"definite; the test statistic cannot be formed.", 459)
	}
	keep = selectindex(lam :> 1e-8 * lmax)
	df   = rows(keep)
	if (df < n_mom) {
		_sgmmj_warn("The covariance of the " + _sgmmj_num(n_mom) + " placebo moment(s) " +
			"is rank " + _sgmmj_num(df) + "; the statistic uses a truncated inverse and " +
			"the degrees of freedom are reduced accordingly.")
	}
	Vk   = X[., keep]
	proj = Vk' * p_hat
	stat = sum((proj :^ 2) :/ lam[keep])
	pval = chi2tail(df, stat)

	mom = (mom_g, mom_t, mom_t - mom_g, p_hat)
	st_matrix(momn, mom)
	st_matrixcolstripe(momn, (J(4, 1, ""), ("g" \ "t" \ "event_time" \ "moment")))
	st_numscalar(Jn, stat)
	st_numscalar(dfn, df)
	st_numscalar(pn, pval)
	st_numscalar(basen, t0)
	st_numscalar(nmomn, n_mom)
	st_numscalar(ndropn, n_dropped)
}

end
