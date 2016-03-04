/** 
 * LSLOGIT - Prediction routine
 * 
 * @package lslogit
 */


/**
 * Wrapper program for prediction after lslogit estimation
 *
 * @param `pc1'      Predict choice probabilities
 * @param `xb'       Predict systematic utility
 * @param `dudes'    Predict approx. probability of negative marginal utility
 * @param `wages'    Predict wage rates (only for joint estimation)
 * @param `increase' Increase estimated wage rate by x % before prediction
 */
program define lslogit_p, rclass
    version 12
    syntax newvarlist(min=1 max=4) [if] [in] [, pc1 xb DUdes WAGEs ///
                                                INCrease(numlist min=1 max=1)]


    //
    // Initialize
    //

    local finished = 0

    // Predict command works only with lslogit estimates
    if ("`e(cmd)'" != "lslogit") error 301

    // What to predict?
    local n_leisure : word count `e(leisure)'
    local n_wincrea : word count `increase'
    local n_varlist : word count `varlist'
    local n_optlist : word count `pc1' `xb' `dudes' `wages'
    if (("`wages'" != "" | "`increase'" != "") & "`e(wagevars)'" == "") {
        di in r "wage options work only after joint estimation"
        exit 498
    }
    else if ("`wages'" != "" & "`increase'" != "") {
        di in r "No, that's not gonna happen. I can either predict wages " ///
                "or increase wages, it's up to you."
        exit 498
    }
    else if (`n_wincrea' > 0 & `n_wincrea' != `n_leisure') {
        di in r "I've got `n_leisure' guys but `n_wincrea' wage " ///
                "increases!? Seriously?"
        exit 498
    }
    else if (`n_varlist' == `n_optlist') {
        local opt `pc1' `xb' `dudes' `wages'
        di as text "(options '`opt'' chosen)"
    }
    else if (`n_varlist' == 1 & `n_optlist' == 0) {
        local opt pc1
        di as text "(option '`opt'' assumed)"
    }
    else {
        di in r "number of variables does not match number of options"
        exit 498
    }

    // Mark the estimation sample
    marksample touse, novarlist
    markout `touse' `e(depvar)' `e(group)' `e(consum)' `e(leisure)' `e(cx)' ///
                    `e(lx1)' `e(lx2)' `e(c2x)' `e(l2x1)' `e(l2x2)'          ///
                    `e(indeps)' `e(wagepred)' `e(days)' `e(tria1)'          ///
                    `e(tria2)' `e(wagevars)'

    // Sort sample by households / BUGGY / DEBUG
    qui sort `e(group)'

    // Restrict sample to those of interest
    preserve
    qui keep if `touse'
    qui count
    local nobs = r(N)


    //
    // Prepare data
    //

    local n_cxias    : word count `e(cx)'
    local n_lx1ias   : word count `e(lx1)'
    local n_lx2ias   : word count `e(lx2)'
    local n_c2xias   : word count `e(c2x)'
    local n_l2x1ias  : word count `e(l2x1)'
    local n_l2x2ias  : word count `e(l2x2)'
    local n_indeps   : word count `e(indeps)'
    local n_randvars : word count `e(randvars)'
    local n_wagep    : word count `e(wagepred)'
    local n_wagesig  : word count `e(wagesigma)'
    local n_hwage    : word count `e(hwage)'
    local n_taxrias1 : word count `e(taxreg_ia1)'
    local n_taxrias2 : word count `e(taxreg_ia2)'
    local n_wagevars : word count `e(wagevars)'
    local rvars      : word count `e(randvars)'


    //
    // Move data to Mata
    //

    // Coefficient vector
    mata: lsl_B     = st_matrix("e(b)")

    // Dependent variable
    mata: lsl_Y     = st_data(., "`e(depvar)'")

    // Functional form of utility function
    mata: lsl_ufunc = "`e(ufunc)'"

    // Number of individuals in household
    mata: lsl_nlei  = `n_leisure'

    // Number of random coefficients
    mata: lsl_joint = ("`e(joint)'" == "joint")

    // Number of random coefficients
    mata: lsl_rvars = `rvars'

    // Random coefficients
    mata: lsl_Rvars = ("`e(randvars)'" != ""                   ///
                         ? strtoreal(tokens("`e(randvars)'"))' ///
                         : J(0, 1, 0))

    // Correlation?
    mata: lsl_corr  = ("`e(randvars)'" != "" & "`e(corr)'" == "1" ? 1 : 0)

    // Number of draws
    mata: lsl_draws = strtoreal("`e(draws)'")

    // Number of draws to burn
    mata: lsl_burn  = strtoreal("`e(burn)'")

    // Hourly wage rates
    mata: lsl_Hwage = ("`e(hwage)'" != ""                   ///
                         ? st_data(., tokens("`e(hwage)'")) ///
                         : J(`nobs', lsl_nlei, 0))

    // Wage observed?
    mata: lsl_Wobs  = (lsl_Y :* (log(lsl_Hwage) :< .))

    // Run Wage Prediction
    mata: lsl_wagep = `e(wagep)'
    mata: lsl_Sigma = ("`e(wagesigma)'" != ""                  ///
                         ? strtoreal(tokens("`e(wagesigma)'")) ///
                         : J(1, `n_leisure', 0))

    // To round, or not to round?
    mata: lsl_round = ("`e(round)'" == "1")

    // Wage/preference correlation?
    mata: lsl_wagecorr    = ("`e(wagecorr)'" != "" ///
                               ? strtoreal("`e(wagecorr)'") : 0)

    // Wage/preference correlation anchor?
    mata: lsl_residanchor = ("`e(residanchor)'" == "1")

    // Choices per group
    mata: lsl_J      = panelsetup(st_data(., "`e(group)'"), 1)

    // Number of groups
    mata: lsl_groups = rows(lsl_J)

    // Dummies enabling or disabling the wage prediction
    mata: lsl_Wpred  = (lsl_wagep ? (strtrim("`e(wagepred)'") != ""  ///
                                       ? st_data(., "`e(wagepred)'") ///
                                       : J(`nobs', lsl_nlei, 1))     ///
                                  : J(`nobs', lsl_nlei, 0))

    // Days of taxyear
    mata: lsl_Days   = ("`e(days)'" != "" ? st_data(., "`e(days)'") ///
                                          : J(`nobs', 1, 365))

    // Weights?
    mata: lsl_Weight = ("`e(weight)'" != "" ? st_data(., "`e(weight)'") ///
                                            : J(`nobs', 1, 1))

    //
    // Build coefficients vector
    //

    mata: lsl_bwage = (`n_wagevars' != 0) + `n_wagevars'
    mata: lsl_bheck = (lsl_wagecorr > 0) + lsl_wagecorr
    mata: lsl_brnd  = (lsl_corr == 1 ? lsl_rvars * (lsl_rvars + 1) / 2 ///
                                     : lsl_rvars)
    mata: lsl_blam  = (lsl_ufunc == "boxcox" ? 1 + lsl_nlei : 0)
    mata: lsl_bfix  = `n_cxias' + 1 + (`n_c2xias' + 1) *                ///
                      (lsl_ufunc != "boxcox") + lsl_nlei + `n_lx1ias' + ///
                      1 + (`n_l2x1ias' + 1) * (lsl_ufunc != "boxcox") + ///
                      (`n_lx2ias' + 1 + (`n_l2x2ias' + 1) *             ///
                      (lsl_ufunc != "boxcox") + 1) * (lsl_nlei == 2) +  ///
                      `n_indeps'

    mata: lsl_iwage = 1         + lsl_bfix
    mata: lsl_iheck = lsl_iwage + lsl_bwage
    mata: lsl_ilam  = lsl_iheck + lsl_bheck
    mata: lsl_irnd  = lsl_ilam  + lsl_blam

    mata: lsl_iC  = `n_cxias' + 1
    mata: lsl_iL1 = `n_cxias' + 1 + `n_c2xias' + 1 + lsl_nlei + `n_lx1ias' + 1
    mata: lsl_iL2 = lsl_iL1 + `n_l2x1ias' + 1 + `n_lx2ias' + 1

    mata: lsl_Bfix  = lsl_B[|1\lsl_bfix|]

    // Wage coefficients
    mata: lsl_Bwage = (lsl_bwage > 0                                    ///
                         ? lsl_B[|lsl_iwage\lsl_iwage + lsl_bwage - 1|] ///
                         : J(0, 0, 0))

    //   Cholesky factor of wage variance
    mata: lsl_Bsig  = (lsl_bheck > 0 ? lsl_B[lsl_iheck] : J(0, 0, 0))

    // Get auxiliary random coefficients
    mata: lsl_Brnd  = (lsl_rvars > 0                                 ///
                         ? lsl_B[|lsl_irnd\lsl_irnd + lsl_brnd - 1|] ///
                         : J(0, 0, 0))

    // Build variance-(covariance) matrix
    mata: lsl_CholB = (lsl_corr == 1 ? lowertriangle(invvech(lsl_Brnd')) ///
                                     : diag(lsl_Brnd'))


    //
    // Tax regression
    //

    if ("`e(taxreg)'" == "1") {
        // Tax regression estimates
        mata: lsl_TaxregB     = st_matrix("e(taxreg_b)")

        // Interaction variables on Mwage1 and Mwage1^2
        mata: lsl_TaxregIas1  = ("`e(taxreg_ia1)'" != ""                   ///
                                   ? st_data(., tokens("`e(taxreg_ia1)'")) ///
                                   : J(`nobs', 0, 0))

        // Interaction variables on Mwage2 and Mwage2^2
        mata: lsl_TaxregIas2  = ("`e(taxreg_ia2)'" != ""                   ///
                                   ? st_data(., tokens("`e(taxreg_ia2)'")) ///
                                   : J(`nobs', 0, 0))

        // Variables that are independent of m_wage
        mata: lsl_TaxregVars  = st_data(., tokens("`e(taxreg_vars)'"))

        // Root Mean Squared error of tax regression
        mata: lsl_taxreg_rmse = strtoreal("`e(taxreg_rmse)'")
    }

    // Box-Cox normalizing constants
    if ("`e(ufunc)'" == "boxcox" & "`e(boxcc)'" != "" & "`e(boxcl)'" != "") {
        mata: lsl_boxcc = (`e(boxcc)' > 0 & `e(boxcc)' < . ? `e(boxcc)' : 1)
        mata: lsl_boxcl = (`e(boxcl)' > 0 & `e(boxcl)' < . ? `e(boxcl)' : 1)
    }
    else mata: lsl_boxcc = lsl_boxcl = 1

    // Consumption, squared and interactions
    mata: lsl_C   = st_data(., "`e(consum)'") :/ lsl_boxcc
    mata: lsl_CX  = ((`n_cxias'  > 0 ? st_data(., tokens("`e(cx)'"))  ///
                                     : J(`nobs', 0, 0)), J(`nobs', 1, 1))
    mata: lsl_C2X = ((`n_c2xias' > 0 ? st_data(., tokens("`e(c2x)'")) ///
                                     : J(`nobs', 0, 0)),              ///
                     J(`nobs', (lsl_ufunc != "boxcox" ? 1 : 0), 1))

    // Leisure, squared and interactions
    forval i = 1/2 {
        local var : word `i' of `e(leisure)'
        if ("`var'" != "") {
            mata: lsl_L`i'   = st_data(., "`var'") :/ lsl_boxcl
            mata: lsl_LX`i'  = ((`n_lx`i'ias'  > 0                    ///
                                   ? st_data(., tokens("`e(lx`i')'")) ///
                                   : J(`nobs', 0, 0)), J(`nobs', 1, 1))
            mata: lsl_L2X`i' = ((`n_l2x`i'ias' > 0                     ///
                                   ? st_data(., tokens("`e(l2x`i')'")) ///
                                   : J(`nobs', 0, 0)),                 ///
                                J(`nobs', (lsl_ufunc != "boxcox" ? 1 : 0), 1))
        }
        else {
            mata: lsl_L`i'   = J(`nobs', 0, 0)
            mata: lsl_LX`i'  = J(`nobs', 0, 0)
            mata: lsl_L2X`i' = J(`nobs', 0, 0)
        }
    }

    // Hypothetical hours of work (caution: undo Box-Cox normalization!)
    mata: lsl_Hours = `e(totaltime)' :- (lsl_L1, lsl_L2) :* lsl_boxcl

    // Dummy variables
    mata: lsl_Xind = (`n_indeps' > 0 ? st_data(., tokens("`e(indeps)'")) ///
                                     : J(`nobs', 0, 0))

    //
    // Build right hand side
    //

    // Translog utility
    if ("`e(ufunc)'" == "tran") {
        mata: lsl_X = (log(lsl_C)  :* (lsl_CX,  log(lsl_C)  :* lsl_C2X,   ///
                                       log(lsl_L1), log(lsl_L2)),         ///
                       log(lsl_L1) :* (lsl_LX1, log(lsl_L1) :* lsl_L2X1), ///
                       log(lsl_L2) :* (lsl_LX2, log(lsl_L2) :* lsl_L2X2), ///
                       log(lsl_L1) :* log(lsl_L2), lsl_Xind)
    }
    // Quadratic utility
    else if ("`e(ufunc)'" == "quad") {
        mata: lsl_X = (lsl_C  :* (lsl_CX,  lsl_C  :* lsl_C2X,   ///
                                  lsl_L1,  lsl_L2),             ///
                       lsl_L1 :* (lsl_LX1, lsl_L1 :* lsl_L2X1), ///
                       lsl_L2 :* (lsl_LX2, lsl_L2 :* lsl_L2X2), ///
                       lsl_L1 :* lsl_L2, lsl_Xind)
    }
    // Box-Cox utility
    else if ("`e(ufunc)'" == "boxcox") {
        mata: lsl_lC  = lsl_B[lsl_ilam]
        mata: lsl_lL1 = lsl_B[lsl_ilam + 1]
        mata: lsl_lL2 = (lsl_nlei == 2 ? lsl_B[lsl_ilam + 2] : 0)
        mata: lsl_B
        mata: lsl_Bfix
        mata: "lsl_lC", "lsl_lL1", "lsl_lL2"
        mata: lsl_lC, lsl_lL1, lsl_lL2
        mata: lsl_X = (lsl_boxcox(lsl_C, lsl_lC)   :*           ///
                        (lsl_CX, lsl_boxcox(lsl_L1, lsl_lL1),   ///
                                 lsl_boxcox(lsl_L2, lsl_lL2)),  ///
                       lsl_boxcox(lsl_L1, lsl_lL1) :*  lsl_LX1, ///
                       lsl_boxcox(lsl_L2, lsl_lL2) :*  lsl_LX2, ///
                       lsl_boxcox(lsl_L1, lsl_lL1) :*           ///
                        lsl_boxcox(lsl_L2, lsl_lL2), lsl_Xind)
    }

    // Generate Halton sequences
    local R_rvars = `rvars' + `e(wagep)' * `n_leisure' + ///
                    !inlist("`e(taxreg_rmse)'", "", ".")
    mata: lsl_R = (`R_rvars' > 0                                  ///
                     ? invnormal(halton(lsl_groups * lsl_draws,   ///
                                        `R_rvars', 1 + lsl_burn)) ///
                     : J(`nobs', 0, 0))


    //
    // Wage prediction? (Buggy! Works only for single decision maker housholds!)
    //

    mata: lsl_CholBW = J(lsl_wagep * lsl_nlei, lsl_rvars, 0)
    if (inlist("`e(joint)'", "joint", "1")) {
        // Predict wages
        mata: lsl_LnWageHat = cross((st_data(., tokens("`e(wagevars)'")), ///
                                     J(`nobs', 1, 1))', lsl_Bwage')

        // Get wage equation residuals
        mata: lsl_LnWres = lsl_Wobs :* (log(lsl_Hwage) :- lsl_LnWageHat)

        // Root MSE of wage equation if not estimated simultaneously
        if (inlist("`e(wagecorr)'", "", "0")) {
            mata: lsl_Bsig = sqrt(cross(lsl_LnWres, lsl_LnWres) / ///
                                  (colsum(lsl_Wobs) - lsl_bwage))
        }

        // Build lower line of Cholesky matrix
        if ("`e(wcorrvars)'" == "`e(randvars)'") {
            mata: lsl_CholBW = (lsl_wagecorr > 0                                     ///
                                  ? lsl_B[|lsl_iheck + 1\lsl_iheck + lsl_bheck - 1|] ///
                                  : J(1, 0, 0))
        }
        else {
            // Cholesky factors between coefficients and wages
            // Correlation structure different from random coefficients
            mata: lsl_CholBW = J(1, 0, 0)
            local v = 1
            local fullcorrlist `e(randvars)'
            foreach vt in `e(wcorrvars)' {
                foreach wt in `fullcorrlist' {
                    if ("`vt'" == "`wt'") {
                        mata: lsl_CholBW = lsl_CholBW, lsl_B[lsl_iheck + `v']
                        local fullcorrlist = subinstr(" `fullcorrlist' ", " `vt' ", " ", .)
                        continue, break
                    }
                    else mata: lsl_CholBW = lsl_CholBW, 0
                }
                local v = `v' + 1
            }
            mata: lsl_CholBW = lsl_CholBW, J(1, lsl_rvars - cols(lsl_CholBW), 0)
        }
        mata: lsl_CholW = diag(lsl_Bsig)

        // Calculate standard error of wage equation
        mata: lsl_SigmaW = (lsl_wagecorr ? sqrt(rowsum(lsl_CholW:^2))' ///
                                         : lsl_Bsig)

        // Store error
        mata: st_numscalar("r(sigma_w1)", lsl_SigmaW)
        return scalar sigma_w1 = `r(sigma_w1)'

        // Replace hourly wage rates
        mata: lsl_Hwage = exp(lsl_LnWageHat :+ lsl_SigmaW:^2/2) :* ///
                          (lsl_Hours[|1,1\.,1|] :!= 0)
    }
    else mata: lsl_CholW = diag(lsl_Sigma)

    // Increase wages
    if (`n_wincrea' > 0) {
        forval i = 1/`n_wincrea' {
            local fac : word `i' of `increase'
            mata: lsl_Hwage[.,`i'] = `fac' :* lsl_Hwage[.,`i']
        }
    }


    // Restore sample
    restore


    //
    // Predict
    //

    // Generate new variables
    foreach v of local varlist {
        qui gen double `v' = .
    }

    // Wage prediction (Buggy! Works only for single decision maker housholds!)
    if ("`wages'" != "") {
        // Get name of new variable
        mata: st_local("newvar", invtokens((tokens("`opt'") :== "wages") :* ///
                                           tokens("`varlist'")))

        // Store predicted wages
        mata: st_store(., "`newvar'", "`touse'", lsl_Hwage)

        // I'm done, drop wage variable from varlist
        mata: st_local("varlist", invtokens((tokens("`opt'") :!= "wages") :* ///
                                            tokens("`varlist'")))
    }

    // Run evaluator and predict pc1/xb/dudes/...
    if (!inlist(trim("`opt'"), "", "wages")) {
        // Print estimated log-likelihood
        if ("`pc1'" != "" | "`opt'" == "pc1") {
            di as text "Estimated ll=`=round(e(ll), 0.0001)'. Now ... " _c
        }

        // Run prediction
        mata: lslogit_p(tokens("`varlist'"), "`touse'", tokens("`opt'"))

        // Store estimated and predicted likelihood
        if ("`pc1'" != "" | "`opt'" == "pc1") {
            return scalar ll_e = `e(ll)'
            return scalar ll_p = lsl_ll_p
        }
    }


    //
    // Clean up
    //

    /*
    foreach m in round ufunc Weight Y Hwage boxcc boxcl C CX C2X ///
                 L1 LX1 L2X1 L2 LX2 L2X2 Xind X joint            ///
                 WageVars Days Hours wagep Wpred Sigma TaxregB   ///
                 TaxregIas1 TaxregIas2 TaxregVars                ///
                 groups J draws burn Rvars corr R lC lL1 lL2 B   ///
                 Bfix bfix rvars {
        cap mata mata drop lsl_`m'
    }
    */
end


***
