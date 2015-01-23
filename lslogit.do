*! version 0.4, 23jan2015, Max Loeffler <loeffler@zew.de>
/**
 * LSLOGIT - ESTIMATING MIXED LOGIT LABOR SUPPLY MODELS WITH STATA
 * 
 *
 * Copyright (C) 2012-2015 Max Löffler <loeffler@zew.de>
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 */


cap program drop lslogit
/**
 * Wrapper programm
 */
program define lslogit
    version 12
    syntax [varlist] [if] [in] [fweight/] [, cov *]

    if ("`cov'" != "") lslogit_cov `0'
    else if (replay()) {
        if ("`e(cmd)'" != "lslogit") error 301
        lslogit_Replay `0'
    }
    else lslogit_Estimate `0'
end


cap program drop lslogit_Replay
/**
 * Display estimation results and model notes
 *
 * @level `level' Singificance level of confindence interval
 * @quiet `quiet' Suppress additional notes on estimated model
 */
program define lslogit_Replay
    version 12
    syntax [, Level(integer `c(level)') Quiet]

    // Set up auxiliary stuff
    local diparm
    if ("`quiet'" == "") {
        // Auxiliary estimation results
        foreach aux in sigma_w1 sigma_w2 dudes {
            if (e(`aux') != .) {
                local val = e(`aux')
                local diparm `diparm' diparm(__lab__, value(`val') ///
                                                      label("[`aux']"))
            }
        }

        // Model description
        local footnote "Model:"
        if ("`e(ufunc)'" == "quad") {
            local footnote "`footnote' - Quadratic utility function"
        }
        else if ("`e(ufunc)'" == "tran") {
            local footnote "`footnote' - Translog utility function"
        }
        else if ("`e(ufunc)'" == "boxcox") {
            local footnote "`footnote' - Box-Cox utility function"
        }
        if ("`e(randvars)'" != "") {
            local footnote "`footnote'" _newline "       " ///
                           "- Random coefficients (`e(randvars)')"
        }
        if ("`e(corr)'" == "1") {
            local footnote "`footnote' with correlation"
        }
        if ("`e(wagep)'" == "1") {
            local footnote "`footnote'" _newline "       " ///
                           "- Wage prediction error integrated out"
        }
        if ("`e(draws)'" != "1") {
            local footnote "`footnote'" _newline "       " ///
                           "- Approximated using `e(draws)' Halton sequences"
        }
        if ("`e(wagevars)'" != "") {
            local footnote "`footnote'" _newline "       " ///
                           "- Joint labor supply and wage estimation"
        }
        if ("`e(wagecorr)'" != "") {
            local footnote "`footnote' with correlation (`e(wcorrvars)')"
        }
        if ("`e(residanchor)'" == "1") {
            local footnote "`footnote'" _newline "       " ///
                           "  (with residual anchor, using observed wages if possible)"
        }
        if ("`e(density)'" != "") {
            local footnote "`footnote'" _newline "       " ///
                           "- Accounted for choice density"
        }
    }

    // Display output
    ml display, level(`level') `diparm'
    if ("`footnote'" != "") di in gr "`footnote'"
end


cap program drop lslogit_Estimate
/**
 * Set up data and run estimation
 */
program define lslogit_Estimate, eclass
    version 12
    syntax varname(numeric) [if] [in] [fweight/],             ///
            GRoup(varname numeric)                            ///
            Consumption(varname numeric)                      ///
            Leisure(varlist numeric min=1 max=2)              ///
            [BOXCox QUADratic TRANslog                        ///
             cx(varlist numeric fv) c2x(varlist numeric fv)   ///
             lx1(varlist numeric fv) l2x1(varlist numeric fv) ///
             lx2(varlist numeric fv) l2x2(varlist numeric fv) ///
             INDeps(varlist fv)                               ///
             TOTALTime(integer 80) DAYs(varname numeric)      ///
             boxcc(integer 1000) boxcl(integer 80)            ///
             HWage(varlist numeric min=1 max=2)               ///
             TAXReg(name)                                     ///
             tria1(varlist numeric) tria2(varlist numeric)    ///
             WAGEPred(varlist numeric min=1 max=2)            ///
             WAGESIGma(numlist min=1 max=2)                   ///
             RANDvars(numlist ascending) corr                 ///
             DRaws(integer 50) burn(integer 15)               ///
             WAGEVars(varlist numeric fv)                     ///
             wagecorr(numlist ascending)                      ///
             noanchor TECHnique(string)                       ///
             round Quiet Verbose lambda(real 0)               ///
             force(varname numeric)                           ///
             difficult trace search(name)                     ///
             iterate(integer 100) method(name)                ///
             gradient hessian debug                           ///
             Level(integer `c(level)') from(string)           ///
             DENsity(varname numeric)                         ///
             WAGEDraws(varname numeric) RANDSAMple]

    /* INITIALIZE ESTIMATOR
     */

    // Mark the estimation sample
    marksample touse
    markout `touse' `varlist' `group' `consumption' `leisure'  ///
                    `cx' `lx1' `lx2' `c2x' `l2x1' `l2x2'       ///
                    `indeps' `wagepred' `days' `tria1' `tria2' ///
                    `wagevars' `density' `wagedraws'

    // Check for observations
    qui count if `touse'
    local nobs = r(N)
    if (`nobs' <= 0) {
        di in r "no observations"
        exit 2000
    }

    // Check dependent variable
    cap assert inlist(`varlist', 0, 1)
    if (_rc != 0) {
        di in r "choice variable has to equal zero or one"
        exit 498
    }

    // Verbose mode
    if ("`verbose'" == "") local qui qui

    // Validate Maximum Likelihood method
    if ("`method'" == "") local method d2
    if (!inlist("`method'", "d0", "d1", "d2")) {
        di in r "method must be either 'd0', 'd1' or 'd2'"
        exit 498
    }

    // Validate utility function
    local ufunc "quad"
    if ("`boxcox'" != "" & "`translog'" == "" & "`quadratic'" == "") {
        local ufunc "boxcox"
    }
    if ("`boxcox'" == "" & "`translog'" != "" & "`quadratic'" == "") {
        local ufunc "tran"
    }
    if (("`boxcox'" != "") + ("`translog'" != "") + ///
        ("`quadratic'" != "") > 1) {
        di in r "utility function can be either 'quad', 'tran' or 'boxcox'"
        exit 498
    }
    // If translog, set up pre-text
    if ("`ufunc'" == "tran") {
        local ln  "ln"
        local pre "`ln'_"
    }
    // If translog or Box-Cox, check for zeros
    if (inlist("`ufunc'", "tran", "boxcox")) {
        qui count if log(`consumption') == . & `touse'
    }
    else {
        qui count if `consumption' < 0       & `touse'
    }
    // Check for negative values
    if (r(N) > 0) {
        di in r "consumption contains values smaller or equal to zero"
        exit 498
    }
    // If Box-Cox, there is no sense in quadratic interactions
    if ("`ufunc'" == "boxcox" ///
      & ("`c2x'" != "" | "`l2x1'" != "" | "`l2x2'" != "")) {
        di in r "options c2x(), l2x1() and l2x2() not allowed with " ///
                "Box-Cox utility function"
        exit 498
    }

    // Validate wage estimation settings
    if ("`wagevars'" != "") local joint joint
    if ("`wagevars'" != "" & "`hwage'" == "") {
        di in r "option hwage() required when estimating jointly"
        exit 498
    }
    // Check sample selection
    if ("`hwage'" != "") {
        foreach hw of local hwage {
            qui count if log(`hw') == . & `touse'
            if (r(N) > 0 & "`wagevars'" == "") {
                di in r "wage variable `hw' censored, use joint estimation"
                exit 498
            }
        }
    }

    // Check dude force settings
    if (`lambda' != 0 & "`force'" == "") {
        di in r "option force() required when marginal utility is contrained"
        exit 498
    }

    // Check choice density settings
    if ("`density'" != "") {
        qui count if `density' <= 0
        if (r(N) > 0) {
            di in r "Choice density contains values smaller or equal to zero."
            exit 498
        }
    }
    if ("`wagevars'" == "" & "`wagedraws'" != "") {
        di in r "What shall we do with the drunken wagedraws? " ///
                "Only useful for joint estimation."
        exit 498
    }
    if ("`wagedraws'" != "" & ("`wagepred'" != "" | "`wagecorr'" != "")) {
        di in r "You can either use wagepred/wagecorr or wagedraws. " ///
                "It's up to you."
        exit 498
    }

    // Get variable count
    local n_leisure   : word count `leisure'
    local n_cxias     : word count `cx'
    local n_lx1ias    : word count `lx1'
    local n_lx2ias    : word count `lx2'
    local n_c2xias    : word count `c2x'
    local n_l2x1ias   : word count `l2x1'
    local n_l2x2ias   : word count `l2x2'
    local n_indeps    : word count `indeps'
    local n_randvars  : word count `randvars'
    local n_wagep     : word count `wagepred'
    local n_hwage     : word count `hwage'
    local n_wagesig   : word count `wagesigma'
    local n_taxrias1  : word count `tria1'
    local n_taxrias2  : word count `tria2'
    local n_wagevars  : word count `wagevars'
    local n_wcorrvars : word count `wagecorr'

    // Validate wage correlation settings
    if ("`wagecorr'" != "") {
        foreach x of local wagecorr {
            if (strpos(" `randvars' ", " `x' ") == 0) {
                local randvars `randvars' `x'
            }
        }
        local n_randvars : word count `randvars'
    }

    // Validate Wage Prediction Options
    if (`n_wagep' == 0) {
        // Joint estimation with correlation?
        if (`n_wcorrvars' != 0) {
            // Run integration over wage errors
            local wagep = 1

            // Fake number of wage prediction indicators
            local n_wagep = `n_leisure'

            // Run wage error prediction for full sample
            local wagep_force = 1
        }
        // Never mind, no wage prediction
        else {
            local wagep       = 0
            local wagep_force = 0
        }
    }
    else {
        // Wage prediction error settings given, stay calm
        local wagep_force = 0

        // Check wage prediction settings
        foreach wagep of local wagepred {
            cap bys `group': assert `wagep' == `wagep'[_N] if `touse'
            if (_rc != 0) {
                di in r "wage prediction error settings ('`wagep'') need " ///
                        "to be constant within households"
                exit 498
            }

            // Check that wage prediction error settings are in line
            //   with wage correlation settings
            if (`n_wcorrvars' != 0) {
                cap assert `wagep' == 1 if `touse'
                if (_rc != 0) local wagep_force = 1
            }
        }

        // Wage prediction enabled
        if (`n_wagep' == `n_leisure' & `n_wagep' == `n_hwage' ///
          & (`n_wagep' == `n_wagesig' | "`wagevars'" != "")) {
            tempvar preds
            qui egen `preds' = rowtotal(`wagepred') if `touse'
            qui count if `preds' >= 1 & (`touse')
            local wagep = (r(N) > 0)
        }
        // Settings incorrect
        else {
            di in r "number of wage prediction variables does not match " ///
                    "the number of leisure terms, hourly wage rates or "  ///
                    "mean squared errors"
            exit 498
        }
    }

    // Show info that wage prediction settings will be ignored
    if (`wagep_force' == 1) {
        di as text "(you selected option 'wagecorr', so wage prediction " ///
                   "errors will be integrated out anyway)"
    }

    // No need to take random draws
    if (`wagep' == 0 & "`randvars'" == "" & "`randsample'" == "") {
        local draws = 1
    }

    // Tax regression or tax benefit calculator needed
    local taxreg_rmse = 0
    if (`wagep' == 1 | "`wagevars'" != "") {
        if ("`taxreg'" != "" & "`taxben'" == "") {
            tempname taxreg_from

            // Load estimates and estimation statistics
            qui est restore `taxreg'
            mat `taxreg_from' = e(b)
            local taxreg_rmse = e(rmse)

            // Fetch "constant" regressors (starting after wage variables
            //   and wage interactions)
            local taxreg_betas   : colnames `taxreg_from'
            local n_taxreg_betas : word count `taxreg_betas'
            local taxreg_vars
            local start = 1 + 2 * (`n_leisure' + `n_taxrias1' + `n_taxrias2')
            forval x = `start'/`n_taxreg_betas' {
                local var : word `x' of `taxreg_betas'
                if ("`var'" != "_cons") local taxreg_vars `taxreg_vars' `var'
            }
        }
        // Run tax benefit calculator
        else if ("`taxreg'" == "" & "`taxben'" != "") {
            //
        }
        // Either taxben or taxreg have to be specified
        else {
            //di in r "either option taxreg() or option taxben() required"
            di in r "option taxreg() required"
            exit 198
        }
    }

    // Build weight settings
    if ("`weight'" != "") {
        local wgtvar `exp'
        local wgt "[`weight'=`wgtvar']"
        qui sum `wgtvar' if `touse', meanonly
        local wgtnobs = r(sum)
    }
    else local wgtnobs = `nobs'

    // Necessary number of random variables
    local rvars = `n_randvars' + `n_wagep' + (`taxreg_rmse' > 0)


    /* LOOK FOR INITIAL VALUES
     */

    if ("`search'" != "off" & "`from'" == "") {
        // Verbose mode on?
        if ("`verbose'" != "") di as text "Looking for initial values..."

        // Set up consumption and leisure
        tempvar c l1 l2
        if (inlist("`ufunc'", "tran", "boxcox")) {
            qui gen `c' = log(cond("`ufunc'" == "boxcox" & `boxcc' > 0,     ///
                                   `consumption' / `boxcc', `consumption')) ///
                          if `touse'
            foreach var of local leisure {
                if (strpos("`leisure'", "`var'") == 1) local lei l1
                else                                   local lei l2
                qui gen ``lei'' = log(cond("`ufunc'" == "boxcox" &  ///
                                                `boxcl' > 0,        ///
                                           `var' / `boxcl', `var')) ///
                                  if `touse'
            }
        }
        else {
            qui gen `c' = `consumption' if `touse'
            foreach var of local leisure {
                if (strpos("`leisure'", "`var'") == 1) local lei l1
                else                                   local lei l2
                qui gen ``lei'' = `var' if `touse'
            }
        }

        // Build up var list for search of initial values
        local initrhs
        if (`n_leisure' == 2) local leisurelist l1 l2
        else                  local leisurelist l1
        foreach ia in c `leisurelist' {
            local f = substr("`ia'", 1, 1)
            if (strlen("`ia'") == 2) local l = substr("`ia'", 2, 1)
            else local l
            foreach var in ``f'x`l'' 0 {
                if ("`var'" != "0") local initrhs `initrhs' c.``ia''#c.`var'
                else                local initrhs `initrhs' ``ia''
            }
            if ("`ufunc'" != "boxcox") {
                foreach var in ``f'2x`l'' 0 {
                    if ("`var'" != "0") {
                        local initrhs `initrhs' c.``ia''#c.``ia''#c.`var'
                    }
                    else local initrhs `initrhs' c.``ia''#c.``ia''
                }
            }
            if ("`ia'" == "c") {
                foreach lei of local leisurelist {
                    local initrhs `initrhs' c.`c'#c.``lei''
                }
            }
        }
        // Leisure cross term
        if (`n_leisure' == 2) local initrhs `initrhs' c.`l1'#c.`l2'
        // Add independent variables to var list
        local initrhs `initrhs' `indeps'

        // Estimate
        tempname init_from
        mat `init_from' = J(1, 2 + 2 * `n_leisure' + `n_cxias' + `n_lx1ias' ///
                                 + `n_lx2ias' + `n_indeps', 0)

        // Get init values from simple conditional logit estimation
        `qui' clogit `varlist' `initrhs' if `touse' `wgt', ///
                    group(`group') iterate(25)
        if (e(converged) == 1) {
            // Save results
            mat `init_from' = e(b)
            local k         = e(k)
            local ll        = e(ll)

            // Update sample
            qui replace `touse' = e(sample)
        }

        // Wage equation
        if ("`wagevars'" != "") {
            tempname init_wage
            foreach w of local hwage {
                tempvar ln`w'
                qui gen `ln`w'' = ln(`w') if `touse' & `varlist'

                // Fetch and store estimates
                `qui' reg `ln`w'' `wagevars' if `touse' & `varlist' `wgt'
                mat `init_wage' = e(b)

                // Save root mean squared error if wage variance is
                //   to be estimated and add correlation terms if
                //   correlation between wages and preferences enabled
                if ("`wagecorr'" != "") {
                    mat `init_wage' = `init_wage', e(rmse), ///
                                      J(1, `n_wcorrvars', 0.0001)
                }
            }
        }

        // Save init options
        local initopt init(`init_from', copy) lf0(`k' `ll')
    }

    // Initial values given
    if ("`initopt'" == "" & "`from'" != "") {
        local initopt init(`from')
    }


    /* PREPARING DATA
     */

    if ("`verbose'" != "") di as text "Preparing data..."

    // Drop missing data
    preserve
    qui keep if `touse'
    sort `group' //`leisure'


    //
    // Setup data
    //

    // To round, or not to round?
    mata: lsl_round  = ("`round'" != "")

    // Utility function
    mata: lsl_ufunc  = st_local("ufunc")

    // Weights?
    mata: lsl_Weight = ("`weight'" != "" ? st_data(., "`wgtvar'") ///
                                         : J(`nobs', 1, 1))

    // Left hand side / Choice indicator
    mata: lsl_Y      = st_data(., "`varlist'")

    // Hourly wage rates
    mata: lsl_Hwage  = (`n_hwage' > 0 ? st_data(., tokens("`hwage'")) ///
                                      : J(`nobs', `n_leisure', 0))

    // Normalizing constant for Box-Cox consumption
    mata: lsl_boxcc  = (lsl_ufunc == "boxcox" & `boxcc' > 0 ? `boxcc' : 1)

    // Normalizing constant for Box-Cox leisure
    mata: lsl_boxcl  = (lsl_ufunc == "boxcox" & `boxcl' > 0 ? `boxcl' : 1)


    //
    // Right hand side
    //

    // Consumption, squared and interactions
    mata: lsl_C   = st_data(., "`consumption'") :/ lsl_boxcc
    mata: lsl_CX  = ((`n_cxias'  > 0 ? st_data(., tokens("`cx'"))  ///
                                     : J(`nobs', 0, 0)),           ///
                     J(`nobs', 1, 1))
    mata: lsl_C2X = ((`n_c2xias' > 0 ? st_data(., tokens("`c2x'")) ///
                                     : J(`nobs', 0, 0)),           ///
                     J(`nobs', (lsl_ufunc != "boxcox" ? 1 : 0), 1))

    // Leisure, squared and interactions
    forval i = 1/2 {
        local var : word `i' of `leisure'
        if ("`var'" != "") {
            mata: lsl_L`i'   = st_data(., "`var'") :/ lsl_boxcl
            mata: lsl_LX`i'  = ((`n_lx`i'ias'  >  0                 ///
                                   ? st_data(., tokens("`lx`i''"))  ///
                                   : J(`nobs', 0, 0)),              ///
                                J(`nobs', 1, 1))
            mata: lsl_L2X`i' = ((`n_l2x`i'ias' >  0                 ///
                                   ? st_data(., tokens("`l2x`i''")) ///
                                   : J(`nobs', 0, 0)),              ///
                                J(`nobs', (lsl_ufunc != "boxcox" ? 1 : 0), 1))
        }
        else {
            mata: lsl_L`i'   = J(`nobs', 0, 0)
            mata: lsl_LX`i'  = J(`nobs', 0, 0)
            mata: lsl_L2X`i' = J(`nobs', 0, 0)
        }
    }

    // Dummy variables
    mata: lsl_Xind = (`n_indeps' > 0 ? st_data(., tokens("`indeps'")) ///
                                     : J(`nobs', 0, 0))

    //
    // Build right hand side
    //

    // Translog utility
    if ("`ufunc'" == "tran") {
        mata: lsl_X = (log(lsl_C)  :* (lsl_CX,  log(lsl_C)  :* lsl_C2X,   ///
                                       log(lsl_L1), log(lsl_L2)),         ///
                       log(lsl_L1) :* (lsl_LX1, log(lsl_L1) :* lsl_L2X1), ///
                       log(lsl_L2) :* (lsl_LX2, log(lsl_L2) :* lsl_L2X2), ///
                       log(lsl_L1) :* log(lsl_L2), lsl_Xind)
    }
    // Quadratic utility
    else if ("`ufunc'" == "quad") {
        mata: lsl_X = (lsl_C  :* (lsl_CX,  lsl_C  :* lsl_C2X,   ///
                                  lsl_L1, lsl_L2),              ///
                       lsl_L1 :* (lsl_LX1, lsl_L1 :* lsl_L2X1), ///
                       lsl_L2 :* (lsl_LX2, lsl_L2 :* lsl_L2X2), ///
                       lsl_L1 :* lsl_L2, lsl_Xind)
    }
    else mata: lsl_X = J(`nobs', 0, 0)


    //
    // Joint wage estimation
    //

    // ML joint estimation?
    mata: lsl_joint       = ("`joint'" != "")

    // Wage observed?
    mata: lsl_Wobs        = (lsl_Y :* (log(lsl_Hwage) :< .))

    // Wage variables
    mata: lsl_WageVars    = (lsl_joint == 1                         ///
                               ? (st_data(., tokens("`wagevars'")), ///
                                  J(`nobs', 1, 1))                  ///
                               : J(`nobs', 0, 0))

    // Days of taxyear
    mata: lsl_Days        = ("`days'" != "" ? st_data(., "`days'") ///
                                            : J(`nobs', 1, 365))

    // Hypothetical hours (caution: Box-Cox normalization!)
    mata: lsl_Hours       = `totaltime' :- (lsl_L1, lsl_L2) :* lsl_boxcl

    // Correlate wage rates and preferences?
    mata: lsl_wagecorr    = `n_wcorrvars'
    mata: st_local("wagecorr", invtokens(strofreal(sort(strtoreal( ///
                                         tokens("`wagecorr'"))', 1)')))

    // Wage correlation coefficients
    mata: lsl_Wcorrvars   = ("`wagecorr'" != ""                   ///
                               ? strtoreal(tokens("`wagecorr'"))' ///
                               : J(0, 1, 0))

    // Use actual residiuals for folks with observed wages
    mata: lsl_residanchor = ("`anchor'" != "noanchor") & (`wagep' == 1) ///
                                                       & lsl_joint

    // Wagedraws? For estimation on random choice set
    mata: lsl_Wagedraws   = ("`wagedraws'" != ""           ///
                               ? st_data(., "`wagedraws'") ///
                               : J(`nobs', 0, 0))


    //
    // Wage Prediction Stuff
    //

    // Run Wage Prediction
    mata: lsl_wagep = `wagep'

    // Dummies enabling or disabling the wage prediction for individuals
    mata: lsl_Wpred = (lsl_wagep ? (`wagep_force'                 ///
                                      ? J(`nobs', `n_wagep', 1)   ///
                                      : st_data(., "`wagepred'")) ///
                                 : J(`nobs', 1 + cols(lsl_L2), 0))

    // Estimated variance of wage equation
    mata: lsl_Sigma = ("`wagesigma'" != ""                  ///
                         ? strtoreal(tokens("`wagesigma'")) ///
                         : J(1, `n_wagesig', 0))

    //
    // Tax regression
    //

    if ("`taxreg'" != "") {
        // Tax regression estimates
        mata: lsl_TaxregB     = st_matrix("`taxreg_from'")

        // Root Mean Squared error of tax regression
        mata: lsl_taxreg_rmse = `taxreg_rmse'

        // Interaction variables on Mwage1 and Mwage1^2
        mata: lsl_TaxregIas1  = ("`tria1'" != ""                   ///
                                   ? st_data(., tokens("`tria1'")) ///
                                   : J(`nobs', 0, 0))

        // Interaction variables on Mwage2 and Mwage2^2
        mata: lsl_TaxregIas2  = ("`tria2'" != ""                   ///
                                   ? st_data(., tokens("`tria2'")) ///
                                   : J(`nobs', 0, 0))

        // Variables that are independent of m_wage
        mata: lsl_TaxregVars  = st_data(., tokens("`taxreg_vars'"))
    }


    //
    // Group level stuff
    //

    // Separate households
    mata: lsl_J = panelsetup(st_data(., "`group'"), 1)

    // Add column with number of choices
    mata: lsl_J = (lsl_J, 1 :+ lsl_J[.,2] :- lsl_J[.,1])

    // Number of groups
    mata: lsl_groups = rows(lsl_J)
    mata: st_numscalar("n_groups", lsl_groups)

    // Density of offered choices
    mata: lsl_Dens = ("`density'" != "" ? st_data(., "`density'") ///
                                        : J(`nobs', 1, 1))


    //
    // Dude force (lagrange constraints)
    //

    // Lagrange multiplier
    mata: lsl_lambda = `lambda'

    // Which observations to include in constraint?
    mata: lsl_force  = ("`force'" != "" ? st_data(., "`force'") ///
                                        : J(`nobs', 1, 0))


    //
    // Random draws
    //

    // Number of draws
    mata: lsl_draws = strtoreal(st_local("draws"))

    // Number of draws to burn
    mata: lsl_burn  = strtoreal(st_local("burn"))
    mata: st_local("randvars", invtokens(strofreal(sort(strtoreal( ///
                                         tokens("`randvars'"))', 1)')))

    // Random coefficients
    mata: lsl_Rvars = ("`randvars'" != ""                   ///
                         ? strtoreal(tokens("`randvars'"))' ///
                         : J(0, 1, 0))

    // Random coefficients correlated?
    mata: lsl_corr  = ("`corr'" != "")

    // Halton sequences
    mata: lsl_R     = (`rvars' > 0                                  ///
                         ? invnormal(halton(lsl_groups * lsl_draws, ///
                                            `rvars', 1 + lsl_burn)) ///
                         : J(`nobs', 0, 0))

    /// Estimate on random choice set?
    mata: lsl_randsample = ("`randsample'" != "")
    if ("`randsample'" != "") {
        mata: lsl_WDRW = invnormal(runiform(`nobs', lsl_draws))
        mata: lsl_WPDF = normalden(lsl_WDRW)
    }
    else mata: lsl_WPDF = J(`nobs', lsl_draws, 1)


    //
    // Calculate log-likelihood of null model
    //

    mata: st_local("SigmaW", strofreal(                                     ///
                   (lsl_joint ? sqrt(cross(lsl_Wobs :* log(lsl_Hwage),      ///
                                           log(lsl_Hwage)) /                ///
                                     (colsum(lsl_Wobs) - `n_wagevars' - 1)) ///
                              : 0)))
    mata: st_local("ll_0", strofreal(                                       ///
                   - colsum(cross(lsl_Weight[lsl_J[.,1]], log(lsl_J[.,3]))) ///
                   + (lsl_joint ? cross(lsl_Weight :* lsl_Wobs,             ///
                                        log(normalden(log(lsl_Hwage) :/     ///
                                        `SigmaW')) :- log(`SigmaW'))        ///
                                : 0)))
    if ("`density'" != "") {
        di as text "(Note: Pseudo-R2 buggy when using choice densities)"
    }

    // Restore data
    restore


    /* RUN ESTIMATION
     */

    if ("`verbose'" != "") di as text "Run estimation..."

    //
    // Set up equations
    //

    // Consumption
    local eq_consum (Cx: `varlist' = `cx')

    // Consumption^2
    if ("`ufunc'" != "boxcox") local eq_consum `eq_consum' (CxC: `c2x')

    // Leisure terms
    local eq_leisure
    foreach var of local leisure {
        local i = 1 + (strpos("`leisure'", "`var'") > 1)

        // Leisure
        local eq_leisure `eq_leisure' (L`i'x: `lx`i'')

        // Leisure^2
        if ("`ufunc'" != "boxcox") local eq_leisure `eq_leisure' ///
                                                    (L`i'xL`i': `l2x`i'')

        // Consumption X leisure interaction
        local eq_consum  `eq_consum' /CxL`i'
    }

    // Leisure term interaction
    if (`n_leisure' == 2) local eq_leisure `eq_leisure' /L1xL2

    // Independent variables / dummies
    if (`n_indeps'  >  0) local eq_indeps  (IND: `indeps', noconst)

    // Joint wage estimation?
    if ("`wagevars'" != "") {
        // Wage equation
        local eq_wages (lnW1: `wagevars')

        // Correlation with preferences?
        if ("`wagecorr'" != "") {
            // Cholesky for Var(W1) and Wagecorr
            local eq_wages `eq_wages' /l_W1_W1
            foreach x of local wagecorr {
                local eq_wages `eq_wages' /l_`x'_W1
            }
        }

        // Add wage estimates to initial values
        if ("`initopt'" != "" & "`init_from'" != "") {
            mat `init_from' = (`init_from', `init_wage')
        }
    }

    // Box-Cox utility function?
    if ("`ufunc'" == "boxcox") {
        local eq_boxcox /l_C /l_L1
        if (`n_leisure' == 2) local eq_boxcox `eq_boxcox' /l_L2

        // Add initial values (assume Stone-Geary model)
        if ("`initopt'" != "" & "`init_from'" != "") {
            mat `init_from' = (`init_from', 0.0, J(1, `n_leisure', 0.0))
        }
    }

    // Random coefficients?
    if (`n_randvars' > 0) {
        local eq_rands
        if ("`corr'" == "") {
            foreach coef of local randvars {
                local eq_rands `eq_rands' /sd_`coef'
            }
            if ("`initopt'" != "" & "`init_from'" != "") {
                mat `init_from' = (`init_from', J(1, `n_randvars', 0.0001))
            }
        }
        else {
            forval i = 1/`n_randvars' {
                local a : word `i' of `randvars'
                forval k = `i'/`n_randvars' {
                    local b : word `k' of `randvars'
                    local eq_rands `eq_rands' /l_`a'_`b'
                }
            }
            if ("`initopt'" != "" & "`init_from'" != "") {
                mat `init_from' = (`init_from',        ///
                                   J(1, `n_randvars' * ///
                                        (`n_randvars' + 1) / 2, 0.0002))
            }
        }
    }


    //
    // Estimate
    //

    scalar lsl_dudes = .

    local callcmd ml model `method'`debug' lslogit_d2() `eq_consum' ///
                            `eq_leisure' `eq_indeps' `eq_wages'     ///
                            `eq_boxcox' `eq_rands'                  ///
                        if `touse' `wgt', group(`group') `initopt'  ///
                            obs(`wgtnobs') search(off) nopreserve   ///
                            iterate(`iterate') max `difficult'      ///
                            `trace' `gradient' `hessian'            ///
                            technique(`technique')
    if ("`debug'" != "" | "`verbose'" != "") di as text "`callcmd'"
    `callcmd'


    /* SAVE RESULTS
     */

    // Save model setup
    ereturn local title         "Mixed Logit Labor Supply Model"
    ereturn local cmd           "lslogit"
    ereturn local predict       "lslogit_p"
    ereturn local ufunc          `ufunc'
    ereturn local wagep          `wagep'
    ereturn local joint          `joint'
    ereturn local residanchor = ("`anchor'" != "noanchor") & (`wagep' == 1) ///
                                                           & ("`joint'" != "")
    ereturn local draws          `draws'
    ereturn local burn           `burn'
    ereturn local totaltime      `totaltime'
    ereturn local round       = ("`round'" != "")
    ereturn local taxreg      = ("`taxreg'" != "")
    ereturn local group          `group'
    ereturn local groups         `n_groups'
    ereturn local depvar         `varlist'
    ereturn local consum         `consumption'
    ereturn local leisure        `leisure'

    // Save varlists
    if ("`cx'"       != "") ereturn local cx       `cx'
    if ("`c2x'"      != "") ereturn local c2x      `c2x'
    if ("`lx1'"      != "") ereturn local lx1      `lx1'
    if ("`l2x1'"     != "") ereturn local l2x1     `l2x1'
    if ("`lx2'"      != "") ereturn local lx2      `lx2'
    if ("`l2x2'"     != "") ereturn local l2x2     `l2x2'
    if ("`indeps'"   != "") ereturn local indeps   `indeps'
    if ("`wagepred'" != "") ereturn local wagepred `wagepred'
    if ("`hwage'"    != "") ereturn local hwage    `hwage'
    if ("`days'"     != "") ereturn local days     `days'
    if ("`wgtvar'"   != "") ereturn local weight   `wgtvar'

    // Box-Cox
    if ("`ufunc'" == "boxcox") {
        if ("`boxcc'" != "") ereturn local boxcc   `boxcc'
        if ("`boxcl'" != "") ereturn local boxcl   `boxcl'
    }

    // Random coefficients
    if (`n_randvars' > 0) {
        ereturn local randvars `randvars'
        ereturn local corr = cond("`corr'" != "", 1, 0)
    }

    // Random coefficients correlated with wages
    if (`n_wcorrvars' > 0) {
        ereturn local wcorrvars `wagecorr'
        ereturn local wagecorr = `n_wcorrvars'
    }

    // Tax regression
    if ("`taxreg'" != "") {
        ereturn matrix taxreg_b = `taxreg_from', copy
        if ("`taxreg_vars'" != "") ereturn local taxreg_vars `taxreg_vars'
        if ("`tria1'" != "") ereturn local taxreg_ia1 `tria1'
        if ("`tria2'" != "") ereturn local taxreg_ia2 `tria2'
        ereturn scalar taxreg_rmse = `taxreg_rmse'
    }
    if ("`wagesigma'" != "") ereturn local wagesigma `wagesigma'
    if ("`wagevars'"  != "") ereturn local wagevars  `wagevars'
    if ("`wagedraws'" != "") ereturn local wagedraws `wagedraws'
    if ("`density'"   != "") ereturn local density   `density'

    // Display some coefficients as auxiliary
    ereturn scalar k_aux = ("`ufunc'" == "boxcox") * (1 + `n_leisure') +    ///
                           `n_randvars' * cond("`corr'" != "",              ///
                                               (`n_randvars' + 1) / 2, 1) + ///
                           ("`wagecorr'" != "") + `n_wcorrvars'

    // Pseudo R2 (may be misleading as it refers to the null-model,
    //            LR and p value refer to init values...)
    if ("`ll_0'" != "") ereturn scalar r2_p = 1 - e(ll)/`ll_0'

    // Additional output
    if (lsl_dudes != .) ereturn scalar dudes = lsl_dudes
    /*
    foreach aux in sigma_w1 sigma_w2 {
        di as error "`aux' = " `aux'
        if (`aux' != "") ereturn scalar `aux' = ``aux''
    }
    */


    /* CLEAN UP
     */

    /*
    foreach m in round ufunc Weight Y Hwage boxcc boxcl C CX C2X ///
                 L1 LX1 L2X1 L2 LX2 L2X2 Xind X joint            ///
                 WageVars Days Hours wagep Wpred Sigma           ///
                 TaxregB TaxregIas1 TaxregIas2 TaxregVars        ///
                 groups J draws burn Rvars corr R {
        cap mata mata drop lsl_`m'
    }
    */


    //
    // Show results
    //

    lslogit_Replay, level(`level') `quiet'
end


// Drop mata functions if they exist
foreach fct in lslogit_d2 lsl_boxcox lsl_boxcox_g lsl_boxcox_h lslogit_p {
    cap mata mata drop `fct'()
}

mata
mata set matastrict on
/**
 * Evaluator
 */
void lslogit_d2(transmorphic scalar ML, real scalar todo, real rowvector B,
                real scalar lnf, real rowvector G, real matrix H) {

    // Functional form
    external string scalar  lsl_ufunc

    // Number of groups
    external real scalar    lsl_groups

    // Left hand side variable
    external real colvector lsl_Y

    // Panel setup data
    external real matrix    lsl_J

    // Right hand side variables
    external real matrix    lsl_X

    // Group weights
    external real colvector lsl_Weight

    // Density of choice
    external real colvector lsl_Dens

    // Number of random draws
    external real scalar    lsl_draws

    //   Halton sequences
    external real matrix    lsl_R

    // Random coefficients
    external real colvector lsl_Rvars

    //   Enable correlation?
    external real scalar    lsl_corr

    // Joint ML estimation?
    external real scalar    lsl_joint

    //   Wage rate observed?
    external real matrix    lsl_Wobs

    //   Right hand side variables
    external real matrix    lsl_WageVars

    //   Number of covariances between wages and leisure preferences?
    external real scalar    lsl_wagecorr

    //   Number of covariances between wages and leisure preferences?
    external real colvector lsl_Wcorrvars

    external real scalar    lsl_residanchor
    external real matrix    lsl_Wagedraws
    external real matrix    lsl_randsample
    external real matrix    lsl_WPDF
    external real matrix    lsl_WDRW

    // Wage Prediction Error?
    external real scalar    lsl_wagep

    //   Prediction dummies
    external real matrix    lsl_Wpred

    //   Number of days per tax year
    external real colvector lsl_Days

    //   Hourly wage rates
    external real matrix    lsl_Hwage

    //   Variance of the wage regression
    external real rowvector lsl_Sigma

    //   Hours of work
    external real matrix    lsl_Hours

    // Tax Regression
    external real rowvector lsl_TaxregB

    //   Root Mean Squared Error
    external real scalar    lsl_taxreg_rmse

    //   Wage independent variables of tax regression
    external real matrix    lsl_TaxregVars

    //   Wage interaction variables of tax regression
    external real matrix    lsl_TaxregIas1

    //   Wage interaction variables of tax regression
    external real matrix    lsl_TaxregIas2

    // To round, or not to round.
    external real scalar    lsl_round

    // Normalizing constant for Box-Cox consumption
    external real scalar    lsl_boxcc

    // Lagrange multiplier for dude correction
    external real scalar    lsl_lambda

    // Which observations to include in constraint?
    external real colvector lsl_force

    // Right hand side and interactions
    external real colvector lsl_C
    external real matrix    lsl_CX
    external real matrix    lsl_C2X
    external real colvector lsl_L1
    external real matrix    lsl_LX1
    external real matrix    lsl_L2X1
    external real matrix    lsl_LX2
    external real matrix    lsl_L2X2
    external real matrix    lsl_Xind

    // Buggy, L2 is colvector in fact, but colvector needs one column...
    external real matrix    lsl_L2


    //
    // Declarations
    //

    real scalar    n, i, c, e, nobs, nlei, ncons, wagep
    real scalar    b, bfix, bwage, bheck, blam, brnd
    real rowvector Bfix, Bwage, Bsig, Brnd, Zeta, Beta

    real scalar    nRV, rvars, iRV, r, iC, iL1, iL2,
                   iwage, iheck, ilam, irnd, ilC, ilL1, ilL2
    real matrix    CholB, CholBW, CholW

    real matrix    Hwage, LnWresPur, Wn
    real matrix    Mwage, TaxregX1, TaxregX2, TaxregX, DCdM, D2CdM2

    real matrix    DUdx, DUdB, DUdlam, DUdBr, DWdBw, Dude,
                   DWdBsig, YmPn_D2UdB2, YmPn_D2UdBr2, YmPn_D2UdBdBr,
                   YmPn_D2Udx2, DUdlC, DUdlL1, DUdlL2
    real colvector DUdC, D2UdC2, DMdH, D2MdH2

    real scalar    ncx, nc2x, nlx1, nl2x1, nlx2, nl2x2, nxind
    real colvector Yn, C, L1, L2
    real matrix    Xnr, CX, C2X, LX1, L2X1, LX2, L2X2, Xind
    real colvector Unr, Enr, Pnr, YmPn

    real scalar    lsum, pni
    real rowvector Gsum
    real matrix    H1sum, H2sum

    real scalar    lC, lL1, lL2
    real colvector BcC, BcL1
    real matrix    BcL2, BcCx, BcL1x, BcL2x
    real scalar    SigmaW

    //lsl_Y = moptimize_util_depvar(ML, 1)     // Left hand side variable


    /* Setup */

    //
    // Definitions
    //

    // Indicates first observation of active group
    i = 1

    // Number of random variables
    rvars = rows(lsl_Rvars)

    // Indicates first random variable to use (column of lsl_R) for wage pred
    nRV   = rvars + 1

    // Number of observations
    nobs  = rows(lsl_Y)

    // Number of right hand side terms
    nlei  = 1 + cols(lsl_L2)    // Leisure terms
    ncx   = cols(lsl_CX)        // Consumption interactions
    nc2x  = cols(lsl_C2X)       //   squared
    nlx1  = cols(lsl_LX1)       // Leisure 1 interactions
    nl2x1 = cols(lsl_L2X1)      //   squared
    nlx2  = cols(lsl_LX2)       // Leisure 2 interactions
    nl2x2 = cols(lsl_L2X2)      //   squared
    nxind = cols(lsl_Xind)      // Independent terms
    ncons = ncx + nc2x + nlei   // Number of terms including consumption


    //
    // Number of coefficients
    //

    // Total number
    b     = cols(B)

    // Number of wage regression coefficients
    bwage = cols(lsl_WageVars)

    // Number of additional wage related coefficients
    //   (correlation with preferences)
    bheck = (lsl_wagecorr > 0) + rows(lsl_Wcorrvars)

    // Number of variance and covariance terms for random coefficients
    brnd  = (lsl_corr  == 1 ? rvars * (rvars + 1) / 2 : rvars)

    // Number of Box-Cox transformation coefficients
    blam  = (lsl_ufunc == "boxcox" ? 1 + nlei : 0)

    // Number of fixed preference coefficients
    bfix  = ncons + nlx1 + nl2x1 + nlx2 + nl2x2 + (nlei == 2) + nxind

    // Maximum Likelihood Parameter
    lnf = 0             // Log-likelihood
    G   = J(1, b, 0)    // Gradient
    H   = J(b, b, 0)    // Hessian matrix

    // Initialize individual dude probability
    Dude = J(nobs, lsl_draws, 0)


    //
    // Build coefficient vector
    //

    // Get coefficient indices
    iwage = 1     + bfix
    iheck = iwage + bwage
    ilam  = iheck + bheck
    irnd  = ilam  + blam

    iC  = ncx
    iL1 = ncons + nlx1
    iL2 = iL1 + nl2x1 + nlx2

    // Split up vector

    // Get fixed coefficients
    Bfix  = B[|1\bfix|]

    // Wage coefficients
    Bwage = (bwage > 0 ? B[|iwage\iwage + bwage - 1|] : J(0, 0, 0))

    //   Cholesky factor of wage variance
    Bsig  = (bheck > 0 ? B[iheck]                     : J(0, 0, 0))

    // Get auxiliary random coefficients
    Brnd  = (rvars > 0 ? B[|irnd\irnd + brnd - 1|]    : J(0, 0, 0))

    // Build Cholesky matrix

    // Cholesky factors of coefficients vector
    CholB  = (lsl_corr ? lowertriangle(invvech(Brnd')) : diag(Brnd))

    // Cholesky factors between coefficients and wages
    CholBW = (lsl_wagecorr > 0 ? B[|iheck + 1\iheck + bheck - 1|]
                               : J(nlei, rvars, 0))

    // Cholesky factors of wage variances
    CholW  = diag(Bsig)

    // DEBUG (works only for singles!)
    if (lsl_wagecorr > 0 & lsl_Wcorrvars != lsl_Rvars) {
        if (nlei == 1) {
            // Cholesky factors between coefficients and wages
            CholBW = J(1, 0, 0)
            real scalar v, w, s
            s = 1
            for (v = 1; v <= rows(lsl_Wcorrvars); v++) {
                for (w = s; w <= rows(lsl_Rvars); w++) {
                    if (lsl_Wcorrvars[v] == lsl_Rvars[w]) {
                        CholBW = CholBW, B[iheck + v]
                        s = w + 1
                        break
                    }
                    else CholBW = CholBW, 0
                }
            }
            CholBW = CholBW, J(1, rows(lsl_Rvars) - cols(CholBW), 0)
        }
        else "DEBUG!"
    }
    // DEBUG

    // Calculate wage variance (by Cholesky or just the coefficient)
    SigmaW = (lsl_wagecorr ? sqrt(rowsum(CholW:^2))' : Bsig)

    // Initialize matrix with random coefficients, every row is a draw
    if (brnd > 0) Zeta = J(rows(lsl_R), bfix, 0)
    // DEBUG: This makes no sense anymore since we fill Zeta for each n/r...
    //        Correct in the future and drop this line here.


    //
    // Box-Cox utility function
    //

    if (lsl_ufunc == "boxcox") {
        // Get Box-Cox-Lambdas
        ilC  = ilam
        ilL1 = ilam + 1
        ilL2 = (nlei == 2 ? ilam + 2 : 0)
        lC   = B[ilC]
        lL1  = B[ilL1]
        lL2  = (nlei == 2 ? B[ilL2] : 0)

        // Transform consumption and leisure
        BcC  = lsl_boxcox(lsl_C, lC)
        BcL1 = lsl_boxcox(lsl_L1, lL1)
        BcL2 = lsl_boxcox(lsl_L2, lL2)

        // Replace lsl_X
        lsl_X = ((lsl_CX, BcL1, BcL2) :* BcC, lsl_LX1 :* BcL1,
                 lsl_LX2 :* BcL2, BcL1 :* BcL2, lsl_Xind)
    }


    //
    // Joint wage estimation
    //

    if (lsl_joint) {
        // Predict log-wages
        Hwage = cross(lsl_WageVars', Bwage')

        // Wage observed?
        //Wobs = lsl_Y :* (log(lsl_Hwage) :< .)

        // Redisuals
        LnWresPur = lsl_Wobs :* (log(lsl_Hwage) :- Hwage)

        // No correlation stuff? Calculate wage variance
        if (!lsl_wagecorr) {
            // Root MSE
            Bsig   = sqrt(cross(LnWresPur, LnWresPur) /
                          (colsum(lsl_Wobs) - bwage))

            // Overwrite wage variance
            SigmaW = Bsig

            // Overwrite Cholesky factors
            CholW  = diag(Bsig)

            // Save RMSE
            for (c = 1; c <= cols(Bsig); c++) {
                st_numscalar("r(sigma_w" + strofreal(c) + ")", Bsig[1,c])
            }
        }

        // Predict wages
        Hwage = exp(Hwage :+ (SigmaW^2 / 2)) :* (lsl_Hours[|1,1\.,1|] :!= 0)

        // Use random wage draws
        if (cols(lsl_Wagedraws) > 0) {
            real colvector Wagedraws
            Wagedraws = invnormal(runiform(nobs, 1))
            Hwage = Hwage :* exp(SigmaW :* lsl_Wagedraws :* (lsl_Wobs :== 0)
                                 :+ LnWresPur)
            lsl_Dens = normalden(lsl_Wagedraws) :* (lsl_Wobs :== 0) :+
                       lsl_Wobs :* normalden(LnWresPur, SigmaW)
        }
    }
    // No joint estimation, get initial data
    else {
        Hwage  = lsl_Hwage
        Bsig   = lsl_Sigma
        SigmaW = Bsig           // Overwrite wage variance
        CholW  = diag(Bsig)     // Overwrite Cholesky factors
    }

    // Round wage rates?
    if (lsl_round & cols(Hwage) > 0) Hwage = round(Hwage, 0.01)


    /* Loop over households */

    for (n = 1; n <= lsl_groups; n++) {
        // Last observation of group n
        i = lsl_J[n,1]
        e = lsl_J[n,2]
        c = e - i + 1

        // Fetch relevant variables
        Yn   = lsl_Y[|i\e|]
        Xnr  = lsl_X[|i,1\e,.|]
        C    = lsl_C[|i\e|]
        CX   = (cols(lsl_CX)   > 0 ?   lsl_CX[|i,1\e,.|] : J(c, 0, 0))
        C2X  = (cols(lsl_C2X)  > 0 ?  lsl_C2X[|i,1\e,.|] : J(c, 0, 0))
        L1   = lsl_L1[|i\e|]
        LX1  = (cols(lsl_LX1)  > 0 ?  lsl_LX1[|i,1\e,.|] : J(c, 0, 0))
        L2X1 = (cols(lsl_L2X1) > 0 ? lsl_L2X1[|i,1\e,.|] : J(c, 0, 0))
        L2   = (cols(lsl_L2)   > 0 ?   lsl_L2[|i\e|]     : J(c, 0, 0))
        LX2  = (cols(lsl_LX2)  > 0 ?  lsl_LX2[|i,1\e,.|] : J(c, 0, 0))
        L2X2 = (cols(lsl_L2X2) > 0 ? lsl_L2X2[|i,1\e,.|] : J(c, 0, 0))
        Xind = (cols(lsl_Xind) > 0 ? lsl_Xind[|i,1\e,.|] : J(c, 0, 0))
        Wn   = Hwage[|i,1\e,.|]

        // Transform consumption and leisure
        if (lsl_ufunc == "boxcox") {
            BcC  = lsl_boxcox(C, lC)
            BcL1 = lsl_boxcox(L1, lL1)
            BcL2 = lsl_boxcox(L2, lL2)
        }

        // Sum over draws
        lsum  = 0
        Gsum  = J(1, b, 0)
        H1sum = J(1, b, 0)
        H2sum = J(b, b, 0)

        // Check if wage prediction needed
        wagep = (lsl_wagep & sum(lsl_Wpred[|i,1\e,.|]) > 0)

        // Wage residual anchor? Use actual wage equation residuals
        //   instead of random draws
        if (lsl_residanchor & lsl_joint & wagep
          & colsum(lsl_Wobs[|i,1\e,1|]) == 1) {
            lsl_R[|lsl_draws * (n - 1) + 1,1\lsl_draws * (n - 1) + lsl_draws,1|] =
                    J(lsl_draws, 1, colsum(LnWresPur[|i,1\e,1|]) / CholW)
        }
        if (lsl_randsample & lsl_joint) {
            lsl_WDRW[|i,1\e,.|] = lsl_WDRW[|i,1\e,.|] :*
                                   (lsl_Wobs[|i,1\e,.|] :== 0) :+
                                  LnWresPur[|i,1\e,1|] / SigmaW
            lsl_WPDF[|i,1\e,.|] = normalden(lsl_WDRW[|i,1\e,.|])
        }

        // Run by random draw
        for (r = 1; r <= lsl_draws; r++) {
            // Indicates the active Halton sequence
            iRV  = lsl_draws * (n - 1) + r

            // Build (random?) coefficients matrix (DEBUG HERE!)
            if (brnd > 0) {
                Zeta[iRV,lsl_Rvars] =
                    cross(lsl_R[|iRV,1\iRV,lsl_wagep * nlei + rvars|]',
                          (CholBW', CholB)')
            }
            Beta = Bfix :+ (brnd > 0 ? Zeta[iRV,.] : 0)


            /* Integrate out wage prediction error */

            if (wagep | lsl_joint) {
                //
                // Calculate monthly earnings
                //

                // Adjust wages with random draws if prediction enabled
                if (lsl_wagep) {
                    Wn = Hwage[|i,1\e,.|] :* exp(cross(CholW', lsl_R[|iRV,1\iRV,nlei|]')' :*
                                                  lsl_Wpred[|i,1\e,.|])
                }
                if (lsl_randsample) Wn = Hwage[|i,1\e,.|] :* exp(cross((CholBW, CholW)',
                                                                 lsl_WDRW[|i,r\e,r|]')')

                // Calculate monthly earnings
                Mwage = (lsl_Days[|i\e|] :/ 12 :/ 7) :*
                         lsl_Hours[|i,1\e,.|] :* Wn

                // Round monthly earnings if enabled
                if (lsl_round) Mwage = round(Mwage, 0.01)


                //
                // Set up tax regression covariates and predict dpi
                //

                // Fill matrix of independent variables for dpi prediction
                TaxregX1 = (Mwage[.,1], Mwage[.,1]:^2,
                            lsl_TaxregIas1[|i,1\e,.|] :* Mwage[.,1],
                            lsl_TaxregIas1[|i,1\e,.|] :* Mwage[.,1]:^2)
                if (nlei == 2) {
                    TaxregX2 = (Mwage[.,2], Mwage[.,2]:^2,
                                lsl_TaxregIas2[|i,1\e,.|] :* Mwage[.,2],
                                lsl_TaxregIas2[|i,1\e,.|] :* Mwage[.,2]:^2)
                } else TaxregX2 = J(c, 0, 0)
                TaxregX = (TaxregX1, TaxregX2,
                           lsl_TaxregVars[|i,1\e,.|], J(c, 1, 1))

                // Predict disposable income
                //   - can't be negative!
                //   - add random draw
                //   - normalize for Box-Cox specification
                C = rowmax((cross(TaxregX', lsl_TaxregB') :+
                            lsl_taxreg_rmse :* lsl_R[iRV,cols(lsl_R)],
                            J(c, 1, 1))) :/ lsl_boxcc

                // Build matrix with independent variables
                if (lsl_ufunc == "tran") {
                    Xnr = (log(C)  :* (CX,  log(C)  :* C2X, log(L1), log(L2)),
                           log(L1) :* (LX1, log(L1) :* L2X1),
                           log(L2) :* (LX2, log(L2) :* L2X2),
                           log(L1) :* log(L2), Xind)
                } else if (lsl_ufunc == "quad") {
                    Xnr = (C  :* (CX,  C  :* C2X,   L1,   L2),
                           L1 :* (LX1, L1 :* L2X1),
                           L2 :* (LX2, L2 :* L2X2), L1 :* L2, Xind)
                } else if (lsl_ufunc == "boxcox") {
                    BcC = lsl_boxcox(C, lC)
                    Xnr = ((CX, BcL1, BcL2) :* BcC, LX1 :* BcL1,
                           LX2 :* BcL2, BcL1 :* BcL2, Xind)
                }
            }


            /* Calculate utility levels */

            // Utility (choices in rows, draws in columns)
            Unr = cross(Xnr', Beta') :- log(lsl_WPDF[|i,r\e,r|])

            // Standardize to avoid missings
            Enr = exp(Unr :+ colmin(-mean(Unr) \ 700 :- colmax(Unr)))

            // Probabilities
            Pnr = Enr :/ colsum(Enr)

            // Simplify
            pni  = cross(Yn, Pnr)   // Probability that choice is chosen
            YmPn = Yn :- Pnr        // Choice minus probabilities

            // Calculate first derivative to consumption (dudes)
            if (todo == 0 | lsl_lambda > 0) {
                if (lsl_ufunc == "boxcox") {
                    Dude[|i,r\e,r|] = cross((CX, lsl_boxcox(L1, lL1),
                                             lsl_boxcox(L2, lL2))',
                                            Beta[|1\ncons|]') :*
                                       (reldif(lC, 0) >= 1e-25 ? C:^(lC - 1)
                                                                  : (1 :/ C))
                } else if (lsl_ufunc == "quad") {
                    Dude[|i,r\e,r|] = cross((CX, 2 :* C2X :* C, L1, L2)',
                                            Beta[|1\ncons|]')
                } else if (lsl_ufunc == "tran") {
                    Dude[|i,r\e,r|] = cross((CX, 2 :* C2X :* log(C), log(L1),
                                            log(L2))', Beta[|1\ncons|]') :/ C
                }
            }


            /* Add to sum over draws */

            // Add to likelihood
            lsum = lsum + pni

            // Calculate gradient vector
            if (todo >= 1) {
                // Calculate gradient of systematic utility
                DUdB = Xnr

                // Box-Cox transformation coefficients
                if (lsl_ufunc == "boxcox") {
                    // Build interaction terms
                    BcCx  = (CX, BcL1, (nlei == 2 ? BcL2 : J(c, 0, 0)),
                             J(c, bfix - ncons, 0))
                    BcL1x = (J(c, ncx, 0), BcC, J(c, (nlei == 2), 0), LX1,
                             J(c, nlx2, 0), (nlei == 2 ? BcL2 : J(c, 0, 0)),
                             J(c, nxind, 0))
                    BcL2x = (nlei == 2 ? (J(c, ncx + 1, 0), BcC, J(c, nlx1, 0),
                                          LX2, BcL1, J(c, nxind, 0))
                                       : J(c, 0, 0))

                    // Calculate gradients
                    DUdlC  = lsl_boxcox_g(C,  lC)  :* cross(BcCx', Beta')
                    DUdlL1 = lsl_boxcox_g(L1, lL1) :* cross(BcL1x', Beta')
                    DUdlL2 = (nlei == 2 ? lsl_boxcox_g(L2, lL2) :*
                                           cross(BcL2x', Beta')
                                        : J(c, 0, 0))
                    DUdlam = (DUdlC, DUdlL1, DUdlL2)
                } else DUdlam = J(c, 0, 0)

                // Full Maximum Likelihood estimation?
                real matrix DUdwage, DWdBwcorr
                if (lsl_joint) {
                    if (lsl_ufunc == "quad") {
                        DUdC = cross((CX, 2 :* C2X :* C, L1, L2)',
                                     Beta[|1\ncons|]')
                    } else if (lsl_ufunc == "boxcox") {
                        DUdC = cross(((CX, BcL1, (nlei == 2 ? BcL2 : J(c, 0, 0)))
                                       :* (reldif(lC, 0) >= 1e-25
                                           ? C:^(lC - 1) :/ lsl_boxcc
                                           : lsl_boxcc :/ C))',
                                     Beta[|1\ncons|]')
                    } else if (lsl_ufunc == "tran") {
                        DUdC = cross(((CX, 2 :* C2X :* log(C), log(L1),
                                       log(L2)) :/ C)', Beta[|1\ncons|]')
                    }
                    DCdM = cross((J(c, 1, 1), 2 :* Mwage,
                                  lsl_TaxregIas1[|i,1\e,.|],
                                  2 :* Mwage :* lsl_TaxregIas1[|i,1\e,.|])',
                                  lsl_TaxregB[|1\2 + 2 * cols(lsl_TaxregIas1)|]')
                    DMdH = (lsl_Days[|i,1\e,1|] :/ 12 :/ 7) :*
                            lsl_Hours[|i,1\e,.|]

                    if (lsl_wagecorr) {
                        DWdBw   = Wn :* (lsl_WageVars[|i,1\e,.|] :-
                                         lsl_residanchor :*
                                          cross(lsl_Wobs[|i,1\e,1|],
                                                lsl_WageVars[|i,1\e,.|]))
                        DWdBsig = Wn :* (Bsig + lsl_R[iRV,1] *
                                          (colsum(lsl_Wobs[|i,1\e,1|]) == 0 |
                                           lsl_residanchor == 0))
                        DWdBwcorr = DUdB[.,lsl_Wcorrvars] :* lsl_R[iRV,1]
                    } else {
                        DWdBw = Wn :* (lsl_WageVars[|i,1\e,.|] :-
                                       cross(LnWresPur, lsl_WageVars) :/
                                        (colsum(lsl_Wobs) - bwage))
                        if (lsl_wagep) {
                            DWdBw = DWdBw :- Wn :*
                                     ((colsum(lsl_Wobs[|i,1\e,1|]) == 0 |
                                       lsl_residanchor == 0) :*
                                      cross((lsl_Wpred[|i,1\e,.|] :*
                                             lsl_R[iRV,1] :/ SigmaW)',
                                            cross(LnWresPur, lsl_WageVars) :/
                                            (colsum(lsl_Wobs) - bwage)) :+
                                      (colsum(lsl_Wobs[|i,1\e,1|]) == 1 &
                                       lsl_residanchor == 1) :*
                                      colsum(lsl_Wpred[|i,1\e,.|] :*
                                       lsl_Wobs[|i,1\e,1|] :*
                                       lsl_WageVars[|i,1\e,.|]))
                        }
                        DWdBsig   = J(c, 0, 0)
                        DWdBwcorr = J(c, 0, 0)
                    }

                    DUdwage = DUdC :* DCdM :* DMdH :* (DWdBw, DWdBsig)

                    if (lsl_residanchor & lsl_wagecorr) {
                        DUdwage[|1,1\c,bwage|] = DUdwage[|1,1\c,bwage|] :-
                                                 cross(cross(DUdB[.,lsl_Wcorrvars]',
                                                             B[|iheck + 1\iheck + bheck - 1|]')',
                                                       cross(lsl_Wobs[|i,1\e,1|],
                                                             lsl_WageVars[|i,1\e,.|])) :/ SigmaW
                        DUdwage[|1,bwage + 1\c,bwage + 1|] = DUdwage[|1,bwage + 1\c,bwage + 1|] :-
                                                             cross(DUdB[.,lsl_Wcorrvars]',
                                                                   B[|iheck + 1\iheck + bheck - 1|]') :*
                                                              colsum(LnWresPur[|i,1\e,1|]) :/ SigmaW:^2
                    }
                } else {
                    DUdwage   = J(c, 0, 0)
                    DWdBwcorr = J(c, 0, 0)
                }

                // Random components
                if (brnd > 0) {
                    DUdBr = (lsl_corr == 1
                               ? cross(DUdB[.,vech(J(1, rvars, lsl_Rvars))]',
                                       diag(vech(J(rvars, 1, lsl_R[|iRV,1 + lsl_wagep * nlei\iRV,lsl_wagep * nlei + rvars|]))))
                               : DUdB[.,lsl_Rvars] :* lsl_R[|iRV,1 + lsl_wagep * nlei\iRV,lsl_wagep * nlei + rvars|])
                } else DUdBr = J(c, 0, 0)

                // Total
                DUdx = (DUdB, DUdwage, DWdBwcorr, DUdlam, DUdBr)
                Gsum = Gsum + pni :* cross(YmPn, DUdx)
            }

            // Calculate Hessian matrix
            if (todo == 2) {
                // Won't work with joint wage estimation!!!

                // Utility
                YmPn_D2UdB2 = J(bfix + blam, bfix + blam, 0)

                // Random components
                YmPn_D2UdBdBr = J(brnd, bfix + blam, 0)
                YmPn_D2UdBr2  = J(brnd, brnd, 0)

                // Box-Cox transformation parameters
                if (lsl_ufunc == "boxcox") {
                    YmPn_D2UdB2[|ilC,1\ilC,ncons|] = cross(YmPn :* lsl_boxcox_g(C, lC), (CX, BcL1, (nlei == 2 ? BcL2 : J(c, 0, 0))))
                    YmPn_D2UdB2[ilC,ilC] = cross(YmPn :* lsl_boxcox_h(C, lC), cross((CX, BcL1, (nlei == 2 ? BcL2 : J(c, 0, 0)))', Beta[|1\ncons|]'))
                    YmPn_D2UdB2[ilL1,ncx + 1] = cross(YmPn :* lsl_boxcox_g(L1, lL1), BcC)
                    YmPn_D2UdB2[|ilL1,ncons + 1\ilL1,ncons + nlx1|] = cross(YmPn :* lsl_boxcox_g(L1, lL1), LX1)
                    YmPn_D2UdB2[ilL1,ilC]  = cross(YmPn :* lsl_boxcox_g(L1, lL1), Beta[ncons - nlei + 1] :* lsl_boxcox_g(C, lC))
                    YmPn_D2UdB2[ilL1,ilL1] = cross(YmPn :* lsl_boxcox_h(L1, lL1), cross((BcC, LX1, (nlei == 2 ? BcL2 : J(c, 0, 0)))',
                                                                                        (Beta[ncx + 1], Beta[|ncons + 1\ncons + nlx1|], (nlei == 2 ? Beta[ncons + nlx1 + nlx2 + 1] : J(1, 0, 0)))'))
                    if (nlei == 2) {
                        YmPn_D2UdB2[ilL1,ncons + nlx1 + nlx2 + 1] = cross(YmPn :* lsl_boxcox_g(L1, lL1), BcL2)
                        YmPn_D2UdB2[ilL2,ncons + nlx1 + nlx2 + 1] = cross(YmPn :* lsl_boxcox_g(L2, lL2), BcL1)
                        YmPn_D2UdB2[ilL2,ncons] = cross(YmPn :* lsl_boxcox_g(L2, lL2), BcC)
                        YmPn_D2UdB2[|ilL2,ncons + nlx1 + 1\ilL2,ncons + nlx1 + nlx2|] = cross(YmPn :* lsl_boxcox_g(L2, lL2), LX2)
                        YmPn_D2UdB2[ilL2,ilC]  = cross(YmPn :* lsl_boxcox_g(L2, lL2), Beta[ncons] :* lsl_boxcox_g(C, lC))
                        YmPn_D2UdB2[ilL2,ilL1] = cross(YmPn :* lsl_boxcox_g(L2, lL2), Beta[ncons + nlx1 + nlx2 + 1] :* lsl_boxcox_g(L1, lL1))
                        YmPn_D2UdB2[ilL2,ilL2] = cross(YmPn :* lsl_boxcox_h(L2, lL2), cross((BcC, LX2, BcL1)', (Beta[ncons], Beta[|ncons + nlx1 + 1\ncons + nlx1 + nlx2 + 1|])'))
                    }

                    // Random coefficients?
                    if (brnd > 0) {
                        if (lsl_corr == 1) {
                            YmPn_D2UdBdBr[.,ilC]  = cross(cross(YmPn :* lsl_boxcox_g(C, lC), BcCx)[.,vech(J(1, rvars, lsl_Rvars))]', diag(vech(J(rvars, 1, lsl_R[|iRV,1\iRV,rvars|]))))'
                            YmPn_D2UdBdBr[.,ilL1] = cross(cross(YmPn :* lsl_boxcox_g(L1, lL1), BcL1x)[.,vech(J(1, rvars, lsl_Rvars))]', diag(vech(J(rvars, 1, lsl_R[|iRV,1\iRV,rvars|]))))'
                            if (nlei == 2) YmPn_D2UdBdBr[.,ilL2] = cross(cross(YmPn :* lsl_boxcox_g(L2, lL2), BcL2x)[.,vech(J(1, rvars, lsl_Rvars))]', diag(vech(J(rvars, 1, lsl_R[|iRV,1\iRV,rvars|]))))'
                        } else {
                            YmPn_D2UdBdBr[.,ilC]  = (cross(YmPn :* lsl_boxcox_g(C, lC), BcCx)[.,lsl_Rvars] :* lsl_R[|iRV,1\iRV,rvars|])'
                            YmPn_D2UdBdBr[.,ilL1] = (cross(YmPn :* lsl_boxcox_g(L1, lL1), BcL1x)[.,lsl_Rvars] :* lsl_R[|iRV,1\iRV,rvars|])'
                            if (nlei == 2) YmPn_D2UdBdBr[.,ilL2] = (cross(YmPn :* lsl_boxcox_g(L2, lL2), BcL2x)[.,lsl_Rvars] :* lsl_R[|iRV,1\iRV,rvars|])'
                        }
                    }
                }

                // Partial second derivatives
                YmPn_D2Udx2 = makesymmetric((YmPn_D2UdB2  , YmPn_D2UdBdBr' \
                                             YmPn_D2UdBdBr, YmPn_D2UdBr2   ))

                // Full Maximum Likelihood estimation?   // BUGGY
                real matrix YmPn_D2UdBdwage, YmPn_D2UdBrdwage, YmPn_D2Udwage2
                YmPn_D2UdBrdwage = J(brnd, bwage + lsl_wagecorr + (lsl_wagecorr > 0), 0)
                YmPn_D2UdBdwage = J(bwage + lsl_wagecorr + (lsl_wagecorr > 0), bfix, 0)
                YmPn_D2Udwage2  = J(bwage + lsl_wagecorr + (lsl_wagecorr > 0), bwage + lsl_wagecorr + (lsl_wagecorr > 0), 0)

                if (lsl_joint == 1) {   // BUGGY
                    // Prepare stubs
                    if (lsl_ufunc == "quad") {
                        YmPn_D2UdBdwage = cross(YmPn :* DCdM :* DMdH :* (DWdBw, DWdBsig, DWdBwcorr), (CX, 2 :* C2X :* C, L1, L2, J(c, cols(Xnr) - ncons, 0)))
                        D2UdC2          = 2 :* cross(C2X', Beta[|ncx + 1\ncx + nc2x|]')
                    } else if (lsl_ufunc == "tran") {
                        YmPn_D2UdBdwage = cross(YmPn :* DCdM :* DMdH :* (DWdBw, DWdBsig, DWdBwcorr), ((CX, 2 :* C2X :* log(C), log(L1), log(L2)) :/ C, J(c, cols(Xnr) - ncons, 0)))
                        D2UdC2          = - cross(((CX, 2 :* C2X :* (log(C) :- 1), log(L1), log(L2)) :/ C:^2)', Beta[|1\ncons|]')
                    }
                    D2CdM2 = cross((J(c, 1, 2), 2 :* lsl_TaxregIas1[|i,1\e,.|])', (lsl_TaxregB[2], lsl_TaxregB[|2 + cols(lsl_TaxregIas1) + 1\2 + 2 * cols(lsl_TaxregIas1)|])')
                    D2MdH2 = 0

                    // Calculate second and cross derivatives
                    if (lsl_wagecorr) {
                        // d2U / dBw2
                        YmPn_D2Udwage2[|1,1\bwage,bwage|] = cross(YmPn :* DWdBw :* (DMdH:^2 :* (D2UdC2 :* DCdM:^2 :+ DUdC :* D2CdM2) + DUdC :* DCdM :* D2MdH2), DWdBw) :+
                                                            cross(YmPn :* DWdBw :* DUdC :* DCdM :* DMdH, (lsl_WageVars[|i,1\e,.|] :- lsl_residanchor :* Bsig :* colsum(lsl_Wobs[|i,1\e,1|] :* lsl_WageVars[|i,1\e,.|])))
                        // d2U / dBsig2
                        YmPn_D2Udwage2[bwage + 1,bwage + 1] = cross(YmPn :* DWdBsig :* (DMdH:^2 :* (D2UdC2 :* DCdM:^2 :+ DUdC :* D2CdM2) + DUdC :* DCdM :* D2MdH2), DWdBsig) :+
                                                              cross(YmPn :* DUdC :* DCdM :* DMdH, DWdBsig :* (lsl_R[iRV,1] :+ Bsig) :+ Wn :* lsl_R[iRV,1])
                        /*
                        if      (lsl_wagecorr == 1) DWdBwcorr = Wn :* (rowsum(lsl_R[|iRV,1\iRV,rvars|] :* ((lsl_Rvars :== iC) :+ (lsl_Rvars :== iL1))') :+ B[iheck + 1])
                        else if (lsl_wagecorr == 2) DWdBwcorr = cross(Wn', ((rowsum(lsl_R[|iRV,1\iRV,rvars|] :* (lsl_Rvars :== iC )') :+ B[iheck + 1]),
                                                    (rowsum(lsl_R[|iRV,1\iRV,rvars|] :* (lsl_Rvars :== iL1)') :+ B[iheck + 2])))
                        YmPn_D2Udwage2[|bwage + 1,bwage + 1\bwage + lsl_wagecorr,bwage + lsl_wagecorr|] = cross(YmPn :* DWdBsig :* (DMdH:^2 :* (D2UdC2 :* DCdM:^2 :+ DUdC :* D2CdM2) + DUdC :* DCdM :* D2MdH2), DWdBsig) :+
                                                                                                          cross(YmPn :* DUdC :* DCdM :* DMdH, DWdBsig :* Bsig :+ Wn)
                        */
                        // d2U / dBw dBsig
                        YmPn_D2Udwage2[|bwage + 1,1\bwage + 1,bwage|] = cross(YmPn :* DWdBsig :* (DMdH:^2 :* (D2UdC2 :* DCdM:^2 :+ DUdC :* D2CdM2) + DUdC :* DCdM :* D2MdH2), DWdBw) :+
                                                                        cross(YmPn :* DUdC :* DCdM :* DMdH, DWdBsig :* (lsl_WageVars[|i,1\e,.|] :- lsl_residanchor :* Bsig :* colsum(lsl_Wobs[|i,1\e,1|] :* lsl_WageVars[|i,1\e,.|])) :- lsl_residanchor :* cross(Wn', colsum(lsl_Wobs[|i,1\e,1|] :* lsl_WageVars[|i,1\e,.|])))
                        // d2U / dBw dBwcorr
                        YmPn_D2Udwage2[|bwage + 2,1\bwage + 1 + lsl_wagecorr,bwage|] = cross(YmPn :* DWdBsig :* (DMdH:^2 :* (D2UdC2 :* DCdM:^2 :+ DUdC :* D2CdM2) + DUdC :* DCdM :* D2MdH2), DWdBw) :+
                                                                                       cross(YmPn :* DUdC :* DCdM :* DMdH, DWdBwcorr :* (lsl_WageVars[|i,1\e,.|] :- lsl_residanchor :* Bsig :* colsum(lsl_Wobs[|i,1\e,1|] :* lsl_WageVars[|i,1\e,.|])))
                    } else {
                        // Buggy
                    }

                    // Partial second derivatives
                    YmPn_D2Udx2 = makesymmetric((YmPn_D2UdB2    , YmPn_D2UdBdwage', YmPn_D2UdBdBr'    \
                                                 YmPn_D2UdBdwage, YmPn_D2Udwage2  , YmPn_D2UdBrdwage' \
                                                 YmPn_D2UdBdBr  , YmPn_D2UdBrdwage, YmPn_D2UdBr2))
                }

                // Total second derivatives
                H1sum = H1sum :+ pni :* (cross(Yn, DUdx) :- cross(Pnr, DUdx))
                H2sum = H2sum :+ pni :* (cross(cross(Yn, DUdx) :- cross(Pnr, DUdx), cross(YmPn, DUdx)) :-
                                         cross(Pnr :* DUdx, DUdx :- cross(Pnr, DUdx)) :+ YmPn_D2Udx2)
            }

            // Dude correction
            if (lsl_lambda > 0 & lsl_joint == 0) {
                lnf = lnf + lsl_lambda :* lsl_Weight[i] :* cross(lsl_force[|i\e|], Dude[|i,r\e,r|])
                if (todo >= 1) {
                    if    (lsl_ufunc == "boxcox") G[|1\ncons|] = G[|1\ncons|]   // BUGGY!!!
                    else if (lsl_ufunc == "quad") G[|1\ncons|] = G[|1\ncons|] :+ lsl_lambda :* colsum(lsl_Weight[i] :* lsl_force[|i\e|] :* (CX, 2 :* C2X :* C, L1, L2))
                    else if (lsl_ufunc == "tran") G[|1\ncons|] = G[|1\ncons|] :+ lsl_lambda :* colsum(lsl_Weight[i] :* lsl_force[|i\e|] :* (CX, 2 :* C2X :* log(C), log(L1), log(L2)) :/ C)
                }
                if (todo == 2 & lsl_ufunc == "boxcox") {
                    // BUGGY!!!
                }
            }
        }

        // Prevent likelihood from becoming exactly zero
        lsum = max((lsum, 1e-25))

        // Add to overall statistics
        lnf = lnf + lsl_Weight[i] * log(lsum / lsl_draws)
        if (todo >= 1) G = G + lsl_Weight[i] * (lsum > 1e-25 ? Gsum / lsum
                                                             : J(1, b, 0))
        if (todo == 2) H = H + lsl_Weight[i] * (lsum > 1e-25
                                                  ? H2sum / lsum
                                                    - cross(Gsum, H1sum)
                                                       / lsum^2
                                                  : J(b, b, 0))
    }

    // Add likelihood of wage equation?
    // *
    if (lsl_joint) {
        lnf = lnf + cross(lsl_Wobs, log(normalden(LnWresPur :/ SigmaW))
                                    :- log(SigmaW))
        if (todo >= 1) {
            if (lsl_wagecorr) {
                G[|iwage\iwage + bwage - 1|] = G[|iwage\iwage + bwage - 1|] :+
                                               cross(LnWresPur, lsl_WageVars
                                                                :/ SigmaW:^2)
                G[iheck] = G[iheck] :+ cross(LnWresPur, LnWresPur) :/ SigmaW:^3 :-
                           colsum(lsl_Wobs :/ SigmaW)
            } else {
                G[|iwage\iwage + bwage - 1|] = G[|iwage\iwage + bwage - 1|] :+
                                               cross(LnWresPur, lsl_WageVars) :/
                                                SigmaW:^2 :-
                                               cross(LnWresPur, LnWresPur) :*
                                                cross(LnWresPur, lsl_WageVars) :/
                                                (colsum(lsl_Wobs) - bwage) :/
                                                SigmaW:^4 :+
                                               colsum(lsl_Wobs) :*
                                                cross(LnWresPur, lsl_WageVars) :/
                                                cross(LnWresPur, LnWresPur)
            }
        }
    }
    // *

    // Calculate dude share
    if (todo == 0 | lsl_lambda > 0) {
        st_numscalar("lsl_dudes", sum(Dude :< 0) / (nobs * lsl_draws))
    }
}


/**
 * Returns the transformed variable x^(l).
 *
 * @param  real matrix Var Variable x
 * @param  real scalar lam Power parameter l
 * @return real matrix Transformed variable
 */
real matrix lsl_boxcox(real matrix Var, real scalar lam) {
    return (reldif(lam, 0) >= 1e-25 ? (Var:^lam :- 1) :/ lam : log(Var))
}


/**
 * Returns the first derivative of x^(l) with respect to l, i.e., dx^(l)/dl.
 *
 * @param  real matrix Var Variable x
 * @param  real scalar lam Power parameter l
 * @return real matrix First derivative
 */
real matrix lsl_boxcox_g(real matrix Var, real scalar lam) {
    return (reldif(lam, 0) >= 1e-25 ? (Var:^lam :* (lam :* log(Var) :- 1) :+ 1) :/ lam^2 : 0.5 :* log(Var):^2)
}


/**
 * Returns the second derivative of x^(l) with respect to l, i.e., d2x^(l)/dl2.
 *
 * @param  real matrix Var Variable x
 * @param  real scalar lam Power parameter l
 * @return real matrix Second derivative
 */
real matrix lsl_boxcox_h(real matrix Var, real scalar lam) {
    return (reldif(lam, 0) >= 1e-25 ? (Var:^lam :* (lam^2 :* log(Var):^2 :- 2 :* lam :* log(Var) :+ 2) :- 2) :/ lam^3 : (1/3) :* log(Var):^3)
}


/**
 * Predicts systematic utilities or choice probabilities.
 *
 * @param string rowvector newvar Name(s) of the new variable(s)
 * @param string scalar touse Sample selection variable
 * @param string rowvector opt Prediction options (xb and/or pc1)
 */
void lslogit_p(string rowvector newvar, string scalar touse, string rowvector opt) {
    //external real rowvector lsl_B
    external real rowvector lsl_Bfix
    external real matrix    lsl_CholBW
    external real matrix    lsl_CholW
    external real matrix    lsl_CholB
    external real matrix    lsl_Brnd
    external real scalar    lsl_bfix
    external real scalar    lsl_nlei
    external real scalar    lsl_boxcc
    external real scalar    lsl_lC, lsl_lL1, lsl_lL2
    external real matrix    lsl_LnWres

    // Right hand side variables
    external real matrix    lsl_X

    // Number of groups
    external real scalar    lsl_groups

    // Number of random draws
    external real scalar    lsl_draws

    // Left hand side variable
    external real colvector lsl_Y

    // Panel setup information
    external real matrix    lsl_J

    //   Halton sequences
    external real matrix    lsl_R

    // Random coefficients
    external real colvector lsl_Rvars

    // Number of random coefficients
    external real scalar    lsl_rvars

    //   Enable correlation?
    external real scalar    lsl_corr

    // Functional form
    external string scalar  lsl_ufunc

    // Wage Prediction Error?
    external real scalar    lsl_wagep

    //   Prediction dummies
    external real matrix    lsl_Wpred

    //   Number of days per tax year
    external real colvector lsl_Days

    //   Hourly wage rates
    external real matrix    lsl_Hwage

    //   Variance of the wage regression
    external real rowvector lsl_Sigma

    //   Variance of the wage regression (buggy?)
    external real rowvector lsl_SigmaW

    //   Hours of work
    external real matrix    lsl_Hours

    // Tax Regression
    external real rowvector lsl_TaxregB

    //   Wage independent variables of tax regression
    external real matrix    lsl_TaxregVars

    //   Wage interaction variables of tax regression
    external real matrix    lsl_TaxregIas1

    //   Wage interaction variables of tax regression
    external real matrix    lsl_TaxregIas2

    external real matrix    lsl_taxreg_rmse

    // To round, or not to round.
    external real scalar    lsl_round
    external real scalar    lsl_joint
    external real scalar    lsl_residanchor
    external real matrix    lsl_LnWres
    external real matrix    lsl_Wobs
    external real colvector lsl_Weight

    external real colvector lsl_C
    external real matrix    lsl_CX
    external real matrix    lsl_C2X
    external real colvector lsl_L1
    external real matrix    lsl_LX1
    external real matrix    lsl_L2X1
    external real matrix    lsl_LX2
    external real matrix    lsl_L2X2
    external real matrix    lsl_Xind

    // Buggy, L2 is colvector in fact
    external real matrix    lsl_L2

    real scalar    nobs, n, r, i, e, c, iRV, nRV, wp, ncons, lsum, lnf,
                   getutils, getprobs, getdudes
    real colvector U, P, D, Un, Pn, Unr, Enr, Pnr, C, L1, BcC, BcL1, Yn
    real rowvector Beta, Zeta, Bsig
    real matrix    Xnr, CX, C2X, LX1, L2X1, L2, BcL2, LX2, L2X2, Xind,
                   Wn, Mwage, TaxregX1, TaxregX2, TaxregX, Dude

    // Indicates first observation of active group
    i = 1

    // Number of observations
    nobs  = rows(lsl_Y)
    nRV   = lsl_rvars + 1
    ncons = cols(lsl_CX) + cols(lsl_C2X) + 1 + cols(lsl_L2)

    // Initialize
    lnf = 0
    U   = J(nobs, 1, 0)
    P   = J(nobs, 1, 0)
    D   = J(nobs, 1, 0)

    // Which options have been selected?
    getutils = (sum(opt :== "xb")    == 1)
    getprobs = (sum(opt :== "pc1")   == 1)
    getdudes = (sum(opt :== "dudes") == 1)

    // Need to calculate Dude shares? Initialize
    if (getdudes == 1) Dude = J(nobs, lsl_draws, 0)

    // Initialize random coefficients vector
    if (lsl_rvars > 0) {
        Zeta = J(rows(lsl_R), lsl_bfix, 0)
        //Zeta[.,lsl_Rvars] = cross(lsl_R[|1,1\.,lsl_rvars|]', lsl_CholB')
    }

    // Loop over individuals
    for (n = 1; n <= lsl_groups; n++) {
        i   = lsl_J[n,1]
        e   = lsl_J[n,2]
        c   = e - i + 1
        Yn  = lsl_Y[|i,1\e,.|]
        Xnr = lsl_X[|i,1\e,.|]

        // Fetch needed right hand side parts
        if (lsl_wagep == 1 | getdudes == 1 | lsl_joint) {
            C    =  lsl_C[|i\e|]             // Get consumption from data
            CX   = (cols(lsl_CX)   > 0 ?   lsl_CX[|i,1\e,.|] : J(c, 0, 0))
            C2X  = (cols(lsl_C2X)  > 0 ?  lsl_C2X[|i,1\e,.|] : J(c, 0, 0))
            L1   = lsl_L1[|i\e|]
            LX1  = (cols(lsl_LX1)  > 0 ?  lsl_LX1[|i,1\e,.|] : J(c, 0, 0))
            L2X1 = (cols(lsl_L2X1) > 0 ? lsl_L2X1[|i,1\e,.|] : J(c, 0, 0))
            L2   = (cols(lsl_L2)   > 0 ?   lsl_L2[|i\e|]     : J(c, 0, 0))
            LX2  = (cols(lsl_LX2)  > 0 ?  lsl_LX2[|i,1\e,.|] : J(c, 0, 0))
            L2X2 = (cols(lsl_L2X2) > 0 ? lsl_L2X2[|i,1\e,.|] : J(c, 0, 0))
            Xind = (cols(lsl_Xind) > 0 ? lsl_Xind[|i,1\e,.|] : J(c, 0, 0))
            Wn   = lsl_Hwage[|i,1\e,.|]

            // Transform consumption and leisure
            if (lsl_ufunc == "boxcox") {
                BcC  = lsl_boxcox(C, lsl_lC)
                BcL1 = lsl_boxcox(L1, lsl_lL1)
                BcL2 = lsl_boxcox(L2, lsl_lL2)
            }
        }
        wp = (lsl_wagep & sum(lsl_Wpred[|i,1\e,.|]) > 0)

        // Wage residual anchor? Use actual wage equation residuals
        //   instead of random draws
        if (lsl_residanchor & lsl_joint & wp
          & colsum(lsl_Wobs[|i,1\e,1|]) == 1) {
            lsl_R[|lsl_draws * (n - 1) + 1,1\lsl_draws * (n - 1) + lsl_draws,1|] =
                    J(lsl_draws, 1, colsum(lsl_LnWres[|i,1\e,1|]) / lsl_CholW)
        }

        lsum = 0
        Un   = J(c, 1, 0)
        Pn   = J(c, 1, 0)

        // Loop over draws
        for (r = 1; r <= lsl_draws; r++) {
            // Indicates the active Halton sequence
            iRV = lsl_draws * (n - 1) + r

            // Build (random?) coefficients matrix (DEBUG HERE!)
            if (lsl_rvars > 0) {
                Zeta[iRV,lsl_Rvars] =
                    cross(lsl_R[|iRV,1\iRV,lsl_wagep * lsl_nlei + lsl_rvars|]',
                          (lsl_CholBW', lsl_CholB)')
            }

            // Build (random?) coefficients vector
            Beta = lsl_Bfix :+ (lsl_rvars > 0 ? Zeta[iRV,.] : 0)

            if (wp | lsl_joint) {

                //
                // Calculate monthly earnings
                //

                // Adjust wages with random draws if prediction enabled
                if (lsl_wagep) {
                    Wn = lsl_Hwage[|i,1\e,.|] :* exp(cross(lsl_CholW', lsl_R[|iRV,1\iRV,lsl_nlei|]')' :*
                                lsl_Wpred[|i,1\e,.|])
                }

                // Calculate monthly earnings
                Mwage = (lsl_Days[|i\e|] :/ 12 :/ 7) :* lsl_Hours[|i,1\e,.|] :* Wn

                // Round monthly earnings if enabled
                if (lsl_round) Mwage = round(Mwage, 0.01)

                //
                // Predict disposable income
                //

                // Fill matrix of independent variables for dpi prediction
                TaxregX1 = (Mwage[.,1], Mwage[.,1]:^2, lsl_TaxregIas1[|i,1\e,.|] :* Mwage[.,1],
                                                       lsl_TaxregIas1[|i,1\e,.|] :* Mwage[.,1]:^2)
                if (cols(L2) == 1) TaxregX2 = (Mwage[.,2], Mwage[.,2]:^2, lsl_TaxregIas2[|i,1\e,.|] :* Mwage[.,2],
                                                                          lsl_TaxregIas2[|i,1\e,.|] :* Mwage[.,2]:^2)
                else           TaxregX2 = J(c, 0, 0)
                TaxregX = (TaxregX1, TaxregX2, lsl_TaxregVars[|i,1\e,.|], J(c, 1, 1))

                // Predict disposable income (can't be negative!)
                C = rowmax((cross(TaxregX', lsl_TaxregB') :+ (lsl_taxreg_rmse :* lsl_R[iRV,cols(lsl_R)]), J(c, 1, 1))) :/ lsl_boxcc

                // Build matrix with independent variables
                if      (lsl_ufunc == "tran") Xnr = (log(C)  :* (CX,  log(C)  :* C2X,   log(L1),   log(L2)),
                                                     log(L1) :* (LX1, log(L1) :* L2X1),
                                                     log(L2) :* (LX2, log(L2) :* L2X2), log(L1) :* log(L2), Xind)
                else if (lsl_ufunc == "quad") Xnr = (C  :* (CX,  C  :* C2X,   L1,   L2),
                                                     L1 :* (LX1, L1 :* L2X1),
                                                     L2 :* (LX2, L2 :* L2X2), L1 :* L2, Xind)
                else if (lsl_ufunc == "boxcox") {
                    BcC = lsl_boxcox(C, lsl_lC)
                    Xnr = ((CX, BcL1, BcL2) :* BcC, LX1 :* BcL1, LX2 :* BcL2, BcL1 :* BcL2, Xind)
                }

            }

            // Calculate systematic utility and choice probability
            Unr = cross(Xnr', Beta')
            if (getprobs == 1) {
                // Standardize to avoid missings
                Enr = exp(Unr :+ colmin(-mean(Unr) \ 700 :- colmax(Unr)))

                // Calculate choice probability
                Pnr = Enr :/ colsum(Enr)

                // Recall "predicted" probability that choice is chosen
                lsum = lsum + cross(Yn, Pnr)
            }

            // Calculate dudes
            if (getdudes == 1) {
                if (lsl_ufunc == "boxcox") {
                    Dude[|i,r\e,r|] = cross((CX, lsl_boxcox(L1, lsl_lL1),
                                             lsl_boxcox(L2, lsl_lL2))',
                                            Beta[|1\ncons|]') :*
                                       (reldif(lsl_lC, 0) >= 1e-25
                                          ? C:^(lsl_lC - 1) : (1 :/ C))
                } else if (lsl_ufunc == "quad") {
                    Dude[|i,r\e,r|] = cross((CX, 2 :* C2X :* C, L1, L2)',
                                            Beta[|1\ncons|]')
                } else if (lsl_ufunc == "tran") {
                    Dude[|i,r\e,r|] = cross((CX, 2 :* C2X :* log(C), log(L1),
                                             log(L2))', Beta[|1\ncons|]') :/ C
                }
            }

            // Sum up
            if (getutils == 1) Un = Un :+ Unr
            if (getprobs == 1) Pn = Pn :+ Pnr
        }

        // Add up to "predicted" log-likelihood
        lnf = lnf + lsl_Weight[i] * log(max((lsum, 1e-25)) / lsl_draws)

        // Next household
        if (getutils == 1) U[|i\e|] = Un :/ lsl_draws
        if (getprobs == 1) P[|i\e|] = Pn :/ lsl_draws
        if (getdudes == 1) D[|i\e|] = rowsum(Dude[|i,1\e,.|] :< 0) :/ lsl_draws
    }

    // Print calculated ("predicted") log-likelihood
    if (getprobs == 1) {
        if (lsl_joint) lnf = lnf + cross(lsl_Wobs, log(normalden(lsl_LnWres :/ lsl_SigmaW))
                                                   :- log(lsl_SigmaW))
        st_numscalar("lsl_ll_p", lnf)
        round(lnf, .0001)
    }

    // Store prediction
    real matrix result
    real scalar a
    result = J(nobs, cols(opt), 0)
    for (a = 1; a <= cols(opt); a++) {
        if      (opt[a] == "pc1")   result[.,a] = P
        else if (opt[a] == "xb")    result[.,a] = U
        else if (opt[a] == "dudes") result[.,a] = D
    }
    if (opt[cols(opt)] == "wages") result = result[|1,1\.,cols(opt)-1|]
    st_store(., newvar, touse, result)
}

mata set matastrict off
end


cap program drop lslogit_p
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

    mata: lsl_CholBW = J(lsl_nlei, lsl_rvars, 0)
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


cap program drop lslogit_cov
/**
 * Conditional Logit but integrating out wage prediction errors (Wrapper programm)
 */
program define lslogit_cov
    version 12
    syntax [, cov NOHEADer NOWAGE NOPREF post]

    // Predict command works only with lslogit estimates
    if ("`e(cmd)'" != "lslogit") error 301

    // Cov calculation appropriate?
    if ("`e(randvars)'" == "") {
        di in r "no random terms in model setup"
        exit 498
    }
    if (`e(corr)' == 0 & `e(wagecorr)' == 0) {
        di in r "no covariance terms estimated"
        exit 498
    }


    //
    // Build estimation command
    //

    local expline

    // Add terms of variance-covariance matrix
    if ("`nopref'" == "") {
        local n_coef : word count `e(randvars)'
        forval c = 1/`n_coef' {
            local coef_c : word `c' of `e(randvars)'
            forval r = `c'/`n_coef' {
                local coef_r : word `r' of `e(randvars)'
                local label = cond(`r' == `c', "var_`coef_c'", "cov_`coef_c'_`coef_r'")
                if (`c' == 1) local expline `expline' (`label': _b[/l_`coef_c'_`coef_r'] * _b[/l_`coef_c'_`coef_c'])
                else {
                    forval i = 1/`c' {
                        local coef_i : word `i' of `e(randvars)'
                        if      (`i' == 1)   local expline `expline' (`label': _b[/l_`coef_i'_`coef_c'] * _b[/l_`coef_i'_`coef_r'] +
                        else if (`i' == `c') local expline `expline'           _b[/l_`coef_c'_`coef_r'] * _b[/l_`coef_c'_`coef_c'])
                        else                 local expline `expline'           _b[/l_`coef_i'_`coef_c'] * _b[/l_`coef_i'_`coef_r'] +
                    }
                }
            }
        }
    }

    // Add wage correlation terms?
    if (`e(wagecorr)' & "`nowage'" == "") {
        local iC  : word count `e(cx)' C
        local iL1 : word count `e(cx)' C C2 `e(leisure)' `e(lx1)' L1
        if (strpos(" `e(randvars)' ", " `iC' ") > 0 & strpos(" `e(randvars)' ", " `iL1' ") > 0) {
            local expline `expline' (cov_C_W1 : _b[/l_`iC'_`iC'] * _b[/l_C_W1]) ///
                                    (cov_L1_W1: _b[/l_`iC'_`iL1'] * _b[/l_C_W1] + _b[/l_`iL1'_`iL1'] * _b[/l_L1_W1]) ///
                                    (var_W1   : _b[/l_C_W1]^2 + _b[/l_L1_W1]^2 + _b[/l_W1_W1]^2)
        }
        else if (strpos(" `e(randvars)' ", " `iC' ") == 0 & strpos(" `e(randvars)' ", " `iL1' ") > 0) {
            local expline `expline' (cov_L1_W1: _b[/l_`iL1'_`iL1'] * _b[/l_L1_W1]) (var_W1: _b[/l_L1_W1]^2 + _b[/l_W1_W1]^2)
        }
        else if (strpos(" `e(randvars)' ", " `iC' ") > 0 & strpos(" `e(randvars)' ", " `iL1' ") == 0) {
            local expline `expline' (cov_C_W1 : _b[/l_`iC'_`iC'] * _b[/l_C_W1]) (var_W1: _b[/l_C_W1]^2 + _b[/l_W1_W1]^2)
        }
    }

    nlcom `expline', `noheader' `post'
end

***
