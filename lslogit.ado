*! version 0.4-3, 04mar2016, Max Loeffler <loeffler@zew.de>
/**
 * LSLOGIT - ESTIMATING MIXED LOGIT LABOR SUPPLY MODELS WITH STATA
 * 
 *
 * Copyright (C) 2012-2016 Max Löffler <loeffler@zew.de>
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


// Load lslogit Mata library (from .mlib-file)
cap mata: mata drop lslogit_d2()
cap mata: mata drop lslogit_p()
cap mata: mata drop lsl_boxcox()
cap mata: mata drop lsl_boxcox_g()
cap mata: mata drop lsl_boxcox_h()
cap mata: mata mlib index


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
                        if `touse', group(`group') `initopt'  ///
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
