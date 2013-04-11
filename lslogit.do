/*****************************************************************************
 *
 * lslogit -- ESTIMATING MIXED LOGIT LABOR SUPPLY MODELS WITH STATA
 * 
 * (c) 2012-2013 Max Löffler
 *
 *****************************************************************************/


cap program drop lslogit
/**
 * Wrapper programm
 */
program define lslogit
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
    syntax [, Level(integer `c(level)') Quiet]
    
    // Set up auxiliary stuff
    local diparm
    if ("`quiet'" == "") {
        // Auxiliary estimation results
        foreach aux in sigma_w1 sigma_w2 dudes {
            if (e(`aux') != .) {
                local val = e(`aux')
                local diparm `diparm' diparm(__lab__, value(`val') label("[`aux']"))
            }
        }
        
        // Model description
        local footnote "Model:"
        if      ("`e(ufunc)'" == "quad")   local footnote "`footnote' - Quadratic utility function"
        else if ("`e(ufunc)'" == "tran")   local footnote "`footnote' - Translog utility function"
        else if ("`e(ufunc)'" == "boxcox") local footnote "`footnote' - Box-Cox utility function"
        if ("`e(randvars)'" != "")         local footnote "`footnote'" _newline "       - Random coefficients (`e(randvars)')"
        if ("`e(corr)'" == "1")            local footnote "`footnote' with correlation"
        if ("`e(wagep)'" == "1")           local footnote "`footnote'" _newline "       - Wage prediction error integrated out"
        if ("`e(draws)'" != "1")           local footnote "`footnote'" _newline "       - Approximated using `e(draws)' Halton sequences"
        if ("`e(heckman)'" != "")          local footnote "`footnote'" _newline "       - Joint labor supply and wage estimation"
        if ("`e(wagecorr)'" == "1")        local footnote "`footnote' with correlation"
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
    syntax varname(numeric) [if] [in] [fweight/], GRoup(varname numeric)                                                    ///
                                                  Consumption(varname numeric) Leisure(varlist numeric min=1 max=2)         ///
                                                  [BOXCox QUADratic TRANslog                                                ///
                                                   cx(varlist numeric)  lx1(varlist numeric)  lx2(varlist numeric)          ///
                                                   c2x(varlist numeric) l2x1(varlist numeric) l2x2(varlist numeric)         ///
                                                   INDeps(varlist) TOTALTime(integer 80) DAYs(varname numeric)              ///
                                                   boxcc(integer 1000) boxcl(integer 80) HWage(varlist numeric min=1 max=2) ///
                                                   TAXReg(name) tria1(varlist numeric) tria2(varlist numeric)               ///
                                                   WAGEPred(varlist numeric min=1 max=2) HECKSIGma(numlist min=1 max=2)     ///
                                                   RANDvars(string) corr DRaws(integer 50) burn(integer 15)                 ///
                                                   HECKMan(varlist) SELect(varlist) wagecorr                                ///
                                                   noround Quiet Verbose lambda(real 0) force(varname numeric)              ///
                                                   difficult trace search(name) iterate(integer 100) method(name)           ///
                                                   gradient hessian debug Level(integer `c(level)') from(string)]
    
    /* INITIALIZE ESTIMATOR
     */
    
    // Mark the estimation sample
    marksample touse
    markout `touse' `varlist' `group' `consumption' `leisure' `cx' `lx1' `lx2' `c2x' `l2x1' `l2x2'  ///
                    `indeps' `wagepred' `days' `tria1' `tria2' `heckman' `select'
    
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
    if ("`boxcox'" != "" & "`translog'" == "" & "`quadratic'" == "") local ufunc "boxcox"
    if ("`boxcox'" == "" & "`translog'" != "" & "`quadratic'" == "") local ufunc "tran"
    if (("`boxcox'" != "") + ("`translog'" != "") + ("`quadratic'" != "") > 1) {
        di in r "utility function can be either 'quad', 'tran' or 'boxcox'"
        exit 498
    }
    // If translog, set up pre-text
    if ("`ufunc'" == "tran") {
        local ln  "ln"
        local pre "`ln'_"
    }
    // If translog or Box-Cox, check for zeros
    if (inlist("`ufunc'", "tran", "boxcox")) qui count if log(`consumption') == . & `touse'
    else                                     qui count if `consumption' < 0       & `touse'
    // Check for negative values
    if (r(N) > 0) {
        di in r "consumption contains values smaller or equal to zero"
        exit 498
    }
    // If Box-Cox, there is no sense in quadratic interactions
    if ("`ufunc'" == "boxcox" & ("`c2x'" != "" | "`l2x1'" != "" | "`l2x2'" != "")) {
        di in r "options c2x(), l2x1() and l2x2() not allowed with Box-Cox utility function"
        exit 498
    }
    // If Box-Cox, take it easy
    if ("`ufunc'" == "boxcox" & "`heckman'" != "") {
        di in r "all work and no play makes lslogit a dull command" _newline   ///
                "all work and no play makes lslogit a dull command" _newline   ///
                "all work and no play makes lslogit a dull command" _newline   ///
                "all work and no play makes lslogit a dull command" _newline   ///
                "all work and no play makes lslogit a dull command" _newline   ///
                "all work and no play makes lslogit a dull command" _newline   ///
                "all work and no play makes lslogit a dull command" _newline   ///
                "all work and no play makes lslogit a dull command" _newline   ///
                "all work and no play makes lslogit a dull command" _newline _newline   ///
                "seriously, take it easy, dude" _newline
        exit 498
    }
    
    // Validate wage estimation settings
    if ("`heckman'" != "") local joint joint
    if ("`heckman'" == "" & "`select'" != "") {
        di in r "option heckman() required when estimating jointly"
        exit 498
    }
    if ("`heckman'" != "" & "`hwage'" == "") {
        di in r "option hwage() required when estimating jointly"
        exit 498
    }
    // Check sample selection
    if ("`hwage'" != "") {
        foreach hw of local hwage {
            qui count if log(`hw') == . & `touse'
            if ("`heckman'" != "") {
                if ("`select'" != "" & r(N) == 0) {
                    di in r "wage variable `hw' never censored because of selection"
                    exit 498
                }
            }
            if (r(N) > 0 & "`heckman'" == "") {
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
    
    // Get variable count
    local n_leisure  : word count `leisure'
    local n_cxias    : word count `cx'
    local n_lx1ias   : word count `lx1'
    local n_lx2ias   : word count `lx2'
    local n_c2xias   : word count `c2x'
    local n_l2x1ias  : word count `l2x1'
    local n_l2x2ias  : word count `l2x2'
    local n_indeps   : word count `indeps'
    local n_randvars : word count `randvars'
    local n_wagep    : word count `wagepred'
    local n_hwage    : word count `hwage'
    local n_hecksig  : word count `hecksigma'
    local n_taxrias1 : word count `tria1'
    local n_taxrias2 : word count `tria2'
    local n_heckvars : word count `heckman'
    local n_selvars  : word count `select'
    
    // Validate wage correlation settings
    if ("`wagecorr'" != "") {
        // Check random variables
        local iC  : word count `cx' _cons
        local iL1 : word count `cx' _cons `c2x' `leisure' `lx1' _cons
        local iL2 : word count `cx' _cons `c2x' `leisure' `lx1' _cons `l2x1' `lx2' _cons // BUGGY for two guys?
        if ("`ufunc'" != "boxcox") {
            local iL1 = `iL1' + 1
            local iL2 = `iL2' + 2
        }
        local wagecorrC  = cond(strpos(" `randvars' ", " `iC' ")  > 0, 1, 0)
        local wagecorrL1 = cond(strpos(" `randvars' ", " `iL1' ") > 0, 1, 0)
        local wagecorrL2 = cond(strpos(" `randvars' ", " `iL2' ") > 0 & `n_leisure' == 2, 1, 0)
        local corrBW     = `wagecorrC' + `wagecorrL1' + `wagecorrL2'
        
        // At least one has to be enabled
        if (`corrBW' == 0) {
            di in r "either preferences for consumption or leisure have to be random when estimating wage correlation"
            exit 498
        }
        else local corr = "corr"                // Turn on preference correlation when wage correlation enabled
    }
    else {
        local corrBW     = 0
        local wagecorrC  = 0
        local wagecorrL1 = 0
        local wagecorrL2 = 0
    }
    
    // Validate Wage Prediction Options
    if (`n_wagep' == 0) local wagep = 0     // No wage prediction
    else {
        // Wage prediction enabled
        if (`n_wagep' == `n_leisure' & `n_wagep' == `n_hwage' & (`n_wagep' == `n_hecksig' | "`heckman'" != "")) {
            tempvar preds
            qui egen `preds' = rowtotal(`wagepred') if `touse'
            qui count if `preds' >= 1 & (`touse')
            local wagep = (r(N) > 0)
        }
        // Settings incorrect
        else {
            di in r "number of wage prediction variables does not match the number of leisure terms, hourly wage rates or mean squared errors"
            exit 498
        }
    }
    
    // No need to take random draws
    if (`wagep' == 0 & "`randvars'" == "") local draws = 1
    
    // Tax regression or tax benefit calculator needed
    local taxreg_rmse = 0
    if (`wagep' == 1 | "`heckman'" != "") {
        if ("`taxreg'" != "" & "`taxben'" == "") {
            tempname taxreg_from
            
            // Load estimates and estimation statistics
            qui est restore `taxreg'
            mat `taxreg_from' = e(b)
            local taxreg_rmse = e(rmse)
            
            // Fetch "constant" regressors (starting after wage variables and wage interactions)
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
    if ("`weight'" != "") local wgt "[`weight'=`exp']"
    
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
            qui gen `c' = log(cond("`ufunc'" == "boxcox" & `boxcc' > 0, `consumption' / `boxcc', `consumption')) if `touse'
            foreach var of local leisure {
                if (strpos("`leisure'", "`var'") == 1) local lei l1
                else                                   local lei l2
                qui gen ``lei'' = log(cond("`ufunc'" == "boxcox" & `boxcl' > 0, `var' / `boxcl', `var')) if `touse'
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
                    if ("`var'" != "0") local initrhs `initrhs' c.``ia''#c.``ia''#c.`var'
                    else                local initrhs `initrhs' c.``ia''#c.``ia''
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
        mat `init_from' = J(1, 2 + 2 * `n_leisure' + `n_cxias' + `n_lx1ias' + `n_lx2ias' + `n_indeps', 0)

        // Get init values from simple conditional logit estimation
        `qui' clogit `varlist' `initrhs' if `touse' `wgt', group(`group') iterate(25)
        if (e(converged) == 1) {
            // Save results
            mat `init_from' = e(b)
            local nobs      = e(N)
            local k         = e(k)
            local ll        = e(ll)
            //local ll_0      = e(ll_0) // calculate null model on our own
            
            // Update sample
            qui replace `touse' = e(sample)
        }
        
        // Wage equation
        if ("`heckman'" != "") {
            tempname init_wage
            foreach w of local hwage {
                tempvar ln`w'
                qui gen `ln`w'' = ln(`w') if `touse' & `varlist'

                // Run Heckman selection model
                if ("`select'" != "") {
                    // Fetch and store estimates
                    `qui' heckman `ln`w'' `heckman' if `touse' & `varlist' `wgt', select(`select')
                    mat `init_wage' = e(b)

                    // We estimate rho and sigma directly...
                    mat `init_wage'[1,colsof(`init_wage')-1] = e(rho)
                    mat `init_wage'[1,colsof(`init_wage')]   = e(sigma)
                }
                // Run simple OLS
                else {
                    // Fetch and store estimates
                    `qui' reg `ln`w'' `heckman' if `touse' & `varlist' `wgt'
                    mat `init_wage' = e(b)

                    // Save root mean squared error if wage variance is to be estimated
                    if ("`wagecorr'" != "") mat `init_wage' = `init_wage', e(rmse)
                }

                // If correlation between wages and preferences enabled
                if ("`wagecorr'" != "") mat `init_wage' = `init_wage', J(1, `corrBW', 0.0001)
            }
        }
        
        // Save init options
        local initopt init(`init_from', copy) obs(`nobs') lf0(`k' `ll')
    }
    else {
        qui count if `touse'
        local nobs = r(N)
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
    
    mata: lsl_round  = ("`round'" != "noround")                             // To round, or not to round?
    mata: lsl_ufunc  = st_local("ufunc")                                    // Utility function
    mata: lsl_Weight = ("`exp'" != "" ? st_data(., "`exp'") ///
                                      : J(`nobs', 1, 1))                    // Weight
    mata: lsl_Y      = st_data(., "`varlist'")                              // Left hand side
    mata: lsl_Hwage  = (`n_hwage' > 0 ? st_data(., tokens("`hwage'"))   ///
                                      : J(`nobs', `n_leisure', 0))          // Hourly wage rates
    mata: lsl_boxcc  = (lsl_ufunc == "boxcox" & `boxcc' > 0 ? `boxcc' : 1)   // Normalizing constant for Box-Cox consumption
    mata: lsl_boxcl  = (lsl_ufunc == "boxcox" & `boxcl' > 0 ? `boxcl' : 1)   // Normalizing constant for Box-Cox leisure
    
    //
    // Right hand side
    //
    
    // Consumption, squared and interactions
    mata: lsl_C   = st_data(., "`consumption'") :/ lsl_boxcc
    mata: lsl_CX  = ((`n_cxias'  > 0 ? st_data(., tokens("`cx'"))   ///
                                     : J(`nobs', 0, 0)), J(`nobs', 1, 1))
    mata: lsl_C2X = ((`n_c2xias' > 0 ? st_data(., tokens("`c2x'"))  ///
                                     : J(`nobs', 0, 0)), J(`nobs', (lsl_ufunc != "boxcox" ? 1 : 0), 1))
    
    // Leisure, squared and interactions
    forval i = 1/2 {
        local var : word `i' of `leisure'
        if ("`var'" != "") {
            mata: lsl_L`i'   = st_data(., "`var'") :/ lsl_boxcl
            mata: lsl_LX`i'  = ((`n_lx`i'ias'  >  0 ? st_data(., tokens("`lx`i''"))     ///
                                                    : J(`nobs', 0, 0)), J(`nobs', 1, 1))
            mata: lsl_L2X`i' = ((`n_l2x`i'ias' >  0 ? st_data(., tokens("`l2x`i''"))    ///
                                                    : J(`nobs', 0, 0)), J(`nobs', (lsl_ufunc != "boxcox" ? 1 : 0), 1))
        }
        else {
            mata: lsl_L`i'   = J(`nobs', 0, 0)
            mata: lsl_LX`i'  = J(`nobs', 0, 0)
            mata: lsl_L2X`i' = J(`nobs', 0, 0)
        }
    }
    
    // Dummy variables
    mata: lsl_Xind = (`n_indeps' > 0 ? st_data(., tokens("`indeps'")) : J(`nobs', 0, 0))
    
    // Build right hand side
    if ("`ufunc'" == "tran") {          // Translog utility
        mata: lsl_X = (log(lsl_C)  :* (lsl_CX,  log(lsl_C)  :* lsl_C2X,   log(lsl_L1),   log(lsl_L2)),         ///
                       log(lsl_L1) :* (lsl_LX1, log(lsl_L1) :* lsl_L2X1),                                    ///
                       log(lsl_L2) :* (lsl_LX2, log(lsl_L2) :* lsl_L2X2), log(lsl_L1) :* log(lsl_L2), lsl_Xind)
    }
    else if ("`ufunc'" == "quad") {     // Quadratic utility
        mata: lsl_X = (lsl_C  :* (lsl_CX,  lsl_C  :* lsl_C2X,   lsl_L1,   lsl_L2),         ///
                       lsl_L1 :* (lsl_LX1, lsl_L1 :* lsl_L2X1),                          ///
                       lsl_L2 :* (lsl_LX2, lsl_L2 :* lsl_L2X2), lsl_L1 :* lsl_L2, lsl_Xind)
    }
    else mata: lsl_X = J(`nobs', 0, 0)
    
    
    //
    // Joint wage estimation
    //
    
    mata: lsl_joint      = ("`joint'" != "")                                                // ML joint estimation?
    mata: lsl_Wobs       = (lsl_Y :* (log(lsl_Hwage) :< .))                                 // Wage observed?
    mata: lsl_HeckmVars  = (lsl_joint == 1 ? (st_data(., tokens("`heckman'")), J(`nobs', 1, 1))    ///
                                                        : J(`nobs', 0, 0))                  // Wage variables
    mata: lsl_SelectVars = (lsl_joint == 1 & "`select'" != "" ? (st_data(., tokens("`select'")), J(`nobs', 1, 1))    ///
                                                                             : J(`nobs', 0, 0))  // Wage variables
    mata: lsl_Days       = ("`days'" != "" ? st_data(., "`days'") : J(`nobs', 1, 365))      // Days of taxyear
    mata: lsl_Hours      = `totaltime' :- (lsl_L1, lsl_L2) :* lsl_boxcl                     // Hypothetical hours (caution: Box-Cox normalization!)
    mata: lsl_wagecorr   = `corrBW'                                                         // Correlate wage rates and preferences?
    
    
    //
    // Wage Prediction Stuff
    //
    
    mata: lsl_wagep = `wagep'                                                   // Run Wage Prediction
    mata: lsl_Wpred = (lsl_wagep == 1 ? st_data(., "`wagepred'")                ///
                                      : J(`nobs', 1 + cols(lsl_L2), 0))         // Dummies enabling or disabling the wage prediction
    mata: lsl_Sigma = ("`hecksigma'" != "" ? strtoreal(tokens("`hecksigma'"))   ///
                                           : J(1, `n_hecksig', 0))              // Estimated variance of Heckman correction
    
    //
    // Tax regression
    //
    
    if ("`taxreg'" != "") {
        mata: lsl_TaxregB     = st_matrix("`taxreg_from'")                                          // Tax regression estimates
        mata: lsl_taxreg_rmse = `taxreg_rmse'                                                       // Root Mean Squared error of tax regression
        mata: lsl_TaxregIas1  = ("`tria1'" != "" ? st_data(., tokens("`tria1'")) : J(`nobs', 0, 0)) // Interaction variables on Mwage1 and Mwage1^2
        mata: lsl_TaxregIas2  = ("`tria2'" != "" ? st_data(., tokens("`tria2'")) : J(`nobs', 0, 0)) // Interaction variables on Mwage2 and Mwage2^2
        mata: lsl_TaxregVars  = st_data(., tokens("`taxreg_vars'"))                                 // Variables that are independent of m_wage
    }
    
    
    //
    // Group level stuff
    //
    
    mata: lsl_J = panelsetup(st_data(., "`group'"), 1)
    mata: lsl_J = (lsl_J, 1 :+ lsl_J[.,2] :- lsl_J[.,1])
    mata: lsl_groups = rows(lsl_J)

    // Calculate log-likelihood of null model
    mata: st_local("ll_0", strofreal(-colsum(log(lsl_J[.,3]))))
    
    
    //
    // Dude force (lagrange constraints)
    //
    
    mata: lsl_lambda = `lambda'                                                         // Lagrange multiplier
    mata: lsl_force  = ("`force'" != "" ? st_data(., "`force'") : J(`nobs', 1, 0))      // Which observations to include in constraint?
    
    
    //
    // Random draws
    //
    
    mata: lsl_draws = strtoreal(st_local("draws"))                                                                      // Number of draws
    mata: lsl_burn  = strtoreal(st_local("burn"))                                                                       // Number of draws to burn
    mata: st_local("randvars", invtokens(strofreal(sort(strtoreal(tokens("`randvars'"))', 1))'))                        // Sort random coefficients
    mata: lsl_Rvars = ("`randvars'" != "" ? strtoreal(tokens("`randvars'"))' : J(0, 1, 0))                              // Random coefficients
    mata: lsl_corr  = ("`corr'" != "")                                                                                  // Random coefficients correlated?
    mata: lsl_R     = (`rvars' > 0 ? invnormal(halton(lsl_groups*lsl_draws, `rvars', 1+lsl_burn)) : J(`nobs', 0, 0))    // Halton sequences
    
    // Restore data
    restore
    
    
    /* RUN ESTIMATION
     */
    
    if ("`verbose'" != "") di as text "Run estimation..."
    
    // Set up equations
    local eq_consum (Cx: `varlist' = `cx')                              // Consumption
    if ("`ufunc'" != "boxcox") local eq_consum `eq_consum' (CxC: `c2x') // Consumption^2
    local eq_leisure
    foreach var of local leisure {
        local i = 1 + (strpos("`leisure'", "`var'") > 1)
        local eq_leisure `eq_leisure' (L`i'x: `lx`i'')                                  // Leisure
        if ("`ufunc'" != "boxcox") local eq_leisure `eq_leisure' (L`i'xL`i': `l2x`i'')  // Leisure^2
        local eq_consum  `eq_consum' /CxL`i'                                            // Consumption X leisure interaction
    }
    if (`n_leisure' == 2) local eq_leisure `eq_leisure' /L1xL2          // Leisure term interaction
    if (`n_indeps'  >  0) local eq_indeps  (IND: `indeps', noconst)     // Independent variables / dummies
    
    // Joint wage estimation?
    if ("`heckman'" != "") {
        // Wage equation
        local eq_wages (lnW1: `heckman')
        
        // Selection equation?
        if ("`select'" != "") {
            local eq_wages `eq_wages' (S1: `select') /rho               // Selection equation and rho needed when Heckman model
            if ("`wagecorr'" == "") local eq_wages `eq_wages' /sigma    // No correlation, estimate variance of wages
        }
        
        // Correlation with preferences?
        if ("`wagecorr'" != "") {
            local eq_wages `eq_wages' /l_W1_W1      // Cholesky for Var(W1)
            if (`wagecorrC'  == 1) local eq_wages `eq_wages' /l_C_W1    // Cholesky for Cov(W1, C)
            if (`wagecorrL1' == 1) local eq_wages `eq_wages' /l_L1_W1   // Cholesky for Cov(W1, L1)
            if (`wagecorrL2' == 1) local eq_wages `eq_wages' /l_L2_W2   // Cholesky for Cov(W2, L2) -- BUGGY!!
        }

        // Add wage estimates to initial values
        if ("`initopt'" != "" & "`init_from'" != "") mat `init_from' = (`init_from', `init_wage')
    }
    
    // Box-Cox utility function?
    if ("`ufunc'" == "boxcox") {
        local eq_boxcox /l_C /l_L1
        if (`n_leisure' == 2) local eq_boxcox `eq_boxcox' /l_L2

        // Add initial values (assume Stone-Geary model)
        if ("`initopt'" != "" & "`init_from'" != "") mat `init_from' = (`init_from', 0.0, J(1, `n_leisure', 0.0))
    }
    
    // Random coefficients?
    if (`n_randvars' > 0) {
        local eq_rands
        if ("`corr'" == "" & "`wagecorr'" == "") {
            forval i = 1/`n_randvars' {
                local coef : word `i' of `randvars'
                local eq_rands `eq_rands' /sd_`coef'
            }
            if ("`initopt'" != "" & "`init_from'" != "") mat `init_from' = (`init_from', J(1, `n_randvars', 0.0001))
        }
        else {
            forval i = 1/`n_randvars' {
                local a : word `i' of `randvars'
                forval k = `i'/`n_randvars' {
                    local b : word `k' of `randvars'
                    local eq_rands `eq_rands' /l_`a'_`b'
                }
            }
            if ("`initopt'" != "" & "`init_from'" != "") mat `init_from' = (`init_from', J(1, `n_randvars' * (`n_randvars' + 1) / 2, 0.0002))
        }
    }

    
    //
    // Estimate
    //
    
    scalar lsl_dudes = .
    
    if ("`debug'" != "") di as text "ml `method'`debug' lslogit_d2() `eq_consum' `eq_leisure' `eq_indeps' `eq_wages' `eq_boxcox' `eq_rands'"
    ml model `method'`debug' lslogit_d2() `eq_consum' `eq_leisure' `eq_indeps' `eq_wages' `eq_boxcox' `eq_rands' ///
            if `touse' `wgt', group(`group') `initopt' search(off) iterate(`iterate') nopreserve max `difficult' `trace' `gradient' `hessian' technique(`technique')
    
    //
    // Save results
    //
    
    // Save model setup
    ereturn local title   "Mixed Logit Labor Supply Model"
    ereturn local cmd     "lslogit"
    ereturn local predict "lslogit_p"
    ereturn local ufunc    `ufunc'
    ereturn local wagep    `wagep'
    ereturn local joint    `joint'
    ereturn local wagecorr = ("`wagecorr'" != "")
    ereturn local draws    `draws'
    ereturn local burn     `burn'
    ereturn local totaltime `totaltime'
    ereturn local round  = ("`round'" != "noround")
    ereturn local taxreg = ("`taxreg'" != "")
    ereturn local group    `group'
    ereturn local groups   `n_groups'
    ereturn local depvar   `varlist'
    ereturn local consum   `consumption'
    ereturn local leisure  `leisure'
    
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
    if ("`ufunc'" == "boxcox") {    // Box-Cox
        if ("`boxcc'" != "") ereturn local boxcc   `boxcc'
        if ("`boxcl'" != "") ereturn local boxcl   `boxcl'
    }
    if (`n_randvars' > 0) {         // Random coefficients
        ereturn local randvars `randvars'
        ereturn local corr = cond("`corr'" != "", 1, 0)
    }
    if ("`taxreg'" != "") {         // Tax regression
        ereturn matrix taxreg_b = `taxreg_from', copy
        if ("`taxreg_vars'" != "") ereturn local taxreg_vars `taxreg_vars'
        if ("`tria1'" != "") ereturn local taxreg_ia1 `tria1'
        if ("`tria2'" != "") ereturn local taxreg_ia2 `tria2'
    }
    if ("`heckman'" != "") ereturn local heckman `heckman'
    if ("`select'"  != "") ereturn local select  `select'
    
    // Display some coefficients as auxiliary
    ereturn scalar k_aux = ("`ufunc'" == "boxcox") * (1 + `n_leisure') + ("`select'" != "") +   ///
                           `n_randvars' * cond("`corr'" != "", (`n_randvars' + 1) / 2, 1)   +   ///
                           ("`wagecorr'" != "" | "`select'" != "") + `corrBW'
    
    // Pseudo R2 (may be misleading as it refers to the null-model, LR and p value refer to init values...)
    if ("`ll_0'" != "") ereturn scalar r2_p = 1 - e(ll)/`ll_0'
    
    // Additional output
    if (lsl_dudes != .) ereturn scalar dudes = lsl_dudes
    /*
    foreach aux in sigma_w1 sigma_w2 dudes {
        di as error "`aux' = " `aux'
        if (`aux' != "") ereturn scalar `aux' = ``aux''
    }
    */
    
    //
    // Clean up
    //
    
    foreach m in round ufunc Weight Y Hwage boxcc boxcl C CX C2X L1 LX1 L2X1 L2 LX2 L2X2 Xind X joint ///
                 HeckmVars SelectVars Days Hours wagep Wpred Sigma TaxregB TaxregIas1 TaxregIas2 TaxregVars ///
                 groups J draws burn Rvars corr R {
        cap mata mata drop lsl_`m'
    }
    
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
    
    external string scalar  lsl_ufunc           // Functional form
    external real scalar    lsl_groups          // Number of groups
    external real colvector lsl_Y               // Left hand side variable
    external real matrix    lsl_J               // Panel setup data
    external real matrix    lsl_X               // Right hand side variables
    external real colvector lsl_Weight          // Group weights
    
    external real scalar    lsl_draws           // Number of random draws
    external real matrix    lsl_R               //   Halton sequences
    
    external real colvector lsl_Rvars           // Random coefficients
    external real scalar    lsl_corr            //   Enable correlation?
    
    external real scalar    lsl_joint           // Joint ML estimation?
    external real matrix    lsl_Wobs            //   Wage rate observed?
    external real matrix    lsl_HeckmVars       //   Right hand side variables
    external real matrix    lsl_SelectVars      //   Selection variables
    external real scalar    lsl_wagecorr        //   Number of covariances between wages and leisure preferences?
    
    external real scalar    lsl_wagep           // Wage Prediction Error?
    external real matrix    lsl_Wpred           //   Prediction dummies
    external real colvector lsl_Days            //   Number of days per tax year
    external real matrix    lsl_Hwage           //   Hourly wage rates
    external real rowvector lsl_Sigma           //   Variance of the wage regression
    external real matrix    lsl_Hours           //   Hours of work
    
    external real rowvector lsl_TaxregB         // Tax Regression
    external real scalar    lsl_taxreg_rmse     //   Root Mean Squared Error
    external real matrix    lsl_TaxregVars      //   Wage independent variables of tax regression
    external real matrix    lsl_TaxregIas1      //   Wage interaction variables of tax regression
    external real matrix    lsl_TaxregIas2      //   Wage interaction variables of tax regression
    
    external real scalar    lsl_round           // To round, or not to round.
    external real scalar    lsl_boxcc           // Normalizing constant for Box-Cox consumption
    
    external real scalar    lsl_lambda          // Lagrange multiplier for dude correction
    external real colvector lsl_force           // Which observations to include in constraint?
    
    external real colvector lsl_C
    external real matrix    lsl_CX
    external real matrix    lsl_C2X
    external real colvector lsl_L1
    external real matrix    lsl_LX1
    external real matrix    lsl_L2X1
    external real matrix    lsl_L2       // Buggy, L2 is colvector in fact, but colvector needs one column...
    external real matrix    lsl_LX2
    external real matrix    lsl_L2X2
    external real matrix    lsl_Xind
    
    
    //
    // Declarations
    //
    
    real scalar    n, i, c, e, nobs, nlei, ncons, wagep
    real scalar    b, bfix, bwage, bsel, bheck, blam, brnd
    real rowvector Bfix, Bwage, Bsel, Brho, Bsig, Brnd, Zeta, Beta
    
    real scalar    nRV, rvars, iRV, r, iC, iL1, iL2,
                   iwage, isel, iheck, ilam, irnd, ilC, ilL1, ilL2
    real matrix    CholB, CholBW, CholW
    
    real matrix    Hwage, LnWres, LnWresPur, Select, Lambda, Wn
    real matrix    Mwage, TaxregX1, TaxregX2, TaxregX, DCdM, D2CdM2
    
    real matrix    DUdx, DUdB, DUdlam, DUdBr, DWdBw, DWdBs, DWdBrho, Dude,
                   DWdBsig, DXdH, YmPn_D2UdB2, YmPn_D2UdBr2, YmPn_D2UdBdBr, YmPn_D2Udx2, DUdlC, DUdlL1, DUdlL2
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
    
    // Definitions
    i     = 1                   // Indicates first observation of active group
    rvars = rows(lsl_Rvars)     // Number of random variables
    nRV   = rvars + 1           // Indicates first random variable to use (column of lsl_R) for wage pred
    nobs  = rows(lsl_Y)         // Number of observations
    nlei  = 1 + cols(lsl_L2)    // Number of leisure terms
    ncx   = cols(lsl_CX)
    nc2x  = cols(lsl_C2X)
    nlx1  = cols(lsl_LX1)
    nl2x1 = cols(lsl_L2X1)
    nlx2  = cols(lsl_LX2)
    nl2x2 = cols(lsl_L2X2)
    nxind = cols(lsl_Xind)
    ncons = ncx + nc2x + nlei   // Number of variables including consumption
    Dude  = J(nobs, lsl_draws, 0)
    
    // Number of coefficients
    b     = cols(B)                                                     // Total number
    bwage = cols(lsl_HeckmVars)                                         // Number of wage regression coefficients
    bsel  = cols(lsl_SelectVars)                                        // Number of selection equation coefficients
    bheck = (bsel > 0) + (bsel > 0 | lsl_wagecorr) + lsl_wagecorr       // Number of additional wage related coefficients
    brnd  = (lsl_corr  == 1 ? rvars * (rvars + 1) / 2 : rvars)          // Number of variance and covariance terms for random coefficients
    blam  = (lsl_ufunc == "boxcox" ? 1 + nlei : 0)                      // Number of Box-Cox transformation coefficients
    bfix  = ncons + nlx1 + nl2x1 + nlx2 + nl2x2 + (nlei == 2) + nxind   // Number of fixed preference coefficients
    
    // Maximum Likelihood Parameter
    lnf = 0             // Log-likelihood
    G   = J(1, b, 0)    // Gradient
    H   = J(b, b, 0)    // Hessian matrix
    
    
    //
    // Build coefficient vector
    //
    
    // Get coefficient indices
    iwage = 1     + bfix
    isel  = iwage + bwage
    iheck = isel  + bsel
    ilam  = iheck + bheck
    irnd  = ilam  + blam
    
    iC  = ncx
    iL1 = ncons + nlx1
    iL2 = iL1 + nl2x1 + nlx2
    
    // Split up vector
    Bfix  = B[|1\bfix|]                                                 // Get fixed coefficients
    Bwage = (bwage > 0 ? B[|iwage\iwage + bwage - 1|] : J(0, 0, 0))     // Wage coefficients
    Bsel  = (bsel  > 0 ? B[|isel\isel + bsel - 1|]    : J(0, 0, 0))     //   Selection coefficients
    Brho  = (bsel  > 0 ? B[iheck]                     : J(0, 0, 0))     //   Heckman rho
    Bsig  = (bheck > 0 ? B[iheck + (bsel > 0)]        : J(0, 0, 0))     //   Cholesky factor of wage variance
    Brnd  = (rvars > 0 ? B[|irnd\irnd + brnd - 1|]    : J(0, 0, 0))     // Get auxiliary random coefficients
    
    // Build Cholesky matrix
    CholB  = (lsl_corr ? lowertriangle(invvech(Brnd')) : diag(Brnd))    // Cholesky factors of coefficients vector
    CholBW = J(nlei, brnd, 0)                                           // Cholesky factors between coefficients and wages
    CholW  = diag(Bsig)                                                 // Cholesky factors of wage variances
    
    // DEBUG
    if (nlei == 1) {
        if      (lsl_wagecorr == 1) CholBW = B[iheck + (bsel > 0) + 1] :* ((lsl_Rvars :== iC) :+ (lsl_Rvars :== iL1))'
        else if (lsl_wagecorr == 2) CholBW = B[iheck + (bsel > 0) + 1] :* (lsl_Rvars :== iC)' :+
                                             B[iheck + (bsel > 0) + 2] :* (lsl_Rvars :== iL1)'
    }
    else "DEBUG!"
    // DEBUG
    
    /*
    B
    "iC=" + strofreal(iC) + ", iL1=" + strofreal(iL1) + ", iL2=" + strofreal(iL2)
    "iwage=" + strofreal(iwage) + ", isel=" + strofreal(isel) + ", iheck=" + strofreal(iheck) + ", ilam=" + strofreal(ilam) + ", irnd=" + strofreal(irnd)
    "bwage=" + strofreal(bwage) + ", bsel=" + strofreal(bsel) + ", bheck=" + strofreal(bheck) + ", blam=" + strofreal(blam) + ", brnd=" + strofreal(brnd)
    Brnd
    Bsig
    CholB
    CholBW
    */
        
    // Calculate wage variance (by Cholesky or just the coefficient)
    SigmaW = (lsl_wagecorr ? sqrt(rowsum((CholBW, CholW):^2))' : Bsig)
    
    // Build matrix with random coefficients (mean zero), every row is a draw
    if (brnd > 0) {
        Zeta = J(rows(lsl_R), bfix, 0)
        Zeta[.,lsl_Rvars] = cross(lsl_R[|1,1\.,rvars|]', CholB')
    }
    // From now on: Beta[rows=lsl_R,cols=Bfix] = Bfix :+ Zeta
    
    
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
        // Selection equation
        if (bsel > 0) {
            Select = cross(lsl_SelectVars', Bsel')          // Selection prediction
            Lambda = normalden(Select) :/ normal(Select)    //   Heckman lambda
            
            // Predict log-wages with sample correction
            Hwage  = cross((lsl_HeckmVars, Lambda)', (Bwage, (rows(Brho :* SigmaW) > 0 ? Brho :* SigmaW : J(1, 0, 0)))')
        }
        // Predict log-wages without sample correction
        else Hwage = cross(lsl_HeckmVars', Bwage')
        
        // Wage observed?
        //Wobs = lsl_Y :* (log(lsl_Hwage) :< .)
        
        // Redisuals
        LnWres    = lsl_Wobs :* (log(lsl_Hwage) :- Hwage)                           // Including Heckman correction
        LnWresPur = lsl_Wobs :* (log(lsl_Hwage) :- cross(lsl_HeckmVars', Bwage'))   // Just wage equation residuals
        
        // Neither Heckman nor correlation, calculate wage variance
        if (!bsel & !lsl_wagecorr) {
            Bsig   = sqrt(cross(LnWresPur, LnWresPur) / (colsum(lsl_Wobs) - bwage))     // Root MSE
            SigmaW = Bsig                                                           // Overwrite wage variance
            CholW  = diag(Bsig)                                                     // Overwrite Cholesky factors
            for (c = 1; c <= cols(Bsig); c++) {
                st_numscalar("r(sigma_w" + strofreal(c) + ")", Bsig[1,c])           // Save RMSE
            }
        }
        
        // Predict wages
        Hwage = exp(Hwage :+ SigmaW^2/2)
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

        // Run by random draw
        for (r = 1; r <= lsl_draws; r++) {
            // Init
            iRV  = lsl_draws * (n - 1) + r       // Indicates the active Halton sequence
            
            // Build (random?) coefficients matrix
            Beta = Bfix :+ (brnd > 0 ? Zeta[iRV,.] : 0)
            

            /* Integrate out wage prediction error */
            
            if (wagep | lsl_joint) {
                //
                // Calculate monthly earnings
                //

                // Adjust wages with random draws if prediction enabled
                if (lsl_wagep) Wn = Hwage[|i,1\e,.|] :* exp(cross((CholBW, CholW)', lsl_R[|iRV,1\iRV,cols(lsl_R) - 1|]')' :* lsl_Wpred[|i,1\e,.|])
                
                // Calculate monthly earnings
                Mwage = (lsl_Days[|i\e|] :/ 12 :/ 7) :* lsl_Hours[|i,1\e,.|] :* Wn

                // Round monthly earnings if enabled
                if (lsl_round) Mwage = round(Mwage, 0.01)

                
                //
                // Set up tax regression covariates and predict dpi
                //
                
                // Fill matrix of independent variables for dpi prediction
                TaxregX1 = (Mwage[.,1], Mwage[.,1]:^2, lsl_TaxregIas1[|i,1\e,.|] :* Mwage[.,1],
                                                       lsl_TaxregIas1[|i,1\e,.|] :* Mwage[.,1]:^2)
                if (nlei == 2) TaxregX2 = (Mwage[.,2], Mwage[.,2]:^2, lsl_TaxregIas2[|i,1\e,.|] :* Mwage[.,2],
                                                                      lsl_TaxregIas2[|i,1\e,.|] :* Mwage[.,2]:^2)
                else           TaxregX2 = J(c, 0, 0)
                TaxregX = (TaxregX1, TaxregX2, lsl_TaxregVars[|i,1\e,.|], J(c, 1, 1))
                
                // Predict disposable income
                //   - can't be negative!
                //   - add random draw
                //   - normalize for Box-Cox specification
                C = rowmax((cross(TaxregX', lsl_TaxregB') :+ (lsl_taxreg_rmse :* lsl_R[iRV,cols(lsl_R)]), J(c, 1, 1))) :/ lsl_boxcc
                
                // Build matrix with independent variables
                if      (lsl_ufunc == "tran") Xnr = (log(C)  :* (CX,  log(C)  :* C2X,   log(L1),   log(L2)),
                                                     log(L1) :* (LX1, log(L1) :* L2X1),
                                                     log(L2) :* (LX2, log(L2) :* L2X2), log(L1) :* log(L2), Xind)
                else if (lsl_ufunc == "quad") Xnr = (C  :* (CX,  C  :* C2X,   L1,   L2),
                                                     L1 :* (LX1, L1 :* L2X1),
                                                     L2 :* (LX2, L2 :* L2X2), L1 :* L2, Xind)
                else if (lsl_ufunc == "boxcox") {
                    BcC = lsl_boxcox(C, lC)
                    Xnr = ((CX, BcL1, BcL2) :* BcC, LX1 :* BcL1, LX2 :* BcL2, BcL1 :* BcL2, Xind)
                }
            }
            

            /* Calculate utility levels */

            // Calculate choice probabilities
            Unr = cross(Xnr', Beta')                                    // Utility (choices in rows, draws in columns)
            Enr = exp(Unr :+ colmin(-mean(Unr) \ 700 :- colmax(Unr)))   // Standardize to avoid missings
            Pnr = Enr :/ colsum(Enr)                                    // Probabilities
            
            // Simplify
            pni  = cross(Yn, Pnr)   // Probability that choice is chosen
            YmPn = Yn :- Pnr        // Choice minus probabilities
            
            // Calculate first derivative to consumption (dudes)
            if (todo == 0 | lsl_lambda > 0) {
                if    (lsl_ufunc == "boxcox") Dude[|i,r\e,r|] = cross((CX, lsl_boxcox(L1, lL1), lsl_boxcox(L2, lL2))', Beta[|1\ncons|]') :*
                                                                (reldif(lC, 0) >= 1e-25 ? C:^(lC - 1) : (1 :/ C))
                else if (lsl_ufunc == "quad") Dude[|i,r\e,r|] = cross((CX, 2 :* C2X :* C, L1, L2)', Beta[|1\ncons|]')
                else if (lsl_ufunc == "tran") Dude[|i,r\e,r|] = cross((CX, 2 :* C2X :* log(C), log(L1), log(L2))', Beta[|1\ncons|]') :/ C
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
                    BcCx  = (CX, BcL1, (nlei == 2 ? BcL2 : J(c, 0, 0)), J(c, bfix - ncons, 0))
                    BcL1x = (J(c, ncx, 0), BcC, J(c, (nlei == 2), 0), LX1, J(c, nlx2, 0), (nlei == 2 ? BcL2 : J(c, 0, 0)), J(c, nxind, 0))
                    BcL2x = (nlei == 2 ? (J(c, ncx + 1, 0), BcC, J(c, nlx1, 0), LX2, BcL1, J(c, nxind, 0)) : J(c, 0, 0))
                    
                    // Calculate gradients
                    DUdlC  = lsl_boxcox_g(C,  lC)  :* cross(BcCx', Beta')
                    DUdlL1 = lsl_boxcox_g(L1, lL1) :* cross(BcL1x', Beta')
                    DUdlL2 = (nlei == 2 ? lsl_boxcox_g(L2, lL2) :* cross(BcL2x', Beta') : J(c, 0, 0))
                    DUdlam = (DUdlC, DUdlL1, DUdlL2)
                } else DUdlam = J(c, 0, 0)
                
                // Full Maximum Likelihood estimation?
                real matrix DUdwage, DWdBwcorr
                if (lsl_joint) {
                    if      (lsl_ufunc == "quad") DUdC = cross((CX, 2 :* C2X :* C, L1, L2)', Beta[|1\ncons|]')
                    else if (lsl_ufunc == "tran") DUdC = cross(((CX, 2 :* C2X :* log(C), L1, L2) :/ C)', Beta[|1\ncons|]')
                    DCdM = cross((J(c, 1, 1), 2 :* Mwage, lsl_TaxregIas1[|i,1\e,.|], 2 :* Mwage :* lsl_TaxregIas1[|i,1\e,.|])', lsl_TaxregB[|1\2 + 2 * cols(lsl_TaxregIas1)|]')
                    DMdH = (lsl_Days[|i,1\e,1|] :/ 12 :/ 7) :* lsl_Hours[|i,1\e,.|]
                    if (bsel > 0) {
                        /*
                        DWdBw   =   Wn :* lsl_HeckmVars[|i,1\e,.|]
                        DWdBs   = - Wn :* Bsig :* Brho :* normalden(Select[|i,1\e,.|]) :* lsl_SelectVars[|i,1\e,.|] :*
                                    (normal(Select[|i,1\e,.|]) :* Select[|i,1\e,.|] :+ normalden(Select[|i,1\e,.|])) :/ normal(Select[|i,1\e,.|]):^2
                        DWdBrho =   Wn :*  Bsig :* Lambda[|i,1\e,.|]
                        DWdBsig =   Wn :* (Brho :* Lambda[|i,1\e,.|] :+ Bsig)
                        DWdBwcorr = J(c, 0, 0)
                        */
                    } else {
                        if (lsl_wagecorr) {
                            DWdBw     = Wn :* lsl_HeckmVars[|i,1\e,.|]
                            DWdBsig   = Wn :* (lsl_R[iRV,nRV] :+ Bsig)
                            if      (lsl_wagecorr == 1) DWdBwcorr = Wn :* (rowsum(lsl_R[|iRV,1\iRV,rvars|] :* ((lsl_Rvars :== iC) :+ (lsl_Rvars :== iL1))') :+ B[iheck + (bsel > 0) + 1])
                            else if (lsl_wagecorr == 2) DWdBwcorr = cross(Wn', ((rowsum(lsl_R[|iRV,1\iRV,rvars|] :* (lsl_Rvars :== iC )') :+ B[iheck + (bsel > 0) + 1]),
                                                                                (rowsum(lsl_R[|iRV,1\iRV,rvars|] :* (lsl_Rvars :== iL1)') :+ B[iheck + (bsel > 0) + 2])))
                            DWdBs     = DWdBrho = J(c, 0, 0)
                        } else {
                            DWdBw =   Wn :* (lsl_HeckmVars[|i,1\e,.|] :- (Bsig + lsl_R[iRV,nRV]) :* colsum(lsl_HeckmVars :* LnWresPur) / (colsum(lsl_Wobs) - bwage))
                            DWdBs = DWdBrho = DWdBsig = DWdBwcorr = J(c, 0, 0)
                        }
                    }
                    DUdwage = DUdC :* DCdM :* DMdH :* (DWdBw, DWdBs, DWdBrho, DWdBsig, DWdBwcorr)
                } else DUdwage = J(c, 0, 0)
                
                // Random components
                if (brnd > 0) {
                    DUdBr = (lsl_corr == 1 ? cross(DUdB[.,vech(J(1, rvars, lsl_Rvars))]',
                                                   diag(vech(J(rvars, 1, lsl_R[|iRV,1\iRV,rvars|]))))
                                           : DUdB[.,lsl_Rvars] :* lsl_R[|iRV,1\iRV,rvars|])
                } else DUdBr = J(c, 0, 0)
                
                // Total
                DUdx = (DUdB, DUdwage, DUdlam, DUdBr)
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
                
                // Full Maximum Likelihood estimation?
                real matrix YmPn_D2UdBdwage, YmPn_D2Udwage2
                YmPn_D2UdBdwage = J(bwage + bsel + 2 * (bsel > 0), bfix, 0)
                YmPn_D2Udwage2  = J(bwage + bsel + 2 * (bsel > 0), bwage + bsel + 2 * (bsel > 0), 0)
                if (lsl_joint == 1) {
                    if (lsl_ufunc == "quad") YmPn_D2UdBdwage = cross(YmPn :* DCdM :* DMdH :* (DWdBw, DWdBs, DWdBrho, DWdBsig), (CX, 2 :* C2X :* C, L1, J(c, cols(Xnr) - ncons, 0)))
                    if (lsl_ufunc == "quad") D2UdC2 = 2 :* cross(C2X', Beta[|ncx + 1\ncx + nc2x|]')
                    D2CdM2 = cross((J(c, 1, 2), 2 :* lsl_TaxregIas1[|i,1\e,.|])', (lsl_TaxregB[2], lsl_TaxregB[|2 + cols(lsl_TaxregIas1) + 1\2 + 2 * cols(lsl_TaxregIas1)|])')
                    D2MdH2 = 0
                    // d2U / dBw2
                    YmPn_D2Udwage2[|1,1\bwage,bwage|] = cross(YmPn :* DWdBw :* (DMdH:^2 :* (D2UdC2 :* DCdM:^2 :+ DUdC :* D2CdM2) + DUdC :* DCdM :* D2MdH2), DWdBw) :+
                                                        cross(YmPn :* DWdBw :* DUdC :* DCdM :* DMdH, lsl_HeckmVars[|i,1\e,.|])
                    // d2U / dBs2
                    YmPn_D2Udwage2[|bwage+1,bwage+1\bwage+bsel,bwage+bsel|] = cross(YmPn :* DWdBs :* (DMdH:^2 :* (D2UdC2 :* DCdM:^2 :+ DUdC :* D2CdM2) + DUdC :* DCdM :* D2MdH2), DWdBs) :+
                                                                              cross(YmPn :* DUdC :* DCdM :* DMdH :* Brho :* Bsig :* normalden(Select[|i,1\e,.|]) :* lsl_SelectVars[|i,1\e,.|],
                                                                                    - ((DWdBs :- Wn :* Select[|i,1\e,.|] :* lsl_SelectVars[|i,1\e,.|]) :* (normal(Select[|i,1\e,.|]) :* Select[|i,1\e,.|] :+ normalden(Select[|i,1\e,.|])) :/ normal(Select[|i,1\e,.|]):^2
                                                                                       :+ Wn :* lsl_SelectVars[|i,1\e,.|] :* (normal(Select[|i,1\e,.|]):^2 :- 2 :* normalden(Select[|i,1\e,.|]):^2 :- 2 :* normal(Select[|i,1\e,.|]) :* normalden(Select[|i,1\e,.|]) :* Select[|i,1\e,.|]) :/ normal(Select[|i,1\e,.|]):^3))
                    // d2U / dBrho2
                    YmPn_D2Udwage2[bwage+bsel+1,bwage+bsel+1] = cross(YmPn :* DWdBrho :* (DMdH:^2 :* (D2UdC2 :* DCdM:^2 :+ DUdC :* D2CdM2) + DUdC :* DCdM :* D2MdH2), DWdBrho) :+
                                                                cross(YmPn :* DWdBrho :* DUdC :* DCdM :* DMdH, Bsig :* Lambda[|i,1\e,.|])
                    // d2U / dBsig2
                    YmPn_D2Udwage2[bwage+bsel+2,bwage+bsel+2] = cross(YmPn :* DWdBsig :* (DMdH:^2 :* (D2UdC2 :* DCdM:^2 :+ DUdC :* D2CdM2) + DUdC :* DCdM :* D2MdH2), DWdBsig) :+
                                                                cross(YmPn :* DUdC :* DCdM :* DMdH, DWdBsig :* (Brho :* Lambda[|i,1\e,.|] :+ Bsig) :+ Wn)
                    // d2U / dBrho dBsig
                    YmPn_D2Udwage2[bwage+bsel+2,bwage+bsel+1] = cross(YmPn :* DWdBrho :* (DMdH:^2 :* (D2UdC2 :* DCdM:^2 :+ DUdC :* D2CdM2) + DUdC :* DCdM :* D2MdH2), DWdBsig) :+
                                                                cross(YmPn :* DUdC :* DCdM :* DMdH, DWdBsig :* Bsig :* Lambda[|i,1\e,.|] :+ Wn :* Lambda[|i,1\e,.|])
                    // d2U / dBrho dBw
                    YmPn_D2Udwage2[|bwage+bsel+1,1\bwage+bsel+1,bwage|] = cross(YmPn :* DWdBrho :* (DMdH:^2 :* (D2UdC2 :* DCdM:^2 :+ DUdC :* D2CdM2) + DUdC :* DCdM :* D2MdH2), DWdBw) :+
                                                                          cross(YmPn :* DUdC :* DCdM :* DMdH, DWdBrho :* lsl_HeckmVars[|i,1\e,.|])
                    // d2U / dBsig dBw
                    YmPn_D2Udwage2[|bwage+bsel+2,1\bwage+bsel+2,bwage|] = cross(YmPn :* DWdBsig :* (DMdH:^2 :* (D2UdC2 :* DCdM:^2 :+ DUdC :* D2CdM2) + DUdC :* DCdM :* D2MdH2), DWdBw) :+
                                                                          cross(YmPn :* DUdC :* DCdM :* DMdH, DWdBsig :* lsl_HeckmVars[|i,1\e,.|])
                    // d2U / dBs dBw
                    YmPn_D2Udwage2[|bwage+1,1\bwage+bsel,bwage|] = cross(YmPn :* DWdBs :* (DMdH:^2 :* (D2UdC2 :* DCdM:^2 :+ DUdC :* D2CdM2) + DUdC :* DCdM :* D2MdH2), DWdBw) :+
                                                                   cross(YmPn :* DWdBs :* DUdC :* DCdM :* DMdH, lsl_HeckmVars[|i,1\e,.|])
                    // d2U / dBs dBrho
                    YmPn_D2Udwage2[|bwage+bsel+1,bwage+1\bwage+bsel+1,bwage+bsel|] = cross(YmPn :* DWdBrho :* (DMdH:^2 :* (D2UdC2 :* DCdM:^2 :+ DUdC :* D2CdM2) + DUdC :* DCdM :* D2MdH2), DWdBs) :+
                                                                                     cross(YmPn :* DUdC :* DCdM :* DMdH, DWdBs :* (1 :/ Brho :+ Bsig :* Lambda[|i,1\e,.|]))
                    // d2U / dBs dBsig
                    YmPn_D2Udwage2[|bwage+bsel+2,bwage+1\bwage+bsel+2,bwage+bsel|] = cross(YmPn :* DWdBsig :* (DMdH:^2 :* (D2UdC2 :* DCdM:^2 :+ DUdC :* D2CdM2) + DUdC :* DCdM :* D2MdH2), DWdBs) :+
                                                                                     cross(YmPn :* DUdC :* DCdM :* DMdH, DWdBs :* (1 :/ Bsig :+ Brho :* Lambda[|i,1\e,.|] :+ Bsig))
                    
                    // Partial second derivatives
                    YmPn_D2Udx2 = makesymmetric((YmPn_D2UdB2    , YmPn_D2UdBdwage' \
                                                 YmPn_D2UdBdwage, YmPn_D2Udwage2  ))
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
                    if    (lsl_ufunc == "boxcox") G[|1\ncons|] = G[|1\ncons|]
                    else if (lsl_ufunc == "quad") G[|1\ncons|] = G[|1\ncons|] :+ lsl_lambda :* colsum(lsl_Weight[i] :* lsl_force[|i\e|] :* (CX, 2 :* C2X :* C, L1, L2))
                    else if (lsl_ufunc == "tran") G[|1\ncons|] = G[|1\ncons|] :+ lsl_lambda :* colsum(lsl_Weight[i] :* lsl_force[|i\e|] :* (CX, 2 :* C2X :* log(C), log(L1), log(L2)) :/ C)
                }
                if (todo == 2 & lsl_ufunc == "boxcox") {

                }
            }
        }

        // Prevent likelihood from becoming exactly zero
        lsum = max((lsum, 1e-25))
        
        // Add to overall statistics
        lnf = lnf + lsl_Weight[i] * log(lsum / lsl_draws)
        if (todo >= 1) G = G + lsl_Weight[i] * (lsum > 1e-25 ? Gsum / lsum : J(1, b, 0))
        if (todo == 2) H = H + lsl_Weight[i] * (lsum > 1e-25 ? H2sum / lsum - cross(Gsum, H1sum) / lsum^2: J(b, b, 0))
    }
    
    // Add likelihood of wage equation?
    if (lsl_joint) {
        lnf = lnf + cross(lsl_Wobs, log(normalden(LnWresPur :/ SigmaW)) :- log(SigmaW))
        /*if (bsel > 0) {
            lnf = lnf + cross(lsl_Wobs, log(normal((Select :+ LnWresPur :* Brho :/ SigmaW) :/ sqrt(1 - Brho^2)))) +
                        cross((lsl_Wobs :== 0) :* lsl_Y, log(normal(-Select)))
            if (todo >= 1) {
                G[|iwage\iwage + bwage + bsel - 1|] = G[|iwage\iwage + bwage + bsel - 1|] :+ (cross(lsl_Wobs :* lsl_HeckmVars, LnWresPur :/ SigmaW^2 :- (normalden((Select :+ LnWresPur :* Brho :/ SigmaW) :/ sqrt(1 - Brho^2)) :/ normal((Select :+ LnWresPur :* Brho :/ SigmaW) :/ sqrt(1 - Brho^2))) :* (Brho / SigmaW) :/ sqrt(1 - Brho^2))',
                                                                                              cross(lsl_SelectVars, lsl_Wobs :* (normalden((Select :+ LnWresPur :* Brho :/ SigmaW) :/ sqrt(1 - Brho^2)) :/ normal((Select :+ LnWresPur :* Brho :/ SigmaW) :/ sqrt(1 - Brho^2))) :/ sqrt(1 - Brho^2) :- (lsl_Wobs :== 0) :* lsl_Y :* normalden(-Select) :/ normal(-Select))')
            }
        } else {*/
            if (todo >= 1) {
                G[|iwage\iwage + bwage - 1|] = G[|iwage\iwage + bwage - 1|] :+ cross(lsl_Wobs, lsl_HeckmVars :* LnWresPur :/ SigmaW:^2)
                if (lsl_wagecorr) G[|iheck\iheck + lsl_wagecorr|] = G[|iheck\iheck + lsl_wagecorr|] :+ cross(lsl_Wobs, LnWresPur:^2 :/ SigmaW:^4 :- 1 :/ SigmaW:^2) :* B[|iheck\iheck + lsl_wagecorr|]
            }
        //}
    }
    
    // Calculate dude share
    if (todo == 0 | lsl_lambda > 0) st_numscalar("lsl_dudes", sum(Dude :< 0) / (nobs * lsl_draws))
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
    external real matrix    lsl_Brnd
    external real scalar    lsl_bfix
    external real scalar    lsl_boxcc
    external real scalar    lsl_lC, lsl_lL1, lsl_lL2
    
    external real matrix    lsl_X               // Right hand side variables
    external real scalar    lsl_groups          // Number of groups
    external real scalar    lsl_draws           // Number of random draws
    external real colvector lsl_Y               // Left hand side variable
    external real matrix    lsl_J               // Panel setup information
    external real matrix    lsl_R               //   Halton sequences    
    external real colvector lsl_Rvars           // Random coefficients
    external real scalar    lsl_rvars           // Number of random coefficients
    external real scalar    lsl_corr            //   Enable correlation?
    external string scalar  lsl_ufunc           // Functional form
    external real scalar    lsl_wagep           // Wage Prediction Error?
    external real matrix    lsl_Wpred           //   Prediction dummies
    external real colvector lsl_Days            //   Number of days per tax year
    external real matrix    lsl_Hwage           //   Hourly wage rates
    external real rowvector lsl_Sigma           //   Variance of the wage regression
    external real matrix    lsl_Hours           //   Hours of work
    external real rowvector lsl_TaxregB         // Tax Regression
    external real matrix    lsl_TaxregVars      //   Wage independent variables of tax regression
    external real matrix    lsl_TaxregIas1      //   Wage interaction variables of tax regression
    external real matrix    lsl_TaxregIas2      //   Wage interaction variables of tax regression
    external real scalar    lsl_round           // To round, or not to round.

    external real colvector lsl_C
    external real matrix    lsl_CX
    external real matrix    lsl_C2X
    external real colvector lsl_L1
    external real matrix    lsl_LX1
    external real matrix    lsl_L2X1
    external real matrix    lsl_L2       // Buggy, L2 is colvector in fact
    external real matrix    lsl_LX2
    external real matrix    lsl_L2X2
    external real matrix    lsl_Xind
    
    real scalar nobs, n, r, i, e, c, iRV, nRV, nrep, wp, ncons, getutils, getprobs, getdudes
    real colvector U, P, D, Un, Pn, Unr, Enr, Pnr, C, L1, BcC, BcL1
    real rowvector Beta, Zeta, Bsig
    real matrix CholB, Xnr, CX, C2X, LX1, L2X1, L2, BcL2, LX2, L2X2, Xind, Wn, Mwage, TaxregX1, TaxregX2, TaxregX, Dude
    
    i     = 1            // Indicates first observation of active group
    nobs  = rows(lsl_Y)  // Number of observations
    nRV   = lsl_rvars + 1
    ncons = cols(lsl_CX) + cols(lsl_C2X) + 1 + cols(lsl_L2)
    
    U = J(nobs, 1, 0)
    P = J(nobs, 1, 0)
    D = J(nobs, 1, 0)
    
    getutils = (sum(opt :== "xb") == 1)
    getprobs = (sum(opt :== "pc1") == 1)
    getdudes = (sum(opt :== "dudes") == 1)
    if (getdudes == 1) Dude = J(nobs, lsl_draws, 0)
    
    Bsig  = lsl_Sigma
    CholB = (lsl_corr == 1 ? lowertriangle(invvech(lsl_Brnd')) : diag(lsl_Brnd'))   // Build variance-(covariance) matrix
    
    // Initialize random coefficients vector
    if (lsl_rvars > 0) {
        Zeta = J(rows(lsl_R), lsl_bfix, 0)
        Zeta[.,lsl_Rvars] = cross(lsl_R[|1,1\.,lsl_rvars|]', CholB')
    }
    
    for (n = 1; n <= lsl_groups; n++) {
        i   = lsl_J[n,1]
        e   = lsl_J[n,2]
        c   = e - i + 1
        Xnr = lsl_X[|i,1\e,.|]
        
        // Fetch needed right hand side parts
        if (lsl_wagep == 1 | getdudes == 1) {
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
        
        wp   = (lsl_wagep == 1 & sum(lsl_Wpred[|i,1\e,.|]) > 0)
        nrep = (lsl_rvars > 0 | wp == 1 ? lsl_draws : 1)
        
        Un = J(c, 1, 0)
        Pn = J(c, 1, 0)
        for (r = 1; r <= nrep; r++) {
            // Indicates the active Halton sequence
            iRV  = lsl_draws * (n - 1) + r
            
            if (wp == 1) {
                //
                // Calculate monthly earnings
                //

                // Adjust wages with random draws if prediction enabled
                if (lsl_wagep == 1) Wn = Wn :* exp(cross(Bsig' :* lsl_R[|iRV,nRV\iRV,.|]', lsl_Wpred[|i,1\e,.|]'))'

                // Calculate monthly earnings
                Mwage = (lsl_Days[|i\e|] :/ 12 :/ 7) :* lsl_Hours[|i,1\e,.|] :* Wn

                // Round monthly earnings if enabled
                if (lsl_round == 1) Mwage = round(Mwage, 0.01)

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
                C = rowmax((cross(TaxregX', lsl_TaxregB') :/ lsl_boxcc, J(c, 1, 1)))

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
            
            // Build (random?) coefficients vector
            Beta = lsl_Bfix :+ (lsl_rvars > 0 ? Zeta[iRV,.] : 0)
            
            // Calculate systematic utility and choice probability
            Unr = cross(Xnr', Beta')
            if (getprobs == 1) {
                Enr = exp(Unr :+ colmin(-mean(Unr) \ 700 :- colmax(Unr)))   // Standardize to avoid missings
                Pnr = Enr :/ colsum(Enr)
            }
            
            // Calculate dudes
            if (getdudes == 1) {
                if    (lsl_ufunc == "boxcox") Dude[|i,r\e,(nrep == 1 ? . : r)|] = cross((CX, lsl_boxcox(L1, lsl_lL1), lsl_boxcox(L2, lsl_lL2))', Beta[|1\ncons|]') :*
                                                                                  (reldif(lsl_lC, 0) >= 1e-25 ? C:^(lsl_lC - 1) : (1 :/ C))
                else if (lsl_ufunc == "quad") Dude[|i,r\e,(nrep == 1 ? . : r)|] = cross((CX, 2 :* C2X :* C, L1, L2)', Beta[|1\ncons|]')
                else if (lsl_ufunc == "tran") Dude[|i,r\e,(nrep == 1 ? . : r)|] = cross((CX, 2 :* C2X :* log(C), log(L1), log(L2))', Beta[|1\ncons|]') :/ C
            }
            
            // Sum up
            if (getutils == 1) Un = Un :+ Unr
            if (getprobs == 1) Pn = Pn :+ Pnr
        }

        // Next household
        if (getutils == 1) U[|i\e|] = Un :/ nrep
        if (getprobs == 1) P[|i\e|] = Pn :/ nrep
        if (getdudes == 1) D[|i\e|] = rowsum(Dude[|i,1\e,.|] :< 0) :/ nrep
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
    st_store(., newvar, touse, result)
}

mata set matastrict off
end


cap program drop lslogit_p
/**
 * Conditional Logit but integrating out wage prediction errors (Wrapper programm)
 * 
 * @param `xb'  Predict systematic utility
 * @param `pc1' Predict choice probabilities
 */
program define lslogit_p
    syntax newvarlist(min=1 max=3) [if] [in] [, xb pc1 dudes]
    
    
    //
    // Initialize
    //
    
    // Predict command works only with lslogit estimates
    if ("`e(cmd)'" != "lslogit") error 301
    
    // What to predict?
    local n_varlist : word count `varlist'
    if ("`xb'" != "" & "`pc1'" != "") {
        if (`n_varlist' != 2) {
            di in r "either option 'xb' or 'pc1' allowed"
            exit 498
        }
        else local opt pc1 xb
    }
    else if ("`xb'" == "" & "`pc1'" == "" & "`dudes'" == "") {
        if      (`n_varlist' == 3) local opt pc1 xb dudes
        else if (`n_varlist' == 2) local opt pc1 xb
        else {
            di as text "(option pc1 assumed; probability of success given one success within group)"
            local opt pc1    // Default: pc1
        }
    }
    else local opt `xb'`pc1'`dudes'
    
    // Mark the estimation sample
    marksample touse, novarlist
    markout `touse' `e(depvar)' `e(group)' `e(consum)' `e(leisure)' `e(cx)' `e(lx1)' `e(lx2)' `e(c2x)'  ///
                    `e(l2x1)' `e(l2x2)' `e(indeps)' `e(wagepred)' `e(days)' `e(tria1)' `e(tria2)'
    
    // Restrict sample to those of interest
    qui sort `e(group)' // BUGGY
    preserve
    qui keep if `touse'
    qui count
    local nobs = r(N)
    
    
    //
    // Prepare data
    //
    
    local n_leisure  : word count `e(leisure)'
    local n_cxias    : word count `e(cx)'
    local n_lx1ias   : word count `e(lx1)'
    local n_lx2ias   : word count `e(lx2)'
    local n_c2xias   : word count `e(c2x)'
    local n_l2x1ias  : word count `e(l2x1)'
    local n_l2x2ias  : word count `e(l2x2)'
    local n_indeps   : word count `e(indeps)'
    local n_randvars : word count `e(randvars)'
    local n_wagep    : word count `e(wagepred)'
    local n_hecksig  : word count `e(hecksigma)'
    local n_hwage    : word count `e(hwage)'
    local n_taxrias1 : word count `e(taxreg_ia1)'
    local n_taxrias2 : word count `e(taxreg_ia2)'
    local n_heckvars : word count `e(heckman)'      // Not yet implemented
    local n_selvars  : word count `e(select)'       // Not yet implemented
    
    local rvars      : word count `e(randvars)'
    
    
    //
    // Move data to Mata
    //
    
    mata: lsl_J = panelsetup(st_data(., "`e(group)'"), 1)   // Choices per group
    mata: lsl_groups = rows(lsl_J)                          // Number of groups
    
    mata: lsl_ufunc = "`e(ufunc)'"              // Utility function
    mata: lsl_nlei  = `n_leisure'               // Utility function
    mata: lsl_draws = strtoreal("`e(draws)'")   // Number of draws
    mata: lsl_burn  = strtoreal("`e(burn)'")    // Number of draws to burn
    mata: lsl_round = ("`e(round)'" == "1")     // To round, or not to round?
    mata: lsl_rvars = `rvars'                                                                           // Random coefficients
    mata: lsl_Rvars = ("`e(randvars)'" != "" ? strtoreal(tokens("`e(randvars)'"))' : J(0, 1, 0))        // Random coefficients?
    mata: lsl_corr  = ("`e(randvars)'" != "" & "`e(corr)'" == "1" ? 1 : 0)                              // Correlation?
    mata: lsl_Sigma = ("`e(hecksigma)'" != "" ? strtoreal(tokens("`e(hecksigma)'")) : J(1, `n_leisure', 0))
    mata: lsl_wagep = `e(wagep)'                                                                        // Run Wage Prediction
    mata: lsl_Wpred = (lsl_wagep == 1 ? st_data(., "`e(wagepred)'") : J(`nobs', lsl_nlei, 0))           // Dummies enabling or disabling the wage prediction
    mata: lsl_Hwage = ("`e(hwage)'" != "" ? st_data(., tokens("`e(hwage)'")) : J(`nobs', lsl_nlei, 0))  // Hourly wage rates
    mata: lsl_Days  = ("`e(days)'" != "" ? st_data(., "`e(days)'") : J(`nobs', 1, 365))                 // Days of taxyear
    
    
    //
    // Tax regression
    //
    
    if ("`e(taxreg)'" == "1") {
        mata: lsl_TaxregB    = st_matrix("e(taxreg_b)")                                            // Tax regression estimates
        mata: lsl_TaxregIas1 = ("`e(taxreg_ia1)'" != "" ? st_data(., tokens("`e(taxreg_ia1)'")) : J(`nobs', 0, 0))   // Interaction variables on Mwage1 and Mwage1^2
        mata: lsl_TaxregIas2 = ("`e(taxreg_ia2)'" != "" ? st_data(., tokens("`e(taxreg_ia2)'")) : J(`nobs', 0, 0))   // Interaction variables on Mwage2 and Mwage2^2
        mata: lsl_TaxregVars = st_data(., tokens("`e(taxreg_vars)'"))                                   // Variables that are independent of m_wage
    }
    
    // Box-Cox normalizing constants
    if ("`e(ufunc)'" == "boxcox" & "`e(boxcc)'" != "" & "`e(boxcl)'" != "") {
        mata: lsl_boxcc = (`e(boxcc)' > 0 & `e(boxcc)' < . ? `e(boxcc)' : 1)
        mata: lsl_boxcl = (`e(boxcl)' > 0 & `e(boxcl)' < . ? `e(boxcl)' : 1)
    }
    else mata: lsl_boxcc = lsl_boxcl = 1
    
    mata: lsl_Y     = st_data(., "`e(depvar)'")                                                             // Dependent variable
    
    // Consumption, squared and interactions
    mata: lsl_C   = st_data(., "`e(consum)'") :/ lsl_boxcc
    mata: lsl_CX  = ((`n_cxias'  > 0 ? st_data(., tokens("`e(cx)'"))   ///
                                     : J(`nobs', 0, 0)), J(`nobs', 1, 1))
    mata: lsl_C2X = ((`n_c2xias' > 0 ? st_data(., tokens("`e(c2x)'"))  ///
                                     : J(`nobs', 0, 0)), J(`nobs', (lsl_ufunc != "boxcox" ? 1 : 0), 1))

    // Leisure, squared and interactions
    forval i = 1/2 {
        local var : word `i' of `e(leisure)'
        if ("`var'" != "") {
            mata: lsl_L`i'   = st_data(., "`var'") :/ lsl_boxcl
            mata: lsl_LX`i'  = ((`n_lx`i'ias'  >  0 ? st_data(., tokens("`e(lx`i')'"))     ///
                                                    : J(`nobs', 0, 0)), J(`nobs', 1, 1))
            mata: lsl_L2X`i' = ((`n_l2x`i'ias' >  0 ? st_data(., tokens("`e(l2x`i')'"))    ///
                                                    : J(`nobs', 0, 0)), J(`nobs', (lsl_ufunc != "boxcox" ? 1 : 0), 1))
        }
        else {
            mata: lsl_L`i'   = J(`nobs', 0, 0)
            mata: lsl_LX`i'  = J(`nobs', 0, 0)
            mata: lsl_L2X`i' = J(`nobs', 0, 0)
        }
    }
    // Hypothetical hours of work
    mata: lsl_Hours = `e(totaltime)' :- (lsl_L1, lsl_L2) :* lsl_boxcl  // caution: Box-Cox normalization!

    // Dummy variables
    mata: lsl_Xind = (`n_indeps' > 0 ? st_data(., tokens("`e(indeps)'")) : J(`nobs', 0, 0))
    
    // Build coefficients vector
    mata: lsl_bfix = cols(lsl_CX)  + cols(lsl_C2X) + lsl_nlei + ///
                     cols(lsl_LX1) + cols(lsl_L2X1) + cols(lsl_LX2) + cols(lsl_L2X2) + ///
                     (lsl_nlei == 2) + cols(lsl_Xind)
    mata: lsl_B    = st_matrix("e(b)")
    mata: lsl_Bfix = lsl_B[|1\lsl_bfix|]
    mata: lsl_Brnd = (lsl_rvars > 0 ? lsl_B[|lsl_bfix + (lsl_ufunc == "boxcox" ? 1 + lsl_nlei : 0) + 1\.|] : J(0, 0, 0))     // Get auxiliary random coefficients
    
    // Build right hand side
    if ("`e(ufunc)'" == "tran") {           // Translog utility
        mata: lsl_X = (log(lsl_C)  :* (lsl_CX,  log(lsl_C)  :* lsl_C2X,   log(lsl_L1),   log(lsl_L2)),      ///
                       log(lsl_L1) :* (lsl_LX1, log(lsl_L1) :* lsl_L2X1),                                   ///
                       log(lsl_L2) :* (lsl_LX2, log(lsl_L2) :* lsl_L2X2), log(lsl_L1) :* log(lsl_L2), lsl_Xind)
    }
    else if ("`e(ufunc)'" == "quad") {      // Quadratic utility
        mata: lsl_X = (lsl_C  :* (lsl_CX,  lsl_C  :* lsl_C2X,   lsl_L1,   lsl_L2),      ///
                       lsl_L1 :* (lsl_LX1, lsl_L1 :* lsl_L2X1),                         ///
                       lsl_L2 :* (lsl_LX2, lsl_L2 :* lsl_L2X2), lsl_L1 :* lsl_L2, lsl_Xind)
    }
    else if ("`e(ufunc)'" == "boxcox") {    // Box-Cox utility
        mata: lsl_lC  = lsl_B[lsl_bfix + 1]
        mata: lsl_lL1 = lsl_B[lsl_bfix + 2]
        mata: lsl_lL2 = (lsl_nlei == 2 ? lsl_B[lsl_bfix + 3] : 0)
        mata: lsl_X   = (lsl_boxcox(lsl_C, lsl_lC)   :* (lsl_CX,  lsl_boxcox(lsl_L1, lsl_lL1),   lsl_boxcox(lsl_L2, lsl_lL2)),    ///
                         lsl_boxcox(lsl_L1, lsl_lL1) :*  lsl_LX1,                                                                 ///
                         lsl_boxcox(lsl_L2, lsl_lL2) :*  lsl_LX2, lsl_boxcox(lsl_L1, lsl_lL1) :* lsl_boxcox(lsl_L2, lsl_lL2), lsl_Xind)
    }
    
    // Generate Halton sequences
    mata: lsl_R = (`rvars' + `n_wagep' > 0 ? invnormal(halton(lsl_groups * lsl_draws, `rvars' + `n_wagep', 1 + lsl_burn)) : J(`nobs', 0, 0))
    
    // Restore sample
    restore
    
    
    //
    // Predict
    //
    
    foreach v of local varlist {
        qui gen double `v' = .
    }
    mata: lslogit_p(tokens("`varlist'"), "`touse'", tokens("`opt'"))
    
    
    //
    // Clean up
    //
    
    foreach m in round ufunc Weight Y Hwage boxcc boxcl C CX C2X L1 LX1 L2X1 L2 LX2 L2X2 Xind X joint ///
                 HeckmVars SelectVars Days Hours wagep Wpred Sigma TaxregB TaxregIas1 TaxregIas2 TaxregVars ///
                 groups J draws burn Rvars corr R lC lL1 lL2 B Bfix bfix rvars {
        cap mata mata drop lsl_`m'
    }
end


cap program drop lslogit_cov
/**
 * Conditional Logit but integrating out wage prediction errors (Wrapper programm)
 */
program define lslogit_cov
    syntax [, cov NOHEADer NOWAGE NOPREF]
    
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
    if (`e(wagecorr)' == 1 & "`nowage'" == "") {
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