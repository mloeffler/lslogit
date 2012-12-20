/*****************************************************************************
 *
 * lslogit -- ESTIMATING MIXED LOGIT LABOR SUPPLY MODELS WITH STATA
 * 
 * (c) 2012 - Max Löffler
 *
 *****************************************************************************/

cap program drop lslogit
/**
 * Conditional Logit but integrating out wage prediction errors (Wrapper programm)
 * 
 * @param `group'  Group identifier variable
 * @param `taxreg' Stored estimates of the tax regression
 */
program define lslogit
    if (replay()) {
        if (`"`e(cmd)'"' != "lslogit")   error 301
        lslogit_Replay `0'
    }
    else lslogit_Estimate `0'
end

cap program drop lslogit_Replay
/**
 * Conditional Logit but integrating out wage prediction errors (Wrapper programm)
 * 
 * @param `group'  Group identifier variable
 * @param `taxreg' Stored estimates of the tax regression
 */
program define lslogit_Replay
    syntax [, Level(integer `c(level)') Quiet]
    
    // Set up auxiliary stuff
    local diparm
    if ("`quiet'" == "") {
        foreach aux in sigma_w1 sigma_w2 dudes {
            if (e(`aux') != .) {
                local val = e(`aux')
                local diparm `diparm' diparm(__lab__, value(`val') label("[`aux']"))
            }
        }
    }
    
    // Display output
    ml display, level(`level') `diparm'
                               /*diparm(ln_consum, f(`dudes') d(0) label("% dU/dc>=0"))*/
end

cap program drop lslogit_Estimate
/**
 * Conditional Logit but integrating out wage prediction errors (Wrapper programm)
 * 
 * @param `group' varname  Group identifier variable
 * @param `taxreg' name Stored estimates of the tax regression
 * @param `burn' integer Number of initial Halton draws to burn
 */
program define lslogit_Estimate, eclass
    syntax varname(numeric) [if] [in] [fweight/], GRoup(varname numeric) Ufunc(name)                                        ///
                                                  Consumption(varname numeric) Leisure(varlist numeric min=1 max=2)         ///
                                                  [cx(varlist numeric)  lx1(varlist numeric)  lx2(varlist numeric)          ///
                                                   c2x(varlist numeric) l2x1(varlist numeric) l2x2(varlist numeric)         ///
                                                   INDeps(varlist) TOTALTime(integer 80) DAYs(varname numeric)              ///
                                                   boxc(integer 1000) boxl(integer 80) HWage(varlist numeric min=1 max=2)   ///
                                                   TAXReg(name) tria1(varlist numeric) tria2(varlist numeric)               ///
                                                   WAGEPred(varlist numeric min=1 max=2) HECKSIGma(numlist min=1 max=2)     ///
                                                   RANDvars(string) corr DRaws(integer 50) burn(integer 15)                 ///
                                                   HECKMan(varlist) SELect(varlist)                                         ///
                                                   noround Quiet Verbose                                                    ///
                                                   difficult trace search(name) iterate(integer 100) method(name)           ///
                                                   gradient hessian debug Level(integer `c(level)') from(string)            ///
                                                   technique(string)]
    
    /* INITIALIZE ESTIMATOR
     */
    
    // Mark the estimation sample
    marksample touse
    markout `touse' `varlist' `group' `consumption' `leisure' `cx' `lx1' `lx2' `c2x' `l2x1' `l2x2'  ///
                    `indeps'`wagepred' `days' `tria1' `tria2' `heckman' `select'
    
    // Verbose mode
    if ("`verbose'" == "") local qui qui
    
    // Validate Maximum Likelihood method
    if ("`method'" == "") local method d2
    if (!inlist("`method'", "d0", "d1", "d2")) {
        di in r "method must be either 'd0', 'd1' or 'd2'"
        exit 498
    }
    
    // Validate utility function
    if (!inlist("`ufunc'", "quad", "tran", "boxcox")) {
        di in r "utility function must be either 'quad', 'tran' or 'boxcox'"
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
        qui count if `consumption' < 0 & `touse'
    }
    // Check for negative values
    if (r(N) > 0) {
        di in r "consumption contains values smaller or equal to zero"
        exit 498
    }
    // If Box-Cox, make it easy.
    if ("`ufunc'" == "boxcox" & ("`c2x'" != "" | "`l2x1'" != "" | "`l2x2'" != "")) {
        di in r "options c2x(), l2x1() and l2x2() not allowed with Box-Cox utility function"
        exit 498
    }
    // If Box-Cox, make it easy.
    if ("`ufunc'" == "boxcox" & ("`heckman'" != "" | "`select'" != "" | "`randvars'" != "" | "`wagepred'" != "")) {
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
    
    // Validate joint estimation settings
    if ("`select'" != "" & "`heckman'" == "") {
        di in r "option heckman() required when estimating jointly"
        exit 498
    }
    if ("`heckman'" != "" & "`hwage'" == "") {
        di in r "option hwage() required when estimating jointly"
        exit 498
    }
    if ("`heckman'" != "") {
        qui count if `hwage' == .
        if ("`select'" != "" & r(N) == 0) {
            di in r "wage variable never censored because of selection"
            exit 498
        }
        else if ("`select'" == "" & r(N) > 0) {
            di in r "wage variable censored, use option select()"
            exit 498
        }
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
    
    // Validate Wage Prediction Options
    if (`n_wagep' == 0) {
        local wagep = 0     // No wage prediction
    }
    else {
        // Wage prediction enabled
        if (`n_wagep' == `n_leisure' & `n_wagep' == `n_hwage' & `n_wagep' == `n_hecksig') {
            tempvar preds
            qui egen `preds' = rowtotal(`wagepred') if `touse'
            qui count if inlist(`preds', 1, 2) & `touse'
            local wagep = (r(N) > 0)
        }
        // Settings incorrect
        else {
            di in r "number of wage prediction variables does not match the number of leisure terms, hourly wage rates or mean squared errors"
            exit 498
        }
    }
    
    // No need for random variables
    if (`wagep' == 0 & "`randvars'" == "") local draws = 1
    
    // Tax regression or tax benefit calculator needed
    if (`wagep' == 1 | "`heckman'" != "") {
        if ("`taxreg'" != "" & "`taxben'" == "") {
            tempname taxreg_from
            // Load tax regression estimates
            qui est restore `taxreg'
            mat `taxreg_from' = e(b)
            local taxreg_betas : colnames `taxreg_from'
            local n_taxreg_betas : word count `taxreg_betas'
            local taxreg_vars
            local start = 1 + 2 * (`n_leisure' + `n_taxrias1' + `n_taxrias2')         // + `n_leisure' * 5
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
            di in r "either option taxreg() or option taxben() required"
            exit 198
        }
    }
    
    // Build weight settings
    if ("`weight'" != "")   local wgt "[`weight'=`exp']"
    
    // Select random variables
    local rvars = `n_randvars' + `n_wagep'
    
    
    /* LOOK FOR INITIAL VALUES
     */
    
    if ("`search'" != "off" & "`from'" == "") {
        // Verbose mode on?
        if ("`verbose'" != "") di as text "Looking for initial values..."
        
        // Set up consumption and leisure
        tempvar c l1 l2
        if (inlist("`ufunc'", "tran", "boxcox")) {
            qui gen `c' = log(cond("`ufunc'" == "boxcox" & `boxc' > 0, `consumption' / `boxc', `consumption')) if `touse'
            foreach var of local leisure {
                if (strpos("`leisure'", "`var'") == 1) local lei l1
                else                                   local lei l2
                qui gen ``lei'' = log(cond("`ufunc'" == "boxcox" & `boxl' > 0, `var' / `boxl', `var')) if `touse'
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
        `qui' clogit `varlist' `initrhs' if `touse' `wgt', group(`group') iterate(25)
        if (e(converged) == 1) {
            // Save results
            mat `init_from' = e(b)
            local nobs      = e(N)
            local k         = e(k)
            local ll        = e(ll)
            // Update sample
            qui replace `touse' = e(sample)
        }
        
        // Wage equation
        if ("`heckman'" != "") {
            tempname init_wage init_w
            if (`n_hwage' == 0) mat `init_wage' = J(1, 1 + `n_heckvars' * `n_leisure', 0)
            else {
                foreach w of local hwage {
                    tempvar ln`w'
                    qui gen `ln`w'' = ln(`w') if `touse' & `varlist'
                    if ("`select'" != "") {
                        `qui' heckman `ln`w'' `heckman' if `touse' & `varlist' `wgt', select(`select')
                        mat `init_w' = e(b)
                        mat `init_w'[1,colsof(`init_w')-1] = e(rho)
                        mat `init_w'[1,colsof(`init_w')] = e(sigma)
                        //mat `init_w' = `init_w'[1,1..colsof(`init_w')-1]
                    }
                    else {
                        `qui' reg `ln`w'' `heckman' if `touse' & `varlist' `wgt'
                        mat `init_w' = e(b)
                    }
                    mat `init_wage' = (nullmat(`init_wage'), `init_w')
                }
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
    
    mata: lsl_round  = ("`round'" != "noround")                                                     // To round, or not to round?
    mata: lsl_ufunc  = st_local("ufunc")                                                            // Utility function
    mata: lsl_Weight = ("`exp'" != "" ? st_data(., "`exp'") : J(`nobs', 1, 1))                      // Weight
    mata: lsl_Y      = st_data(., "`varlist'")                                                      // Left hand side
    mata: lsl_Hwage  = (`n_hwage' > 0 ? st_data(., tokens("`hwage'")) : J(`nobs', `n_leisure', 0))  // Hourly wage rates
    
    //
    // Right hand side
    //
    
    // Consumption, squared and interactions
    mata: lsl_C   = st_data(., "`consumption'") :/ (lsl_ufunc == "boxcox" & `boxc' > 0 ? `boxc' : 1)
    mata: lsl_CX  = ((`n_cxias'  > 0 ? st_data(., tokens("`cx'"))  : J(`nobs', 0, 0)), J(`nobs', 1, 1))
    mata: lsl_C2X = ((`n_c2xias' > 0 ? st_data(., tokens("`c2x'")) : J(`nobs', 0, 0)), J(`nobs', ("`ufunc'" != "boxcox" ? 1 : 0), 1))
    
    // Leisure, squared and interactions
    forval i = 1/2 {
        local var : word `i' of `leisure'
        if ("`var'" != "") {
            mata: lsl_L`i'   = st_data(., "`var'") :/ (lsl_ufunc == "boxcox" & `boxl' > 0 ? `boxl' : 1)
            mata: lsl_LX`i'  = ((`n_lx`i'ias'  >  0 ? st_data(., tokens("`lx`i''"))  : J(`nobs', 0, 0)), J(`nobs', 1, 1))
            mata: lsl_L2X`i' = ((`n_l2x`i'ias' >  0 ? st_data(., tokens("`l2x`i''")) : J(`nobs', 0, 0)), J(`nobs', ("`ufunc'" != "boxcox" ? 1 : 0), 1))
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
    
    mata: lsl_heckm      = ("`heckman'" != "")                                                                                       // Run joint estimation?
    mata: lsl_HeckmVars  = (lsl_heckm == 1 ? (st_data(., tokens("`heckman'")), J(`nobs', 1, 1)) : J(`nobs', 0, 0))                    // Wage variables
    mata: lsl_SelectVars = (lsl_heckm == 1 & "`select'" != "" ? (st_data(., tokens("`select'")), J(`nobs', 1, 1)) : J(`nobs', 0, 0))  // Wage variables
    mata: lsl_Days       = ("`days'" != "" ?  st_data(., st_local("days")) : J(`nobs', 1, 365))                                      // Days of taxyear
    mata: lsl_Hours      = `totaltime' :- (lsl_L1, lsl_L2)                                                                             // Hypothetical hours
    
    //
    // Wage Prediction Stuff
    //
    
    mata: lsl_wagep = `wagep'                                                            // Run Wage Prediction
    if (`wagep' == 1) {
        mata: lsl_Wpred = st_data(., ("`wagepred'"))                                     // Dummies enabling or disabling the wage prediction
        mata: lsl_Sigma = J(1, `n_hecksig', 0)                                           // Estimated variance of Heckman correction
        forval i = 1/`n_hecksig' {
            local sig : word `i' of `hecksigma'
            mata: lsl_Sigma[1,`i'] = `sig'
        }
    }
    
    //
    // Tax regression
    //
    if ("`taxreg'" != "") {
        mata: lsl_TaxregB    = st_matrix("`taxreg_from'")                                            // Tax regression estimates
        mata: lsl_TaxregIas1 = ("`tria1'" != "" ? st_data(., tokens("`tria1'")) : J(`nobs', 0, 0))   // Interaction variables on Mwage1 and Mwage1^2
        mata: lsl_TaxregIas2 = ("`tria2'" != "" ? st_data(., tokens("`tria2'")) : J(`nobs', 0, 0))   // Interaction variables on Mwage2 and Mwage2^2
        mata: lsl_TaxregVars = st_data(., tokens("`taxreg_vars'"))                                   // Variables that are independent of m_wage
    }
    
    //
    // Group level stuff
    //
    qui duplicates report `group'
    mata: lsl_groups = st_numscalar("r(unique_value)")   // Number of groups
    tempvar choices
    by `group': gen `choices' = _N
    mata: lsl_J = st_data(., st_local("choices"))        // Choices per group
    
    //
    // Random draws
    //
    mata: lsl_draws = strtoreal(st_local("draws"))                                                                   // Number of draws
    mata: lsl_burn  = strtoreal(st_local("burn"))                                                                    // Number of draws to burn
    mata: st_local("randvars", invtokens(strofreal(sort(strtoreal(tokens("`randvars'"))', 1))'))                    // Sort random coefficients
    mata: lsl_Rvars = ("`randvars'" != "" ? strtoreal(tokens("`randvars'"))' : J(0, 1, 0))                           // Random coefficients
    mata: lsl_corr  = ("`corr'" != "")                                                                               // Random coefficients correlated?
    mata: lsl_R     = (`rvars' > 0 ? invnormal(halton(lsl_groups*lsl_draws, `rvars', 1+lsl_burn)) : J(`nobs', 0, 0))    // Halton sequences
    
    // Restore data
    restore
    
    
    /* RUN ESTIMATION
     */
    
    if ("`verbose'" != "") di as text "Run estimation..."
    
    // Set up equations
    local eq_consum (C: `varlist' = `cx')                               // Consumption
    if ("`ufunc'" != "boxcox") local eq_consum `eq_consum' (CC: `c2x')  // Consumption^2
    local eq_leisure
    foreach var of local leisure {
        local i = 1 + (strpos("`leisure'", "`var'") > 1)
        local eq_leisure `eq_leisure' (L`i': `lx`i'')                                   // Leisure
        if ("`ufunc'" != "boxcox") local eq_leisure `eq_leisure' (L`i'L`i': `l2x`i'')   // Leisure^2
        local eq_consum  `eq_consum' /CXL`i'                                            // Consumption X leisure interaction
    }
    if (`n_leisure' == 2) local eq_leisure  `eq_leisure' /L1XL2         // Leisure term interaction
    if (`n_indeps'  >  0) local eq_indeps   (IND: `indeps', noconst)    // Independent variables / dummies
    // Box-Cox utility function?
    if ("`ufunc'" == "boxcox") {
        local eq_boxcox /l_C /l_L1
        if (`n_leisure' == 2) local eq_boxcox `eq_boxcox' /l_L2
        if ("`initopt'" != "" & "`init_from'" != "") mat `init_from' = (`init_from', 0.0, J(1, `n_leisure', 0.0))
    }
    // Random coefficients?
    if (`n_randvars' > 0) {
        local eq_rands
        if ("`corr'" == "") {
            forval i = 1/`n_randvars' {
                local sd : word `i' of `randvars'
                local eq_rands `eq_rands' /sd_`sd'
            }
            if ("`initopt'" != "" & "`init_from'" != "") mat `init_from' = (`init_from', J(1, `n_randvars', 0.0001))
        }
        else {
            forval i = 1/`n_randvars' {
                local a : word `i' of `randvars'
                forval k = `i'/`n_randvars' {
                    local b : word `k' of `randvars'
                    if (`a' == `b') local lab sd_`a'
                    else local lab s_`a'_`b'
                    local eq_rands `eq_rands' /`lab'
                }
            }
            if ("`initopt'" != "" & "`init_from'" != "") mat `init_from' = (`init_from', J(1, `n_randvars' * (`n_randvars' + 1) / 2, 0.0001))
        }
    }
    // Joint wage estimation?
    if ("`heckman'" != "") {
        local eq_heckm (lnW: `heckman')
        if ("`select'"  != "") local eq_heckm `eq_heckm' (S: `select') /rho /sigma
        if ("`initopt'" != "" & "`init_from'" != "") mat `init_from' = (`init_from', `init_wage')
    }
    
    //
    // Estimate
    //
    
    ml model `method'`debug' lslogit_d2() `eq_consum' `eq_leisure' `eq_indeps' `eq_boxcox' `eq_rands' `eq_heckm' ///
            if `touse' `wgt', group(`group') `initopt' search(off) iterate(`iterate') nopreserve max `difficult' `trace' `gradient' `hessian' technique(`technique')
    
    // Save results
    foreach aux in sigma_w1 sigma_w2 dudes {
        if (r(`aux') != .) ereturn scalar `aux' = r(`aux')
    }
    ereturn local  title    "Mixed Logit Labor Supply Model"
    ereturn local  cmd      "lslogit"
    ereturn local  predict  "lslpred"
    ereturn local  depvar    `varlist'
    ereturn local  group     `group'
    ereturn local  ufunc    "`ufunc'"
    ereturn scalar draws   = `draws'
    if ("`ufunc'" == "boxcox") ereturn scalar k_aux = 1 + `n_leisure'
    else if ("`select'" != "") ereturn scalar k_aux = 2
    else {
        if ("`corr'" != "") {
            ereturn scalar corr  = 1
            ereturn scalar k_aux = `n_randvars' * (`n_randvars' + 1) / 2
        }
        else {
            ereturn scalar corr = 0
            ereturn scalar k_aux = `n_randvars'
        }
    }
    
    // Show results
    lslogit_Replay, level(`level') `quiet'
end

cap mata mata drop lslogit_d2()
cap mata mata drop lsl_boxcox()
cap mata mata drop lsl_boxcox_g()
cap mata mata drop lsl_boxcox_h()
mata:
mata set matastrict on
/**
 * Standard Conditional Logit but integrating out wage prediction errors (Evaluator)
 * 
 * @param B_s Stata matrix of coefficients
 */
void lslogit_d2(transmorphic scalar ML, real scalar todo, real rowvector B,
                real scalar lnf, real rowvector G, real matrix H) {
    
    external string scalar  lsl_ufunc           // Functional form
    external real scalar    lsl_groups          // Number of groups
    external real colvector lsl_Y               // Left hand side variable
    external real colvector lsl_J               // Number of choices per group
    external real matrix    lsl_X               // Right hand side variables
    external real colvector lsl_Weight          // Group weights
    
    external real scalar    lsl_draws           // Number of random draws
    external real matrix    lsl_R               //   Halton sequences
    
    external real colvector lsl_Rvars           // Random coefficients
    external real scalar    lsl_corr            //   Enable correlation?
    
    external real scalar    lsl_heckm           // Joint wage estimation?
    external real matrix    lsl_HeckmVars       //   Right hand side variables
    external real matrix    lsl_SelectVars      //   Selection variables
    
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
    
    //
    // Declarations
    //
    
    real scalar    n, i, c, e, nobs, nlei, ncons, dudes
    real scalar    b, brnd, bheck, bsel, bfix
    real rowvector Bfix, Brnd, Bheck, Bsel, Brho, Bsig, Beta
    
    real scalar    nRV, rvars, iRV, r, rv, rv2, iRow, iCol, nRows, nCols
    real matrix    Sigm
    
    real scalar    hv
    real matrix    Hwage, Hwres, Select, SelRes, Lambda, Wn
    real matrix    Mwage, TaxregX1, TaxregX2, TaxregX, DCdM, D2CdM2
    
    real matrix    DUdB, DUdBr, DUdBw, DUdBs, DUdBrho, DUdBsig, DWdBw, DWdBs, DWdBrho, DWdBsig, DXdH, YmPn_D2UdB2, YmPnD2UdBw2, DUdlC, DUdlL1, DUdlL2
    real colvector DUdC, D2UdC2, DMdH, D2MdH2
    
    real scalar    ncx, nc2x, nlx1, nl2x1, nlx2, nl2x2, nxind
    real colvector Yn, C, L1, L2
    real matrix    Xnr, CX, C2X, LX1, L2X1, LX2, L2X2, Xind, PXn, YXn
    real colvector Unr, Enr, Pnr, YmPn
    
    real scalar    lsum, pni
    real rowvector Gsum, Gnr, H1, S1, W1
    real matrix    H1sum, H2sum, H2, Svar, S2xx, S2xy, W2xx, W2xy
    
    real scalar    lC, lL1, lL2
    real colvector BcC, BcL1
    real matrix    BcL2
    
    //lsl_Y = moptimize_util_depvar(ML, 1)     // Left hand side variable
    
    
    /* Setup */
    
    // Definitions
    i     = 1                                                               // Indicates first observation of active group
    nRV   = 1                                                               // Indicates next random variable to use (column of lsl_R)
    rvars = rows(lsl_Rvars)                                                  // Number of random variables
    nobs  = rows(lsl_Y)                                                      // Number of observations
    nlei  = 1 + cols(lsl_L2)                                                 // Number of leisure terms
    ncx   = cols(lsl_CX)
    nc2x  = cols(lsl_C2X)
    nlx1  = cols(lsl_LX1)
    nl2x1 = cols(lsl_L2X1)
    nlx2  = cols(lsl_LX2)
    nl2x2 = cols(lsl_L2X2)
    nxind = cols(lsl_Xind)
    ncons = ncx + nc2x + nlei  // Number of variables including consumption
    
    // Number of coefficients
    b     = cols(B)                                             // Total number
    brnd  = (lsl_corr  == 1 ? rvars * (rvars + 1) / 2 : rvars)  // Number of variance and covariance terms for random coefficients
    bheck = cols(lsl_HeckmVars)                                 // Number of Heckman wage regression coefficients
    bsel  = cols(lsl_SelectVars)                                // Number of wage selection coefficients (+ rho)
    bfix  = b - brnd - bheck - bsel - 2 * (bsel > 0)            // Number of fix preference coefficients
    
    // Maximum Likelihood Parameter
    lnf = 0             // Log-likelihood
    G   = J(1, b, 0)    // Gradient
    H   = J(b, b, 0)    // Hessian matrix
    
    //
    // Build coefficient vector
    //
    
    Bfix  = B[|1\bfix|]                                                                             // Get fixed coefficients
    Beta  = (lsl_ufunc == "boxcox" ? Bfix[|1\bfix - nlei - 1|]     : Bfix)                           //   Separate Betas from Lambdas for Box-Cox
    Brnd  = (rvars > 0 ? B[|bfix + 1\bfix + brnd|]                : J(0, 0, 0))                     // Get auxiliary random coefficients
    Sigm  = (lsl_corr == 1 ? lowertriangle(invvech(Brnd'))         : diag(Brnd'))                    //   Build variance-(covariance) matrix
    Bheck = (bheck > 0 ? B[|bfix + brnd + 1\bfix + brnd + bheck|] : J(0, 0, 0))                     // Wage coefficients
    Bsel  = (bsel  > 0 ? B[|bfix + brnd + bheck + 1\b - 2|]       : J(0, 0, 0))                     //   Selection coefficients
    //Brho = (bsel > 0 ? (exp(2 :* B[1,b - 1]) :- 1) / (exp(2 :* B[1,b - 1]) :+ 1) : J(0, 0, 0))      //   Heckman rho
    //Bsig = (bsel > 0 ? exp(B[1,b]) : J(0, 0, 0))                                                    //   Heckman sigma (lnSig)
    Brho = (bsel > 0 ? B[b - 1] : J(0, 0, 0))                                                       //   Heckman rho
    Bsig = (bsel > 0 ? B[b]     : J(0, 0, 0))                                                       //   Heckman sigma
    
    // Build matrix with random coefficients (mean zero), every row is a draw
    if (brnd > 0) {
        Brnd = J(rows(lsl_R), bfix, 0)
        for (rv = 1; rv <= rvars; rv++) {
            if (lsl_corr == 0) Brnd[.,lsl_Rvars[rv]] = lsl_R[.,rv] :* Sigm[rv,rv]
            else {
                for (rv2 = rv; rv2 <= rvars; rv2++) {
                    Brnd[.,lsl_Rvars[rv2]] = Brnd[.,lsl_Rvars[rv2]] + lsl_R[.,rv] :* Sigm[rv2,rv]
                }
            }
            nRV = nRV + 1
        }
    }
    // From now on: Beta[rows=lsl_R,cols=Bfix] = Bfix :+ Brnd
    
    //
    // Box-Cox utility function
    //
    
    if (lsl_ufunc == "boxcox") {
        real scalar ilC, ilL1, ilL2
        
        // Get Box-Cox-Lambdas
        ilC  = bfix - nlei
        ilL1 = bfix - nlei + 1
        ilL2 = bfix
        lC  = Bfix[ilC]
        lL1 = Bfix[ilL1]
        lL2 = (nlei == 2 ? Bfix[ilL2] : 0)
        
        // Transform consumption and leisure
        BcC  = lsl_boxcox(lsl_C, lC)
        BcL1 = lsl_boxcox(lsl_L1, lL1)
        BcL2 = lsl_boxcox(lsl_L2, lL2)
        
        // Replace lsl_X
        lsl_X = ((lsl_CX, BcL1, BcL2) :* BcC, lsl_LX1 :* BcL1, lsl_LX2 :* BcL2, BcL1 :* BcL2, lsl_Xind)
    }
    
    //
    // Dude share
    //
    
    if (lsl_heckm == 0 & lsl_wagep == 0 & lsl_corr == 0) {
        if (lsl_ufunc == "boxcox") DUdC = cross((lsl_CX, BcL1, BcL2)', Bfix[|1\ncons|]') :* (reldif(lC, 0) >= 1e-25 ? lsl_C:^(lC-1) : (1 :/ lsl_C))
        else                       DUdC = cross((lsl_CX, 2 :* lsl_C2X :* lsl_C, lsl_L1, lsl_L2)', Bfix[|1\ncons|]')
        if (brnd > 0) dudes = 1 - colsum(normal(DUdC :/ sqrt(colsum((lsl_Rvars :<= ncons) :* diagonal(Sigm):^2)))) :/ nobs
        else          dudes = colsum(DUdC :< 0) / nobs
        st_numscalar("r(dudes)", dudes)
    }
    
    //
    // Joint wage estimation
    //
    
    if (lsl_heckm == 1) {
        // Selection equation
        if (bsel > 0) {
            Select = cross(lsl_SelectVars', Bsel')                               // Selection prediction
            SelRes = lsl_Y :* ((lsl_Hwage :< .) :- Select)                        //   Selection residuals
            Lambda = normalden(Select) :/ normal(Select)    //   Heckman lambda
        } else Lambda = J(nobs, 0, 0)
        
        // Predict log-wages
        Hwage = cross((lsl_HeckmVars, Lambda)', (Bheck, Brho :* Bsig)')
        
        if (bsel == 0) {
            Hwres = lsl_Y :* (log(lsl_Hwage) :- Hwage)                  // Residuals
            Bsig  = sqrt(cross(Hwres, Hwres) / (lsl_groups - bheck))   // RMSE
            for (c = 1; c <= cols(Bsig); c++) {                      // Save RMSE
                st_numscalar("r(sigma_w" + strofreal(c) + ")", Bsig[1,c])
            }
        }
        
        // Predict wages
        Hwage = exp(Hwage :+ Bsig^2/2)
        
        //printf("SelSig = " + strofreal(SelSig) + ", Brho = " + strofreal(B[1,b - 1]) + "/" + strofreal(Brho) + ", Bsig = " + strofreal(B[1,b]) + "/" + strofreal(Bsig) + "\n");
    } else {
        Hwage = lsl_Hwage
        Bsig  = lsl_Sigma
    }
    // Round wage rates?
    if (lsl_round == 1 & cols(Hwage) > 0) Hwage = round(Hwage, 0.01)
    
    
    /* Loop over households */
    
    for (n = 1; n <= lsl_groups; n++) {
        // Last observation of group n
        c   = lsl_J[i]
        e   = i + c - 1
        Yn  = lsl_Y[|i\e|]
        
        // Fetch needed right hand side parts
        if (lsl_wagep == 1 | lsl_heckm == 1) {
            C    =  lsl_C[|i\e|]             // Get consumption from data
            CX   = (cols(lsl_CX)   > 0 ?   lsl_CX[|i\e|]     : J(c, 0, 0))
            C2X  = (cols(lsl_C2X)  > 0 ?  lsl_C2X[|i,1\e,.|] : J(c, 0, 0))
            L1   = lsl_L1[|i\e|]
            LX1  = (cols(lsl_LX1)  > 0 ?  lsl_LX1[|i,1\e,.|] : J(c, 0, 0))
            L2X1 = (cols(lsl_L2X1) > 0 ? lsl_L2X1[|i,1\e,.|] : J(c, 0, 0))
            L2   = (cols(lsl_L2)   > 0 ?   lsl_L2[|i\e|]     : J(c, 0, 0))
            LX2  = (cols(lsl_LX2)  > 0 ?  lsl_LX2[|i,1\e,.|] : J(c, 0, 0))
            L2X2 = (cols(lsl_L2X2) > 0 ? lsl_L2X2[|i,1\e,.|] : J(c, 0, 0))
            Xind = (cols(lsl_Xind) > 0 ? lsl_Xind[|i,1\e,.|] : J(c, 0, 0))
            Wn   = Hwage[|i,1\e,.|]
        } else Xnr = lsl_X[|i,1\e,.|]

        // Sum over draws
        lsum  = 0
        Gsum  = J(1, b, 0)
        H1sum = J(1, b, 0)
        H2sum = J(b, b, 0)

        // Run by random draw
        for (r = 1; r <= lsl_draws; r++) {
            // Init
            iRV  = lsl_draws * (n - 1) + r       // Indicates the active Halton sequence
            
            // Build (random?) coefficients matrix
            if (brnd > 0) Beta = Bfix :+ (brnd > 0 ? Brnd[iRV,.] : 0)


            /* Integrate out wage prediction error */

            if (lsl_wagep == 1 | lsl_heckm == 1) {
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
                TaxregX1 = (Mwage[.,1], Mwage[.,1]:^2, Mwage[.,1] :* lsl_TaxregIas1[|i,1\e,.|], Mwage[.,1]:^2 :* lsl_TaxregIas1[|i,1\e,.|])
                if (nlei == 2) TaxregX2 = (Mwage[.,2], Mwage[.,2]:^2, Mwage[.,2] :* lsl_TaxregIas2[|i,1\e,.|], Mwage[.,2]:^2 :* lsl_TaxregIas2[|i,1\e,.|])
                else           TaxregX2 = J(c, 0, 0)
                TaxregX = (TaxregX1, TaxregX2, lsl_TaxregVars[|i,1\e,.|], J(c, 1, 1))

                // Predict disposable income (can't be negative!)
                C = rowmax((cross(TaxregX', lsl_TaxregB'), J(c, 1, 1)))

                // Build matrix with independent variables
                if      (lsl_ufunc == "tran") Xnr = (log(C)  :* (CX,  log(C)  :* C2X,   log(L1),   log(L2)),
                                                     log(L1) :* (LX1, log(L1) :* L2X1),
                                                     log(L2) :* (LX2, log(L2) :* L2X2), log(L1) :* log(L2), Xind)
                else if (lsl_ufunc == "quad") Xnr = (C  :* (CX,  C  :* C2X,   L1,   L2),
                                                     L1 :* (LX1, L1 :* L2X1),
                                                     L2 :* (LX2, L2 :* L2X2), L1 :* L2, Xind)
            }


            /* Calculate utility levels */
            
            // Calculate choice probabilities
            Unr = cross(Xnr', Beta')                                    // Utility (choices in rows, draws in columns)
            Enr = exp(Unr :+ colmin(-mean(Unr) \ 700 :- colmax(Unr)))   // Standardize to avoid missings
            Pnr = Enr :/ colsum(Enr)                                    // Probabilities
            
            // Simplify
            pni  = cross(Yn, Pnr)   // Probability that choice is chosen
            YmPn = Yn :- Pnr        // Choice minus probabilities
            PXn  = cross(Pnr, Xnr)  // Right hand side cross by probs
            YXn  = cross(Yn, Xnr)   // Right hand side cross by choice


            /* Add to sum over draws */

            // Add to likelihood
            lsum = lsum + pni

            // Calculate gradient vector
            if (todo >= 1) {
                // Calculate gradient of systematic utility
                if (lsl_ufunc == "tran" | lsl_ufunc == "quad") DUdB = Xnr
                else {  // Box-Cox utility
                    DUdlC  = lsl_boxcox_g(lsl_C[|i\e|],  lC)  :* cross((lsl_CX[|i,1\e,.|], BcL1[|i\e|], (nlei == 2 ? BcL2[|i\e|] : J(c, 0, 0)))', B[|1\ncons|]')
                    DUdlL1 = lsl_boxcox_g(lsl_L1[|i\e|], lL1) :* cross((BcC[|i\e|], lsl_LX1[|i,1\e,.|], (nlei == 2 ? BcL2[|i\e|] : J(c, 0, 0)))',
                                                                       (B[ncx + 1], B[|ncons + 1\ncons + nlx1|],
                                                                        (nlei == 2 ? B[ncons + nlx1 + nlx2 + 1] : J(1, 0, 0)))')
                    DUdlL2 = (nlei == 2 ? lsl_boxcox_g(lsl_L2[|i\e|], lL2) :* cross((BcC[|i\e|], lsl_LX2[|i,1\e,.|], BcL1[|i\e|])',
                                                                                    (B[ncons], B[|ncons + nlx1 + 1\ncons + nlx1 + nlx2 + 1|])') : J(c, 0, 0))
                    DUdB   = Xnr, DUdlC, DUdlL1, DUdlL2
                }
                
                // Random components
                DUdBr = J(c, 0, 0)
                if (brnd > 0) {
                    if (lsl_corr == 1) {
                        for (rv = 1; rv <= rvars; rv++) {
                            DUdBr = (DUdBr, DUdB[.,lsl_Rvars[|rv\rvars|]] :* lsl_R[iRV,rv])
                        }
                    } else DUdBr = DUdB[.,lsl_Rvars] :* lsl_R[iRV,.]
                }
                
                // Heckman?
                DUdBw   = J(c, 0, 0)
                DUdBs   = J(c, 0, 0)
                DUdBrho = J(c, 0, 0)
                DUdBsig = J(c, 0, 0)
                if (lsl_heckm == 1) {
                    if      (lsl_ufunc == "quad") DUdC = cross((CX, 2 :* C2X :* C, L1, L2)', Beta[|1\ncons|]')
                    else if (lsl_ufunc == "tran") DUdC = cross(((CX, 2 :* C2X :* log(C), L1, L2) :/ C)', Beta[|1\ncons|]')
                    DCdM = cross((J(c, 1, 1), 2 :* Mwage, lsl_TaxregIas1[|i,1\e,.|], 2 :* Mwage :* lsl_TaxregIas1[|i,1\e,.|])', lsl_TaxregB[|1\2 + 2 * cols(lsl_TaxregIas1)|]')
                    DMdH = (lsl_Days[|i,1\e,1|] :/ 12 :/ 7) :* lsl_Hours[|i,1\e,.|]
                    if (bsel > 0) {
                        DWdBw   = Wn :* lsl_HeckmVars[|i,1\e,.|]
                        DWdBs   = - Wn :* Bsig :* Brho :* normalden(Select[|i,1\e,.|]) :* lsl_SelectVars[|i,1\e,.|] :*
                                   (normal(Select[|i,1\e,.|]) :* Select[|i,1\e,.|] :+ normalden(Select[|i,1\e,.|])) :/ normal(Select[|i,1\e,.|]):^2
                        DWdBrho = Wn :*  Bsig :* Lambda[|i,1\e,.|]
                        DWdBsig = Wn :* (Brho :* Lambda[|i,1\e,.|] :+ Bsig)
                        Gnr     = (Gnr, DUdC :* DCdM :* DMdH :* (DWdBw, DWdBs, DWdBrho, DWdBsig))
                    } else {
                        DWdBw = Wn :* (lsl_HeckmVars[|i,1\e,.|] :- colsum(lsl_HeckmVars :* Hwres) / (lsl_groups - bheck))
                        DUdBw = DUdC :* DCdM :* DMdH :* DWdBw
                        Gnr   = (Gnr, DUdBw)
                    }
                }

                // Total
                Gsum = Gsum + pni :* cross(YmPn, (DUdB, DUdBr, DUdBw, DUdBs, DUdBrho, DUdBsig))
            }

            // Calculate Hessian matrix
            if (todo == 2) {
                if (lsl_ufunc == "tran" | lsl_ufunc == "quad") YmPn_D2UdB2 = 0
                else {
                    YmPn_D2UdB2 = J(bfix, bfix, 0)
                    YmPn_D2UdB2[|ilC,1\ilC,ncons|] = cross(YmPn :* lsl_boxcox_g(lsl_C[|i\e|], lC), (lsl_CX[|i,1\e,.|], BcL1[|i\e|], (nlei == 2 ? BcL2[|i\e|] : J(c, 0, 0))))
                    YmPn_D2UdB2[ilC,ilC] = cross(YmPn :* lsl_boxcox_h(lsl_C[|i\e|], lC), cross((lsl_CX[|i,1\e,.|], BcL1[|i\e|], (nlei == 2 ? BcL2[|i\e|] : J(c, 0, 0)))', B[|1\ncons|]'))
                    
                    YmPn_D2UdB2[ilL1,ncx + 1] = cross(YmPn :* lsl_boxcox_g(lsl_L1[|i\e|], lL1), BcC[|i\e|])
                    YmPn_D2UdB2[|ilL1,ncons + 1\ilL1,ncons + nlx1|] = cross(YmPn :* lsl_boxcox_g(lsl_L1[|i\e|], lL1), lsl_LX1[|i,1\e,.|])
                    YmPn_D2UdB2[ilL1,ilC] = cross(YmPn :* lsl_boxcox_g(lsl_L1[|i\e|], lL1), B[ncons - nlei + 1] :* lsl_boxcox_g(lsl_C[|i\e|], lC))
                    YmPn_D2UdB2[ilL1,ilL1] = cross(YmPn :* lsl_boxcox_h(lsl_L1[|i\e|], lL1), cross((BcC[|i\e|], lsl_LX1[|i,1\e,.|], (nlei == 2 ? BcL2[|i\e|] : J(c, 0, 0)))',
                                                                                                   (B[ncx + 1], B[|ncons + 1\ncons + nlx1|],
                                                                                                   (nlei == 2 ? B[ncons + nlx1 + nlx2 + 1] : J(1, 0, 0)))'))
                    if (nlei == 2) {
                        YmPn_D2UdB2[ilL1,ncons + nlx1 + nlx2 + 1] = cross(YmPn :* lsl_boxcox_g(lsl_L1[|i\e|], lL1), BcL2[|i\e|])
                        YmPn_D2UdB2[ilL2,ncons + nlx1 + nlx2 + 1] = cross(YmPn :* lsl_boxcox_g(lsl_L2[|i\e|], lL2), BcL1[|i\e|])
                        YmPn_D2UdB2[ilL2,ncons] = cross(YmPn :* lsl_boxcox_g(lsl_L2[|i\e|], lL2), BcC[|i\e|])
                        YmPn_D2UdB2[|ilL2,ncons + nlx1 + 1\ilL2,ncons + nlx1 + nlx2|] = cross(YmPn :* lsl_boxcox_g(lsl_L2[|i\e|], lL2), lsl_LX2[|i,1\e,.|])
                        YmPn_D2UdB2[ilL2,ilC]  = cross(YmPn :* lsl_boxcox_g(lsl_L2[|i\e|], lL2), B[ncons] :* lsl_boxcox_g(lsl_C[|i\e|], lC))
                        YmPn_D2UdB2[ilL2,ilL1] = cross(YmPn :* lsl_boxcox_g(lsl_L2[|i\e|], lL2), B[ncons + nlx1 + nlx2 + 1] :* lsl_boxcox_g(lsl_L1[|i\e|], lL1))
                        YmPn_D2UdB2[ilL2,ilL2] = cross(YmPn :* lsl_boxcox_h(lsl_L2[|i\e|], lL2), cross((BcC[|i\e|], lsl_LX2[|i,1\e,.|], BcL1[|i\e|])',
                                                                                                       (B[ncons], B[|ncons + nlx1 + 1\ncons + nlx1 + nlx2 + 1|])'))
                    }
                    YmPn_D2UdB2 = makesymmetric(YmPn_D2UdB2)
                }
                
                // Utility
                H1 = pni :* (cross(Yn, DUdB) :- cross(Pnr, DUdB))
                H2 = pni :* (cross(cross(Yn, DUdB) :- cross(Pnr, DUdB), cross(YmPn, DUdB)) :-
                             cross(Pnr :* DUdB, DUdB :- cross(Pnr, DUdB)) :+ YmPn_D2UdB2)

                // Random components
                S1   = J(1, 0, 0)
                S2xx = (brnd > 0 ? J((lsl_corr == 1 ? brnd : 0), (lsl_corr == 1 ? brnd : 1), 0) : J(0, 0, 0))
                S2xy = J(0, bfix, 0)
                if (brnd > 0) {
                    if (lsl_corr == 1) {
                        iCol = 1
                        for (rv = 1; rv <= rvars; rv++) {
                            nCols = rvars - rv + 1
                            iRow = iCol
                            for (rv2 = rv; rv2 <= rvars; rv2++) {
                                nRows = rvars - rv2 + 1
                                S1   = S1,     pni :*  (YXn[.,lsl_Rvars[rv2,1]] - PXn[.,lsl_Rvars[rv2,1]])  * lsl_R[iRV,rv]
                                S2xy = S2xy \  pni :* (cross(YXn[.,lsl_Rvars[rv2,1]]                  - PXn[.,lsl_Rvars[rv2,1]],                 cross(YmPn, Xnr)) -
                                                       cross(Xnr[.,lsl_Rvars[rv2,1]]                 :- PXn[.,lsl_Rvars[rv2,1]], (Pnr :* Xnr))) * lsl_R[iRV,rv]
                                Svar =         pni :* (cross(YXn[.,lsl_Rvars[|rv2,1\rv2+nRows-1,1|]] :- PXn[.,lsl_Rvars[|rv2,1\rv2+nRows-1,1|]], cross(YmPn, Xnr[.,lsl_Rvars[|rv,1\rv+nCols-1,1|]])) -
                                                       cross(Xnr[.,lsl_Rvars[|rv2,1\rv2+nRows-1,1|]] :- PXn[.,lsl_Rvars[|rv2,1\rv2+nRows-1,1|]],      Pnr :* Xnr[.,lsl_Rvars[|rv,1\rv+nCols-1,1|]])) *
                                                    lsl_R[iRV,rv] * lsl_R[iRV,rv2]
                                S2xx[|iRow,iCol\iRow+nRows-1,iCol+nCols-1|] = Svar
                                if (iRow != iCol) S2xx[|iCol,iRow\iCol+nCols-1,iRow+nRows-1|] = Svar'
                                iRow = iRow + nRows
                            }
                            iCol = iCol + nCols
                        }
                    } else {
                        for (rv = 1; rv <= rvars; rv++) {
                            S1   = S1,        pni :*       (YXn[.,lsl_Rvars[rv,1]] - PXn[.,lsl_Rvars[rv,1]]) * lsl_R[iRV,rv]
                            S2xy = S2xy \     pni :* (cross(YXn[.,lsl_Rvars[rv,1]] - PXn[.,lsl_Rvars[rv,1]], cross(YmPn, Xnr))                    - cross(Xnr[.,lsl_Rvars[rv,1]] :- PXn[.,lsl_Rvars[rv,1]], Pnr :* Xnr)) * lsl_R[iRV,rv]
                            S2xx = S2xx \     pni :* (cross(YXn[.,lsl_Rvars[rv,1]] - PXn[.,lsl_Rvars[rv,1]], cross(YmPn, Xnr[.,lsl_Rvars[rv,1]]))  - cross(Xnr[.,lsl_Rvars[rv,1]] :- PXn[.,lsl_Rvars[rv,1]], Pnr :* Xnr[.,lsl_Rvars[rv,1]])) * lsl_R[iRV,rv]:^2
                            for (rv2 = rv + 1; rv2 <= rvars; rv2++) {
                                S2xx = S2xx \ pni :* (cross(YXn[.,lsl_Rvars[rv,1]] - PXn[.,lsl_Rvars[rv,1]], cross(YmPn, Xnr[.,lsl_Rvars[rv2,1]])) - cross(Xnr[.,lsl_Rvars[rv,1]] :- PXn[.,lsl_Rvars[rv,1]], Pnr :* Xnr[.,lsl_Rvars[rv2,1]])) * lsl_R[iRV,rv] * lsl_R[iRV,rv2]
                            }
                        }
                        S2xx = invvech(S2xx)
                    }

                    H2sum = H2sum + (H2, S2xy' \ S2xy, S2xx)
                }

                // Heckman
                W1   = J(1, 0, 0)
                W2xy = J(0, bfix, 0)
                W2xx = J(0, 0, 0)
                if (lsl_heckm == 1) {
                    DXdH   = (CX, C2X :* 2 :* C, L1, J(c, cols(Xnr) - (cols(CX) + cols(C2X) + 2 + nlei), 0)) :*
                             ((J(c, 1, 1), 2 :* Mwage, lsl_TaxregIas1[|i,1\e,.|], 2 :* Mwage :* lsl_TaxregIas1[|i,1\e,.|]) * lsl_TaxregB[|1,1\1,2 + 2 * cols(lsl_TaxregIas1)|]') :*
                             (lsl_Days[|i,1\e,1|] :/ 12 :/ 7) :* lsl_Hours[|i,1\e,.|]
                    D2UdC2 = 2 * Beta[|1,cols(CX) + 2\1,cols(CX) + cols(C2X) + 2|]
                    D2CdM2 = cross((J(c, 1, 2), 2 :* lsl_TaxregIas1[|i,1\e,.|])', (lsl_TaxregB[1,2], lsl_TaxregB[|1,2 + cols(lsl_TaxregIas1) + 1\1,2 + 2 * cols(lsl_TaxregIas1)|])')
                    D2MdH2 = 0

                    YmPnD2UdBw2 = cross(YmPn :* DWdBw :* (DMdH:^2 :* (D2UdC2 :* DCdM:^2 :+ DUdC :* D2CdM2) + DUdC :* DCdM :* D2MdH2), DWdBw)
                    for (hv = 1; hv <= bheck; hv++) {
                        YmPnD2UdBw2[hv,.] = YmPnD2UdBw2[hv,.] + cross(YmPn :* DUdC :* DCdM :* DMdH, DWdBw :* (lsl_HeckmVars[|i,hv\e,hv|] :- colsum(lsl_HeckmVars[.,hv] :* Hwres) / (lsl_groups - bheck)) :+
                                                                                                    cross(Wn', cross(lsl_Y :* lsl_HeckmVars[.,hv], lsl_HeckmVars) :/ (lsl_groups - bheck)))
                    }

                    W1   =   pni :* (cross(Yn, DUdBw) - cross(Pnr, DUdBw))
                    W2xy =   pni :* (cross(cross(Yn, DUdBw) - cross(Pnr, DUdBw), cross(YmPn, Xnr)) -
                                     cross(DUdBw :- cross(Pnr, DUdBw), Pnr :* Xnr) +
                                     cross(YmPn :* DXdH, /*Wn :* lsl_HeckmVars[|i,1\e,.|]*/DWdBw)')
                    W2xx =   pni :* (cross(cross(Yn, DUdBw) - cross(Pnr, DUdBw), cross(YmPn, DUdBw)) -
                                     cross(DUdBw :- cross(Pnr, DUdBw), Pnr :* DUdBw) + YmPnD2UdBw2)

                    H2sum = H2sum + (H2, W2xy' \ W2xy, W2xx)
                }

                // Total
                H1sum = H1sum + (H1, S1, W1)
                if (brnd == 0 & lsl_heckm == 0) H2sum = H2sum + H2
            }

        }

        // Prevent likelihood from becoming exactly zero
        lsum = max((lsum, 1e-25))

        // Add to overall statistics
        lnf = lnf + lsl_Weight[i] * log(lsum / lsl_draws)
        if (todo >= 1) G = G + lsl_Weight[i] * (lsum > 1e-25 ? Gsum / lsum : J(1, b, 0))
        if (todo == 2) H = H + lsl_Weight[i] * (lsum > 1e-25 ? H2sum / lsum - cross(Gsum, H1sum) / lsum^2: J(b, b, 0))

        // Next household
        i = i + c
    }
    //lnf
}

real matrix lsl_boxcox(real matrix Var, real scalar lam) {
    return (reldif(lam, 0) >= 1e-25 ? (Var:^lam :- 1) :/ lam : log(Var))
}

real matrix lsl_boxcox_g(real matrix Var, real scalar lam) {
    return (reldif(lam, 0) >= 1e-25 ? (Var:^lam :* (lam :* log(Var) :- 1) :+ 1) :/ (lam)^2 : 0.5 :* log(Var):^2)
}

real matrix lsl_boxcox_h(real matrix Var, real scalar lam) {
    return (reldif(lam, 0) >= 1e-25 ? (Var:^lam :* (lam:^2 :* log(Var):^2 :- 2 :* lam :* log(Var) :+ 2) :- 2) :/ (lam)^3 : (1/3) :* log(Var):^3)
}

mata set matastrict off
end

***
