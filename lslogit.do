/*****************************************************************************
 *
 * izaML_clogitI: CONDITIONAL LOGIT MODEL,
 *                INTEGRATING WAGE PREDICTION ERRORS OUT
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
    syntax [, Level(integer `c(level)')]
    local dudes = 1 - e(dudes)
    ml display, level(`level') /*diparm(ln_consum, f(`dudes') d(0) label("% dU/dc>=0"))*/
end

cap program drop lslogit_Estimate
/**
 * Conditional Logit but integrating out wage prediction errors (Wrapper programm)
 * 
 * @param `group'  Group identifier variable
 * @param `taxreg' Stored estimates of the tax regression
 */
program define lslogit_Estimate, eclass
    syntax varname(numeric) [if] [in] [fweight/], GRoup(varname numeric) TYP(name) Ufunc(name)                          ///
                                                  Consumption(varname numeric) Leisure(varlist numeric min=1 max=2)          ///
                                                  [cx(varlist numeric) l1x(varlist numeric) l2x(varlist numeric) INDeps(varlist)          ///
                                                   /*Wages(varlist numeric min=1 max=2) Hours(varlist numeric min=1 max=2)*/ ///
                                                   difficult trace search(name) iterate(integer 100)                         ///
                                                   Level(integer `c(level)') gradient hessian debug burn(integer 15)        ///
                                                   DRaws(integer 50) Verbose WAGEPred TAXReg(name) RANDvars(string) method(name)]
    
    /* INITIALIZE ESTIMATOR
     */
    
    // Mark the estimation sample
    marksample touse
    markout `touse' `varlist' `group' `consumption' `leisure' `cx' `l1x' `l2x' `indeps'
    
    // Validate Maximum Likelihood method
    if ("`method'" == "") local method d2
    if (!inlist("`method'", "d0", "d1", "d2")) {
        di in r "method must be either 'd0', 'd1' or 'd2'"
        exit 498
    }
    
    // Validate utility function
    if (!inlist("`ufunc'", "quad", "tran")) {
        di in r "utility function must be either 'quad' or 'tran'"
        exit 498
    }
    
    // If translog, set up pre-text and check for zeros
    if ("`ufunc'" == "tran") {
        local ln  "ln"
        local pre "`ln'_"
        qui count if log(`consumption') == . & `touse'
        if (r(N) > 0) {
            di in r "consumption variable has to be positive with translog utility"
            exit 498
        }
    }
    
    // Validate Wage Prediction Options
    if ("`taxreg'" != "") local wagepred = "wagepred"
    if ("`wagepred'" == "" & "`randvars'" == "") local draws = 1
    else {
        if ("`taxreg'" != "") {
            tempname taxreg_from
            // Load tax regression estimates
            qui est restore `taxreg'
            mat `taxreg_from' = e(b)
            local taxreg_betas : colnames `taxreg_from'
            local taxreg_const
            local taxreg_wXias
            foreach beta of local taxreg_betas {
                if (strpos("`beta'", "m_wage")  == 0 & "`beta'" != "_cons")   local taxreg_const `taxreg_const' `beta'
                else if (strpos("`beta'", "Xm_wage") != 0) {
                    local beta = substr("`beta'", 1, strpos("`beta'", "Xm_wage") - 1)
                    local taxreg_wXias `taxreg_wXias' `beta'
                }
            }
        }
        /*
        else {
            di in r "option taxreg() required"
            exit 198
        }
        */
    }
    
    // Build weight settings
    if ("`weight'" != "")   local wgt "[`weight'=`exp']"
    
    // Get variable count
    local n_leisure  : word count `leisure'
    local n_cxias    : word count `cx'
    local n_l1xias   : word count `l1x'
    local n_l2xias   : word count `l2x'
    local n_indeps   : word count `indeps'
    local n_randvars : word count `randvars'
    
    // Adjust gender
    local genderc = cond(inlist("`typ'", "sg_m", "sg_f"), 1, 2)
    if (`genderc' == 2)   local genders   m f
    else                  local genders = substr("`typ'", -1, 1)
    
    // Select random variables
    local max_rvars = `n_randvars' + `genderc' * ("`wagepred'" != "")
    forval i = 1/`max_rvars' {
        local rvars `rvars' r`i'_1-r`i'_`draws'
    }
    
    
    /* LOOK FOR INITIAL VALUES
     */
    
    if ("`search'" != "off") {
        // Verbose mode on?
        if ("`verbose'" == "") local qui qui
        else                   di as text "Looking for initial values..."
        
        // Set up consumption and leisure
        tempvar c l1 l2
        if ("`ufunc'" == "tran") {
            qui gen `c' = log(`consumption') if `touse'
            foreach var of local leisure {
                if (strpos("`leisure'", "`var'") == 1) local lei l1
                else                                   local lei l2
                qui gen ``lei'' = log(`var') if `touse'
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
            if ("`ia'" != "c") local cxl `c'
            else               local cxl
            foreach var in ``ia'x' ``ia'' 0 `cxl' {
                if ("`var'" != "0") {
                    tempvar `ia'X`var'
                    qui gen ``ia'X`var'' = ``ia'' * `var' if `touse'
                    local initrhs `initrhs' ``ia'X`var''
                }
                else local initrhs `initrhs' ``ia''
            }
        }
        // Leisure cross term
        if (`n_leisure' == 2) {
            tempvar l1Xl2
            qui gen `l1Xl2' = `l1' * `l2' if `touse'
            local initrhs `initrhs' `l1Xl2'
        }
        // Add independent variables to var list
        local initrhs `initrhs' `indeps'
        
        // Estimate
        tempname init_from
        `qui' clogit `varlist' `initrhs' if `touse' `wgt', group(`group') iterate(25)
        if (e(converged) == 1) {
            // Save results
            mat `init_from' = e(b)
            local nobs      = e(N)
            local k         = e(k)
            local ll        = e(ll)
            local initopt     init(`init_from', copy) obs(`nobs') lf0(`k' `ll')
            // Update sample
            qui replace `touse' = e(sample)
        }
    }
    else {
        qui count if `touse'
        local nobs = r(N)
    }
    
    
    /* PREPARING DATA
     */
    
    if ("`verbose'" != "")   di as text "Preparing data..."
    
    // Drop missing data
    preserve
    qui keep if `touse'
    sort `group' //`leisure'
    
    // Setup data
    mata: ml_ufunc  = st_local("ufunc")                      // Utility function
    mata: ml_typ    = st_local("typ")                        // Labor supply type
    mata: ml_Weight = st_data(., st_local("exp"))            // Weight
    mata: ml_Y      = st_data(., st_local("varlist"))        // Left hand side
    
    // Set up right hand side
    mata: ml_C  = st_data(., st_local("consumption"))                                           // Consumption
    mata: ml_CX = (`n_cxias' > 0 ? st_data(., tokens(st_local("cx"))) : J(`nobs', 0, 0))         // Interactions
    forval i = 1/2 {                                                                            // Right hand side independent of consumption
        local var : word `i' of `leisure'
        mata: ml_L`i'  = ("`var'"     != "" ? st_data(., st_local("var")) : J(`nobs', 0, 0))
        mata: ml_L`i'X = (`n_l`i'xias' >  0 ? st_data(., tokens(st_local("l`i'x"))) : J(`nobs', 0, 0))
    }
    mata: ml_Xind = (`n_indeps' > 0 ? st_data(., tokens(st_local("indeps"))) : J(`nobs', 0, 0)) // Dummy variables
    if ("`ufunc'" == "tran") {                                                                  // Right hand side (translog)
        mata: ml_X = (ml_CX :*log(ml_C),  log(ml_C) :^2, log(ml_C),                         ///
                      ml_L1X:*log(ml_L1), log(ml_L1):^2, log(ml_L1), log(ml_L1):*log(ml_C), ///
                      ml_L2X:*log(ml_L2), log(ml_L2):^2, log(ml_L2), log(ml_L2):*log(ml_C), ///
                                                                     log(ml_L1):*log(ml_L2), ml_Xind)
    }
    else if ("`ufunc'" == "quad") {                                                             // Right hand side (quad)
        mata: ml_X = (ml_CX :*ml_C,  ml_C :^2, ml_C,                ///
                      ml_L1X:*ml_L1, ml_L1:^2, ml_L1, ml_L1:*ml_C,  ///
                      ml_L2X:*ml_L2, ml_L2:^2, ml_L2, ml_L2:*ml_C,  ///
                                                      ml_L1:*ml_L2, ml_Xind)
    }
    else mata: ml_X = J(`nobs', 0, 0)
    
    mata: ml_wagep = ("`wagepred'" != "")                   // Run Wage Prediction

    mata: ml_Days  = st_data(., ("days"))                   // Days of taxyear
    mata: ml_Hwage = J(`nobs', `genderc', 0)                // Hourly wage rates
    mata: ml_Mwage = J(`nobs', `genderc', 0)                // Monthly wages
    mata: ml_Hours = J(`nobs', `genderc', 0)                // Hypothetical hours
    mata: ml_Sigma = J(1, `genderc', 0)                     // Estimated variance of Heckman correction
    
    // Add data by gender
    forval i = 1/`genderc' {
        local s : word `i' of `genders'
        mata: ml_Hours[.,`i'] = st_data(., tokens("hyp_hour_`s'"))              // Hypothetical hours
        mata: ml_Hwage[.,`i'] = st_data(., tokens("h_wage_`s'"))                // Hourly wage rates
        mata: ml_Mwage[.,`i'] = st_data(., tokens("m_wage_`s'"))                // Monthly wages
        mata: ml_Sigma[1,`i'] = strtoreal(st_global("log_heckman_sigma_`s'"))   // Estimated variance of Heckman correction
    }
    
    // Tax regression
    if ("`wagepred'" != "" & "`taxreg'" != "") {
        mata: ml_TaxregB     = st_matrix("`taxreg_from'")                      // Tax regression estimates
        mata: ml_TaxregIas   = st_data(., tokens(st_local("taxreg_wXias")))    // Interaction variables on m_wage
        mata: ml_TaxregConst = st_data(., tokens(st_local("taxreg_const")))    // Variables that are independent of m_wage
    }
    
    // Fetch individual id
    qui duplicates report `group'
    mata: ml_groups = st_numscalar("r(unique_value)")   // Number of groups

    // Random draws
    mata: ml_draws = strtoreal(st_local("draws"))                   // Number of draws
    mata: ml_burn  = strtoreal(st_local("burn"))                    // Number of draws to burn
    mata: ml_Rvars = ("`randvars'" != "" ? strtoreal(tokens(st_local("randvars")))' : J(0, 0, 0))   // Random coefficients
    /*mata: ml_R     = st_data(., tokens(st_local("rvars")))          // Random draws*/
    mata: ml_R     = (`max_rvars' > 0 ? invnormal(halton(ml_groups*ml_draws, `max_rvars', 1+`burn')) : J(`nobs', 0, 0))     // Random draws
    
    // Number of choices per group
    tempvar choices
    by `group': gen `choices' = _N
    mata: ml_J = st_data(., st_local("choices"))
    
    // Restore data
    restore
    
    
    /* RUN ESTIMATION
     */
    
    if ("`verbose'" != "")   di as text "Run estimation..."
    
    // Set up equations
    /*
    local eq_consum (`pre'`consumption':  `varlist' = `cx' `consumption')
    local eq_leisure
    foreach var of local leisure {
        local i = 1 + (strpos("`leisure'", "`var'") > 1)
        local eq_leisure `eq_leisure' (`pre'`var': `l`i'x' `var') /`ln'CX`ln'L`i'
    }
    if (`n_leisure' == 2)  local eq_leisure  `eq_leisure' /`ln'L1X`ln'L2
    if (`n_indeps'  > 0)   local eq_indepvar (indeps: `indeps', noconst)
    */
    local eq_consum (C:  `varlist' = `cx' `consumption')
    local eq_leisure
    foreach var of local leisure {
        local i = 1 + (strpos("`leisure'", "`var'") > 1)
        local eq_leisure `eq_leisure' (L`i': `l`i'x' `var') /`ln'CX`ln'L`i'
    }
    if (`n_leisure' == 2)  local eq_leisure  `eq_leisure' /L1XL2
    if (`n_indeps'  > 0)   local eq_indepvar (IND: `indeps', noconst)
    forval i = 1/`n_randvars' {
        local eq_indepvar `eq_indepvar' /SD`i'
        if ("`initopt'" != "") mat `init_from' = (`init_from', 0.1)
    }
    
    // Estimate
    ml model `method'`debug' lslogit_d2() `eq_consum' `eq_leisure' `eq_indepvar' ///
            if `touse' `wgt', group(`group') `initopt' search(off) iterate(`iterate') nopreserve max `difficult' `trace' `gradient' `hessian'
    
    // Calculate marginal utility of consumption
    local dudes = rnormal(1,1)
    
    // Save results
    ereturn local  title    "Mixed Logit Labor Supply Model"
    ereturn local  cmd      "lslogit"
    ereturn local  predict  "izamodP"
    ereturn local  depvar    `varlist'
    ereturn local  group     `group'
    ereturn local  ufunc    "`ufunc'"
    ereturn scalar draws   = `draws'
    ereturn scalar dudes   = `dudes'
    
    // Show results
    lslogit_Replay, level(`level')
end

cap mata mata drop lslogit_d2()
mata:
/**
 * Standard Conditional Logit but integrating out wage prediction errors (Evaluator)
 * 
 * @param B_s Stata matrix of coefficients
 */
void lslogit_d2(transmorphic scalar M, real scalar todo, real rowvector B,
                real scalar lnf, real rowvector G, real matrix H) {
    external ml_wagep
    external ml_ufunc
    external ml_groups
    external ml_draws
    external ml_burn
    external ml_Sigma
    external ml_Weight
    external ml_typ
    external ml_R
    external ml_Rvars
    external ml_J
    external ml_TaxregB         // Coefficients of tax regression
    external ml_TaxregConst     // Wage independent variables of tax regression
    external ml_TaxregIas       // Wage interaction variables of tax regression
    external ml_Days
    external ml_Hwage
    external ml_Mwage
    external ml_Hours
    
    external ml_C
    external ml_CX
    external ml_L1
    external ml_L1X
    external ml_L2
    external ml_L2X
    external ml_Xind
    
    external ml_X
    external ml_Y               // Get dependent variable
    
    
    /* Setup */
    
    /*
    // If not all household members need error integration, which do?
    c_integr = ((ml_typ == "sg_m" | ml_typ == "sg_f") ? 1 : 2)
    if (c_integr == 1) {
        integr_a = 1
        integr_b = 1
    } else {
        integr_a = ((ml_typ == "co_m" | ml_typ == "co_v") ? 1 : 2)
        integr_b = ((ml_typ == "co_f" | ml_typ == "co_v") ? 2 : 1)
    }

    c_taxregIas = cols(ml_TaxregIas)
    */

    // Definitions
    i   = 1         // Indicates first observation of active group
    iRV = 1         // Indicates next random variable to use (column of ml_R)
    
    // Build coefficient vector
    Bfix = B[|1,1\1,cols(B)-rows(ml_Rvars)|]                            // Get fixed coefficients
    Brnd = (rows(ml_Rvars) > 0 ? B[|1,cols(Bfix)+1\1,.|] : J(0, 0, 0))  // Get auxiliary random coefficients
    //Sigm = diag(exp(Brnd))                                              // Build variance-(covariance) matrix
    Sigm = diag(Brnd)                                              // Build variance-(covariance) matrix
    
    // Build matrix with random coefficients (mean zero), every row is a draw
    if (rows(ml_Rvars) > 0) {
        Brnd = J(rows(ml_R), cols(Bfix), 0)
        for (rv = 1; rv <= rows(ml_Rvars); rv++) {
            Brnd[.,ml_Rvars[rv,1]] = ml_R[.,iRV] :* Sigm[iRV,iRV]
            iRV = iRV + 1
        }
    }
    // From now on: Beta[rows=ml_R,cols=Bfix] = Bfix :+ Brnd
    
    // Maximum Likelihood Parameter
    lnf  = 0                                    // Log-likelihood
    G    = J(1, cols(B), 0)                     // Gradient
    H    = J(cols(B), cols(B), 0)               // Hessian matrix
    
    
    /* Loop over households */

    for (n = 1; n <= ml_groups; n++) {
        // Last observation of group n
        c   = ml_J[i,1]
        e   = i + c - 1
        Yn  = ml_Y[|i,1\e,1|]
        Xnr = ml_X[|i,1\e,.|]
        
        /*
        Mwage = ml_Mwage[|i,1\e,.|]
        C     =  ml_C[|i,1\e,1|]   // Get consumption from data
        CX    = (cols(ml_CX)   > 0 ?   ml_CX[|i,1\e,.|] : J(c, 0, 0))
        L1    = ml_L1[|i,1\e,1|]
        L1X   = (cols(ml_L1X)  > 0 ?  ml_L1X[|i,1\e,.|] : J(c, 0, 0))
        L2    = (cols(ml_L2)   > 0 ?   ml_L2[|i,1\e,1|] : J(c, 0, 0))
        L2X   = (cols(ml_L2X)  > 0 ?  ml_L2X[|i,1\e,.|] : J(c, 0, 0))
        Xind  = (cols(ml_Xind) > 0 ? ml_Xind[|i,1\e,.|] : J(c, 0, 0))
        */
        
        // Sum over draws
        lsum = 0
        Gsum = J(1, cols(B), 0)
        if (ml_draws == 1) Hsum = J(cols(B), cols(B), 0)
        else {
            H1sum = J(1, cols(B), 0)
            H2sum = J(cols(B), cols(B), 0)
        }
        
        // Run by random draw
        for (r = 1; r <= ml_draws; r++) {
            // Init
            iRV = 0     // Indicates the active random variable (starting at 0!!)
            
            // If wage prediction enabled
            if (ml_wagep == 1) {
                printf("wagep\n");
                /* Calculate monthly earnings for draw r */
                
                // Adjust wages with random draws
                for (s = integr_a; s <= integr_b; s++) {
                    Mwage[.,s] = round(((ml_Days[|i,1\e,1|] :/ 12) :/ 7) :* ml_Hours[|i,s\e,s|] :* ml_Hwage[|i,s\e,s|] :* exp(ml_Sigma[1,s] * ml_R[|i,ml_draws*iRV+r\e,ml_draws*iRV+r|]))
                    iRV = iRV + 1
                }
                // Monthly wages to the power of...
                Mwage2 = Mwage:^2 :/ 100^2
                Mwage3 = Mwage:^3 :/ 100^3
                Mwage4 = Mwage:^4 :/ 100^4
                Mwage5 = Mwage:^5 :/ 100^5
                
                printf("blaaaa!!!!\n");
                
                
                /* Predict disposable income */
                
                // Container with independent variables for dpi prediction
                TaxregX = (J(c, 12*c_integr, 0), ml_TaxregConst[|i,1\e,.|], J(c, 1, 1))

                // Fill matrix of independent variables for dpi prediction
                for (s = 1; s <= c_integr; s++) {
                    TaxregX[|1,12*(s-1)+1\c,12*s|] = (Mwage[.,s], Mwage2[.,s], Mwage3[.,s], Mwage4[.,s], Mwage5[.,s], Mwage[.,s] :* ml_TaxregIas[|i,1+(s-1)*c_taxregIas/c_integr\e,s*c_taxregIas/c_integr|])
                }

                // Predict disposable income (can't be negative!)
                C = rowmax(((TaxregX * ml_TaxregB'), J(c, 1, 0.01)))
                
                // Build matrix with independent variables
                if      (ml_ufunc == "tran") Xnr = (CX:*log(C), log(C):^2, log(C), L1X:*log(L1), log(L1):^2, log(L1), log(L1):*log(C), L2X:*log(L2), log(L2):^2, log(L2), log(L2):*log(C), log(L1):*log(L2), Xind)
                else if (ml_ufunc == "quad") Xnr = (CX:*C, C:^2, C, L1X:*L1, L1:^2, L1, L1:*C, L2X:*L2, L2:^2, L2, L2:*C, L1:*L2, Xind)
            }
            
            
            /* Calculate utility levels */
            
            // Build (random?) coefficients matrix
            Beta = Bfix :+ (rows(ml_Rvars) > 0 ? Brnd[ml_draws * (n-1) + r,.] : 0)
            
            // Calculate choice probabilities
            Unr = Xnr * Beta'                                           // Utility (choices in rows, draws in columns)
            Enr = exp(Unr :+ colmin(-mean(Unr) \ 700 :- colmax(Unr)))   // Standardize to avoid missings
            Pnr = Enr :/ colsum(Enr)                                    // Probabilities


            /* Add to sum over draws */
            
            // Add to likelihood
            lsum = lsum + Yn' * Pnr

            // Calculate gradient vector
            if (todo >= 1) {
                Gnr = colsum(Yn' * Pnr * (Yn - Pnr) :* Xnr)
                for (rv = 1; rv <= rows(ml_Rvars); rv++) {
                    Gnr = (Gnr, (Pnr' * Yn :* (Yn' :- Pnr') * Xnr[.,ml_Rvars[rv,1]])' * ml_R[ml_draws * (n-1) + r,rv])
                }
                Gsum = Gsum + Gnr
            }
            
            // Calculate Hessian matrix
            if (todo == 2) {
                if (ml_draws == 1) Hsum = Hsum - (Pnr :* Xnr)' * (Xnr :- Pnr' * Xnr)
                else {
                    H1nr = - (Yn' * Pnr) :*  (Yn' * Xnr - Pnr' * Xnr)
                    H2nr =   (Yn' * Pnr) :* ((Yn' * Xnr - Pnr' * Xnr)' * ((Yn - Pnr)' * Xnr) - (Pnr :* Xnr)' * (Xnr :- Pnr' * Xnr))
                    S1xy = J(1, 0, 0)
                    S2xx = J(0, 1, 0)
                    S2xy = J(0, cols(H2nr), 0)
                    for (rv = 1; rv <= rows(ml_Rvars); rv++) {
                        S1xy = S1xy, - (Yn' * Pnr) :*  (Yn' * Xnr[.,ml_Rvars[rv,1]] - Pnr' * Xnr[.,ml_Rvars[rv,1]]) * ml_R[ml_draws * (n-1) + r,rv]
                        S2xy = S2xy \  (Yn' * Pnr) :* ((Yn' * Xnr[.,ml_Rvars[rv,1]] - Pnr' * Xnr[.,ml_Rvars[rv,1]])' * ((Yn - Pnr)' * Xnr) - (Xnr[.,ml_Rvars[rv,1]] :- Pnr' * Xnr[.,ml_Rvars[rv,1]])' * (Pnr :* Xnr)) * ml_R[ml_draws * (n-1) + r,rv]
                        S2xx = S2xx \  (Yn' * Pnr) :* ((Yn' * Xnr[.,ml_Rvars[rv,1]] - Pnr' * Xnr[.,ml_Rvars[rv,1]])' * ((Yn - Pnr)' * Xnr[.,ml_Rvars[rv,1]]) - (Pnr :* Xnr[.,ml_Rvars[rv,1]])' * (Xnr[.,ml_Rvars[rv,1]] :- Pnr' * Xnr[.,ml_Rvars[rv,1]])) * ml_R[ml_draws * (n-1) + r,rv]:^2
                        for (rv2 = rv + 1; rv2 <= rows(ml_Rvars); rv2++) {
                            S2xx = S2xx \ (Yn' * Pnr) :* ((Yn' * Xnr[.,ml_Rvars[rv,1]] - Pnr' * Xnr[.,ml_Rvars[rv,1]])' * ((Yn - Pnr)' * Xnr[.,ml_Rvars[rv2,1]]) - (Pnr :* Xnr[.,ml_Rvars[rv2,1]])' * (Xnr[.,ml_Rvars[rv,1]] :- Pnr' * Xnr[.,ml_Rvars[rv,1]])) * ml_R[ml_draws * (n-1) + r,rv] * ml_R[ml_draws * (n-1) + r,rv2]
                        }
                    }
                    H1sum = H1sum + (H1nr, S1xy)
                    H2sum = H2sum + (H2nr, S2xy' \ S2xy, invvech(S2xx))
                }
            }
        }

        // Prevent likelihood from becoming exactly zero
        lsum = max((lsum, 1e-25))
        
        // Add to overall statistics
        lnf = lnf + ml_Weight[i,1] :* log(lsum / ml_draws)
        if (todo >= 1) G = G + ml_Weight[i,1] :* (lsum > 1e-25 ? Gsum / lsum : J(1, cols(G), 0))
        if (todo == 2) H = H + ml_Weight[i,1] :* (lsum > 1e-25 ? (ml_draws == 1 ? Hsum : Gsum' * H1sum / lsum^2 + H2sum / lsum) : J(rows(H), cols(H), 0))
        
        // Next household
        i = i + c
    }
}
end

***
