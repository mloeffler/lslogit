/** 
 * LSLOGIT - Mata source code (functions are compiled to llslogit.mlib)
 * 
 * @package lslogit
 */


// Initialize
version 13.0
mata:
mata clear
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
                               : J(lsl_wagep * nlei, rvars, 0)) // DEBUG 2016-02-14: Is this correct?

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


// Finish library
mata set matastrict off
mata mlib create llslogit, dir(.) replace
mata mlib add llslogit lslogit_d2() lslogit_p() lsl_boxcox*()
end


***
