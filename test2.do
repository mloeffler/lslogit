di as error _n(5) "###############################################################################"
di as error       "###############################################################################"
di as error       "###############################################################################" _n(5)

local y choice
local group hhnrakt
local typlist sg_m sg_f co_m co_f co_v
local draws 10

local x_sg_m age_m age2_m handc_m
local i_sg_m D_pt?_m
local x_sg_f age_f age2_f handc_f care child2 child3_6 child7_16
local i_sg_f D_pt?_f
local x_co_m `x_sg_m'
local i_co_m `i_sg_m'
local x_co_f `x_sg_f'
local i_co_f `i_sg_f'
local x_co_v `x_sg_m' `x_sg_f'
local i_co_v `i_sg_m' `i_sg_f'

foreach ufunc in tran {
    if ("`ufunc'" == "tran") local ln l
    else                     local ln
    di as error _n(3) "###########################   UFUNC = `ufunc'   ###########################"
    foreach typ of local typlist {
        if ("`typ'" == "co_v") local end m f
        else                   local end = substr("`typ'", 4, 1)
        di as text _n(2) "                ###########     TYP = `typ'   ###########"
        
        
        //
        // Conditional Logit
        //
        
        //clogit `y' ${log_indepvars_`typ'} ${log_dummies_`typ'} if ${cond_`typ'} & !gino, group(`group')

        local l
        local rL
        foreach ending of local end {
            local l `l' freiz_`ending'
            if ("`ending'" == "m")    local lxvars : word count `x_sg_m'
            if ("`ending'" == "f")    local lxvars : word count `x_sg_f'
            local rL `rL' `lxvars'
        }
        if ("`end'" == "m")    local lx l1x(`x_sg_m')
        if ("`end'" == "f")    local lx l1x(`x_sg_f')
        if ("`typ'" == "co_v") local lx l1x(`x_sg_m') l2x(`x_sg_f')
        //lslogit `y' if ${cond_`typ'} & !gino, group(`group') ufunc(`ufunc') c(consum) cx(`x_`typ'') l(`l') `lx' ind(`i_`typ'') search(off) iterate(20)
        
        
        //
        // Random Coefficients
        //
        
        local hecksigma
        local wagep
        local hwage
        local vars " ${log_indepvars_`typ'} ${log_dummies_`typ'} "
        local vars : subinstr local vars " lconsum " " ", all
        local rand lconsum
        foreach ending of local end {
            local vars : subinstr local vars " lfreiz_`ending' " " ", all
            local rand `rand' lfreiz_`ending'
            local hecksigma `hecksigma' ${log_heckman_sigma_`ending'}
            local wagep `wagep' wagep_`ending'
            local hwage `hwage' h_wage_`ending'
        }
        //mixlogit `y' `vars' if ${cond_`typ'} & !gino, group(`group') nrep(`draws') rand(`rand')
        
        local rvC : word count `x_`typ''
        local rvC = `rvC' + 2
        local rvL
        local tria
        foreach r of local rL {
            local x1 : word 1 of `rL'
            local x2 : word 2 of `rL'
            if ("`x2'" == "`r'") {
                local x2 = 1 + `x1' + 2
                local tria `tria' tria2(east married child2 child3_6 child7_13 child14_17 child_num_tot)
            }
            else {
                local x2 = 0
                local tria `tria' tria1(east married child2 child3_6 child7_13 child14_17 child_num_tot)
            }
            local r = `rvC' + `r' + `x2' + 2
            local rvL `rvL' `r'
        }
        //lslogit `y' if ${cond_`typ'} & !gino, group(`group') ufunc(`ufunc') c(consum) cx(`x_`typ'') l(`l') `lx' ind(`i_`typ'') rand(`rvC' `rvL') iterate(20) dr(`draws') corr
        
        
        //
        // Integrating Out Wage Prediction Errors
        //
        
        lslogit `y' if ${cond_`typ'} & !gino [fw=BSweight$BS], group(`group') ufunc(`ufunc')     ///
                                                                  c(consum) cx(`x_`typ'')           ///
                                                                  l(`l') `lx' ind(`i_`typ'')        ///
                                                                  rand(`rvC' `rvL') corr            ///
                                                                  dr(`draws')                       ///
                                                                  taxreg(taxreg_`typ')              ///
                                                                  `tria'                            ///
                                                                  hecksigma(`hecksigma')            ///
                                                                  wagep(`wagep') hwage(`hwage') round iterate(30)
        
        lslogit `y' if ${cond_`typ'} & !gino [fw=BSweight$BS], group(`group') ufunc(`ufunc')     ///
                                                                  c(consum) cx(`x_`typ'')           ///
                                                                  l(`l') `lx' ind(`i_`typ'')        ///
                                                                  rand(`rvC' `rvL') corr            ///
                                                                  dr(`draws')                       ///
                                                                  taxreg(taxreg_`typ')              ///
                                                                  `tria'                            ///
                                                                  hecksigma(`hecksigma')            ///
                                                                  wagep(`wagep') hwage(`hwage') iterate(30)
    }
}

***