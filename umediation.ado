*! umediation 1.0.0 23oct2026
*! Huijie SUN, University of Science and Technology of China
*! shj0512@mail.ustc.edu.cn | https://github.com/shj0512
*! SPDX-License-Identifier: MIT
*! Requires Stata 17 or later.
* MIT License
*
* Copyright (c) 2026 Huijie SUN
*
* Permission is hereby granted, free of charge, to any person obtaining a copy
* of this software and associated documentation files (the "Software"), to deal
* in the Software without restriction, including without limitation the rights
* to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
* copies of the Software, and to permit persons to whom the Software is
* furnished to do so, subject to the following conditions:
*
* The above copyright notice and this permission notice shall be included in all
* copies or substantial portions of the Software.
*
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
* IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
* FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
* AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
* LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
* OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
* SOFTWARE.

program define umediation, eclass
    version 17.0
    if replay() {
        if "`e(cmd)'" != "umediation" error 301
        syntax
        _umediation_display
        exit
    }

    syntax varname(numeric) [if] [in], X(varname numeric) ///
        MEDiator(varname numeric) [CONTrols(varlist numeric fv) ///
        ABSorb(varlist) CLuster(varname) AT(numlist sort) ///
        ATPERcentiles(numlist sort) RANGE(numlist min=2 max=2 sort) ///
        RANGEPERcentiles(numlist min=2 max=2 sort) ///
        DIRection(string) REPS(integer 1999) SEED(string) ///
        Level(cilevel) SAVing(string) REPlace NODOTS ///
        INFERence(string) WILDREPS(integer 0) WILDWeight(string) ///
        WILDVARiant(string) WILDSAVing(string)]

    local y `varlist'
    local inference = lower(strtrim("`inference'"))
    if "`inference'" == "" local inference pairs
    if !inlist("`inference'","pairs","jackknife") {
        di as error "inference() must be pairs or jackknife"
        exit 198
    }
    local wildweight = lower(strtrim("`wildweight'"))
    if "`wildweight'" == "" local wildweight webb
    local wildvariant = lower(strtrim("`wildvariant'"))
    if "`wildvariant'" == "" local wildvariant score
    if !inlist("`wildweight'","webb","rademacher") | !inlist("`wildvariant'","classic","score") {
        di as error "wildweight() must be webb or rademacher; wildvariant() must be classic or score"
        exit 198
    }
    if `wildreps' < 0 | (`wildreps' > 0 & `wildreps' < 99) {
        di as error "wildreps() must be zero (disabled) or at least 99"
        exit 198
    }
    if ("`inference'" == "jackknife" | `wildreps' > 0) & "`cluster'" == "" {
        di as error "jackknife and wild inference require cluster()"
        exit 198
    }
    if `wildreps' == 0 & `"`wildsaving'"' != "" {
        di as error "wildsaving() requires wildreps()"
        exit 198
    }
    if "`y'" == "`x'" | "`y'" == "`mediator'" | "`x'" == "`mediator'" {
        di as error "outcome, x(), and mediator() must be different variables"
        exit 198
    }
    if `reps' < 50 & "`inference'" == "pairs" {
        di as error "reps() must be at least 50; use 1999 or more for reported intervals"
        exit 198
    }
    if "`at'" != "" & "`atpercentiles'" != "" {
        di as error "specify at() or atpercentiles(), not both"
        exit 198
    }
    if "`range'" != "" & "`rangepercentiles'" != "" {
        di as error "specify range() or rangepercentiles(), not both"
        exit 198
    }
    local direction = lower(strtrim("`direction'"))
    if "`direction'" == "" local direction either
    if !inlist("`direction'", "posneg", "negpos", "either") {
        di as error "direction() must be posneg, negpos, or either"
        exit 198
    }
    if "`at'`atpercentiles'" == "" local atpercentiles 5 25 50 75 95
    if "`range'`rangepercentiles'" == "" local rangepercentiles 5 95
    foreach plist in atpercentiles rangepercentiles {
        foreach p of local `plist' {
            if `p' <= 0 | `p' >= 100 {
                di as error "percentile values must be strictly between 0 and 100"
                exit 198
            }
        }
    }
    local kpoints : word count `at' `atpercentiles'
    if `kpoints' > 100 {
        di as error "at most 100 evaluation positions are allowed"
        exit 198
    }
    if "`replace'" != "" & `"`saving'`wildsaving'"' == "" {
        di as error "replace requires saving() or wildsaving()"
        exit 198
    }
    if `"`saving'"' != "" {
        if lower(substr(`"`saving'"', -4, .)) != ".dta" local saving `"`saving'.dta"'
        if "`replace'" == "" confirm new file `"`saving'"'
    }

    if `"`wildsaving'"' != "" {
        if lower(substr(`"`wildsaving'"',-4,.)) != ".dta" local wildsaving `"`wildsaving'.dta"'
        if "`replace'" == "" confirm new file `"`wildsaving'"'
        local same_output 0
        if `"`saving'"' != "" {
            local output_base `"`c(pwd)'"'
            local output_os "`c(os)'"
            mata: st_local("same_output", strofreal(umediation_same_output( ///
                st_local("saving"), st_local("wildsaving"), ///
                st_local("output_base"), st_local("output_os"))))
        }
        if `same_output' {
            di as error "saving() and wildsaving() must use different files"
            exit 198
        }
    }
    local rawcontrols
    if "`controls'" != "" {
        quietly fvrevar `controls', list
        local rawcontrols `r(varlist)'
        local forbidden `y' `x' `mediator'
        local overlap : list rawcontrols & forbidden
        if "`overlap'" != "" {
            di as error "controls() may not contain the outcome, x(), or mediator()"
            exit 198
        }
    }
    local forbidden `y' `x' `mediator'
    local overlap : list absorb & forbidden
    if "`overlap'" != "" {
        di as error "absorb() may not contain the outcome, x(), or mediator()"
        exit 198
    }
    marksample touse
    markout `touse' `x' `mediator' `rawcontrols'
    markout `touse' `absorb' `cluster', strok
    quietly count if `touse'
    if r(N) < 20 {
        di as error "fewer than 20 complete observations in the common estimation sample"
        exit 2001
    }

    local keepvars `y' `x' `mediator' `rawcontrols' `absorb' `cluster'
    local keepvars : list uniq keepvars
    tempname work draw oldest bootstore
    local oldframe `c(frame)'
    local oldrng `c(rngstate)'
    local held 0
    capture _estimates hold `oldest', restore nullok
    if !_rc local held 1
    quietly frame put `keepvars' if `touse', into(`work')

    local opts x(`x') mediator(`mediator') reps(`reps') level(`level') ///
        direction(`direction') work(`work') draw(`draw') matastore(`bootstore') ///
        inference(`inference') wildreps(`wildreps') wildweight(`wildweight') wildvariant(`wildvariant')
    foreach op in controls absorb cluster at atpercentiles range rangepercentiles {
        if "``op''" != "" local opts `opts' `op'(``op'')
    }
    if "`seed'" != "" local opts `opts' seed(`seed')
    if `"`saving'"' != "" local opts `opts' saving(`"`saving'"')
    if `"`wildsaving'"' != "" local opts `opts' wildsaving(`"`wildsaving'"')
    local opts `opts' `replace' `nodots' 

    capture noisily frame `work': _umediation_work `y', `opts'
    local rc = _rc
    capture mata: mata drop `bootstore'
    if !`rc' {
        tempname B V CI TAB GRID ATVAL RANGEVAL ROOT IUT DIAG COMP ROOTSET INFL SIZES WSTAT WREV
        matrix `B' = r(b)
        matrix `V' = r(V)
        matrix `CI' = r(ci)
        matrix `TAB' = r(table)
        matrix `GRID' = r(grid)
        matrix `ATVAL' = r(at)
        matrix `RANGEVAL' = r(range)
        matrix `ROOT' = r(root_info)
        matrix `IUT' = r(iut)
        matrix `DIAG' = r(diagnostics)
        matrix `SIZES' = r(cluster_diagnostics)
        if "`inference'" == "jackknife" {
            matrix `COMP' = r(component_ci)
            matrix `ROOTSET' = r(root_set)
            matrix `INFL' = r(influence)
            local rootsettype `r(root_set_type)'
        }
        if `wildreps' > 0 {
            matrix `WSTAT' = r(wild_stats)
            matrix `WREV' = r(wild_reversal)
            local wildused = r(wild_reps)
            local enumerated = r(wild_enumerated)
        }
        local used = r(reps_used)
        local N = r(N)
        local nc = r(N_clust)
        local nested `r(nested_absorb)'
        local classification `r(classification)'
        local rootstatus `r(root_status)'
    }
    capture frame change `oldframe'
    capture frame drop `draw'
    capture frame drop `work'
    capture quietly set rngstate `oldrng'
    if `rc' {
        if `held' capture _estimates unhold `oldest'
        exit `rc'
    }
    if `held' capture _estimates drop `oldest'

    ereturn post `B' `V', esample(`touse') obs(`N') depname(`y')
    ereturn matrix ci = `CI'
    ereturn matrix table = `TAB'
    ereturn matrix grid = `GRID'
    ereturn matrix at = `ATVAL'
    ereturn matrix range = `RANGEVAL'
    ereturn matrix root_info = `ROOT', copy
    ereturn matrix iut = `IUT'
    ereturn matrix diagnostics = `DIAG', copy
    ereturn matrix cluster_diagnostics = `SIZES'
    if "`inference'" == "jackknife" {
        ereturn matrix component_ci = `COMP'
        ereturn matrix root_set = `ROOTSET'
        ereturn matrix influence = `INFL'
        ereturn local root_set_type "`rootsettype'"
        ereturn scalar joint_components = 4
        ereturn scalar df_r = `nc'-1
    }
    if `wildreps' > 0 {
        ereturn matrix wild_stats = `WSTAT'
        ereturn matrix wild_reversal = `WREV'
        ereturn scalar wild_reps = `wildused'
        ereturn scalar wild_enumerated = `enumerated'
        ereturn local wildweight "`wildweight'"
        ereturn local wildvariant "`wildvariant'"
        ereturn local wildsaving `"`wildsaving'"'
    }
    ereturn local inference "`inference'"
    ereturn scalar reps_requested = `used'
    ereturn scalar reps_valid = `DIAG'[1,1]
    ereturn scalar reps_failed = `DIAG'[1,2]
    ereturn scalar N_clust = `nc'
    ereturn scalar level = `level'
    ereturn scalar root = `ROOT'[1,1]
    ereturn scalar root_inside = `ROOT'[1,6]
    ereturn scalar sup_critical = `DIAG'[1,3]
    ereturn scalar grid_points = 101
    ereturn local title "Umediation: quadratic first-path mediation"
    ereturn local release_status "first_release"
    ereturn local calibration_status "not_uniformly_calibrated"
    ereturn local depvar `y'
    ereturn local x `x'
    ereturn local mediator `mediator'
    ereturn local controls `controls'
    ereturn local absorb `absorb'
    ereturn local cluster `cluster'
    ereturn local nested_absorb `nested'
    ereturn local direction `direction'
    ereturn local classification `classification'
    ereturn local root_status `rootstatus'
    ereturn local seed `seed'
    ereturn local saving `"`saving'"'
    ereturn local quantile_method "linear_type7"
    ereturn local vce "bootstrap"
    ereturn local vcetype "Bootstrap"
    if "`inference'" == "jackknife" {
        ereturn local vce "clusterjackknife"
        ereturn local vcetype "Cluster jackknife CV3J"
        ereturn local quantile_method "bonferroni4_t_projection"
        ereturn scalar component_critical = `DIAG'[1,3]
        ereturn scalar sup_critical = .
    }
    ereturn local cmdline `"umediation `0'"'
    ereturn local version "1.0.0"
    ereturn local cmd "umediation"
    _umediation_display
end

capture program drop _umediation_work
program define _umediation_work, rclass
    version 17.0
    syntax varname, X(varname) MEDiator(varname) WORK(name) DRAW(name) MATAStore(name) ///
        REPS(integer) Level(real) DIRection(string) ///
        [CONTrols(string) ABSorb(varlist) CLuster(varname) ///
        AT(numlist) ATPERcentiles(numlist) RANGE(numlist) ///
        RANGEPERcentiles(numlist) SEED(string) SAVing(string) REPlace NODOTS ///
        INFERence(string) WILDREPS(integer 0) WILDWeight(string) WILDVARiant(string) WILDSAVing(string)]
    local y `varlist'
    local nobs = _N

    tempvar x2 clid newclid original_order
    generate long `original_order' = _n
    generate double `x2' = `x'^2
    quietly count if missing(`x2')
    if r(N) {
        di as error "x() squared overflows; rescale x() before estimation"
        exit 459
    }
    quietly summarize `x', meanonly
    local xmin = r(min)
    local xmax = r(max)
    if `xmin' >= `xmax' {
        di as error "x() has no variation in the common sample"
        exit 459
    }
    local nc .
    if "`cluster'" != "" {
        quietly egen long `clid' = group(`cluster')
        quietly summarize `clid', meanonly
        local nc = r(max)
        if `nc' < 2 {
            di as error "cluster() must identify at least two clusters"
            exit 459
        }
    }
    local fevars
    local nestedvars
    local nestednames
    foreach f of local absorb {
        tempvar fe
        quietly egen long `fe' = group(`f')
        local fevars `fevars' `fe'
        if "`cluster'" != "" {
            tempvar clmin clmax
            quietly bysort `fe': egen long `clmin' = min(`clid')
            quietly bysort `fe': egen long `clmax' = max(`clid')
            quietly count if `clmin' != `clmax'
            if !r(N) {
                local nestedvars `nestedvars' `fe'
                local nestednames `nestednames' `f'
            }
            drop `clmin' `clmax'
        }
    }

    quietly sort `original_order'
    tempname ATVAL RANGEVAL
    if "`atpercentiles'" != "" {
        quietly _pctile `x', percentiles(`atpercentiles')
        local nk : word count `atpercentiles'
        matrix `ATVAL' = J(1,`nk',.)
        forvalues j=1/`nk' {
            matrix `ATVAL'[1,`j'] = r(r`j')
        }
    }
    else {
        local nk : word count `at'
        matrix `ATVAL' = J(1,`nk',.)
        forvalues j=1/`nk' {
            matrix `ATVAL'[1,`j'] = real(word("`at'",`j'))
        }
    }
    if "`rangepercentiles'" != "" {
        quietly _pctile `x', percentiles(`rangepercentiles')
        local L = r(r1)
        local U = r(r2)
    }
    else {
        local L : word 1 of `range'
        local U : word 2 of `range'
    }
    if `L' >= `U' | missing(`L',`U') {
        di as error "the data range must have two distinct finite endpoints"
        exit 459
    }
    if `L' < `xmin' | `U' > `xmax' {
        di as error "range() must stay within the observed x() range"
        exit 459
    }
    forvalues j=1/`nk' {
        if `ATVAL'[1,`j'] < `xmin' | `ATVAL'[1,`j'] > `xmax' {
            di as error "evaluation positions must stay within the observed x() range"
            exit 459
        }
    }
    matrix `RANGEVAL' = (`L',`U')
    matrix colnames `RANGEVAL' = lower upper
    local cnames
    forvalues j=1/`nk' {
        local cnames `cnames' ie`j'
    }
    matrix colnames `ATVAL' = `cnames'
    local fitopts x(`x') x2(`x2') mediator(`mediator') ///
        positions(`ATVAL') lower(`L') upper(`U')
    if "`controls'" != "" local fitopts `fitopts' controls(`controls')
    if "`fevars'" != "" local fitopts `fitopts' absorb(`fevars')
    capture noisily _umediation_fit `y', `fitopts'
    if _rc {
        di as error "the original common-sample model is not identified"
        di as error "check x, x-squared, mediator variation, and collinearity with controls/fixed effects"
        exit 459
    }
    tempname POINT V CI TAB GRID ROOT IUT DIAG
    local BOOT `matastore'
    matrix `POINT' = r(b)
    local k = colsof(`POINT')
    if "`inference'" == "jackknife" {
        if `nc' < 3 {
            di as error "jackknife requires at least three clusters; this is a computational minimum"
            exit 459
        }
        local reps = `nc'
    }
    tempname COMP ROOTSET INFL SIZES WSTAT WREV
    matrix `SIZES' = J(1,5,.)
    if "`cluster'" != "" mata: st_matrix("`SIZES'",umediation_cluster_sizes(st_data(.,"`clid'")))
    matrix colnames `SIZES' = min_size median_size max_size size_effective_G largest_share
    mata: `BOOT' = J(`reps',`k'+3,.)
    if "`seed'" != "" {
        capture quietly set seed `seed'
        if _rc {
            capture mata: mata drop `BOOT'
            di as error "seed() is not a valid Stata random-number seed"
            exit 198
        }
    }
    if "`nodots'" == "" di as text "`inference' (`reps' fits requested): " _continue
    forvalues r=1/`reps' {
        quietly frame copy `work' `draw', replace
        if "`inference'" == "jackknife" {
            quietly frame `draw': drop if `clid' == `r'
        }
        else if "`cluster'" == "" {
            quietly frame `draw': bsample
        }
        else {
            quietly frame `draw': bsample, cluster(`clid') idcluster(`newclid')
            foreach f of local nestedvars {
                tempvar clone
                quietly frame `draw': egen long `clone' = group(`newclid' `f')
                quietly frame `draw': replace `f' = `clone'
                quietly frame `draw': drop `clone'
            }
        }
        capture frame `draw': _umediation_fit `y', `fitopts'
        local fitrc = _rc
        mata: `BOOT'[`r',`k'+3] = `fitrc'
        frame `draw': local drawN = _N
        mata: `BOOT'[`r',`k'+2] = `drawN'
        if !`fitrc' {
            mata: `BOOT'[`r',1..`k'] = st_matrix("r(b)")
            mata: `BOOT'[`r',`k'+1] = st_numscalar("r(root)")
        }
        if "`nodots'" == "" {
            if `fitrc' di as error "x" _continue
            else di as text "." _continue
            if mod(`r',50)==0 di as text " `r'" _continue
        }
    }
    if "`nodots'" == "" di ""
    local outnames `V' `CI' `TAB' `GRID' `ROOT' `IUT' `DIAG'
    if "`inference'" == "pairs" {
        capture noisily mata: umediation_summarize(`BOOT', st_matrix("`POINT'"), ///
            st_matrix("`ATVAL'"),`L',`U',`level',"`direction'",tokens("`outnames'"))
    }
    else {
        local jout `outnames' `COMP' `ROOTSET' `INFL'
        capture noisily mata: umediation_jksummarize(`BOOT', st_matrix("`POINT'"), ///
            st_matrix("`ATVAL'"),`L',`U',`level',"`direction'",tokens("`jout'"))
    }
    local sumrc = _rc
    if `sumrc' {
        capture mata: mata drop `BOOT'
        exit `sumrc'
    }
    local nvalid = `DIAG'[1,1]
    if ("`inference'" == "pairs" & `nvalid' < max(50,ceil(.8*`reps'))) | ("`inference'" == "jackknife" & `nvalid' != `nc') {
        capture mata: mata drop `BOOT'
        di as error "only `nvalid' of `reps' requested replications were valid"
        di as error "pairs requires 50 and 80% valid draws; jackknife requires EVERY deletion fit; no retries or exclusions"
        exit 498
    }
    if `"`saving'"' != "" {
        local statnames : colnames `POINT'
        quietly frame `draw': clear
        quietly frame `draw': set obs `reps'
        quietly frame `draw': generate long replication = _n
        foreach v of local statnames {
            quietly frame `draw': generate double `v' = .
        }
        foreach v in root draw_n return_code {
            quietly frame `draw': generate double `v' = .
        }
        frame `draw': mata: st_store(., tokens("`statnames' root draw_n return_code"), `BOOT')
        quietly frame `draw': generate byte valid = return_code == 0
        frame `draw': char _dta[umediation_version] "1.0.0"
        frame `draw': char _dta[release_status] "first_release"
        frame `draw': char _dta[inference] "`inference'"
        frame `draw': char _dta[replication_meaning] "pairs draw or integer cluster group deleted"
        frame `draw': char _dta[reps_requested] "`reps'"
        frame `draw': char _dta[reps_valid] "`nvalid'"
        frame `draw': char _dta[seed] "`seed'"
        frame `draw': char _dta[cluster] "`cluster'"
        frame `draw': char _dta[quantile_method] "linear_type7"
        capture frame `draw': save `"`saving'"', `replace'
        local saverc = _rc
        if `saverc' {
            capture mata: mata drop `BOOT'
            exit `saverc'
        }
    }
    capture mata: mata drop `BOOT'
    if `wildreps' > 0 {
        if "`seed'" != "" quietly set seed `seed'
        local wopts x(`x') x2(`x2') mediator(`mediator') clid(`clid') ///
            lower(`L') upper(`U') reps(`wildreps') weight(`wildweight') variant(`wildvariant') level(`level')
        if "`controls'" != "" local wopts `wopts' controls(`controls')
        if "`fevars'" != "" local wopts `wopts' absorb(`fevars')
        if "`nestedvars'" != "" local wopts `wopts' nested(`nestedvars')
        if `"`wildsaving'"' != "" local wopts `wopts' saving(`"`wildsaving'"')
        capture noisily _umediation_wild `y', `wopts' `replace'
        if _rc exit _rc
        matrix `WSTAT' = r(stats)
        matrix `WREV' = r(reversal)
        local wildused = r(reps)
        local enumerated = r(enumerated)
    }
    matrix colnames `V' = `: colnames `POINT''
    matrix rownames `V' = `: colnames `POINT''
    matrix colnames `CI' = `: colnames `POINT''
    matrix rownames `CI' = lower upper
    matrix rownames `TAB' = `: colnames `POINT''
    matrix colnames `TAB' = estimate boot_se ci_lower ci_upper
    if "`inference'" == "jackknife" matrix colnames `TAB' = estimate jackknife_se ci_lower ci_upper
    matrix colnames `GRID' = x ie boot_se point_lower point_upper sim_lower sim_upper
    if "`inference'" == "jackknife" matrix colnames `GRID' = x ie jackknife_se point_lower point_upper sim_lower sim_upper
    matrix colnames `ROOT' = estimate ci_lower ci_upper finite_reps finite_share inside denominator_clear path_clear
    matrix colnames `IUT' = left_bound right_bound alpha reject
    matrix rownames `IUT' = posneg negpos
    matrix colnames `DIAG' = valid failed sup_critical all_positive all_negative
    return scalar reps_used = `reps'
    return matrix cluster_diagnostics = `SIZES'
    if "`inference'" == "jackknife" {
        matrix colnames `COMP' = estimate jackknife_se ci_lower ci_upper
        matrix rownames `COMP' = slope_mid curvature_xm bpath kappa_d
        matrix colnames `ROOTSET' = lower upper lower_unbounded upper_unbounded
        matrix rownames `ROOTSET' = interval1 interval2
        matrix colnames `INFL' = cluster_group scaled_max_shift delta_a1 delta_a2 delta_b delta_kappa_m point_direction_changed
        return matrix component_ci = `COMP'
        return matrix root_set = `ROOTSET'
        return matrix influence = `INFL'
        return local root_set_type `rootsettype'
    }
    if `wildreps' > 0 {
        return matrix wild_stats = `WSTAT'
        return matrix wild_reversal = `WREV'
        return scalar wild_reps = `wildused'
        return scalar wild_enumerated = `enumerated'
    }
    return matrix b = `POINT'
    return matrix V = `V'
    return matrix ci = `CI'
    return matrix table = `TAB'
    return matrix grid = `GRID'
    return matrix at = `ATVAL'
    return matrix range = `RANGEVAL'
    return matrix root_info = `ROOT'
    return matrix iut = `IUT'
    return matrix diagnostics = `DIAG'
    return scalar N = `nobs'
    return scalar N_clust = `nc'
    return local nested_absorb `nestednames'
    return local classification `classification'
    return local root_status `rootstatus'
end

capture program drop _umediation_fit
program define _umediation_fit, rclass
    version 17.0
    syntax varname, X(varname) X2(varname) MEDiator(varname) ///
        POSitions(name) LOWER(real) UPPER(real) [CONTrols(string) ABSorb(varlist)]
    local y `varlist'
    local fedummies
    local firstfe : word 1 of `absorb'
    local nfe : word count `absorb'
    if `nfe' > 1 {
        forvalues j=2/`nfe' {
            local f : word `j' of `absorb'
            local fedummies `fedummies' i.`f'
        }
    }
    if "`firstfe'" == "" {
        quietly regress `mediator' `controls' `x' `x2'
    }
    else {
        quietly areg `mediator' `controls' `fedummies' `x' `x2', absorb(`firstfe')
    }
    if e(N) != _N | e(df_r) <= 0 error 459
    tempname A1 A2 BP C1 C2 OMIT MBC
    matrix `MBC' = e(b)
    quietly _ms_omit_info `MBC'
    matrix `OMIT' = r(omit)
    foreach v in `x' `x2' {
        local j = colnumb(`MBC', "`v'")
        if missing(`j') error 459
        if `OMIT'[1,`j'] error 459
    }
    scalar `A1' = _b[`x']
    scalar `A2' = _b[`x2']
    if "`firstfe'" == "" {
        quietly regress `y' `controls' `x' `x2' `mediator'
    }
    else {
        quietly areg `y' `controls' `fedummies' `x' `x2' `mediator', absorb(`firstfe')
    }
    if e(N) != _N | e(df_r) <= 0 error 459
    matrix `MBC' = e(b)
    quietly _ms_omit_info `MBC'
    matrix `OMIT' = r(omit)
    foreach v in `x' `x2' `mediator' {
        local j = colnumb(`MBC', "`v'")
        if missing(`j') error 459
        if `OMIT'[1,`j'] error 459
    }
    scalar `BP' = _b[`mediator']
    scalar `C1' = _b[`x']
    scalar `C2' = _b[`x2']
    tempname B
    matrix `B' = (`A1',`A2',`BP',`C1',`C2',`A1'*`BP', ///
        2*`A2'*`BP',2*`C2',2*`C2'+2*`A2'*`BP', ///
        (`A1'+2*`A2'*`lower')*`BP', ///
        (`A1'+2*`A2'*`upper')*`BP',2*`A2'*`BP'*(`upper'-`lower'))
    local cnames a1 a2 bpath c1 c2 eta kappa_m kappa_d kappa_t ie_l ie_u delta_ie
    local nk = colsof(`positions')
    forvalues j=1/`nk' {
        matrix `B' = `B', ((`A1'+2*`A2'*`positions'[1,`j'])*`BP')
        local cnames `cnames' ie`j'
    }
    matrix colnames `B' = `cnames'
    mata: st_local("badcore", strofreal(hasmissing(st_matrix("`B'"))))
    if `badcore' error 459
    return matrix b = `B'
    return scalar root = cond(`A2'==0 | `BP'==0,.,-`A1'/(2*`A2'))
end

capture program drop _umediation_display
program define _umediation_display
    version 17.0
    tempname T A R Q D
    matrix `T' = e(table)
    matrix `A' = e(at)
    matrix `R' = e(root_info)
    matrix `Q' = e(range)
    matrix `D' = e(diagnostics)
    di ""
    di as text "Umediation 1.0.0"
    di as text "Finite-sample interval coverage and reversal error rates depend on the design."
    if "`e(cluster)'" != "" & "`e(inference)'" == "pairs" {
        di as text "Few or imbalanced clusters can cause undercoverage and excess false reversals."
        di as text "See help umediation for validation results and inference choices."
    }
    if "`e(inference)'" == "jackknife" {
        di as text "CV3J t(G-1), four-component Bonferroni rectangle and projection."
        di as text "Projection may be conservative; small-cluster calibration is not guaranteed."
    }
    di as text "Outcome: " as result "`e(depvar)'" as text "   X: " as result "`e(x)'" ///
        as text "   Mediator: " as result "`e(mediator)'"
    di as text "Common sample: " as result %10.0fc e(N) ///
        as text "   Fits valid/requested: " as result e(reps_valid) "/" e(reps_requested)
    if "`e(absorb)'" != "" di as text "Fixed effects: " as result "`e(absorb)'"
    if "`e(cluster)'" != "" di as text "Resampling clusters: " as result "`e(cluster)' (" e(N_clust) ")"
    di as text "Data range: [" as result %10.5g `Q'[1,1] as text ", " ///
        as result %10.5g `Q'[1,2] as text "]   Direction: " as result "`e(direction)'"
    di ""
    di as text "Inference: " as result "`e(inference)'"
    if "`e(inference)'" == "pairs" di as text "Statistic          Estimate     Boot. SE       Nominal percentile " e(level) "% CI"
    else di as text "Statistic          Estimate     JK SE          Nominal projection " e(level) "% CI"
    di as text "------------------------------------------------------------------------"
    foreach name in a1 a2 bpath kappa_m kappa_d kappa_t ie_l ie_u delta_ie {
        local j = rownumb(`T',"`name'")
        di as text %-15s "`name'" as result %12.6g `T'[`j',1] " " ///
            %12.6g `T'[`j',2] "   [" %12.6g `T'[`j',3] ", " %12.6g `T'[`j',4] "]"
    }
    di ""
    di as text "Local effects at fixed evaluation positions"
    if "`e(inference)'" == "pairs" di as text "Statistic               X          IE(x)              Percentile CI"
    else di as text "Statistic               X          IE(x)              Projection CI"
    forvalues k=1/`=colsof(`A')' {
        local j = rownumb(`T',"ie`k'")
        di as text %-10s "ie`k'" as result %16.6g `A'[1,`k'] " " ///
            %14.6g `T'[`j',1] "   [" %12.6g `T'[`j',3] ", " %12.6g `T'[`j',4] "]"
    }
    di ""
    di as text "Nominal-rule classification: " as result "`e(classification)'"
    if "`e(classification)'" == "uncertain" {
        di as text "A reversal was not established; this does not establish absence of reversal."
    }
    di as text "Candidate IE zero: " as result %12.6g `R'[1,1] ///
        as text "   Status: " as result "`e(root_status)'"
    if `R'[1,4] >= 2 & "`e(inference)'" == "pairs" {
        di as text "Root percentile interval (diagnostic): [" as result %12.6g `R'[1,2] ///
            as text ", " as result %12.6g `R'[1,3] as text "]"
    }
    if `R'[1,7] == 0 di as text "Root denominator interval includes zero; location inference is unstable."
    if `R'[1,8] == 0 di as text "The mediator-to-outcome path interval includes zero."
    if e(reps_failed) > 0 di as text "Failed draws: " as result e(reps_failed) ///
        as text "; counted in requested draws and never retried."
    if e(reps_requested) < 999 & "`e(inference)'" == "pairs" di as text "Exploratory run: at least 1999 draws improve tail-quantile precision, not calibration."
    if "`e(inference)'" == "pairs" {
        di as text "Intervals: nominal linear-interpolated percentiles; 101-point centered max-standardized band."
        di as text "Root estimates are descriptive; reversal decisions use the endpoint tests."
    }
    else {
        di as text "Joint region covers slope_mid, 2*a2, b and 2*c2; c1 has a marginal t interval."
        di as text "Root confidence-set type: " as result "`e(root_set_type)'"
        di as text "See e(root_set); flags identify infinite endpoints. It concerns the first-path zero."
        di as text "A weak b path prevents interpreting a first-path zero as an established mediated reversal."
    }
    if "`e(cluster)'" != "" {
        tempname CS
        matrix `CS' = e(cluster_diagnostics)
        di as text "Cluster size min/median/max: " as result `CS'[1,1] "/" `CS'[1,2] "/" `CS'[1,3]
        di as text "Size-only effective count: " as result %8.3f `CS'[1,4] ///
            as text "; descriptive, not an information-based validity guarantee."
    }
    if "`e(wildvariant)'" != "" {
        tempname WR
        matrix `WR' = e(wild_reversal)
        di as text "Separate experimental WCR-`e(wildvariant)' reversal test (`e(wildweight)' weights):"
        di as text "  posneg p=" as result %7.4f `WR'[1,1] as text "  negpos p=" ///
            as result %7.4f `WR'[2,1] as text "  either p=" as result %7.4f `WR'[3,1]
        di as text "This test does not replace the primary classification or its confidence intervals."
    }
    di as text "Read help umediation for few-cluster, weak-path, boundary, and double-zero limitations."
end

capture mata: mata drop umediation_q()
capture mata: mata drop umediation_summarize()
capture mata: mata drop umediation_output_path()
capture mata: mata drop umediation_same_output()
mata:
string scalar umediation_output_path(string scalar path, string scalar base, string scalar os)
{
    string scalar part, prefix, resolved
    string rowvector parts
    real scalar cut, minimum
    if (os=="Windows") {
        path=subinstr(path,char(92),"/")
        base=subinstr(base,char(92),"/")
    }
    resolved=pathresolve(base,path)
    if (os=="Windows") resolved=subinstr(resolved,char(92),"/")
    if (os=="Windows" & substr(resolved,1,1)=="/" & substr(resolved,1,2)!="//" & substr(base,2,1)==":") resolved=substr(base,1,2)+resolved
    prefix=(substr(resolved,1,2)=="//" ? "//" : (substr(resolved,1,1)=="/" ? "/" : ""))
    minimum=(prefix=="//" ? 2 : (substr(resolved,2,1)==":" ? 1 : 0))
    parts=J(1,0,"")
    while (resolved!="") {
        cut=strpos(resolved,"/")
        if (cut==0) {
            part=resolved
            resolved=""
        }
        else {
            part=substr(resolved,1,cut-1)
            resolved=substr(resolved,cut+1,.)
        }
        if (part=="" | part==".") continue
        if (part=="..") {
            if (cols(parts)>minimum) {
                if (cols(parts)==1) parts=J(1,0,"")
                else parts=parts[1..(cols(parts)-1)]
            }
        }
        else parts=(parts,part)
    }
    resolved=prefix+invtokens(parts,"/")
    return(os=="Windows" ? strlower(resolved) : resolved)
}
real scalar umediation_same_output(string scalar first, string scalar second, string scalar base, string scalar os)
{
    return(umediation_output_path(first,base,os)==umediation_output_path(second,base,os))
}
real scalar umediation_q(real colvector x, real scalar p)
{
    real scalar h, j, g
    x = select(x, x :< .)
    if (!rows(x)) return(.)
    x = sort(x,1)
    h = 1 + (rows(x)-1)*p
    j = floor(h)
    g = h-j
    if (j >= rows(x)) return(x[rows(x)])
    return((1-g)*x[j] + g*x[j+1])
}

void umediation_summarize(real matrix D, real rowvector b, real rowvector at,
    real scalar L, real scalar U, real scalar level, string scalar direction,
    string rowvector out)
{
    real matrix Z, V, CI, T, G, I, R, DG, H, C
    real colvector valid, roots, grid, s, sup, delta
    real rowvector mu
    real scalar k, n, j, alpha, branch, c, rv, root, denomclear, pathclear
    real scalar posneg, negpos, positive, negative
    string scalar classification, rootstatus
    k = cols(b)
    valid = (D[,k+3] :== 0) :& (rowmissing(D[,1..k]) :== 0)
    Z = select(D[,1..k],valid)
    n = rows(Z)
    if (n < 2) _error(2001,"fewer than two valid bootstrap draws")
    mu = mean(Z)
    C = Z :- mu
    V = quadcross(C,C)/(n-1)
    alpha = 1-level/100
    CI = J(2,k,.)
    for (j=1;j<=k;j++) {
        CI[1,j] = umediation_q(Z[,j],alpha/2)
        CI[2,j] = umediation_q(Z[,j],1-alpha/2)
    }
    T = (b',sqrt(diagonal(V)),CI')
    grid = rangen(L,U,101)
    H = Z[,6]*J(1,101,1) + Z[,7]*grid'
    s = sqrt(diagonal(variance(H)))
    G = (grid, (b[6]:+b[7]*grid), s, J(101,4,.))
    for (j=1;j<=101;j++) {
        G[j,4] = umediation_q(H[,j],alpha/2)
        G[j,5] = umediation_q(H[,j],1-alpha/2)
    }
    sup = J(n,1,0)
    for (j=1;j<=101;j++) {
        if (s[j]>0 & s[j]<.) {
            delta = abs((H[,j]:-G[j,2])/s[j])
            sup = rowmax((sup,delta))
        }
    }
    c = umediation_q(sup,1-alpha)
    G[,6] = G[,2] - c*s
    G[,7] = G[,2] + c*s
    branch = direction=="either" ? alpha/2 : alpha
    I = J(2,4,.)
    posneg = 0
    negpos = 0
    if (direction=="posneg" | direction=="either") {
        I[1,] = (umediation_q(Z[,10],branch),umediation_q(Z[,11],1-branch),branch,0)
        posneg = I[1,1]>0 & I[1,2]<0
        I[1,4] = posneg
    }
    if (direction=="negpos" | direction=="either") {
        I[2,] = (umediation_q(Z[,10],1-branch),umediation_q(Z[,11],branch),branch,0)
        negpos = I[2,1]<0 & I[2,2]>0
        I[2,4] = negpos
    }
    positive = min(G[,6])>0
    negative = max(G[,7])<0
    classification = "uncertain"
    if (positive) classification = "support_same_sign_positive"
    if (negative) classification = "support_same_sign_negative"
    if (posneg) classification = "support_reversal_posneg"
    if (negpos) classification = "support_reversal_negpos"
    root = (b[2]!=0 & b[3]!=0) ? -b[1]/(2*b[2]) : .
    roots = select(D[,k+1],valid)
    rv = sum(roots:<.)
    denomclear = (CI[1,2]>0 | CI[2,2]<0)
    pathclear = (CI[1,3]>0 | CI[2,3]<0)
    R = (root,umediation_q(roots,alpha/2),umediation_q(roots,1-alpha/2),
         rv,rv/n,(root<. ? (root>L & root<U) : .),denomclear,pathclear)
    rootstatus = "outside_interval"
    if (root==.) rootstatus = "undefined"
    else if (!denomclear | !pathclear | rv<n) rootstatus = "unstable_or_weak_path"
    else if (root==L | root==U) rootstatus = "at_boundary"
    else if (root>L & root<U) rootstatus = "inside_interval"
    DG = (n,rows(D)-n,c,positive,negative)
    st_matrix(out[1],V)
    st_matrix(out[2],CI)
    st_matrix(out[3],T)
    st_matrix(out[4],G)
    st_matrix(out[5],R)
    st_matrix(out[6],I)
    st_matrix(out[7],DG)
    st_local("classification",classification)
    st_local("rootstatus",rootstatus)
}
end

capture program drop _umediation_wild
program define _umediation_wild, rclass
    version 17.0
    syntax varname, X(varname) X2(varname) MEDiator(varname) CLID(varname) ///
        LOWER(real) UPPER(real) REPS(integer) WEIGHT(string) VARIANT(string) Level(real) ///
        [CONTrols(string) ABSorb(varlist) NESTED(varlist) SAVing(string) REPlace]
    local y `varlist'
    quietly summarize `clid', meanonly
    local nc = r(max)
    tempname W S T O R Z saveframe
    mata: `W' = umediation_wild_weights(`nc',`reps',"`weight'")
    mata: st_local("actual",strofreal(rows(`W')))
    mata: st_local("enumerated",strofreal(umediation_wild_enumerated(`nc',`reps',"`weight'")))
    mata: `S' = J(3,6,.)
    mata: `T' = J(rows(`W'),3,.)
    local partial : word 1 of `nested'
    local opts x(`x') x2(`x2') mediator(`mediator') clid(`clid') ///
        lower(`lower') upper(`upper') weightmat(`W') variant(`variant') enumerated(`enumerated')
    if "`controls'" != "" local opts `opts' controls(`controls')
    if "`absorb'" != "" local opts `opts' absorb(`absorb')
    if "`partial'" != "" local opts `opts' partial(`partial')
    capture noisily _umediation_wildfit `mediator', `opts' equation(m) output(`O')
    local rc = _rc
    if !`rc' {
        mata: `S'[1..2,] = `O'[1..6,]'
        mata: `T'[,1..2] = `O'[7..rows(`O'),]
        capture noisily _umediation_wildfit `y', `opts' equation(y) output(`O')
        local rc = _rc
    }
    if !`rc' {
        mata: `S'[3,] = `O'[1..6,]'
        mata: `T'[,3] = `O'[7..rows(`O'),]
        mata: st_matrix("`Z'",`S')
        mata: st_matrix("`R'",umediation_wild_combine(`S',1-`level'/100))
    }
    if !`rc' & `"`saving'"' != "" {
        quietly frame create `saveframe'
        quietly frame `saveframe': set obs `actual'
        quietly frame `saveframe': generate long replication = _n
        local wn
        forvalues j=1/`nc' {
            quietly frame `saveframe': generate double w`j' = .
            local wn `wn' w`j'
        }
        foreach vv in t_left t_right t_b {
            quietly frame `saveframe': generate double `vv' = .
        }
        frame `saveframe': mata: st_store(.,tokens("`wn' t_left t_right t_b"),(`W',`T'))
        frame `saveframe': char _dta[umediation_version] "1.0.0"
        frame `saveframe': char _dta[wildvariant] "`variant'"
        frame `saveframe': char _dta[wildweight] "`weight'"
        frame `saveframe': char _dta[enumerated] "`enumerated'"
        frame `saveframe': char _dta[cluster_order] "egen group(cluster), sorted source values; w1 is group 1"
        capture frame `saveframe': save `"`saving'"', `replace'
        local rc = _rc
        capture frame drop `saveframe'
    }
    capture mata: mata drop `W' `S' `T' `O'
    if `rc' exit `rc'
    matrix rownames `Z' = left_slope right_slope bpath
    matrix colnames `Z' = estimate cluster_se t p_positive p_negative degenerate_draws
    matrix rownames `R' = posneg negpos either
    matrix colnames `R' = p reject
    return matrix stats = `Z'
    return matrix reversal = `R'
    return scalar reps = `actual'
    return scalar enumerated = `enumerated'
end

capture program drop _umediation_wildfit
program define _umediation_wildfit
    version 17.0
    syntax varname, X(varname) X2(varname) MEDiator(varname) CLID(varname) ///
        LOWER(real) UPPER(real) WEIGHTMAT(name) VARIANT(string) EQUATION(string) OUTPUT(name) ENUMERATED(integer) ///
        [CONTrols(string) ABSorb(varlist) PARTIAL(varname)]
    local y `varlist'
    local rhs `controls' `x' `x2'
    if "`equation'" == "y" local rhs `rhs' `mediator'
    foreach f of local absorb {
        if "`f'" != "`partial'" local rhs `rhs' i.`f'
    }
    if "`partial'" == "" quietly regress `y' `rhs'
    else quietly areg `y' `rhs', absorb(`partial')
    if e(N) != _N | e(df_r) <= 0 error 459
    tempname BB OMIT R
    matrix `BB' = e(b)
    quietly _ms_omit_info `BB'
    matrix `OMIT' = r(omit)
    local names : colnames `BB'
    local terms
    forvalues j=1/`=colsof(`BB')' {
        local term : word `j' of `names'
        if !`OMIT'[1,`j'] & "`term'" != "_cons" local terms `terms' `term'
    }
    * Explicit retained levels must not acquire a new base level in fvrevar.
    local literalterms
    foreach term of local terms {
        local parts : subinstr local term "#" " ", all
        local literal
        foreach part of local parts {
            if regexm("`part'","^([0-9]+)[.]") {
                local lev = regexs(1)
                local part : subinstr local part "`lev'." "`lev'bn.", all
            }
            if "`literal'" == "" local literal `part'
            else local literal `literal'#`part'
        }
        local literalterms `literalterms' `literal'
    }
    quietly fvrevar `literalterms'
    local rv `r(varlist)'
    local n : word count `rv'
    local cons = "`partial'" == ""
    local p = `n'+`cons'
    if `p' > 500 {
        di as error "experimental wild engine permits at most 500 identified design columns after one nested FE is absorbed"
        exit 908
    }
    local ix : list posof "`x'" in terms
    local ix2 : list posof "`x2'" in terms
    local im : list posof "`mediator'" in terms
    if `ix' == 0 | `ix2' == 0 error 459
    if "`equation'" == "m" {
        matrix `R' = J(`p',2,0)
        matrix `R'[`ix',1] = 1
        matrix `R'[`ix2',1] = 2*`lower'
        matrix `R'[`ix',2] = 1
        matrix `R'[`ix2',2] = 2*`upper'
    }
    else {
        if `im' == 0 error 459
        matrix `R' = J(`p',1,0)
        matrix `R'[`im',1] = 1
    }
    mata: `output' = umediation_wild_calc(st_data(.,tokens("`rv'")), ///
        st_data(.,"`y'"),st_data(.,"`clid'"),"`partial'", ///
        st_matrix("`R'"),`weightmat',"`variant'",`enumerated')
end

capture mata: mata drop umediation_interval_product()
capture mata: mata drop umediation_interval_linear()
capture mata: mata drop umediation_root_region()
capture mata: mata drop umediation_jksummarize()
capture mata: mata drop umediation_cluster_sizes()
capture mata: mata drop umediation_wild_enumerated()
capture mata: mata drop umediation_wild_weights()
capture mata: mata drop umediation_wild_calc()
capture mata: mata drop umediation_wild_combine()
mata:
real rowvector umediation_interval_product(real rowvector a, real rowvector b)
{
    real rowvector v
    v=(a[1]*b[1],a[1]*b[2],a[2]*b[1],a[2]*b[2])
    return((min(v),max(v)))
}
real rowvector umediation_interval_linear(real rowvector a, real rowvector b, real scalar c)
{
    return(a+(c>=0 ? c*b : c*b[(2,1)]))
}
real matrix umediation_root_region(real rowvector s, real rowvector k, real scalar mid)
{
    real matrix R
    real rowvector v
    real scalar j
    string scalar kind
    R=J(2,4,.)
    kind="undefined"
    if (k[1]>0 | k[2]<0) {
        v=(mid-s[1]/k[1],mid-s[1]/k[2],mid-s[2]/k[1],mid-s[2]/k[2])
        R[1,]=(min(v),max(v),0,0)
        kind="bounded"
    }
    else if (!(k[1]==0 & k[2]==0)) {
        if (s[1]<=0 & s[2]>=0) {
            if (s[1]==0 & s[2]==0) {
                R[1,]=(mid,mid,0,0)
                kind="bounded"
            }
            else if ((k[1]==0 & s[1]==0) | (k[2]==0 & s[2]==0)) {
                R[1,]=(.,mid,1,0)
                kind="unbounded"
            }
            else if ((k[1]==0 & s[2]==0) | (k[2]==0 & s[1]==0)) {
                R[1,]=(mid,.,0,1)
                kind="unbounded"
            }
            else {
                R[1,]=(.,.,1,1)
                kind="all_real"
            }
        }
        else {
            j=1
            if (k[1]<0) {
                if (s[1]>0) R[j,]=(mid-s[1]/k[1],.,0,1)
                else R[j,]=(.,mid-s[2]/k[1],1,0)
                j++
            }
            if (k[2]>0) {
                if (s[1]>0) R[j,]=(.,mid-s[1]/k[2],1,0)
                else R[j,]=(mid-s[2]/k[2],.,0,1)
                j++
            }
            kind=(j==3 ? "disjoint" : "unbounded")
            if (j==3 & R[1,3]==0) R=R[(2,1),]
        }
    }
    st_local("rootsettype",kind)
    return(R)
}
real rowvector umediation_cluster_sizes(real colvector id)
{
    real colvector sizes
    real scalar g,j,n
    g=max(id); n=rows(id); sizes=J(g,1,0)
    for (j=1;j<=g;j++) sizes[j]=sum(id:==j)
    return((min(sizes),umediation_q(sizes,.5),max(sizes),n^2/sum(sizes:^2),max(sizes)/n))
}
void umediation_jksummarize(real matrix D, real rowvector b, real rowvector at,
    real scalar L, real scalar U, real scalar level, string scalar direction, string rowvector out)
{
    real matrix Z,C,V,F,VC,CC,CI,T,G,I,R,DG,RS,IF,shift
    real rowvector theta,se,lo,hi,s,k,bb,d,a1,km,eta,lc,uc
    real colvector valid,grid,sd,roots
    real scalar n,p,alpha,mid,crit,mcrit,j,x,posneg,negpos,positive,negative,root,rv,denomclear,pathclear
    string scalar classification,rootstatus
    p=cols(b); valid=(D[,p+3]:==0):&(rowmissing(D[,1..p]):==0)
    if (!all(valid)) _error(459,"all cluster deletion fits must identify every target path; no failed clusters may be omitted")
    Z=D[,1..p]; n=rows(Z); C=Z:-mean(Z); V=(n-1)/n*quadcross(C,C)
    alpha=1-level/100; mid=(L+U)/2
    F=(Z[,1]+2*mid*Z[,2],2*Z[,2],Z[,3],2*Z[,5])
    theta=(b[1]+2*mid*b[2],2*b[2],b[3],2*b[5])
    C=F:-mean(F); VC=(n-1)/n*quadcross(C,C); se=sqrt(diagonal(VC))'
    crit=invttail(n-1,alpha/8); mcrit=invttail(n-1,alpha/2)
    lo=theta-crit*se; hi=theta+crit*se; CC=(theta',se',lo',hi')
    s=(lo[1],hi[1]); k=(lo[2],hi[2]); bb=(lo[3],hi[3]); d=(lo[4],hi[4])
    a1=umediation_interval_linear(s,k,-mid)
    km=umediation_interval_product(k,bb); eta=umediation_interval_product(a1,bb)
    lc=umediation_interval_product(umediation_interval_linear(s,k,L-mid),bb)
    uc=umediation_interval_product(umediation_interval_linear(s,k,U-mid),bb)
    CI=J(2,p,.)
    CI[,1]=a1'; CI[,2]=k'/2; CI[,3]=bb'
    CI[,4]=(b[4]-mcrit*sqrt(V[4,4])\b[4]+mcrit*sqrt(V[4,4]))
    CI[,5]=d'/2; CI[,6]=eta'; CI[,7]=km'; CI[,8]=d'; CI[,9]=(km+d)'
    CI[,10]=lc'; CI[,11]=uc'; CI[,12]=(U-L)*km'
    for (j=1;j<=cols(at);j++) {
        CI[,12+j]=umediation_interval_product(umediation_interval_linear(s,k,at[j]-mid),bb)'
    }
    T=(b',sqrt(diagonal(V)),CI')
    grid=rangen(L,U,101)
    sd=sqrt(rowmax((J(101,1,0),V[6,6]:+2*grid*V[6,7]+grid:^2*V[7,7])))
    G=(grid,b[6]:+b[7]*grid,sd,J(101,4,.))
    for (j=1;j<=101;j++) {
        lc=umediation_interval_product(umediation_interval_linear(s,k,grid[j]-mid),bb)
        G[j,4..7]=(lc,lc)
    }
    lc=CI[,10]'; uc=CI[,11]'
    posneg=lc[1]>0 & uc[2]<0 & (direction=="either" | direction=="posneg")
    negpos=lc[2]<0 & uc[1]>0 & (direction=="either" | direction=="negpos")
    positive=min(G[,6])>0; negative=max(G[,7])<0
    I=(lc[1],uc[2],alpha,posneg\lc[2],uc[1],alpha,negpos)
    classification="uncertain"
    if (positive) classification="support_same_sign_positive"
    if (negative) classification="support_same_sign_negative"
    if (posneg) classification="support_reversal_posneg"
    if (negpos) classification="support_reversal_negpos"
    root=(b[2]!=0 & b[3]!=0) ? -b[1]/(2*b[2]) : .
    roots=D[,p+1]; rv=sum(roots:<.)
    denomclear=(k[1]>0 | k[2]<0); pathclear=(bb[1]>0 | bb[2]<0)
    RS=umediation_root_region(s,k,mid)
    R=(root,.,.,rv,rv/n,(root<. ? (root>L & root<U) : .),denomclear,pathclear)
    if (RS[1,3]==0 & RS[1,4]==0 & RS[2,3]==.) R[1,2..3]=RS[1,1..2]
    rootstatus="outside_interval"
    if (root==.) rootstatus="undefined"
    else if (!denomclear | !pathclear) rootstatus="unstable_or_weak_path"
    else if (root==L | root==U) rootstatus="at_boundary"
    else if (root>L & root<U) rootstatus="inside_interval"
    shift=abs(F:-theta)
    for (j=1;j<=4;j++) {
        if (se[j]>0) shift[,j]=shift[,j]/se[j]
        else shift[,j]=J(n,1,0)
    }
    IF=((1::n),rowmax(shift),Z[,1]:-b[1],Z[,2]:-b[2],Z[,3]:-b[3],Z[,7]:-b[7],
        (sign(Z[,10]):!=sign(b[10])):|(sign(Z[,11]):!=sign(b[11])))
    DG=(n,0,crit,positive,negative)
    st_matrix(out[1],V); st_matrix(out[2],CI); st_matrix(out[3],T); st_matrix(out[4],G)
    st_matrix(out[5],R); st_matrix(out[6],I); st_matrix(out[7],DG)
    st_matrix(out[8],CC); st_matrix(out[9],RS); st_matrix(out[10],IF)
    st_local("classification",classification); st_local("rootstatus",rootstatus)
}
real scalar umediation_wild_enumerated(real scalar g, real scalar B, string scalar kind)
{
    real scalar m
    m=(kind=="webb" ? 6 : 2)
    return((g*ln(m)<=ln(B)+1e-12) & (g*ln(m)<=ln(100000)+1e-12))
}
real matrix umediation_wild_weights(real scalar g, real scalar B, string scalar kind)
{
    real rowvector values
    real matrix W
    real scalar m,count,j
    real colvector code
    values=(kind=="webb" ? (-sqrt(1.5),-1,-sqrt(.5),sqrt(.5),1,sqrt(1.5)) : (-1,1))
    m=cols(values)
    if (umediation_wild_enumerated(g,B,kind)) {
        count=round(m^g); W=J(count,g,.); code=(0::count-1)
        for (j=1;j<=g;j++) W[,j]=values[mod(floor(code/m^(j-1)),m):+1]'
    }
    else {
        W=J(B,g,.)
        for (j=1;j<=g;j++) W[,j]=values[floor(m*runiform(B,1)):+1]'
    }
    return(W)
}
real matrix umediation_wild_calc(real matrix X, real colvector y, real colvector id,
    string scalar partial, real matrix R, real matrix W, string scalar variant, real scalar enumerated)
{
    real matrix q,XX,xy,inv,AG,XYG,A,vv,F,score,delta_star,Av,result,betaleave,Fi
    real colvector levels,f,ix,yg,beta,vr,r,v,br,brg,fv,se,ts,invalid
    real rowvector mu
    real scalar n,p,g,j,k,a,b,observed,actualse,estimate,plus,tol,above,below
    if (partial!="") {
        f=st_data(.,partial); levels=uniqrows(sort(f,1))
        for (j=1;j<=rows(levels);j++) {
            ix=selectindex(f:==levels[j]); X[ix,]=X[ix,]:-mean(X[ix,]); y[ix]=y[ix]:-mean(y[ix])
        }
    }
    else X=(X,J(rows(X),1,1))
    n=rows(X); p=cols(X); g=max(id); XX=quadcross(X,X); xy=quadcross(X,y)
    inv=invsym(XX)
    if (any(diagonal(inv):==0) | hasmissing(inv)) _error(459,"wild design is singular after fixed effects")
    beta=inv*xy; AG=J(g*p,p,.); XYG=J(g,p,.)
    for (j=1;j<=g;j++) {
        ix=selectindex(id:==j); AG[(j-1)*p+1..j*p,]=quadcross(X[ix,],X[ix,])
        XYG[j,]=quadcross(X[ix,],y[ix])'
    }
    result=J(6+rows(W),cols(R),.)
    plus=enumerated ? 0 : 1
    for (k=1;k<=cols(R);k++) {
        r=R[,k]; v=inv*r; estimate=(r'*beta)[1,1]
        F=J(g,p,.); Av=J(g,p,.)
        for (j=1;j<=g;j++) {
            A=AG[(j-1)*p+1..j*p,]
            F[j,]=XYG[j,]-(A*beta)'; Av[j,]=(A*v)'
        }
        actualse=sqrt(g/(g-1)*sum((F*v):^2))
        if (actualse<=1e-14 | actualse>=.) _error(459,"wild observed cluster SE is zero/too small; rescale or check identification")
        observed=estimate/actualse
        br=beta-v*(estimate/(r'*v)[1,1])
        for (j=1;j<=g;j++) {
            A=AG[(j-1)*p+1..j*p,]
            if (variant=="classic") F[j,]=XYG[j,]-(A*br)'
            else {
                Fi=invsym(XX-A)
                if (any(diagonal(Fi):==0) | hasmissing(Fi)) _error(459,"WCR-score requires every leave-cluster design to be identified; no fallback")
                brg=Fi*(xy-XYG[j,]'); vr=Fi*r
                brg=brg-vr*((r'*brg)[1,1]/(r'*vr)[1,1])
                F[j,]=XYG[j,]-(A*brg)'
            }
        }
        fv=F*v; delta_star=W*F*inv
        score=W:*(J(rows(W),1,1)*fv')-delta_star*Av'
        se=sqrt(g/(g-1)*rowsum(score:^2)); ts=(W*fv):/se
        invalid=(se:<=1e-14):|(se:>=.):|(ts:>=.)
        ts=editmissing(ts,.)
        tol=1e-10*(1+abs(observed))
        above=sum(invalid:|(ts:>=observed-tol)); below=sum(invalid:|(ts:<=observed+tol))
        result[1..6,k]=(estimate\actualse\observed\(above+plus)/(rows(W)+plus)\(below+plus)/(rows(W)+plus)\sum(invalid))
        result[7..rows(result),k]=ts
    }
    return(result)
}
real matrix umediation_wild_combine(real matrix S, real scalar alpha)
{
    real rowvector b
    real scalar pp,pn,pe
    b=(max((S[3,4],S[1,4],S[2,5])),max((S[3,5],S[1,5],S[2,4])),
       max((S[3,4],S[1,5],S[2,4])),max((S[3,5],S[1,4],S[2,5])))
    pp=min((1,2*min(b[1..2]))); pn=min((1,2*min(b[3..4]))); pe=min((1,4*min(b)))
    return((pp,pp<=alpha\pn,pn<=alpha\pe,pe<=alpha))
}
end
