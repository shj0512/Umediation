*! umediationplot 1.0.0 23oct2026
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

program define umediationplot
    version 17.0
    syntax [, POINTwise *]
    if "`e(cmd)'" != "umediation" error 301
    di as text "Umediation: confidence levels are nominal; see help umediation."
    tempname G
    matrix `G' = e(grid)
    local xname `e(x)'
    local cilevel = e(level)
    local lo 6
    local hi 7
    local band "Simultaneous grid band"
    if "`pointwise'" != "" {
        local lo 4
        local hi 5
        local band "Pointwise percentile intervals"
    }
    if "`e(inference)'" == "jackknife" local band "joint-region projection band"
    preserve
    quietly clear
    quietly set obs 101
    quietly generate double x = .
    quietly generate double ie = .
    quietly generate double lower = .
    quietly generate double upper = .
    mata: st_store(.,("x","ie","lower","upper"),st_matrix("`G'")[,(1,2,`lo',`hi')])
    twoway (rarea lower upper x, color(navy%15) lcolor(none)) ///
        (line lower x, lcolor(gs8) lpattern(shortdash) lwidth(thin)) ///
        (line upper x, lcolor(gs8) lpattern(shortdash) lwidth(thin)) ///
        (line ie x, lcolor(navy) lwidth(medthick)), ///
        yline(0, lcolor(gs8) lpattern(dash)) ///
        xtitle("`xname'") ytitle("Local indirect effect") ///
        legend(order(4 "IE(x)" 2 "Nominal `cilevel'% `band'")) `options'
    restore
end
