version 17.0
* Run from the source root: do validation/check.do
* Huijie SUN, University of Science and Technology of China
clear all
set more off
capture set processors 1
capture mkdir "validation/output"
capture erase "validation/output/PASSED.txt"
capture log close
log using "validation/output/check.log", text replace
adopath ++ "."
quietly run "umediation.ado"

* Independent OLS coefficients on the same sample.
set seed 42
set obs 400
generate long firm=ceil(_n/8)
generate byte year=mod(_n-1,8)+1
generate byte category=mod(firm,3)
generate byte region=ceil(firm/5)
generate double x=runiform(-2,2)
replace x=-1.9 in 1
replace x=1.9 in 400
generate double z=rnormal()
generate double m=.25*x-.25*x^2+.3*z+.15*sin(firm)+rnormal()
generate double y=.15*x+.1*x^2+.6*m+.2*z+.15*cos(firm)+rnormal()
generate double x2=x^2
quietly regress m c.z i.category i.firm i.year x x2
scalar a1=_b[x]
scalar a2=_b[x2]
quietly regress y c.z i.category i.firm i.year x x2 m
scalar bp=_b[m]
scalar c1=_b[x]
scalar c2=_b[x2]
local state `c(rngstate)'
quietly datasignature
local signature `r(datasignature)'
local originalframe `c(frame)'
umediation y, x(x) mediator(m) controls(c.z i.category) absorb(firm year) ///
    cluster(firm) at(-1 0 1) range(-1.5 1.5) reps(99) seed(3201) nodots ///
    saving("validation/output/draws.dta") replace
assert "`e(version)'"=="1.0.0"
assert "`e(release_status)'"=="first_release"
assert "`e(calibration_status)'"=="not_uniformly_calibrated"
assert e(N)==400 & e(N_clust)==50
assert abs(_b[a1]-a1)<1e-10 & abs(_b[a2]-a2)<1e-10
assert abs(_b[bpath]-bp)<1e-10 & abs(_b[c1]-c1)<1e-10 & abs(_b[c2]-c2)<1e-10
assert abs(_b[kappa_m]-2*a2*bp)<1e-10
assert abs(_b[kappa_t]-(2*c2+2*a2*bp))<1e-10
assert abs(_b[delta_ie]-3*_b[kappa_m])<1e-10
assert "`c(rngstate)'"=="`state'"
assert "`c(frame)'"=="`originalframe'"
quietly datasignature
assert "`r(datasignature)'"=="`signature'"
matrix PB=e(b)
matrix PCI=e(ci)
matrix PV=e(V)
local pc `e(classification)'
umediation
umediationplot, name(check_pairs,replace)
umediationplot, pointwise name(check_pointwise,replace)
assert mreldif(PB,e(b))<1e-14 & mreldif(PCI,e(ci))<1e-14
graph export "validation/output/pointwise.png", width(1200) replace
preserve
use "validation/output/draws.dta", clear
assert _N==99 & return_code==0
assert "`: char _dta[umediation_version]'"=="1.0.0"
assert "`: char _dta[release_status]'"=="first_release"
assert replication==_n
restore

* Reproducible draws and optional-test separation.
umediation y, x(x) mediator(m) controls(c.z i.category) absorb(firm year) ///
    cluster(firm) at(-1 0 1) range(-1.5 1.5) reps(99) seed(3201) nodots ///
    wildreps(199) saving("validation/output/paired_draws.dta") ///
    wildsaving("validation/output/wild.dta") replace
assert mreldif(PB,e(b))<1e-14 & mreldif(PCI,e(ci))<1e-14
assert mreldif(PV,e(V))<1e-14 & "`e(classification)'"=="`pc'"
matrix W=e(wild_stats)
matrix WR=e(wild_reversal)
assert e(wild_reps)==199 & e(wild_enumerated)==0
assert W[1,4]>0 & W[1,5]>0 & W[3,4]>0 & W[3,5]>0
umediation y, x(x) mediator(m) controls(c.z i.category) absorb(firm year) ///
    cluster(firm) at(-1 0 1) range(-1.5 1.5) inference(jackknife) ///
    seed(3201) wildreps(199) nodots
assert mreldif(W,e(wild_stats))<1e-10
assert mreldif(WR,e(wild_reversal))<1e-14
assert e(reps_valid)==50 & e(reps_failed)==0 & e(df_r)==49
assert e(joint_components)==4 & missing(e(sup_critical))
matrix J=e(ci)
matrix JB=e(b)
local jc `e(classification)'
umediation y, x(x) mediator(m) controls(c.z i.category) absorb(firm year) ///
    cluster(firm) at(-1 0 1) range(-1.5 1.5) inference(jackknife) nodots
assert mreldif(J,e(ci))<1e-14 & "`e(classification)'"=="`jc'"
umediationplot, name(check_projection,replace)
graph export "validation/output/projection.png", width(1200) replace

* Exposure translation and numeric/string identifiers.
generate double shifted_x=x+3
umediation y, x(shifted_x) mediator(m) controls(c.z i.category) absorb(firm year) ///
    cluster(firm) at(2 3 4) range(1.5 4.5) inference(jackknife) nodots
matrix S=e(ci)
matrix JA=J[1..2,7..15]
matrix SA=S[1..2,7..15]
assert mreldif(JA,SA)<1e-9
tostring firm, generate(firm_text) format(%03.0f)
tostring year, generate(year_text) format(%03.0f)
umediation y, x(x) mediator(m) controls(c.z i.category) absorb(firm_text year_text) ///
    cluster(firm_text) at(-1 0 1) range(-1.5 1.5) inference(jackknife) nodots
assert mreldif(JB,e(b))<1e-10 & mreldif(J,e(ci))<1e-9

* A coarser cluster identifier and deletion failure handling.
umediation y, x(x) mediator(m) controls(z) absorb(firm year) cluster(region) ///
    inference(jackknife) wildreps(99) seed(54) nodots
assert "`e(nested_absorb)'"!="" & e(N_clust)==10
matrix BEFORE=e(b)
local previousclass `e(classification)'
replace x=0 if firm!=1
bysort firm: replace x=-2+4*(_n-1)/(_N-1) if firm==1
quietly datasignature
local signature `r(datasignature)'
local state `c(rngstate)'
capture noisily umediation y, x(x) mediator(m) controls(z) absorb(firm year) ///
    cluster(firm) inference(jackknife) at(-1 0 1) range(-1.5 1.5) seed(10) nodots
assert _rc==459
assert mreldif(BEFORE,e(b))<1e-14 & "`e(classification)'"=="`previousclass'"
assert "`c(rngstate)'"=="`state'"
quietly datasignature
assert "`r(datasignature)'"=="`signature'"

* Complete enumeration and finite support.
clear
set seed 98
set obs 180
generate byte g=ceil(_n/60)
generate double x=runiform(-2,2)
generate double m=.25*x-.25*x^2+rnormal()
generate double y=.6*m+.1*x^2+rnormal()
umediation y, x(x) mediator(m) cluster(g) inference(jackknife) ///
    wildreps(499) wildweight(webb) wildvariant(classic) nodots
assert e(wild_reps)==216 & e(wild_enumerated)==1
matrix W=e(wild_stats)
assert abs(W[1,4]*216-round(W[1,4]*216))<1e-9
umediation y, x(x) mediator(m) cluster(g) inference(jackknife) ///
    wildreps(99) wildweight(rademacher) wildvariant(score) nodots
assert e(wild_reps)==8 & e(wild_enumerated)==1

* Root confidence-set edge cases.
mata: st_matrix("R",umediation_root_region((1,2),(1,2),0))
assert R[1,1]==-2 & R[1,2]==-.5 & R[1,3]==0 & R[1,4]==0
mata: st_matrix("R",umediation_root_region((1,2),(-1,1),0))
assert R[1,3]==1 & R[1,2]==-1 & R[2,1]==1 & R[2,4]==1
mata: st_matrix("R",umediation_root_region((-1,1),(-1,1),0))
assert R[1,3]==1 & R[1,4]==1 & missing(R[2,3])
mata: st_matrix("R",umediation_root_region((0,1),(0,1),0))
assert R[1,3]==1 & R[1,2]==0 & R[1,4]==0
mata: st_matrix("R",umediation_root_region((-1,0),(0,1),0))
assert R[1,1]==0 & R[1,3]==0 & R[1,4]==1

* Invalid requests leave current estimates in place.
capture noisily umediation y, x(x) mediator(m) inference(jackknife)
assert _rc==198
capture noisily umediation y, x(x) mediator(m) cluster(g) wildreps(1)
assert _rc==198
capture noisily umediation y, x(x) mediator(m) cluster(g) inference(unknown)
assert _rc==198
capture noisily umediation y, x(x) mediator(m) cluster(g) wildreps(99) ///
    saving("same.dta") wildsaving("same.dta") replace
assert _rc==198

* Output path aliases are rejected before replacing an existing file.
file open sentinel using "validation/output/guard.dta", write replace
file write sentinel "preserve this file" _n
file close sentinel
local absolute_alias `"`c(pwd)'/validation/output/../output/guard.dta"'
capture noisily umediation y, x(x) mediator(m) cluster(g) inference(jackknife) ///
    wildreps(99) saving("validation/output/guard.dta") ///
    wildsaving(`"`absolute_alias'"') replace nodots
assert _rc==198
if "`c(os)'"=="Windows" {
    capture noisily umediation y, x(x) mediator(m) cluster(g) inference(jackknife) ///
        wildreps(99) saving("validation/output/guard.dta") ///
        wildsaving("VALIDATION/OUTPUT/GUARD.dta") replace nodots
    assert _rc==198
}
file open sentinel using "validation/output/guard.dta", read text
file read sentinel sentinel_line
file close sentinel
assert `"`sentinel_line'"'=="preserve this file"
mata: assert(umediation_same_output("./a b/../draws.dta","C:/work/draws.dta","C:/work","Windows"))
mata: assert(umediation_same_output("C:/work/a/../draws.dta","C:/work/draws.dta","C:/work","Windows"))
mata: assert(umediation_same_output("C:/WORK/draws.dta","c:"+char(92)+"work"+char(92)+"DRAWS.dta","C:/work","Windows"))
mata: assert(umediation_same_output("../draws.dta","/work/draws.dta","/work/sub","Unix"))
mata: assert(!umediation_same_output("/work/DRAWS.dta","/work/draws.dta","/work","Unix"))
mata: assert(umediation_same_output("//server/share/a/../draws.dta","//SERVER/share/DRAWS.dta","C:/work","Windows"))

* Native SMCL rendering and a complete standalone example.
help umediation
help umediationplot
translate "umediation.sthlp" "validation/output/umediation_help.txt", translator(smcl2log) replace
translate "umediationplot.sthlp" "validation/output/umediationplot_help.txt", translator(smcl2log) replace
foreach rendered in umediation_help umediationplot_help {
    file open rendered_help using "validation/output/`rendered'.txt", read text
    file read rendered_help line
    while r(eof)==0 {
        assert strpos(`"`line'"',"{")==0
        assert strpos(`"`line'"',"}")==0
        file read rendered_help line
    }
    file close rendered_help
}
quietly datasignature
local signature `r(datasignature)'
local state `c(rngstate)'
do "umediation_example.do" 199
assert "`c(rngstate)'"=="`state'"
quietly datasignature
assert "`r(datasignature)'"=="`signature'"
file open passed using "validation/output/PASSED.txt", write replace
file write passed "PASS: numerical, inference, state, failure, enumeration, help and example checks." _n
file close passed
log close
