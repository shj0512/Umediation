version 17.0
* Huijie SUN, University of Science and Technology of China
* Self-contained synthetic panel example. Optional argument: resampling draws.
args draws
if "`draws'" == "" local draws 199
capture program drop _umediation_example
program define _umediation_example
    version 17.0
    args draws
    preserve
    clear
    set seed 42
    set obs 400
    generate long firm = ceil(_n/8)
    generate byte year = mod(_n-1,8)+1
    generate double x = runiform(-2,2)
    replace x = -1.9 in 1
    replace x = 1.9 in 400
    generate double z = rnormal()
    generate double u = rnormal()
    generate double v = rnormal()
    bysort firm (year): replace u = .65*u[_n-1] + sqrt(1-.65^2)*u if _n>1
    bysort firm (year): replace v = .65*v[_n-1] + sqrt(1-.65^2)*v if _n>1
    generate double m = .25*x-.25*x^2+.3*z+.1*sin(firm)+.05*cos(year)+u
    generate double y = .15*x+.1*x^2+.6*m+.2*z+.2*cos(firm)+.05*sin(year)+v

    * Primary pairs bootstrap inference.
    umediation y, x(x) mediator(m) controls(z) absorb(firm year) cluster(firm) ///
        at(-1.5 0 1.5) range(-1.5 1.5) reps(`draws') seed(42) nodots
    umediationplot, name(umediation_pairs, replace)

    * Cluster jackknife projection and the separate wild reversal test.
    umediation y, x(x) mediator(m) controls(z) absorb(firm year) cluster(firm) ///
        at(-1.5 0 1.5) range(-1.5 1.5) inference(jackknife) ///
        wildreps(`draws') seed(42) nodots
    matrix list e(wild_reversal)
    matrix list e(root_set)
    umediationplot, name(umediation_projection, replace)
    display as text "Synthetic example. See help umediation for interpretation and validation."
    restore
end
local saved_rng `c(rngstate)'
capture noisily _umediation_example `draws'
local rc = _rc
quietly set rngstate `saved_rng'
capture program drop _umediation_example
if `rc' exit `rc'
