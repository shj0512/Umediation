{smcl}
{* *! version 1.0.0 23oct2026}{...}
{vieweralsosee "umediationplot" "help umediationplot"}{...}
{vieweralsosee "regress" "help regress"}{...}
{vieweralsosee "areg" "help areg"}{...}
{viewerjumpto "Syntax" "umediation##syntax"}{...}
{viewerjumpto "Description" "umediation##description"}{...}
{viewerjumpto "Options" "umediation##options"}{...}
{viewerjumpto "Remarks and examples" "umediation##remarks"}{...}
{viewerjumpto "Stored results" "umediation##results"}{...}
{viewerjumpto "Methods and references" "umediation##methods"}{...}
{viewerjumpto "Author and contact" "umediation##author"}{...}

{title:Title}

{p2colset 5 20 22 2}{...}
{p2col:{bf:umediation} {hline 2}}Mediation with a quadratic exposure-to-mediator relationship{p_end}
{p2colreset}{...}

{marker syntax}{...}
{title:Syntax}

{p 8 16 2}
{cmd:umediation} {it:y} {ifin}{cmd:,}
{opt x(varname)} {opt mediator(varname)}
[{opt controls(fvvarlist)} {opt absorb(varlist)} {opt cluster(varname)}
{opt at(numlist)} {opt atpercentiles(numlist)}
{opt range(numlist)} {opt rangepercentiles(numlist)}
{opt direction(string)} {opt inference(string)}
{opt reps(#)} {opt seed(#)} {opt level(#)}
{opt saving(filename)} {opt wildreps(#)} {opt wildweight(string)}
{opt wildvariant(string)} {opt wildsaving(filename)} {opt replace} {opt nodots}]

{p 8 16 2}
{cmd:umediation}

{pstd}
The second syntax replays the last result. Stata 17 or later is required.
The command uses built-in Stata and Mata routines; no additional packages are required.

{marker description}{...}
{title:Description}

{pstd}
{cmd:umediation} fits two least-squares equations with a quadratic exposure-to-mediator
path and a linear mediator-to-outcome path. It estimates the local indirect effect
as a function of the exposure, decomposes outcome curvature into mediated and direct
components, and assesses whether the indirect effect changes sign over a specified
exposure interval. Both equations use the same complete-case sample, controls,
and fixed-effect variables. Their control and fixed-effect coefficients may differ.

{pstd}
The default inference method is a pairs bootstrap, resampling observations or
whole clusters. Cluster jackknife projection intervals and a separate wild cluster
bootstrap reversal test are also available. {help umediationplot} graphs the
estimated indirect effect and its confidence limits.

{pstd}
The model applies to continuous exposure, mediator, and outcome variables, with one
mediator and one clustering dimension. Estimation is unweighted. Fixed effects
account for group-specific intercepts; a causal mediation interpretation requires
additional identification assumptions about exposure, mediator, and outcome errors.

{marker options}{...}
{title:Options}

{dlgtab:Model}

{phang}
{opt x(varname)} specifies the numeric exposure. The command includes its linear
and squared terms in both equations. This option is required.

{phang}
{opt mediator(varname)} specifies the numeric mediator. The outcome {it:y},
exposure, and mediator must be distinct variables. This option is required.

{phang}
{opt controls(fvvarlist)} specifies common numeric controls and allows
{help fvvarlist:factor-variable notation}. The list may not contain {it:y},
the exposure, or the mediator.

{phang}
{opt absorb(varlist)} specifies categorical fixed-effect identifiers, which may
be numeric or string. The first is absorbed using {help areg}; additional fixed
effects enter as indicators. The list may not contain {it:y}, the exposure,
or the mediator.

{phang}
{opt cluster(varname)} specifies a numeric or string cluster identifier. Pairs
bootstrap draws resample whole clusters. Fixed-effect identifiers nested within
clusters are relabeled for each sampled copy. At least two clusters are required;
jackknife inference requires at least three. These are computational requirements.
Inference assumes independent clusters and allows dependence within clusters.

{dlgtab:Exposure positions and research interval}

{phang}
{opt at(numlist)} specifies literal exposure values at which local indirect effects
are reported. {opt atpercentiles(numlist)} instead specifies sample percentiles;
the default is {cmd:atpercentiles(5 25 50 75 95)}. The two options are mutually
exclusive. At most 100 positions are allowed.

{phang}
{opt range(numlist)} specifies two distinct endpoints {it:L} and {it:U} of the
research interval. {opt rangepercentiles(numlist)} instead specifies two sample
percentiles; the default is {cmd:rangepercentiles(5 95)}. The two options are
mutually exclusive. Percentile arguments must lie strictly between 0 and 100.

{phang}
Exposure positions and endpoints are sorted and must lie within the observed
exposure support in the estimation sample. They are determined once and held
fixed in all bootstrap draws and cluster deletions. The research interval
defines the region assessed for reversal; it does not restrict the estimation
sample. Positions specified in {cmd:at()} may lie outside that interval.

{phang}
{opt direction(string)} selects {cmd:posneg}, {cmd:negpos}, or {cmd:either}
(the default). {cmd:posneg} assesses IE({it:L})>0 and IE({it:U})<0;
{cmd:negpos} assesses the opposite signs. Choose the interval and direction
from the research question before inspecting test results.

{dlgtab:Primary inference}

{phang}
{opt inference(string)} selects {cmd:pairs} (the default) or {cmd:jackknife}.
Pairs inference provides percentile intervals, a simultaneous band on a
101-point grid, and directional endpoint tests. Jackknife inference constructs
a four-component confidence region and projects it to the reported quantities.
The projection and wild procedures remain experimental inference methods;
their finite-sample performance is described under
{help umediation##validation:Validation and interpretation}.

{phang}
{opt reps(#)} requests the number of pairs bootstrap draws. The default is 1999
and the minimum is 50. At least 50 draws and 80% of requested draws must be
valid. Failed fits are counted and are not retried. With {cmd:inference(jackknife)},
this option has no effect: each of the {it:G} clusters is deleted once, and
every deletion must identify all target coefficients. A failed deletion stops
the command without dropping the affected cluster from inference.

{phang}
{opt level(#)} specifies the nominal confidence level and the corresponding test
level. The default is Stata's current confidence level, {cmd:c(level)}, normally 95.

{phang}
{opt seed(#)} sets the seed used internally. The caller's random-number state
is restored on exit. When a wild test is requested, its random-number stream
is reset to the supplied seed independently of the primary inference method.

{dlgtab:Additional wild cluster test}

{phang}
{opt wildreps(#)} requests a separate wild cluster bootstrap test of reversal.
The default is 0, which disables the test; positive values must be at least 99.
The option requires {cmd:cluster()}. Use more draws, such as 1999 or 9999,
when estimating tail probabilities. The test supplies p-values in
{cmd:e(wild_reversal)}; it leaves the primary confidence intervals and
{cmd:e(classification)} unchanged.

{phang}
{opt wildweight(string)} selects {cmd:webb} (the default) or {cmd:rademacher}.
Webb weights take six equally likely values: +/-sqrt(1/2), +/-1, and
+/-sqrt(3/2). Rademacher weights are -1 and 1 with equal probability. The
complete weight support is enumerated when its size is at most the requested
number of draws and at most 100,000. The actual number of schedules can then
be smaller than {cmd:wildreps()}. With few clusters, the finite support
restricts the attainable p-values.

{phang}
{opt wildvariant(string)} selects {cmd:score} (the default, WCR-S) or
{cmd:classic} (WCR-C). WCR-S uses null-restricted scores adjusted with
restricted leave-cluster estimates; WCR-C uses unadjusted restricted scores.
Both use CV1 studentization. The score variant requires identified
leave-cluster designs and stops on failure.

{phang}
The wild engine permits at most 500 identified design columns after one fixed
effect nested within clusters has been absorbed. Other fixed effects enter
as indicators. This limit applies separately to the two equations.

{dlgtab:Saving and display}

{phang}
{opt saving(filename)} saves all pairs draws or cluster-deletion estimates as
a Stata dataset. The file contains {cmd:replication}, the coefficient and
estimand columns, {cmd:root}, {cmd:draw_n}, and {cmd:return_code}. A nonzero
return code records a failed fit. For jackknife inference, {cmd:replication}
identifies the deleted integer group created by {cmd:egen group(cluster)},
ordered by source cluster values.

{phang}
{opt wildsaving(filename)} saves the wild weight schedules {cmd:w1}, ..., {cmd:wG}
and bootstrap statistics {cmd:t_left}, {cmd:t_right}, and {cmd:t_b}, with method
and enumeration metadata. Weight columns follow the same integer cluster order.
This option requires positive {cmd:wildreps()}. The two output paths must
identify different files. Common relative/absolute path aliases and Windows
case-only differences are rejected. Avoid filesystem-link aliases for the
output files. A {cmd:.dta} extension is appended when absent.

{phang}
{opt replace} permits overwriting the files requested by {cmd:saving()} or
{cmd:wildsaving()}. {opt nodots} suppresses replication progress.

{marker remarks}{...}
{title:Remarks and examples}

{pstd}
The fitted model is

{p 8 8 2}M = a0 + a1*X + a2*X^2 + C*gM + FE_M + u{p_end}
{p 8 8 2}Y = c0 + c1*X + c2*X^2 + b*M + C*gY + FE_Y + v.{p_end}

{pstd}
The local indirect effect is IE(x)=b*(a1+2*a2*x). Its slope with respect to x
is the mediated curvature kappa_m=2*a2*b. Direct and total curvature are
kappa_d=2*c2 and kappa_t=kappa_m+kappa_d. The change in the indirect effect
over the interval is delta_ie=(U-L)*kappa_m. The parameter eta=a1*b is
IE(0), even when zero lies outside the observed exposure support.

{pstd}
Observations missing any model variable, control, fixed effect, or cluster
identifier are excluded from both equations. At least 20 complete observations
are required, with positive residual degrees of freedom in each regression.
The command preserves the caller's data, observation order, active frame,
and random-number state. Successful estimation replaces {cmd:e()}; an
estimation failure restores the previous estimates.

{pstd}
Example with a generated dataset:

{phang2}{cmd:clear}{p_end}
{phang2}{cmd:set seed 42}{p_end}
{phang2}{cmd:set obs 400}{p_end}
{phang2}{cmd:generate double x = runiform(-2,2)}{p_end}
{phang2}{cmd:generate double z = rnormal()}{p_end}
{phang2}{cmd:generate double m = .25*x - .25*x^2 + .3*z + rnormal()}{p_end}
{phang2}{cmd:generate double y = .15*x + .1*x^2 + .6*m + .2*z + rnormal()}{p_end}
{phang2}{cmd:umediation y, x(x) mediator(m) controls(z) reps(1999) seed(42)}{p_end}
{phang2}{cmd:umediationplot}{p_end}

{pstd}
With panel data, replace the variable names in the following commands with
the outcome, exposure, mediator, controls, and identifiers in your dataset:

{phang2}{cmd:umediation y, x(x) mediator(m) controls(z) absorb(firm year) cluster(firm) reps(1999) seed(42)}{p_end}
{phang2}{cmd:umediation y, x(x) mediator(m) controls(z) absorb(firm year) cluster(firm) inference(jackknife) wildreps(1999) seed(42)}{p_end}
{phang2}{cmd:matrix list e(wild_reversal)}{p_end}
{phang2}{cmd:matrix list e(root_set)}{p_end}

{pstd}
The accompanying {cmd:umediation_example.do} generates a panel dataset and
demonstrates both primary methods and the wild test. Obtain it with
{cmd:net get umediation, from("installation_source")}, using the same source
as installation, and run {cmd:do umediation_example.do}. Its default 199
wild draws are demonstration settings; {cmd:do umediation_example.do 1999}
requests more draws. The example preserves existing data and random-number state.

{marker validation}{...}
{pstd}{bf:Validation and interpretation}

{pstd}
Finite-sample coverage and reversal error rates depend on the design. In a
clustered simulation with N=400 and G=5 (5,000 outer replications), nominal
95% mediated-curvature coverage was 84.46% for pairs inference and 99.96%
for jackknife projection. In the corresponding b=0 design, false reversal
rates were 16.08%, 0.62%, and 1.08% for pairs, projection, and WCR-S.
Correct-direction power under the regular alternative was 93.28%, 23.02%,
and 40.02%. These results show substantial coverage and power tradeoffs.

{pstd}
An imbalanced G=10, N=400 design produced 91.15% coverage for the nominal
95% four-component region (2,000 outer replications). WCR-S had no detected
reversals under the evaluated high-leverage alternative. Few clusters,
unequal cluster influence, weak paths, heavy tails, and roots near the
research-interval boundary can impair inference. There is no universal safe
cluster-count threshold. Increasing bootstrap repetitions improves simulation
precision without correcting design-related distortions. The validation
materials accompanying the source describe the simulation settings and uncertainty.

{pstd}
A root estimate inside the interval alone does not establish a reversal.
Failure to establish a reversal does not establish its absence. The
same-sign classifications require the reported band to remain strictly on
one side of zero at every grid position. Specify the primary method in
advance and report its uncertainty, including unbounded root sets where applicable.

{pstd}
For X=(rawX-center)/scale, the local indirect effect per unit of rawX equals
IE(X)/scale; mediated curvature equals kappa_m/scale^2; and the root on the
raw scale equals center+scale*root. The command performs no automatic centering,
scaling, trimming, or control selection. Weights, multiple mediators, an M^2
term, an X-by-M interaction, instrumental-variable estimation, and multiway
clustering are outside the implemented model.

{marker results}{...}
{title:Stored results}

{pstd}
{cmd:umediation} stores the following in {cmd:e()}.

{synoptset 27 tabbed}{...}
{p2col 5 27 30 2: Matrices}{p_end}
{synopt:{cmd:e(b)}}a1, a2, bpath, c1, c2, eta, kappa_m, kappa_d, kappa_t, ie_l, ie_u, delta_ie, and ie1, ...{p_end}
{synopt:{cmd:e(V)}}pairs or centered cluster-deletion covariance of the stored estimands{p_end}
{synopt:{cmd:e(ci)}}lower and upper confidence limits, with columns matching e(b){p_end}
{synopt:{cmd:e(table)}}rows matching e(b); estimate, method-specific SE, ci_lower, ci_upper{p_end}
{synopt:{cmd:e(at)}}fixed exposure positions for ie1, ...{p_end}
{synopt:{cmd:e(range)}}lower and upper research-interval endpoints{p_end}
{synopt:{cmd:e(grid)}}101 rows; x, ie, method-specific SE, point_lower, point_upper, sim_lower, sim_upper{p_end}
{synopt:{cmd:e(iut)}}rows posneg and negpos; left_bound, right_bound, alpha, reject{p_end}
{synopt:{cmd:e(root_info)}}estimate, ci_lower, ci_upper, finite_reps, finite_share, inside, denominator_clear, path_clear{p_end}
{synopt:{cmd:e(diagnostics)}}valid, failed, sup_critical, all_positive, all_negative{p_end}
{synopt:{cmd:e(cluster_diagnostics)}}min_size, median_size, max_size, size_effective_G, largest_share; missing without clustering{p_end}
{synopt:{cmd:e(component_ci)}}jackknife only: rows slope_mid, curvature_xm, bpath, kappa_d; estimate, jackknife_se, ci_lower, ci_upper{p_end}
{synopt:{cmd:e(root_set)}}jackknife only: interval1 and interval2; lower, upper, lower_unbounded, upper_unbounded{p_end}
{synopt:{cmd:e(influence)}}jackknife only: cluster_group, scaled_max_shift, delta_a1, delta_a2, delta_b, delta_kappa_m, point_direction_changed{p_end}
{synopt:{cmd:e(wild_stats)}}wild test only: rows left_slope, right_slope, bpath; estimate, cluster_se, t, p_positive, p_negative, degenerate_draws{p_end}
{synopt:{cmd:e(wild_reversal)}}wild test only: rows posneg, negpos, either; p, reject{p_end}
{p2colreset}{...}

{pstd}
The size-only effective cluster count is N^2/sum(n_g^2); it describes cluster
sizes, not regressor-specific information or leverage. The final influence
column records whether deletion changes the sign of an endpoint point
estimate. It is a point-estimate diagnostic.

{pstd}
The root-set unbounded flags identify infinite endpoints stored as numeric
missing values. An entirely missing second row means that no second interval
is present. Use {cmd:e(root_set_type)} to distinguish {cmd:bounded},
{cmd:unbounded}, {cmd:disjoint}, {cmd:all_real}, and {cmd:undefined} sets.
The set concerns the first-path zero; when b is zero or weak, it does not
establish a mediated reversal. For jackknife inference, {cmd:e(root_info)}
gives a finite interval summary only when the projected set is bounded.

{synoptset 27 tabbed}{...}
{p2col 5 27 30 2: Scalars}{p_end}
{synopt:{cmd:e(N)}}common estimation sample size{p_end}
{synopt:{cmd:e(N_clust)}}number of clusters; missing without cluster(){p_end}
{synopt:{cmd:e(level)}}nominal confidence level{p_end}
{synopt:{cmd:e(reps_requested)}}requested pairs draws or number of jackknife deletions{p_end}
{synopt:{cmd:e(reps_valid)}}valid primary fits{p_end}
{synopt:{cmd:e(reps_failed)}}failed primary fits{p_end}
{synopt:{cmd:e(grid_points)}}101{p_end}
{synopt:{cmd:e(root)}}-a1/(2*a2); missing if a2 or bpath is zero{p_end}
{synopt:{cmd:e(root_inside)}}whether the estimated root lies strictly between L and U{p_end}
{synopt:{cmd:e(sup_critical)}}pairs maximum-standardized-deviation critical value; missing for jackknife{p_end}
{synopt:{cmd:e(component_critical)}}jackknife component t critical value{p_end}
{synopt:{cmd:e(joint_components)}}4 for jackknife{p_end}
{synopt:{cmd:e(df_r)}}G-1 for jackknife{p_end}
{synopt:{cmd:e(wild_reps)}}actual wild weight schedules, if requested{p_end}
{synopt:{cmd:e(wild_enumerated)}}1 for complete enumeration and 0 for sampled weights{p_end}
{p2colreset}{...}

{synoptset 27 tabbed}{...}
{p2col 5 27 30 2: Macros and function}{p_end}
{synopt:{cmd:e(cmd)}}umediation{p_end}
{synopt:{cmd:e(cmdline)}}command as entered{p_end}
{synopt:{cmd:e(title)}}model description{p_end}
{synopt:{cmd:e(version)}}1.0.0{p_end}
{synopt:{cmd:e(release_status)}}first_release{p_end}
{synopt:{cmd:e(calibration_status)}}not_uniformly_calibrated{p_end}
{synopt:{cmd:e(depvar)}, {cmd:e(x)}}outcome and exposure names{p_end}
{synopt:{cmd:e(mediator)}}mediator name{p_end}
{synopt:{cmd:e(controls)}}controls specification{p_end}
{synopt:{cmd:e(absorb)}}fixed-effect identifiers{p_end}
{synopt:{cmd:e(cluster)}}cluster identifier{p_end}
{synopt:{cmd:e(nested_absorb)}}fixed-effect identifiers nested within clusters{p_end}
{synopt:{cmd:e(inference)}}pairs or jackknife{p_end}
{synopt:{cmd:e(vce)}}bootstrap or clusterjackknife{p_end}
{synopt:{cmd:e(vcetype)}}Bootstrap or Cluster jackknife CV3J{p_end}
{synopt:{cmd:e(quantile_method)}}linear_type7 or bonferroni4_t_projection{p_end}
{synopt:{cmd:e(direction)}}posneg, negpos, or either{p_end}
{synopt:{cmd:e(classification)}}support_reversal_posneg, support_reversal_negpos, support_same_sign_positive, support_same_sign_negative, or uncertain{p_end}
{synopt:{cmd:e(root_status)}}inside_interval, outside_interval, at_boundary, unstable_or_weak_path, or undefined{p_end}
{synopt:{cmd:e(root_set_type)}}jackknife root-set topology{p_end}
{synopt:{cmd:e(seed)}}supplied seed; empty if omitted{p_end}
{synopt:{cmd:e(saving)}}primary replication filename{p_end}
{synopt:{cmd:e(wildweight)}}Webb or Rademacher specification, stored in lowercase{p_end}
{synopt:{cmd:e(wildvariant)}}score or classic, if the wild test is requested{p_end}
{synopt:{cmd:e(wildsaving)}}wild replication filename{p_end}
{synopt:{cmd:e(sample)}}common estimation-sample indicator{p_end}
{p2colreset}{...}

{marker methods}{...}
{title:Methods and references}

{pstd}{bf:Pairs bootstrap}

{pstd}
Intervals use linearly interpolated empirical quantiles with index
h=1+(B-1)*p (type 7) across valid draws. This convention can differ from
Stata's built-in finite-bootstrap percentile convention. The simultaneous
grid band uses the bootstrap quantile of the maximum absolute deviation
from the original IE(x), divided by the bootstrap SD at each grid position.
It uses common SDs across draws. Directional endpoint tests use one-sided
percentile bounds at alpha=1-level/100 for a prespecified direction and
alpha/2 for each direction under {cmd:direction(either)}. The root
percentile interval is a descriptive ratio diagnostic.

{pstd}{bf:Cluster jackknife projection}

{pstd}
Let m0=(L+U)/2, s=a1+2*a2*m0, k=2*a2, and theta=(s,k,b,kappa_d).
For G delete-one-cluster estimates, the centered CV3J covariance is
(G-1)/G times the sum of outer products of theta_(-g)-mean(theta_(-g)).
The four two-sided component intervals use t_(G-1,1-alpha/8). The
Bonferroni rectangle has nominal joint level 1-alpha, conditional on
the calibration of the marginal t intervals.

{pstd}
Interval arithmetic projects this region through IE(x)=b*(s+k*(x-m0)),
kappa_m=b*k, and kappa_t=kappa_d+b*k. The projection applies to every
x in the interval; the stored graph uses 101 positions. The c1 interval
is an ordinary marginal CV3J t interval outside the four-component family.
Reported nonlinear standard errors come from the deletion covariance;
projected limits are computed from the component region.

{pstd}
The first-path zero m0-s/k is projected over the (s,k) rectangle with
k=0 excluded. Its set can be bounded, one ray, two rays, all real numbers,
or undefined. This procedure is a rectangular-region projection; it does
not use Fieller inversion or invert the wild tests.

{pstd}{bf:Wild reversal test}

{pstd}
The test evaluates the linear contrasts dL=a1+2*a2*L, dU=a1+2*a2*U,
and b. For a positive-to-negative indirect-effect reversal, the two
branches are (b>0,dL>0,dU<0) and (b<0,dL<0,dU>0). Each branch uses
the largest of its three one-sided p-values. The directional p-value is
min(1,2*min(branch p-values)); the {cmd:either} p-value is
min(1,4*min(all four branch p-values)). Calibration of this combination
depends on the calibration of its component tests.

{pstd}
The same cluster weights are used for all three contrasts. Sampled weights
use (1+tail count)/(B+1); complete enumeration uses tail count/B.
Ties are included within 1e-10*(1+abs(observed t)). Schedules with a
standard error at most 1e-14, or an undefined statistic, count against
rejection in both tails. The observed component standard error must exceed
1e-14. The reported wild SE omits the common (N-1)/(N-p) CV1 factor;
this factor cancels in the observed/bootstrap t comparison.

{pstd}
The wild procedure supplies linear-component tests combined for reversal.
Product confidence intervals and an adaptive bootstrap for double-zero
mediation nulls are outside this implementation.

{pstd}{bf:References}

{phang}
MacKinnon, J. G., M. O. Nielsen, and M. D. Webb. 2023.
Fast and reliable jackknife and bootstrap methods for cluster-robust inference.
{browse "https://arxiv.org/abs/2301.04527":arXiv:2301.04527}.

{marker author}{...}
{title:Author and contact}

{pstd}
Huijie SUN{break}
University of Science and Technology of China{break}
{browse "mailto:shj0512@mail.ustc.edu.cn":shj0512@mail.ustc.edu.cn}{break}
{browse "https://github.com/shj0512":github.com/shj0512}

{pstd}
Copyright (c) 2026 Huijie SUN. Distributed under the MIT License.

{title:Also see}

{psee}
{help umediationplot}, {help regress}, {help areg}, {help bootstrap},
{help jackknife}, {help fvvarlist}, {help net}{p_end}
