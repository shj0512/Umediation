{smcl}
{* *! version 1.0.0 23oct2026}{...}
{vieweralsosee "umediation" "help umediation"}{...}
{vieweralsosee "twoway options" "help twoway_options"}{...}

{title:Title}

{p2colset 5 23 25 2}{...}
{p2col:{bf:umediationplot} {hline 2}}Plot local indirect effects after umediation{p_end}
{p2colreset}{...}

{title:Syntax}

{p 8 16 2}
{cmd:umediationplot} [{cmd:,} {opt pointwise} {it:twoway_options}]

{title:Description}

{pstd}
{cmd:umediationplot} plots IE(x) and confidence limits from the last
{help umediation} result at 101 positions over {cmd:e(range)}. The default
limits are the simultaneous grid band for pairs inference and the
joint-region projection band for jackknife inference. Confidence levels
are nominal; see {help umediation##validation:Validation and interpretation}.

{title:Options}

{phang}
{opt pointwise} selects percentile intervals for pairs inference. Under
jackknife inference it displays the same projected limits as the default.

{phang}
{it:twoway_options} are passed to the graph. See {help twoway_options}.

{title:Remarks and examples}

{pstd}
Run {cmd:umediation} first. The graph uses the primary inference method;
the optional wild reversal test does not supply its limits. The command
preserves the data and estimation results and stores no additional
estimation results.

{phang2}{cmd:umediationplot}{p_end}
{phang2}{cmd:umediationplot, pointwise}{p_end}
{phang2}{cmd:umediationplot, name(indirect_effect, replace) title("Local indirect effect")}{p_end}

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
{help umediation}, {help twoway}, {help twoway_options}{p_end}
