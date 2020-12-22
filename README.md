Randomization-based causal inference framework to analyze 16s rRNA gut microbiome data.
=======================================================================================

Framework
---------

![Image of Graphical abstract](misc/Fig1_graphical_abstract.png)

Data access
-----------

The KORA cohort study ([Holle et
al. (2005)](https://pubmed.ncbi.nlm.nih.gov/16032513/)) data is only
accessible after applying to a research project on the
[KORA-passt](https://epi.helmholtz-muenchen.de) platform.

In the [`microbiome_ASV_data`](microbiome_ASV_data) folder are the
details about the data pre-processing and the data files handed-in after
successful application to a KORA project.

Stage 2: Design
---------------

The R code for our pair matching implementation and diagnostic plots
generation can be found in the [`design`](design) file. The matrix of
10,000 possible randomization of the intervention assignment is also
generated directly after matching.

Note 1: the matching functions
[Stephane\_matching.R](misc/Stephane_matching.R) were written in Rcpp by
Dr. Stéphane Shao. Note 2: other matching strategies are valid. The
researcher should take the conceptual hypothetical experiment into
account when choosing its strategy. Note 3: to make the matching easier
we re-formated/coded the original KORA variables. See
[`misc/format_KORA_variables.R`](misc/format_KORA_variables.R) file.

Stage 3: Analysis
-----------------

### Diversity

#### Alpha and richness

See [`1_alpha_diversity`](1_alpha_diversity)

Willis’ R package [`breakaway`](https://github.com/adw96/breakaway) for
richness and [`DivNet`](https://github.com/adw96/DivNet) for shannon
index.

#### Beta

See [`2_beta_diversity`](2_beta_diversity)

Zhao R package documentation:
[`MiRKAT`](https://cran.r-project.org/web/packages/MiRKAT/index.html).

### Compostion

#### Compositional equivalence

See [`3_mean_diff_test`](3_mean_diff_test) for our implementation of
this test statistic.

Code:
[`composition-two-sampe-test`](https://github.com/yuanpeicao/composition-two-sampe-test)
github.

#### Differential abundance

See [`4_differential_abundance`](4_differential_abundance)

Brill R package: [`dacomp`](https://github.com/barakbri/dacomp) github.

#### Correlation structure

See [`5_networks`](5_networks)

Peschel R package: [`NetCoMi`](https://github.com/stefpeschel/NetCoMi)
github.

#### Further analyses (metabolites)

References
----------

\[Holle et al., 2005\] Holle R, Happich M, Löwel H, Wichmann HE (2005);
[MONICA/KORA Study Group. KORA–a research platform for population based
health research.](https://pubmed.ncbi.nlm.nih.gov/16032513/)
*Gesundheitswesen*, 67.
