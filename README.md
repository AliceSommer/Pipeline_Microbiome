Randomization-based causal inference framework to analyze 16s rRNA gut microbiome data.
=======================================================================================

Framework
---------

![Image of Graphical abstract](misc/Fig1_graphical_abstract.png)

Data access
-----------

Blabla on KORA ([Holle et
al. (2005)](https://pubmed.ncbi.nlm.nih.gov/16032513/)) and links to
access.

See [`microbiome_ASV_data`](microbiome_ASV_data) file for details about
the data pre-processing. See Bla for recode.

Stage 2: Design
---------------

Matching explanation. Links to other options more “mainstream”.

See [`design`](design) file for the code of our pair matching
implementation and diagnostic plots.

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

\[Holle et al., 2005\] Holle R, Happich M, Löwel H, Wichmann HE;
[MONICA/KORA Study Group. KORA–a research platform for population based
health research.](https://pubmed.ncbi.nlm.nih.gov/16032513/)
Gesundheitswesen. 2005 Aug;67
