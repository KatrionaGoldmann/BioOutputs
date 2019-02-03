-   [BioOutputs](#biooutputs)
    -   [Gallery](#gallery)
    -   [Required packages](#required-packages)
    -   [Install Package](#install-package)
    -   [bio\_corr](#bio_corr)
    -   [bio\_frequency](#bio_frequency)
    -   [bio\_volcano](#bio_volcano)
    -   [bio\_bires](#bio_bires)
    -   [bio\_mods](#bio_mods)
    -   [bio\_treatmentGroups](#bio_treatmentgroups)

------------------------------------------------------------------------

BioOutputs
==========

This package contains common R scripts I use in my day to day data
analysis of biological data. The scripts are primarily for plotting and
visualisation, with some data organisation thrown in as well.

------------------------------------------------------------------------

Gallery
-------

<table border="1">
  <tr>
    <td align="center" style="vertical-align:top" height="200"><a href="#bires">Bi-Results</a><img src= ./figs/bio_bires.png  height="150" width="330"/></td>
    <td align="center" style="vertical-align:top" height="200"><a href="#corrp">Correlation Plot</a><img src= ./figs/bio_corr.png  height="150" width="330"/></td>
    <td style="vertical-align:top" align="center" height="200"><a href="#freq">Frequency Table</a><img src= ./figs/bio_freq.png  height="50" width="330"/></td>
  </tr>
<tr>
    <td align="center" style="vertical-align:top" height="200"><a href="#mods">Module Plot</a><img src= ./figs/bio_mods.png  height="150" width="330"/></td>
    <td align="center" style="vertical-align:top" height="200"><a href="#tgcomp">Treatment Group Comparisons</a><img src= ./figs/bio_treatmentGroups.png  height="150" width="330"/></td>
    <td align="center" style="vertical-align:top" height="200"><a href="#volc">Volcano Plot</a><img src= ./figs/bio_volcano.png  height="150" width="330"/></td>
  </tr>
</table>

------------------------------------------------------------------------

Required packages
-----------------

[`dplyr`](https://dplyr.tidyverse.org)
[`ggplot2`](https://ggplot2.tidyverse.org)
[`ggrepel`](https://cran.r-project.org/web/packages/ggrepel/vignettes/ggrepel.html)
[`ComplexHeatmap`](http://bioconductor.org/packages/release/bioc/html/ComplexHeatmap.html)
[`RColorBrewer`](https://www.rdocumentation.org/packages/RColorBrewer/versions/1.1-2/topics/RColorBrewer)
[`Rmisc`](https://www.rdocumentation.org/packages/Rmisc/versions/1.5)
[`ggpubr`](https://cran.r-project.org/web/packages/ggpubr/index.html)
[`grid`](https://www.rdocumentation.org/packages/grid/versions/3.5.2)
[`pBrackets`](https://cran.r-project.org/web/packages/pBrackets/index.html)

------------------------------------------------------------------------

Install Package
---------------

First install devtools to allow installation from gitub and any other
required packages.

    install.packages("devtools")
    library("devtools")

    ## Warning: package 'cowplot' was built under R version 3.5.2

    ## Loading required package: ggplot2

    ## 
    ## Attaching package: 'cowplot'

    ## The following object is masked from 'package:ggplot2':
    ## 
    ##     ggsave

Now install the BioOutputs package.

    install_github("KatrionaGoldmann/BioOutputs")
    library("BioOutputs")

------------------------------------------------------------------------

<a id="corrp"></a>

bio\_corr
---------

Create a correlation plot. Taken from kassambara/ggpubr just changed the
default arguments.

So we can use the classic example with the *mtcars* data frames:

    kable(head(mtcars))

<table>
<thead>
<tr class="header">
<th></th>
<th align="right">mpg</th>
<th align="right">cyl</th>
<th align="right">disp</th>
<th align="right">hp</th>
<th align="right">drat</th>
<th align="right">wt</th>
<th align="right">qsec</th>
<th align="right">vs</th>
<th align="right">am</th>
<th align="right">gear</th>
<th align="right">carb</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td>Mazda RX4</td>
<td align="right">21.0</td>
<td align="right">6</td>
<td align="right">160</td>
<td align="right">110</td>
<td align="right">3.90</td>
<td align="right">2.620</td>
<td align="right">16.46</td>
<td align="right">0</td>
<td align="right">1</td>
<td align="right">4</td>
<td align="right">4</td>
</tr>
<tr class="even">
<td>Mazda RX4 Wag</td>
<td align="right">21.0</td>
<td align="right">6</td>
<td align="right">160</td>
<td align="right">110</td>
<td align="right">3.90</td>
<td align="right">2.875</td>
<td align="right">17.02</td>
<td align="right">0</td>
<td align="right">1</td>
<td align="right">4</td>
<td align="right">4</td>
</tr>
<tr class="odd">
<td>Datsun 710</td>
<td align="right">22.8</td>
<td align="right">4</td>
<td align="right">108</td>
<td align="right">93</td>
<td align="right">3.85</td>
<td align="right">2.320</td>
<td align="right">18.61</td>
<td align="right">1</td>
<td align="right">1</td>
<td align="right">4</td>
<td align="right">1</td>
</tr>
<tr class="even">
<td>Hornet 4 Drive</td>
<td align="right">21.4</td>
<td align="right">6</td>
<td align="right">258</td>
<td align="right">110</td>
<td align="right">3.08</td>
<td align="right">3.215</td>
<td align="right">19.44</td>
<td align="right">1</td>
<td align="right">0</td>
<td align="right">3</td>
<td align="right">1</td>
</tr>
<tr class="odd">
<td>Hornet Sportabout</td>
<td align="right">18.7</td>
<td align="right">8</td>
<td align="right">360</td>
<td align="right">175</td>
<td align="right">3.15</td>
<td align="right">3.440</td>
<td align="right">17.02</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">3</td>
<td align="right">2</td>
</tr>
<tr class="even">
<td>Valiant</td>
<td align="right">18.1</td>
<td align="right">6</td>
<td align="right">225</td>
<td align="right">105</td>
<td align="right">2.76</td>
<td align="right">3.460</td>
<td align="right">20.22</td>
<td align="right">1</td>
<td align="right">0</td>
<td align="right">3</td>
<td align="right">1</td>
</tr>
</tbody>
</table>

    bio_corr(mtcars, "qsec", "wt")

    ## Loading required package: bitops

![](README_files/figure-markdown_strict/bio_corr-1.png)

------------------------------------------------------------------------

<a id="freq"></a>

bio\_frequency
--------------

The *bio\_frequency()* function generates a frequency table from factor
or character vector columns in a data frame. This has the following
arguments:

<table style="width:100%;">
<colgroup>
<col width="27%" />
<col width="72%" />
</colgroup>
<thead>
<tr class="header">
<th>Argument</th>
<th></th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td>data</td>
<td>A data frame containing columns to be counted</td>
</tr>
<tr class="even">
<td>columns</td>
<td>Column names or indices to be counted in data</td>
</tr>
<tr class="odd">
<td>freq.percent</td>
<td>Whether the table should include frequency counts, percentages or both (options = c(&quot;freq&quot;, &quot;percent&quot;, &quot;both&quot;)). Default=&quot;both&quot;</td>
</tr>
<tr class="even">
<td>include.na</td>
<td>Include NA values (options are TRUE/FALSE, default=TRUE)</td>
</tr>
<tr class="odd">
<td>remove.vars</td>
<td>Character vector of variables not to be included in the counts (e.g. remove.vars = c(&quot;&quot;) remove blanks from the count)</td>
</tr>
</tbody>
</table>

Then if we want to see the breakdown of, say, the gear column in mtcars
we can apply:

    kable(bio_frequency(mtcars, "gear"))

<table>
<thead>
<tr class="header">
<th></th>
<th align="left">3</th>
<th align="left">4</th>
<th align="left">5</th>
<th align="left">Total</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td>gear</td>
<td align="left">15 (47%)</td>
<td align="left">12 (38%)</td>
<td align="left">5 (16%)</td>
<td align="left">n = 32</td>
</tr>
</tbody>
</table>

And if wanted we can remove one variable from the table. This is useful
if we have unknowns or the likes.

    kable(bio_frequency(mtcars, "gear", remove.vars=c("5")))

<table>
<thead>
<tr class="header">
<th></th>
<th align="left">3</th>
<th align="left">4</th>
<th align="left">Total</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td>gear</td>
<td align="left">15 (56%)</td>
<td align="left">12 (44%)</td>
<td align="left">n = 27</td>
</tr>
</tbody>
</table>

------------------------------------------------------------------------

<a id="volc"></a>

bio\_volcano
------------

This function generates a volcano plot from a top table using ggplot.
The function contains many parameters, use `?bio_volcano` to interogate.

Lets look at the ALS patients carrying the C9ORF72 data set

<!-- ```{r, message=FALSE, warning=FALSE, eval=FALSE} -->
<!-- BiocManager::install("leukemiasEset", version = "3.8") -->
<!-- ``` -->
<!-- ```{r, message=FALSE, warning=FALSE} -->
<!-- library(leukemiasEset) -->
<!-- library(limma) -->
<!-- ``` -->
<!-- ```{r, message=FALSE, warning=FALSE} -->
<!-- data(leukemiasEset) -->
<!-- ourData <- leukemiasEset[, leukemiasEset$LeukemiaType %in% c("ALL", "NoL")] -->
<!-- ourData$LeukemiaType <- factor(ourData$LeukemiaType) -->
<!-- design <- model.matrix(~ ourData$LeukemiaType) -->
<!-- fit <- lmFit(ourData, design) -->
<!-- fit <- eBayes(fit) -->
<!-- toptable <- topTable(fit) -->
<!-- toptable$pvalue = toptable$P.Value -->
<!-- ``` -->
    toptable <- read.table("https://gist.githubusercontent.com/stephenturner/806e31fce55a8b7175af/raw/1a507c4c3f9f1baaa3a69187223ff3d3050628d4/results.txt", sep="", header=T)
    rownames(toptable) = toptable$Gene

    bio_volcano(toptable, fc.col="log2FoldChange", label.row.indices=1:10, main="ALS patients carrying the C9ORF72", add.lines=TRUE)

![](README_files/figure-markdown_strict/bio_volcano-1.png)

------------------------------------------------------------------------

<a id="bires"></a>

bio\_bires
----------

This function colours data according to whether it is below or above a
defined plane. The plane is plotted as a line and data can be output as
either a line or markers/points.

    data(beavers)
    df.plane = beaver1
    df.data = beaver2
    df.plane$temp <- df.plane$temp + 0.5
    bio_bires(x="time", y="temp", df.data, x.plane="time", y.plane="temp", df.plane)

![](README_files/figure-markdown_strict/bires-1.png)

<a id="mods"></a>

bio\_mods
---------

This function splits expression data into customisable modules and
averages over catagories in a given variable. In this example we will
look at two Li modules and a custom one I made up.

    ## Warning: namespace 'DOSE' is not available and has been replaced
    ## by .GlobalEnv when processing object '.DOSEEnv'

    exp = rld.syn
    meta <- rld.metadata.syn
    mean.var = "Pathotype"
    cols = NULL
    cluster.rows=FALSE
    bio_mods(exp=exp, mod.list=mod_list, meta=meta, mean.var = "Pathotype", show.names=TRUE, cluster.rows=FALSE)

![](README_files/figure-markdown_strict/mods-1.png)

<a id="tgcomp"></a>

bio\_treatmentGroups
--------------------

    library(survival)
    data(kidney)

    df <- kidney
    y.col="frail"
    group.col="status"
    x.col="time"
    id.col="id"
    df <- df[df$time < 150, ]
    bio_treatmentGroups(df, group.col="status", y.col="frail", x.col="time")

    ## Warning in checkConv(attr(opt, "derivs"), opt$par, ctrl =
    ## control$checkConv, : Model failed to converge with max|grad| = 0.24405 (tol
    ## = 0.002, component 1)

    ## Warning in checkConv(attr(opt, "derivs"), opt$par, ctrl = control$checkConv, : Model is nearly unidentifiable: very large eigenvalue
    ##  - Rescale variables?

    ## Warning in as_lmerModLT(model, devfun): Model may not have converged with 1
    ## eigenvalue close to zero: 6.1e-11

    ## Warning in qt(conf.interval/2 + 0.5, datac$N - 1): NaNs produced

    ## Warning: Removed 41 rows containing missing values (geom_errorbar).

![](README_files/figure-markdown_strict/unnamed-chunk-11-1.png)
