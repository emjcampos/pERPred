
<!-- README.md is generated from README.Rmd. Please edit that file -->

# pERPred <a href = 'https://github.com/emjcampos/pERPred'><img src = 'man/figures/pERPred.png' align = "right" height = "138.5" /></a>

<!-- badges: start -->

<!-- badges: end -->

The `pERPred` package is a tool for conducting ERP analyses using the
Principle ERP Reduction algorithm put forth in (cite our paper) by
Campos and Hazlett et al (2019). This package contains the functions
necessary to derive the principle ERPs discussed in that publication as
well as the tools for performing “pERP-space analysis”. Let this serve
as a guide to our approach to ERP analysis in `R`.

## Installation

You can install the released version of pERPred from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("pERPred")
```

And the development version from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("emjcampos/pERPred")
```

## Example

### Simulated Data Description

The generation model for the simulated data is described in the
publication as follows. The observed signal is assumed to be a linear
combination of the principal ERP components (pERPs). Let ![\\scriptstyle
Y\_{i,v,e}(t)](https://latex.codecogs.com/png.latex?%5Cscriptstyle%20Y_%7Bi%2Cv%2Ce%7D%28t%29
"\\scriptstyle Y_{i,v,e}(t)") denote the ERP signal observed for subject
![i](https://latex.codecogs.com/png.latex?i "i") at task
![v](https://latex.codecogs.com/png.latex?v "v") and electrode
![e](https://latex.codecogs.com/png.latex?e "e"). The data generation
model in the simulation will be: ![\\scriptstyle Y\_{i,v,e}(t) =
\\sum\_{c=1}^C k\_{c,v,e}\\phi^{\\star}\_c(t) + \\sum\_{c=1}^C
\\xi\_{c,i,v,e}\\phi^{\\star}\_c(t) +
\\zeta\_{i,v,e}(t).](https://latex.codecogs.com/png.latex?%5Cscriptstyle%20Y_%7Bi%2Cv%2Ce%7D%28t%29%20%3D%20%5Csum_%7Bc%3D1%7D%5EC%20k_%7Bc%2Cv%2Ce%7D%5Cphi%5E%7B%5Cstar%7D_c%28t%29%20%2B%20%5Csum_%7Bc%3D1%7D%5EC%20%5Cxi_%7Bc%2Ci%2Cv%2Ce%7D%5Cphi%5E%7B%5Cstar%7D_c%28t%29%20%2B%20%5Czeta_%7Bi%2Cv%2Ce%7D%28t%29.
"\\scriptstyle Y_{i,v,e}(t) = \\sum_{c=1}^C k_{c,v,e}\\phi^{\\star}_c(t) + \\sum_{c=1}^C \\xi_{c,i,v,e}\\phi^{\\star}_c(t) + \\zeta_{i,v,e}(t).")

In the simulation, ![\\scriptstyle N
= 100](https://latex.codecogs.com/png.latex?%5Cscriptstyle%20N%20%3D%20100
"\\scriptstyle N = 100"), ![\\scriptstyle V
= 9](https://latex.codecogs.com/png.latex?%5Cscriptstyle%20V%20%3D%209
"\\scriptstyle V = 9"), ![\\scriptstyle C
= 5](https://latex.codecogs.com/png.latex?%5Cscriptstyle%20C%20%3D%205
"\\scriptstyle C = 5"), ![\\scriptstyle T
= 384](https://latex.codecogs.com/png.latex?%5Cscriptstyle%20T%20%3D%20384
"\\scriptstyle T = 384"), and each subject has a common ![E
= 40](https://latex.codecogs.com/png.latex?E%20%3D%2040 "E = 40"). The
first term reflects that each task and electrode is composed of a
weighted average of pERPs, where ![\\scriptstyle
\\phi^{\\star}\_c(t)](https://latex.codecogs.com/png.latex?%5Cscriptstyle%20%5Cphi%5E%7B%5Cstar%7D_c%28t%29
"\\scriptstyle \\phi^{\\star}_c(t)") are the \`\`true’’ pERPs for
purposes of simulation. To simulate a single ![\\scriptstyle
\\phi^{\\star}\_c(t)](https://latex.codecogs.com/png.latex?%5Cscriptstyle%20%5Cphi%5E%7B%5Cstar%7D_c%28t%29
"\\scriptstyle \\phi^{\\star}_c(t)"), we draw a function that is a
superposition of Gaussian kernels centered at different time points to
produce a smooth, low-frequency function intended to mimic the type of
source signals we would expect to contribute to ERP waveforms. These
simulated signals are then rotated by ICA to form maximally independent
bases and are plotted
here.

<img src="man/figures/README-true_pERPs-1.png" width="60%" style="display: block; margin: auto;" />

The coefficients, ![\\scriptstyle
k\_{c,v,e}](https://latex.codecogs.com/png.latex?%5Cscriptstyle%20k_%7Bc%2Cv%2Ce%7D
"\\scriptstyle k_{c,v,e}"), are drawn from a normal distribution with
mean zero and variance ![\\scriptstyle
\\sigma^2\_{k}=0.25](https://latex.codecogs.com/png.latex?%5Cscriptstyle%20%5Csigma%5E2_%7Bk%7D%3D0.25
"\\scriptstyle \\sigma^2_{k}=0.25").

The second term reflects the structured stochastic canonical component
with dependencies within a subject across tasks and electrodes. The
subject specific scores ![\\scriptstyle
\\xi\_{c,i,v,e}](https://latex.codecogs.com/png.latex?%5Cscriptstyle%20%5Cxi_%7Bc%2Ci%2Cv%2Ce%7D
"\\scriptstyle \\xi_{c,i,v,e}") will be simulated from a Matrix Normal
(to create a depedence structure within tasks and across electrodes)
with mean zero and covariance matrices ![\\scriptstyle
\\Sigma\_{c,v}](https://latex.codecogs.com/png.latex?%5Cscriptstyle%20%5CSigma_%7Bc%2Cv%7D
"\\scriptstyle \\Sigma_{c,v}") and ![\\scriptstyle
\\Sigma\_{c,e}](https://latex.codecogs.com/png.latex?%5Cscriptstyle%20%5CSigma_%7Bc%2Ce%7D
"\\scriptstyle \\Sigma_{c,e}") of dimension ![\\scriptstyle V\\ast
V](https://latex.codecogs.com/png.latex?%5Cscriptstyle%20V%5Cast%20V
"\\scriptstyle V\\ast V") and ![\\scriptstyle E\\ast
E](https://latex.codecogs.com/png.latex?%5Cscriptstyle%20E%5Cast%20E
"\\scriptstyle E\\ast E") (identically and independently drawn over
subjects only). The covariance matrices are created in the same way for
the electrode and task dimensions. A matrix with 0.5 on the diagonal and
0.1 elsewhere created. Each of these is multiplied by a factor (0.1,
0.2, 0.3, 0.4, 0.5), so that each true ERP component is represented
differently. In the results, we expect the component with the highest
multiplier to explain the most variation in the observed signals.

The last term is the independent and identically distributed (iid)
measurement error. The measurement error ![\\scriptstyle
\\zeta\_{i,v,e}(t)](https://latex.codecogs.com/png.latex?%5Cscriptstyle%20%5Czeta_%7Bi%2Cv%2Ce%7D%28t%29
"\\scriptstyle \\zeta_{i,v,e}(t)") will be simulated as independent
draws (over ![i](https://latex.codecogs.com/png.latex?i "i"),
![v](https://latex.codecogs.com/png.latex?v "v"), and
![e](https://latex.codecogs.com/png.latex?e "e")) from a normal
distribution with mean zero and variance
![\\scriptstyle\\sigma\_\\text{error}^2](https://latex.codecogs.com/png.latex?%5Cscriptstyle%5Csigma_%5Ctext%7Berror%7D%5E2
"\\scriptstyle\\sigma_\\text{error}^2") equal to some proportion of the
signal variance. Define the variance of the simulated signal terms
![\\scriptstyle\\sum\_{c=1}^C k\_{c,v,e}\\phi^{\\star}\_c(t) +
\\sum\_{c=1}^C
\\xi\_{c,i,v,e}\\phi^{\\star}\_c(t)](https://latex.codecogs.com/png.latex?%5Cscriptstyle%5Csum_%7Bc%3D1%7D%5EC%20k_%7Bc%2Cv%2Ce%7D%5Cphi%5E%7B%5Cstar%7D_c%28t%29%20%2B%20%5Csum_%7Bc%3D1%7D%5EC%20%5Cxi_%7Bc%2Ci%2Cv%2Ce%7D%5Cphi%5E%7B%5Cstar%7D_c%28t%29
"\\scriptstyle\\sum_{c=1}^C k_{c,v,e}\\phi^{\\star}_c(t) + \\sum_{c=1}^C \\xi_{c,i,v,e}\\phi^{\\star}_c(t)")
to be
![\\scriptstyle\\sigma\_\\text{signal}^2](https://latex.codecogs.com/png.latex?%5Cscriptstyle%5Csigma_%5Ctext%7Bsignal%7D%5E2
"\\scriptstyle\\sigma_\\text{signal}^2"). For the purposes of this
vignette, we will only explore the low noise case in which 1/3 of the
simulated ERP will be noise so ![\\scriptstyle\\sigma^2\_\\text{error} =
\\sqrt{0.5}\\times\\sigma\_\\text{signal}](https://latex.codecogs.com/png.latex?%5Cscriptstyle%5Csigma%5E2_%5Ctext%7Berror%7D%20%3D%20%5Csqrt%7B0.5%7D%5Ctimes%5Csigma_%5Ctext%7Bsignal%7D
"\\scriptstyle\\sigma^2_\\text{error} = \\sqrt{0.5}\\times\\sigma_\\text{signal}").

The data to be used in the `pERPred` function must be averaged over
trials. Then there must be one column for the “Task”, “Subject”, and
“Time”, and the remaining columns are the named electrodes that were
observed. All subjects must have the same number of time points observed
as well as the same number of tasks. However, there may be a different
number of electrodes (use NA for the missing electrodes). For example,
our simulated data looks like this:

``` r
head(simulated_data[, 1:7])
#>     Task     Subject      Time        AF3       AF4        AFZ         C3
#> 1 Task_1 Subject_001 0.0026042 -0.8315047 1.6557474 -1.4637006 -1.2174339
#> 2 Task_1 Subject_001 0.0052083  0.9028223 1.1525489 -2.3697299 -3.4455528
#> 3 Task_1 Subject_001 0.0078125 -0.0770906 0.5224889 -0.8622194 -0.6448022
#> 4 Task_1 Subject_001 0.0104167 -0.2268654 0.0654785 -2.8212200 -3.6817255
#> 5 Task_1 Subject_001 0.0130208 -1.1375707 1.3043730 -1.1123895 -1.5617383
#> 6 Task_1 Subject_001 0.0156250 -1.8166604 0.7649564 -0.9293482 -2.6024852
```

### pERP-RED

Begin by loading the `pERPred` package.

``` r
library(pERPred)
```

Now using the `pERPred` function, you can estimate the pERPs. There are
two choices that are chosen by the user:  
1\. the number of pERPs to estimate; and  
2\. the proportion of variation each of the PCA steps must explain.

By default, the percent of variation is set to 80. Raising this value
means more components will be kept in each of the PCA steps and the
computation time will rise accordingly. Also, the number of pERPs is
chosen based on an
![\\scriptstyleR^2](https://latex.codecogs.com/png.latex?%5CscriptstyleR%5E2
"\\scriptstyleR^2") value defined as: ![\\scriptstyle R\_{test}^2 = 1 -
\\frac{|| Y\_j - \\Phi \\beta\_j ||\_{\\mathcal{F}}^2}{|| Y\_j
||\_{\\mathcal{F}}^2}](https://latex.codecogs.com/png.latex?%5Cscriptstyle%20R_%7Btest%7D%5E2%20%3D%201%20-%20%5Cfrac%7B%7C%7C%20Y_j%20-%20%5CPhi%20%5Cbeta_j%20%7C%7C_%7B%5Cmathcal%7BF%7D%7D%5E2%7D%7B%7C%7C%20Y_j%20%7C%7C_%7B%5Cmathcal%7BF%7D%7D%5E2%7D
"\\scriptstyle R_{test}^2 = 1 - \\frac{|| Y_j - \\Phi \\beta_j ||_{\\mathcal{F}}^2}{|| Y_j ||_{\\mathcal{F}}^2}")
and can be calculated using the `R2_test` function.

When estimating the pERPs, split the data into a training and test set,
then calculate the ![\\scriptstyle
R^2\_{test}](https://latex.codecogs.com/png.latex?%5Cscriptstyle%20R%5E2_%7Btest%7D
"\\scriptstyle R^2_{test}"). After the number of pERPs is chosen, you
can re-estimate the pERPs using the whole dataset. If preferred, it is
possible to parallelize the estimation to speed up this process using
the `future` package (see commented code).

``` r
subject_list <- unique(simulated_data$Subject)
electrode_list <- names(simulated_data)[-c(1:3)]
task_list <- unique(simulated_data$Task)
train_subjects <- 
  sort(sample(subject_list, round(2 * length(subject_list) / 3)))

pERPs <- map(3:7, 
             ~ pERPred(simulated_data[simulated_data$Subject %in%
                                        train_subjects, ],
                       num_pERPs = .x,
                       percent_variation = 80))

# library(future)
# plan(multiprocess) 
# pERPs <- future_map(3:7, 
#                     ~ pERPred(simulated_data[simulated_data$Subject %in% 
#                                                train_subjects, ], 
#                               num_pERPs = .x, 
#                               percent_variation = 80), 
#                     .progress = TRUE)
```

Based on the ![\\scriptstyle
R^2\_{test}](https://latex.codecogs.com/png.latex?%5Cscriptstyle%20R%5E2_%7Btest%7D
"\\scriptstyle R^2_{test}") values, we would choose 5 as the appropriate
number of pERPs since it is the point at which the ![\\scriptstyle
R^2\_{test}](https://latex.codecogs.com/png.latex?%5Cscriptstyle%20R%5E2_%7Btest%7D
"\\scriptstyle R^2_{test}") tapers off.

``` r
map_dfr(pERPs, ~R2_test(simulated_data, .x)) 
```

    #> # A tibble: 5 x 2
    #>      R2 pERPs
    #>   <dbl> <int>
    #> 1 0.468     3
    #> 2 0.575     4
    #> 3 0.671     5
    #> 4 0.672     6
    #> 5 0.673     7

``` r
pERPs5 <- pERPred(simulated_data,
                 num_pERPs = 5,
                 percent_variation = 80)
```

``` r
pERPs5 %>%
  mutate(Time = unique(simulated_data$Time)) %>%
  gather(pERP, Amplitude, -Time) %>%
  ggplot() +
  geom_line(aes(x = Time, y = Amplitude)) +
  facet_wrap(~ pERP, ncol = 2)
```

<img src="man/figures/README-plot_pERPs-1.png" width="60%" style="display: block; margin: auto;" />

### pERP-space Analysis

There are many analyses we can perform using these pERPs. In the paper,
we discuss a few. First, we want to compute the weights ![\\scriptstyle
\\omega\_j](https://latex.codecogs.com/png.latex?%5Cscriptstyle%20%5Comega_j
"\\scriptstyle \\omega_j") for each record using the `pERP_scorer`
function.

``` r
individual_scores <- pERP_scorer(simulated_data, pERPs5)
```

``` r
head(individual_scores)
#> # A tibble: 6 x 8
#> # Groups:   Task, Subject, Electrode [2]
#>   Task   Subject     Electrode term   estimate std.error statistic  p.value
#>   <chr>  <chr>       <chr>     <chr>     <dbl>     <dbl>     <dbl>    <dbl>
#> 1 Task_1 Subject_001 AF3       pERP …   0.562     0.0456    12.3   1.64e-29
#> 2 Task_1 Subject_001 AF3       pERP …  -0.153     0.0456    -3.35  8.83e- 4
#> 3 Task_1 Subject_001 AF3       pERP …  -0.383     0.0456    -8.38  1.02e-15
#> 4 Task_1 Subject_001 AF3       pERP …   0.682     0.0456    15.0   4.69e-40
#> 5 Task_1 Subject_001 AF3       pERP …   0.734     0.0456    16.1   9.10e-45
#> 6 Task_1 Subject_001 AF4       pERP …   0.0114    0.0480     0.238 8.12e- 1
```

Now we can reconstruct each record using their weights and the pERPs.
Here is an example using subject 1, electrode CZ, task 1.

``` r
record <- simulated_data %>% 
  filter(Subject == "Subject_001", 
         Task == "Task_1") %>% 
  select(Task, Time, "CZ")
scores <- individual_scores %>% 
  filter(Subject == "Subject_001", 
         Electrode == "CZ", 
         Task == "Task_1") %>% 
  pull(estimate)
projected <- data.frame("Projected" = as.matrix(pERPs5) %*% 
                          as.matrix(scores)) %>% 
  as.data.frame() %>% 
  mutate(Time = record$Time) %>% 
  mutate(Projected = Projected + mean(record$CZ))

ggplot() + 
  geom_line(data = record, 
            aes(x = Time, y = CZ, color = "Observed")) + 
  geom_line(data = projected, 
            aes(x = Time, y = Projected, color = "Reconstructed")) +
  theme(legend.title = element_blank()) + 
  labs(y = TeX("$\\mu V$"))
```

<img src="man/figures/README-individual_reconstruction-1.png" width="60%" style="display: block; margin: auto;" />

Keeping in mind that this is simulated data, we can visualize how these
weights are represented across the scalp for all subjects using the
`coefficient_headmap` function. Create a dataframe of the average scores
where the column of averages is named `average`. Split the dataframe by
Task to create a list of dataframes. This will allow the scale to remain
the same over all tasks. If you want the scale to stay the same over
only selected Tasks, only include those tasks in the list of dataframes.

``` r
average_scores <- individual_scores %>% 
  group_by(Task, Electrode, term) %>% 
  summarise(average = mean(estimate, na.rm = TRUE)) %>% 
  split(.$Task)

# keep scale the same across all tasks
coefficient_headmap("Task_1", average_scores)
```

<img src="man/figures/README-headmaps-1.png" width="60%" style="display: block; margin: auto;" />

``` r

# keep scale the same across only tasks 1 and 2
coefficient_headmap("Task_1", average_scores[c("Task_1", "Task_2")])
```

<img src="man/figures/README-headmaps-2.png" width="60%" style="display: block; margin: auto;" />

Here is an example of how one would conduct tests for differences among
groups. We will test if there is a difference in scores for Task\_1 and
Task\_2 at CZ.

``` r
pERP_difference(scores = individual_scores,
                electrode = "CZ",
                task1 = "Task_1",
                task2 = "Task_2")
#>      pERP Task_1 Overall Mean Task_1 Overall APSD Task_1 Overall SE
#> 1 pERP 01           0.3421022           0.3522178        0.03522178
#> 2 pERP 02           0.2765889           0.2360526        0.02360526
#> 3 pERP 03          -0.2785253           0.2056983        0.02056983
#> 4 pERP 04           0.2168410           0.3839734        0.03839734
#> 5 pERP 05          -0.7411287           0.3180623        0.03180623
#>   Task_1 Overall t Task_1 Overall N Task_2 Overall Mean
#> 1         9.712802              100         -0.01720551
#> 2        11.717256              100         -0.74303513
#> 3       -13.540475              100         -0.53359646
#> 4         5.647292              100         -0.26786325
#> 5       -23.301369              100          0.38377263
#>   Task_2 Overall APSD Task_2 Overall SE Task_2 Overall t Task_2 Overall N
#> 1           0.3367229        0.03367229       -0.5109693              100
#> 2           0.2540876        0.02540876      -29.2432679              100
#> 3           0.2073085        0.02073085      -25.7392477              100
#> 4           0.3917388        0.03917388       -6.8378031              100
#> 5           0.3675706        0.03675706       10.4407882              100
#>   Mean Difference SE Difference          t signif
#> 1       0.3593077    0.04872778   7.373774      *
#> 2       1.0196240    0.03468160  29.399566      *
#> 3       0.2550712    0.02920421   8.734053      *
#> 4       0.4847042    0.05485388   8.836280      *
#> 5      -1.1249013    0.04860779 -23.142408      *
```

If we had another grouping variable, such as diagnosis, we could add
that to the individual scores dataframe.

``` r
set.seed(1234)
group <- data.frame(
  Subject = distinct(simulated_data, Subject),
  group_member = sample(
    c("Control", "Treatment"),
    100,
    replace = TRUE,
    prob = c(.5, .5)
  )
)

individual_scores_groups <- full_join(individual_scores, group, by = "Subject")

pERP_difference(scores = individual_scores_groups,
                electrode = "CZ",
                task1 = "Task_1",
                group1 = "Control", 
                group2 = "Treatment")
#>      pERP Task_1 Control Mean Task_1 Control APSD Task_1 Control SE
#> 1 pERP 01           0.3023001           0.3446887        0.05138316
#> 2 pERP 02           0.2997375           0.2228006        0.03321315
#> 3 pERP 03          -0.2922279           0.2063993        0.03076820
#> 4 pERP 04           0.2160037           0.3539055        0.05275711
#> 5 pERP 05          -0.7609320           0.3456268        0.05152300
#>   Task_1 Control t Task_1 Control N Task_1 Treatment Mean
#> 1         5.883251               45             0.3746675
#> 2         9.024661               45             0.2576492
#> 3        -9.497727               45            -0.2673141
#> 4         4.094306               45             0.2175260
#> 5       -14.768784               45            -0.7249260
#>   Task_1 Treatment APSD Task_1 Treatment SE Task_1 Treatment t
#> 1             0.3580919          0.04828510           7.759486
#> 2             0.2467701          0.03327447           7.743148
#> 3             0.2063382          0.02782264          -9.607788
#> 4             0.4101740          0.05530784           3.933005
#> 5             0.2958655          0.03989450         -18.171076
#>   Task_1 Treatment N Mean Difference SE Difference           t signif
#> 1                 55    -0.072367483    0.07051014 -1.02634145       
#> 2                 55     0.042088294    0.04701387  0.89523136       
#> 3                 55    -0.024913886    0.04148230 -0.60059070       
#> 4                 55    -0.001522269    0.07643475 -0.01991593       
#> 5                 55    -0.036006055    0.06516280 -0.55255539
```

The user can insert these tables directly into their manuscript using
`kable`, `xtable` or any other of a number of tools in `R`.
