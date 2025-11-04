# BOBgmms-examples

This repository contains source code to reproduce the results from the article "[BOB: Bayesian optimized bootstrap for approximate posterior sampling in Gaussian mixture models](https://doi.org/10.1007/s11222-025-10763-y)" (Marin et al., 2026).

The `BOBgmms` R package is available at the Github repository: [https://github.com/marinsantiago/BOBgmms](https://github.com/marinsantiago/BOBgmms)

## Contents

- [Folder Structure ](#folders)
- [Software Requirements](#software)
- [System Requirements](#system)
- [Data](#data)
- [Running the Scripts](#run)
- [Citation](#cite)
- [References](#refs)

## <a name="folders"></a> Folder Structure 
```BOBgmms-examples```

```├──``` `R/`: R scripts with helper functions required throughout the analysis.

```├──``` `simulations/`: Source code to reproduce the results in the *Simulations* section from the article.

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; ```├──``` `sim_ar_neg_05/`: Settings with covariance matrices of the form $({\Sigma}_{k})_{j,j'} = \{-0.5^{|j-j'|}\}_{j,j'=1}^{d}$.
    
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; ```├──``` `sim1/`: First simulation setting.

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; ```├──``` `sim2/`: Second simulation setting.

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; ```├──``` `sim3/`: Third simulation setting.

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; ```├──``` `sim4/`: Fourth simulation setting.

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; ```├──``` `sim5/`: Fifth simulation setting.

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; ```├──``` `sim6/`: Sixth simulation setting.

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; ```├──``` `sim7/`: Seventh simulation setting.

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; ```├──``` `sim8/`: Eighth simulation setting.

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; ```├──``` `sim9/`: Ninth simulation setting.

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; ```└──``` `sim-results/`: Code to reproduce the plots and tables with the results from the nine experimental settings.
    
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; ```├──``` `sim-illustrative/`: Code to reproduce the simulation results used as an illustrative example.

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; ```├──``` `sim-varying-n/`: Code to reproduce the simulation results with a varying sample size.

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; ```├──``` `sim-trivial-init/`: Code to reproduce the simulation results with a trivial initialization.

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; ```└──``` `sim-informative-alpha/`: Simulation settings with a strong informative prior on the mixture proportions.

```├──``` `data-analyses/`: Source code to reproduce the results from the real-world data analyses.

  &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; ```├──``` `wine/`: Source code to reproduce the results from the analysis of the *wines* data.

  &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; ```└──``` `kernels/`: Source code to reproduce the results from the analysis of the *seeds* data.
    
```├──``` `stan/`: Stan code used to fit a $K$ component Bayesian Gaussian mixture, assuming conjugate priors.

```└──``` `one-dimensional-densities/`: Source code to generate one-dimensional density plots. 

## <a name="software"></a> Software Requirements

We run all our simulations and data analyses in R (version 4.3.2) via Rstudio (version 2024.12.0.467). Additionally, we make use of the following R packages (along with their dependencies):

  - `BOBgmms` v0.1.0 (available at: [https://github.com/marinsantiago/BOBgmms](https://github.com/marinsantiago/BOBgmms))
  - `clusterHD` v1.0.2
  - `LaplacesDemon` v16.1.6
  - `MASS` v7.3-60.0.1
  - `mclust` v6.1
  - `mvnfast` v0.2.8
  - `parallel` v4.3.2
  - `pbmcapply` v1.5.1
  - `RColorBrewer` v 1.1-3
  - `rstan` v 2.32.6 (for instructions on how to download and install `rstan`, please see [RStan Quick Start Guide](https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started))
  - `sparcl` v1.0.4
  
Before downloading and installing `rstan`, you need to set up your R installation to be able to compile C++ code. For additional details and instructions, see [RStan Getting Started](https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started) under **Configuring C++ Toolchain**. 

## <a name="system"></a> System Requirements

Parallelization over multiple CPU workers is conducted via *forking* 
(rather than *sockets*), so it is only available on POSIX systems (e.g., macOS, 
Linux, Unix, BSD), not on Windows. On non-POSIX platforms (such as Windows), 
the scripts will still run, but the number of CPU workers will be automatically 
set to one. For further details, see 
[Package '`parallel`'](https://stat.ethz.ch/R-manual/R-devel/library/parallel/doc/parallel.pdf). 

One can verify the OS type by running the following R code:

``` r
.Platform$OS.type
```

We run all our code on macOS (version 15.3.1). 

## <a name="run"></a> Running the Scripts

Step-by-step instructions on how to reproduce the results from the article:

#### Simulation Results

  - To reproduce the simulation results from the article, start by downloading the 
  folder `simulations/sim_ar_neg_05/setting1/` and execute the scripts 
  `sim1_data.R`, `sim1_rw.R`, `sim1_stan.R`, and `sim1_results.R` (in that order). 
  It is recommended to restart the R session each time you execute a script in 
  order to terminate any remaining "zombie processes" from the parallelization 
  exercise. Then, repeat these previous steps, but with the folders 
  `simulations/sim_ar_neg_05/setting2/`-to-`simulations/sim_ar_neg_05/setting9/`. 
  To obtain plots and tables summarizing the results from the nine experimental 
  settings, execute the scripts `boxplots.R`, `optimal_x_vals.R`, and  `table_medians.R` 
  from the `simulations/sim_ar_neg_05/sim-results/` folder.

  - To reproduce the results from our illustrative example, download the folder 
  `simulations/sim-illustrative/` and execute the script `sim_illustrative.R`.

  - To reproduce the simulation results with a varying sample sizes, download 
  the folder `simulations/sim-varying-n/` and execute the scripts `sim_n50.R`, 
  `sim_n125.R`, `sim_n250.R`, `sim_n375.R`, `sim_n500.R`, and  `sim_n_results.R` 
  (in that order). Once again, it is recommended to restart the R session each 
  time you execute a script to terminate any remaining "zombie processes" from 
  the parallelization exercise. 

  - To reproduce the simulation results with a trivial initialization, download 
  the folder `simulations/sim-trivial-init/` and execute the scripts 
  `sim_d5_trivial_init.R`, `sim_d10_trivial_init.R`, `sim_d15_trivial_init.R`, 
  and `simtrivial_init_results.R` (in that order). Once again, it is recommended 
  to restart the R session each time you execute a script to terminate any 
  remaining "zombie processes" from the parallelization exercise. 

  - To reproduce the simulation results with a strong informative prior on the 
  mixture proportions, download the folder `sim-informative-alpha/` and execute 
  the scripts `sim_alphas_d5.R`, `sim_alphas_d10.R`, `sim_alphas_d15.R`, and 
  `sim_alphas_results.R` (in that order). Once again, it is recommended to 
  restart the R session each time you execute a script to terminate any remaining 
  "zombie processes" from the parallelization exercise. 

#### Results from Benchmark Data Analysis

  - To reproduce the results from the analysis of benchmark data, start by 
  downloading the folder `data-analyses/wine/`  and execute the script `wines.R`. 
  Then, execute the script `kernels.R` from the folder `data-analyses/kernels/`.

## <a name="data"></a> Data

We make use of the [*wine*](https://archive.ics.uci.edu/dataset/109/wine) (Aeberhard and Forina, 1991) and [*seeds*](https://archive.ics.uci.edu/dataset/236/seeds) (Charytanowicz et al., 2012) data, which are publicly available from the [UC Irvine Machine Learning Repository](https://archive.ics.uci.edu/).

## <a name="cite"></a> Citation

If you use any part of this code in your work, please consider citing our *Statistics and Computing* paper:

```
@article{marin_bob,
  title   = {BOB: Bayesian optimized bootstrap for approximate posterior sampling in Gaussian mixture models},
  author  = {Santiago Marin and Bronwyn Loong and Anton H. Westveld},
  journal = {Statistics and Computing},
  volume  = {36},
  pages   = {14},
  year    = {2026},
  doi     = {10.1007/s11222-025-10763-y}
}
```

## <a name="refs"></a> References

Aeberhard, S., and Forina, M. (1991), "Wine." *UCI Machine Learning Repository*. https://doi.org/10.24432/C5PC7J.

Charytanowicz, M., Niewczas, J., Kulczycki, P., Kowalski, P., and Lukasik, S. (2012), "Seeds." *UCI Machine Learning Repository*. https://doi.org/10.24432/C5H30K.

Marin, S., Loong, B., and Westveld, A. H. (2026), "BOB: Bayesian optimized bootstrap for approximate posterior sampling in Gaussian mixture models." *Statistics and Computing*, **36**, 14. [doi:10.1080/10618600.2025.2572327](https://doi.org/10.1007/s11222-025-10763-y)

