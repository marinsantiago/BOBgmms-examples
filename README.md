# BOBgmms-examples

This repository contains source code to reproduce the results from the article "[BOB: Bayesian Optimized Bootstrap for Approximate Posterior Sampling in Gaussian Mixture Models](https://arxiv.org/abs/2311.03644)" (Marin, Loong and Westveld, 2025+).

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

- `R/`: R scripts with functions required throughout the analysis.
- `simulations/`: Source code to reproduce the results in the *Simulations* section from the article.
  - `sim_ar_neg_05/`: Source code to reproduce the simulation results, assuming covariance matrices of the form $`(\boldsymbol{\Sigma}_{k})_{j,j'} = \{-0.5^{|j-j'|}\}_{j,j'=1}^{d}`$, for all $k\in\{1,\dots,K\}$.
    - `sim1/`: Source code to reproduce the results for the 1st simulation setting.
    - `sim2/`: Source code to reproduce the results for the 2nd simulation setting.
    - `sim3/`: Source code to reproduce the results for the 3rd simulation setting.
    - `sim4/`: Source code to reproduce the results for the 4th simulation setting.
    - `sim5/`: Source code to reproduce the results for the 5th simulation setting.
    - `sim6/`: Source code to reproduce the results for the 6th simulation setting.
    - `sim7/`: Source code to reproduce the results for the 7th simulation setting.
    - `sim8/`: Source code to reproduce the results for the 8th simulation setting.
    - `sim9/`: Source code to reproduce the results for the 9th simulation setting.
    - `sim-results/`: Source code to generate plots and tables with results from the above 9 simulation studies.
  - `sim-illustrative/`: Source code to reproduce the simulation results used as an illustrative example.
  - `sim-varying-n/`: Source code to reproduce the simulation results with a varying sample size.
  - `sim-trivial-init/`: Source code to reproduce the simulation results with a trivial initialization.
  - `sim-informative-alpha/`: Source code to reproduce the simulation results with a strong informative prior on the mixture proportions.
- `data-analyses/`: Source code to reproduce the results from the real-world data analyses.
  - `wine/`: Source code to reproduce the results from the analysis of the *wines* data.
  - `kernels/`: Source code to reproduce the results from the analysis of the *seeds* data.
- `stan/`: Stan code used to fit a $K$ component Bayesian Gaussian Mixture, assuming conjugate priors.
- `one-dimensional-densities/`: Source code to generate one-dimensional density plots. 

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

Parallelization over multiple CPUs is conducted via *forking* (rather than *sockets*), so it only works on POSIX systems (i.e., macOS, Linux, Unix, BSD), not on Windows. To run the scripts on Windows, one would need to set the number of CPU workers to one. For further details, see [Package '`parallel`'](https://stat.ethz.ch/R-manual/R-devel/library/parallel/doc/parallel.pdf). We run all our code on macOS (version 15.3.1). 

## <a name="run"></a> Running the Scripts

Step-by-step instructions on how to reproduce the results from the article:

#### Simulation Results

To reproduce the simulation results from the article, start by downloading the folder `simulations/sim_ar_neg_05/setting1/` and execute the scripts `sim1_data.R`, `sim1_rw.R`, `sim1_stan.R`, and `sim1_results.R` (in that order). It is recommended to restart the R session each time you execute a script in order to terminate any remaining "zombie processes" from the parallelization exercise. Then, repeat the same steps for folders `simulations/sim_ar_neg_05/setting2/`-to-`simulations/sim_ar_neg_05/setting9/`. To obtain box plots and tables summarizing the results for all nine settings, execute the scripts `boxplots.R`, `optimal_x_vals.R`, and  `table_medians.R` from the `simulations/sim_ar_neg_05/sim-results/` folder.

To reproduce the results from our illustrative example, download the folder `simulations/sim-illustrative/` and execute the script `sim_illustrative.R`.

To reproduce the simulation results with a varying sample sizes, download the folder `simulations/sim-varying-n/` and execute the scripts `sim_n50.R`, `sim_n125.R`, `sim_n250.R`, `sim_n375.R`, `sim_n500.R`, and  `sim_n_results.R` (in that order). Once again, it is recommended to restart the R session each time you execute a script to terminate any remaining "zombie processes" from the parallelization exercise. 

To reproduce the simulation results with a trivial initialization, download the folder `simulations/sim-trivial-init/` and execute the scripts `sim_d5_trivial_init.R`, `sim_d10_trivial_init.R`, `sim_d15_trivial_init.R`, and `simtrivial_init_results.R` (in that order). Once again, it is recommended to restart the R session each time you execute a script to terminate any remaining "zombie processes" from the parallelization exercise. 

To reproduce the simulation results with a strong informative prior on the mixture proportions, download the folder `sim-informative-alpha/` and execute the scripts `sim_alphas_d5.R`, `sim_alphas_d10.R`, `sim_alphas_d15.R`, and `sim_alphas_results.R` (in that order). Once again, it is recommended to restart the R session each time you execute a script to terminate any remaining "zombie processes" from the parallelization exercise. 

#### Results from Benchmark Data Analysis

To reproduce the results from the analysis of benchmark data, start by downloading the folder `data-analyses/wine/`  and execute the script `wines.R`. Then, execute the script `kernels.R` from the folder `data-analyses/kernels/`.

## <a name="data"></a> Data

We make use of the [*wine*](https://archive.ics.uci.edu/dataset/109/wine) (Aeberhard and Forina, 1991) and [*seeds*](https://archive.ics.uci.edu/dataset/236/seeds) (Charytanowicz et al., 2012) data sets. Both data sets are publicly available from the [UC Irvine Machine Learning Repository](https://archive.ics.uci.edu/).

## <a name="cite"></a> Citation

If you use any part of this code in your work, please consider citing our paper:

```
@misc{marin_bob,
  title         = {BOB: Bayesian Optimized Bootstrap for Approximate Posterior Sampling in Gaussian Mixture Models}, 
  author        = {Santiago Marin and Bronwyn Loong and Anton H. Westveld},
  year          = {2024},
  eprint        = {2311.03644},
  archivePrefix = {arXiv},
  primaryClass  = {stat.ME}
}
```

## <a name="refs"></a> References

Aeberhard, S., and Forina, M. (1991), "Wine." *UCI Machine Learning Repository*. https://doi.org/10.24432/C5PC7J.

Charytanowicz, M., Niewczas, J., Kulczycki, P., Kowalski, P., and Lukasik, S. (2012), "Seeds." *UCI Machine Learning Repository*. https://doi.org/10.24432/C5H30K.

Marin, S., Loong, B., and Westveld, A. H. (2025+), "BOB: Bayesian Optimized Bootstrap for Uncertainty Quantification in Gaussian Mixture Models."
