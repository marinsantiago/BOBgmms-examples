# BOBgmms-examples

<!-- badges: start -->

[![arXiv](https://img.shields.io/badge/arXiv-2311.03644-blue.svg)](https://arxiv.org/abs/2311.03644)

<!-- badges: end -->

This repository contains source code to reproduce the results from the article "[BOB: Bayesian Optimized Bootstrap for Uncertainty Quantification in Gaussian Mixture Models](https://arxiv.org/abs/2311.03644)" (Marin, Loong and Westveld, 2024).

The `BOBgmms` R package is available at the Github repository: [https://github.com/marinsantiago/BOBgmms](https://github.com/marinsantiago/BOBgmms)

## Contents

- [Folders Structure ](#folders)
- [Software Requirements](#software)
- [System Requirements](#system)
- [Data](#data)
- [Running the Scripts](#run)
- [Citation](#cite)
- [References](#refs)

## <a name="folders"></a> Folders Structure 

- `R`: Helper functions required throughout the analysis.

- `simulations`: Source code to reproduce the results in the *Simulations* section from the article.
  - `sim1`: Source code to reproduce the results for the 1st simulation setting.
  - `sim2`: Source code to reproduce the results for the 2nd simulation setting.
  - `sim3`: Source code to reproduce the results for the 3rd simulation setting.
  - `sim4`: Source code to reproduce the results for the 4th simulation setting.
  - `sim5`: Source code to reproduce the results for the 5th simulation setting.
  - `sim6`: Source code to reproduce the results for the 6th simulation setting.
  - `sim7`: Source code to reproduce the results for the 7th simulation setting.
  - `sim8`: Source code to reproduce the results for the 8th simulation setting.
  - `sim9`: Source code to reproduce the results for the 9th simulation setting.
  - `sim-illustrative`: Source code to reproduce the results used as an illustrative example.
  - `sim-varying-n`: Source code to reproduce the results with a growing sample size.
    - `sim_n1`: Case with $n=50$.
    - `sim_n2`: Case with $n=125$.
    - `sim_n3`: Case with $n=250$.
    - `sim_n4`: Case with $n=375$.
    - `sim_n5`: Case with $n=500$.
  - `sim_results`: Source code to generate plots and tables with results from the simulation experiments.

- `stan`: Stan code to fit a Bayesian mixture of Gaussians assuming conjugate priors.

- `wheat-kernels`: Source code to reproduce the results from the analysis of the *Seeds* data.

- `wines`: Source code to reproduce the results from the analysis of the *Wine* data.

## <a name="software"></a> Software Requirements

We run all our simulations and data analyses in R (version 4.3.2) via Rstudio (version 2023.12.1+402). Additionally, we make use of the following R packages (along with their dependencies):

  - `clusterHD` v1.0.2
  
  - `coda` v0.19-4.1
  
  - `mclust` v6.1
  
  - `mvnfast`v0.2.8
  
  - `rstan` v2.32.6 (for instructions on how to download and install `rstan`, please see [RStan Quick Start Guide](https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started))
  
  - `sparcl` v1.0.4
  
  - `BOBgmms` v1.0.0
  
  - `LaplacesDemon` v16.1.6
  
  - `MASS` v7.3-60.0.1
  
  - `RColorBrewer` v1.1-3

Please be aware that before downloading and installing `rstan`, you need to set up your R installation to be able to compile C++ code. For additional details and instructions, see [RStan Getting Started](https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started) under **Configuring C++ Toolchain**. 

## <a name="system"></a> System Requirements

The parallelization is conducted via *forking* (rather than *sockets*), so it only works on POSIX systems (i.e., Mac, Linux, Unix, BSD), not on Windows. For further details, see [Package '`parallel`'](https://stat.ethz.ch/R-manual/R-devel/library/parallel/doc/parallel.pdf). We run our code on Mac (macOS version 14.4.1). 

## <a name="data"></a> Data

We make use of the [*Wine*](https://archive.ics.uci.edu/dataset/109/wine) (Aeberhard and Forina, 1991) and [*Seeds*](https://archive.ics.uci.edu/dataset/236/seeds) (Charytanowicz et al., 2012) data sets. Both data sets are publicly available from the [UC Irvine Machine Learning Repository](https://archive.ics.uci.edu/).

## <a name="run"></a> Running the Scripts

Step-by-step instructions on how to reproduce the results from the article:

#### Simulation Results

To reproduce the simulation results from the article, start by downloading the folder `sim1` and execute the scripts `sim1_data.R`, `sim1_wbb.R`, `sim1_bob.R`, `sim1_nuts.R`, and `sim1_results.R`, in that order. **It is strongly recommended to restart the R session** each time you execute a script in order to terminate any remaining "zombie processes" from the parallelization exercise. Then, repeat the same steps for folders `sim2` -- `sim9`. To obtain box plots and tables summarizing the results for all nine settings, execute the scripts `boxplots.R` and `medians.R` from the `sim_results` folder.

To reproduce the results from our illustrative example, download the folder `sim-illustrative` and execute the scripts `sim_illustrative_rw.R`, `sim_illustrative_nuts.R`, and `sim_illustrative_results.R`, in that order. Again, **it is strongly recommended to restart the R session** each time you execute a script in order to terminate any remaining "zombie processes" from the parallelization exercise.

To reproduce the simulation results for varying sample sizes, download the folder `sim-varying-n` and execute the `sim_n1.R` script from the `sim_n1` folder. Repeat this process for the folders `sim_n2` -- `sim_n5`. Once again, **it is strongly recommended to restart the R session** each time you execute a script to terminate any remaining "zombie processes" from the parallelization exercise. To generate line plots summarizing the results, execute the script `lineplots.R` located in the `sim_results` folder.

#### Results from Benchmark Data Analysis

To reproduce the results from our benchmark data analysis, start by downloading the folder `wines` and execute the scripts `wines_rw.R`, `wines_nuts.R`, and `wines_results.R`, in that order. **It is strongly recommended to restart the R session** each time you execute a script in order to terminate any remaining "zombie processes" from the parallelization exercise. The data is available in the file `wine.Rdata` (Aeberhard and Forina, 1991). Then, repeat this process for the `wheat-kernels` folder. In such a case, the data is available in the file `seeds.Rdata` (Charytanowicz et al., 2012).

## <a name="cite"></a> Citation

If you use any part of this code in your work, please consider citing our paper:

```
@misc{marin2023bob,
      title         = {BOB: Bayesian Optimized Bootstrap for Uncertainty Quantification in Gaussian Mixture Models}, 
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

Marin, S., Loong, B., and Westveld, A. H. (2024), "BOB: Bayesian Optimized Bootstrap for Uncertainty Quantification in Gaussian Mixture Models."
