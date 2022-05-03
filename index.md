This repository contains the R package [conformalInference.fd](https://cran.r-project.org/web/packages/conformalInference.fd/index.html) (now available also on CRAN), which can produce valid prediction regions at levels 1-α or 1-2α under the basic assumption of _i.i.d._ regression data. 

The package was developed as part of my MSc. final thesis in Mathematical Engineering at Politecnico di Milano, as a functional extension of the main methods for Conformal Prediction for regression in the univariate response case.

### Code Structure

There are three main famililies of functions:

- Prediction methods
- Regression methods
- Plot methods


The central idea upon which the package is designed is the following: regression methods **should not** be included into the prediction methods themselves. Final users can pass as input to the prediction methods custom-coded regression algorithms, which may be more suitable for the prediction task at hand. Anyways the most common regression methods are implemented in the package.

### Main Functions

<br/>
<div align="center">


| Syntax      | Description |
| ----------- | ---------------- |
|concurrent| Build concurrent regression model|
|conformal.fun.jackplus | Computes Jackknife+ prediction sets|
|conformal.fun.split| Computes Split Conformal prediction sets|
|conformal.fun.msplit| Computes Multi Split Conformal prediction sets|
|mean_lists |Build regression method with mean|
|plot_fun |Plot the output prediction methods|
  
  </div>


### Detailed description

A complete description of the theory underpinning the package, an analysis of all the main functions as well as a case study is presented in my final MSc. thesis paper, availble at the following [link]().

### Acknownledgments

Prof. Simone Vantini - _Politecnico di Milano_

Dr. Jacopo Diquigiovanni

Dr. Matteo Fontana

Prof. Aldo Solari - _Università Bicocca di Milano_
