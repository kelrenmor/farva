# FActor Regression for Verbal Autopsy (FARVA)

## Model summary

FARVA is a probabilistic model for VA data in which some decedents have known COD (i.e., there exists some labeled training data) and there is interest in learning about individual COD, population CSMF, and/or the mean and covariance structure of responses in the symptom questionnaire. The explicit goals of the FARVA model are to:
- Capture dependence of symptoms given a cause.
- Share information across causes via hierarchical modeling to improve estimates associated with causes having few observed deaths.
- Allow both the conditional prevalence and the conditional association between symptoms to vary with covariates (e.g., age of patient, time of year, geographic region).
- Probabilistically predict cause of death for a new individual given their symptoms.
- Improve on COD and CSMF estimation relative to current state-of-the-art VA algorithms.

More details on the algorithm and the problem background are available in the manuscript associated with this code base. The submission draft is available on arXiv at [https://arxiv.org/abs/1908.07632](https://arxiv.org/abs/1908.07632). 

## Installation and use

To install the farva library, use the `devtools` package.

```
library(devtools)
devtools::install_github('kelrenmor/farva', quiet=T, upgrade=F, dependencies=T)
library(farva)
```

An example using the FARVA model to analyze the PHMRC data from Mexico City is provided in the [examples/example_phmrc_run.R](examples/example_phmrc_run.R) file. This code estimates the individual probabilities and the CSMF for a hold-out 'test' set of the data.

-----------------------
Original work Copyright (c) 2019 Kelly R. Moran

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.
