# smesir: Hierarchically modeling disease transmission rates

## Getting Started

Installation Instructions
```
library(devtools)
install_github("davidbuch/smesir")
```

## Example

```
# Load a COVID-19 dataset
data()

# Fit a single region and inspect the model fit
# for single-region analyses, the "data" argument is expected to be a
sr_fit <- smesir(y ~ x1 + x2, data  = , etc)

# Fit a hierarchical model to a panel of regions
# note here that the "data" argument must be a list of matrices
# rather than a data.frame object, which is a list of vectors
mr_fit <- smesir(y ~ x1 + x2, data = , etc)

# Inspect the model fit

# Forecast upcoming deaths across the regions
Jf <- 20
smesir_predict(Jf,mr_fit,new_x) 
# new_x is covariate data
# it is not required for predictions under the null model

```


