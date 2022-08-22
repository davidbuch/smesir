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
data(smesir_simstudy)
epi_params <- list(region_populations = c(2e5,6e5,5e5,3e5,2e5),
                  outbreak_times = c(1,3,2,2,1),
                  mean_removal_time = 1.5, 
                  incidence_probabilities = 0.01*c(0.25,0.5,0.25))

# Fit a single region and inspect the model fit
# for single-region analyses, the "data" argument is expected to be a
sr_fit <- smesir(y ~ x1 + x2, data  = , etc)

# Fit a hierarchical model to a panel of regions
# note here that the "data" argument must be a list of matrices
# rather than a data.frame object, which is a list of vectors
sfit <- smesir(deaths ~ X1 + X2 + X3, data = smesir_simstudy, epi_params  = epi_params, region_names = paste0("R", 1:K), quiet = FALSE, seed = 101)

sfor <- smesir_forecast(Jf,sfit,new_x = list(X1 = X1[(Jo + 1):J,], X2 = X2[(Jo + 1):J,], X3 = X3[(Jo + 1):J,]))

par(mfrow = c(2,3))
for(k in 1:K){
  plot(1:Jo,Y[1:Jo,k], xlab = "time", ylab = NA, main = paste("Region",k),
       xlim = c(1,J),ylim=c(0,max(sfor$confints[,,k])), type = "l",lwd=2)
  lines(Jo:J,Y[Jo:J,k],lty="dashed",lwd=2)
  plot_confint(sfor$confints[(Jo + 1):J,,k],x=(Jo + 1):J,density=15,col="blue")
}
par(mfrow = c(1,1))

```


