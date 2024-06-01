## `app.R`

Shiny app for the interactive model of evolutionary failure.

## `Fig2.R` and `Fig3.R`

Perform multiple simulations and plot the results to create raw figures.

## `master-burden-simulator.R`

Script for performing stand-alone simulations

## Dependencies

This code has been tested with R 4.2.2 on a machine running MacOSX 14.5 in an environment with these R packages/versions installed.

```text
tidyverse/2.0.0
shiny/1.7.4.1
adaptivetau/2.3-1 
deSolve/1.36
```

Running the `Fig2.R` simulations takes a couple of minutes on an M3 MacBook Pro. Running the `Fig2.R` simulations takes roughly two hours 