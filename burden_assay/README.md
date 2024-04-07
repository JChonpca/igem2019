## igem2019 burden analysis

Contains the scripts and input data needed for performing the burden analysis.

`commands.sh` is a master script that runs the entire analysis.

Steps are:

1. `burden_fit.R`

Fit the growth rates and GFP production rates for each plate.

2. `burden_merge.R`

Merge the results from all plates. Also removes GFP.rate from plasmids expressing GFP.

3. `burden_normalize.R`

Normalizes readings within each plate to one another based on no-burden peak and internal control strains.
