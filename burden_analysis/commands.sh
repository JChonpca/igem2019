#!/usr/bin/env sh

## One-time setup

MAINDIR=$PWD
SCRIPTDIR=$MAINDIR/scripts

## Setup for certain settings

MINOD=0.03
MAXMETHOD=1
FITTWOPOINTS=F
FPPERODATTIME=T
TIMEPOINTDELTA=4
SETTINGS=minOD_${MINOD}_maxmethod_${MAXMETHOD}_fit2pts_${FITTWOPOINTS}_timeptdelta_${TIMEPOINTDELTA}_fpperofattimept_${FPPERODATTIME}
VECTOR=pSB1C3,pSB1A2



cd $MAINDIR/input-plate-data
rm -r $MAINDIR/01-output-plate-fits
batch_run.pl "mkdir -p $MAINDIR/01-output-plate-fits/#d && $SCRIPTDIR/burden_fit.R -p -i #d -o $MAINDIR/01-output-plate-fits/#d/#d --min-OD ${MINOD} --max-method ${MAXMETHOD} --two-point-fit ${FITTWOPOINTS} --FP-per-OD-at-time ${FPPERODATTIME} --time-point-delta ${TIMEPOINTDELTA}"


##For graphing to spot check
#batch_run.pl -p exp018 "mkdir -p $MAINDIR/01-output-plate-fits/#d && $SCRIPTDIR/burden_fit.R -i #d -o $MAINDIR/01-output-plate-fits/#d/#d --min-OD ${MINOD} --max-method ${MAXMETHOD} --two-point-fit ${FITTWOPOINTS} --time-point-delta ${TIMEPOINTDELTA}"


## Create merged summary filemkdir -p $MAINDIR/02-output-summary-merged/$SETTINGS
mkdir -p $MAINDIR/02-output-summary-merged/$SETTINGS
cd $MAINDIR/01-output-plate-fits
batch_run.pl "cp #d.rates.summary.csv $MAINDIR/02-output-summary-merged/$SETTINGS"
$SCRIPTDIR/burden_merge.R -i $MAINDIR/02-output-summary-merged/$SETTINGS -o $MAINDIR/02-output-summary-merged/$SETTINGS.csv -m $MAINDIR/igem2019_strain_metadata.csv
rm -r $MAINDIR/02-output-summary-merged/$SETTINGS

## Create merged all results file
mkdir -p $MAINDIR/02-output-all-merged/$SETTINGS
cd $MAINDIR/01-output-plate-fits
batch_run.pl "cp #d.rates.all.csv $MAINDIR/02-output-all-merged/$SETTINGS"
$SCRIPTDIR/burden_merge.R -i $MAINDIR/02-output-all-merged/$SETTINGS -o $MAINDIR/02-output-all-merged/$SETTINGS.csv -m $MAINDIR/igem2019_strain_metadata.csv
rm -r $MAINDIR/02-output-all-merged/$SETTINGS

mkdir -p $MAINDIR/03-plate-variation
cd $MAINDIR/02-output-summary-merged
$SCRIPTDIR/igem2019_plate_variation.R -i $SETTINGS.csv -o $MAINDIR/03-plate-variation/$SETTINGS.csv

mkdir -p $MAINDIR/03-setting-variation
$SCRIPTDIR/igem2019_setting_variation.R -i $MAINDIR/02-output-summary-merged -o $MAINDIR/03-setting-variation



mkdir -p $MAINDIR/04-normalization
$SCRIPTDIR/burden_normalize.R \
  -i $MAINDIR/02-output-all-merged/$SETTINGS.csv \
  -o $MAINDIR/04-normalization/${SETTINGS}_${VECTOR} \
  -m $MAINDIR/igem2019_plate_metadata.csv \
  -s $MAINDIR/igem2019_strain_metadata.csv \
  -v $VECTOR \
  -c "JEB1204,JEB1205,JEB1206,JEB1207,JEB1208"
