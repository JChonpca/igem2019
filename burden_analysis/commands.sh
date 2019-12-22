MINOD=0.08
MAXMETHOD=2
FITTWOPOINTS=F
TIMEPOINTDELTA=2
SETTINGS=minOD_${MINOD}_maxmethod_${MAXMETHOD}_fit2pts_${FITTWOPOINTS}_timeptdelta_${TIMEPOINTDELTA}

MAINDIR=$PWD
SCRIPTDIR=$MAINDIR/scripts

cd $MAINDIR/input-plate-data
rm -r $MAINDIR/01-output-plate-fits
batch_run.pl "mkdir -p $MAINDIR/01-output-plate-fits/#d && $SCRIPTDIR/burden_fit.R -i #d -o $MAINDIR/01-output-plate-fits/#d/#d -p --min-OD ${MINOD} --max-method ${MAXMETHOD} --two-point-fit ${FITTWOPOINTS} --time-point-delta ${TIMEPOINTDELTA}"


##For graphing to spot check
#batch_run.pl -p exp018 "mkdir -p $MAINDIR/01-output-plate-fits/#d && $SCRIPTDIR/burden_fit.R -i #d -o $MAINDIR/01-output-plate-fits/#d/#d --min-OD ${MINOD} --max-method ${MAXMETHOD} --two-point-fit ${FITTWOPOINTS} --time-point-delta ${TIMEPOINTDELTA}"


## Create merged summary file
mkdir -p $MAINDIR/02-output-summary-merged/$SETTINGS
cd $MAINDIR/01-output-plate-fits
batch_run.pl "cp #d.rates.summary.csv $MAINDIR/02-output-summary-merged/$SETTINGS"
$SCRIPTDIR/burden_merge.R -i $MAINDIR/02-output-summary-merged/$SETTINGS -o $MAINDIR/02-output-summary-merged/$SETTINGS.csv
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
  -o $MAINDIR/04-normalization/$SETTINGS \
  -m $MAINDIR/igem2019_plate_metadata.csv \
  -c "JEB1204,JEB1205,JEB1206,JEB1207,JEB1208"
