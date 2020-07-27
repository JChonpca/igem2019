#!/usr/bin/env sh

## One-time setup

MAINDIR=$PWD
SCRIPTDIR=$MAINDIR/scripts

## Settings for fitting OD/FP curves

BACKGROUNDTIME=60
MINOD=0.03
MAXMETHOD=2
FITTWOPOINTS=F
TIMEPOINTDELTA=4
VECTOR=pSB1C3,pSB1A2
GROWTHRATEMETHOD="exponential"
FLUORESCENCERATEMETHOD="exponential-with-production"

#For testing multiple options
#SETTINGS=minOD_${MINOD}_maxmethod_${MAXMETHOD}_bgtime_${BACKGROUNDTIME}_timeptdelta_${TIMEPOINTDELTA}_gmethod_${GROWTHRATEMETHOD}_fmethod_${FLUORESCENCERATEMETHOD}_${VECTOR}

#For just running with final options
SETTINGS=output

rm -r $MAINDIR/01-output-plate-fits
cd $MAINDIR/input-plate-data
batch_run.pl "mkdir -p $MAINDIR/01-output-plate-fits/#d && $SCRIPTDIR/burden_fit.R -p -i #d -o $MAINDIR/01-output-plate-fits/#d/#d  --avg-end-time ${BACKGROUNDTIME} --min-OD ${MINOD} --max-method ${MAXMETHOD} --two-point-fit ${FITTWOPOINTS} --growth-rate-method ${GROWTHRATEMETHOD} --fluorescence-rate-method ${FLUORESCENCERATEMETHOD} --time-point-delta ${TIMEPOINTDELTA}"

#Use this command instead if you want to produce plots
#batch_run.pl "mkdir -p $MAINDIR/01-output-plate-fits/#d && $SCRIPTDIR/burden_fit.R -i #d -o $MAINDIR/01-output-plate-fits/#d/#d --avg-end-time ${BACKGROUNDTIME}--min-OD ${MINOD} --max-method ${MAXMETHOD} --two-point-fit ${FITTWOPOINTS} --growth-rate-method ${GROWTHRATEMETHOD} --fluorescence-rate-method ${FLUORESCENCERATEMETHOD} --time-point-delta ${TIMEPOINTDELTA}"


## Test of including data measured at MSU
## Need extended time window for fitting MSU data.
# $SCRIPTDIR/burden_fit.R -i $MAINDIR/input-plate-data/exp057/exp057 -o $MAINDIR/01-output-plate-fits/exp057/exp057 --avg-end-time ${BACKGROUNDTIME} --min-OD ${MINOD} --max-method ${MAXMETHOD} --growth-rate-method ${GROWTHRATEMETHOD} --fluorescence-rate-method ${FLUORESCENCERATEMETHOD} --time-point-delta ${TIMEPOINTDELTA} --max-time=1200

##For graphing one experiment to spot check
#batch_run.pl -p exp050 "mkdir -p $MAINDIR/01-output-plate-fits/#d && $SCRIPTDIR/burden_fit.R -i #d -o $MAINDIR/01-output-plate-fits/#d/#d --min-OD ${MINOD} --max-method ${MAXMETHOD} --two-point-fit ${FITTWOPOINTS} --time-point-delta ${TIMEPOINTDELTA}"


## Create merged summary file
mkdir -p $MAINDIR/02a-output-summary-merged/$SETTINGS
cd $MAINDIR/01-output-plate-fits
batch_run.pl "cp #d.rates.summary.csv $MAINDIR/02a-output-summary-merged/$SETTINGS"
$SCRIPTDIR/burden_merge.R -i $MAINDIR/02a-output-summary-merged/$SETTINGS -o $MAINDIR/02a-output-summary-merged/$SETTINGS.csv -m $MAINDIR/igem2019_strain_metadata.csv
rm -r $MAINDIR/02a-output-summary-merged/$SETTINGS

## Create merged all results file
mkdir -p $MAINDIR/02b-output-all-merged/$SETTINGS
cd $MAINDIR/01-output-plate-fits
batch_run.pl "cp #d.rates.all.csv $MAINDIR/02b-output-all-merged/$SETTINGS"
$SCRIPTDIR/burden_merge.R -i $MAINDIR/02b-output-all-merged/$SETTINGS -o $MAINDIR/02b-output-all-merged/$SETTINGS.csv -m $MAINDIR/igem2019_strain_metadata.csv
rm -r $MAINDIR/02b-output-all-merged/$SETTINGS

cd $MAINDIR
mkdir -p $MAINDIR/03a-plate-variation
cd $MAINDIR/02a-output-summary-merged
$SCRIPTDIR/igem2019_plate_variation.R -i $SETTINGS.csv -o $MAINDIR/03a-plate-variation/$SETTINGS.csv

mkdir -p $MAINDIR/03b-setting-variation
$SCRIPTDIR/igem2019_setting_variation.R -i $MAINDIR/02a-output-summary-merged -o $MAINDIR/03b-setting-variation

mkdir -p $MAINDIR/04-normalization
$SCRIPTDIR/burden_normalize.R \
  -i $MAINDIR/02b-output-all-merged/$SETTINGS.csv \
  -o $MAINDIR/04-normalization/${SETTINGS} \
  -m $MAINDIR/igem2019_plate_metadata.csv \
  -s $MAINDIR/igem2019_strain_metadata.csv \
  -v $VECTOR \
  -c "JEB1204,JEB1205,JEB1206,JEB1207,JEB1208" \
  --no-control-normalization \
  --growth-rate-bw 0.014 \
  --GFP-rate-bw 300


# Add these options for comparing MSU data
#  -x exp057 \
#  -n exp057
mkdir -p $MAINDIR/05-burden-final-output
$SCRIPTDIR/igem2019_graph_normalized.R \
  -i $MAINDIR/04-normalization/$SETTINGS.part.burden.csv \
  -o $MAINDIR/05-burden-final-output/${SETTINGS} \
  -m $MAINDIR/igem2019_strain_metadata.csv

#### Summary statistics

$SCRIPTDIR/igem2019_calculate_statistics.R

##### RFP-series
mkdir -p $MAINDIR/10-RFP-series-output
$SCRIPTDIR/burden_fit.R -i $MAINDIR/input-plate-data-RFP-series/exp057/exp057 -o $MAINDIR/10-RFP-series-output/exp057 --min-OD ${MINOD} --max-method ${MAXMETHOD} --two-point-fit ${FITTWOPOINTS} --growth-rate-method ${GROWTHRATEMETHOD} --fluorescence-rate-method ${FLUORESCENCERATEMETHOD} --time-point-delta ${TIMEPOINTDELTA}
$SCRIPTDIR/burden_summary.R -i $MAINDIR/10-RFP-series-output/exp057.rates.summary.csv

##### BFP-series
mkdir -p $MAINDIR/11-BFP-series-output
$SCRIPTDIR/burden_fit.R -i $MAINDIR/input-plate-data-BFP-series/exp030/exp030 -o $MAINDIR/11-BFP-series-output/exp030 --min-OD ${MINOD} --max-method ${MAXMETHOD} --two-point-fit ${FITTWOPOINTS} --growth-rate-method ${GROWTHRATEMETHOD} --fluorescence-rate-method ${FLUORESCENCERATEMETHOD} --time-point-delta ${TIMEPOINTDELTA}
$SCRIPTDIR/burden_summary.R -i $MAINDIR/11-BFP-series-output/exp030.rates.summary.csv

##### Graphs
cd $MAINDIR
$SCRIPTDIR/igem2019_BFP_RFP_series_graphs.R

##### MSU comparison
## Not included in the paper. You need Need to re-enable MSU plate data and run previous commands for this to work

## create a file for graphing just this information
#head -n 1 $MAINDIR/04-normalization/output.plate-strain.burden.csv > $MAINDIR/20-MSU-compare-output/rates.summary.csv
#grep 'I13001\|I13602\|I3452\|J23070\|J23101\|J23104\|J23107\|J23113\|J23115\|J23117\|J52022\|J61000\|K541501\|K617003\|K634001\|K731721\|K733013\|K819010\|K819017\|S05050' $MAINDIR/04-normalization/output.plate-strain.burden.csv >> $MAINDIR/20-MSU-compare-output/rates.summary.csv
##change accession/strain for graphing script to work?
#sed -i .bak 's/strain/_strain/g' $MAINDIR/20-MSU-compare-output/rates.summary.csv
#sed -i .bak 's/accession/strain/g' $MAINDIR/20-MSU-compare-output/rates.summary.csv
#sed -i .bak 's/plate/isolate/g' $MAINDIR/20-MSU-compare-output/rates.summary.csv
#sed -i .bak 's/normalized.growth.rate.mean/growth.rate/g' $MAINDIR/20-MSU-compare-output/rates.summary.csv
#sed -i .bak 's/normalized.growth.rate.sd/growth.rate.sd/g' $MAINDIR/20-MSU-compare-output/rates.summary.csv
#sed -i .bak 's/normalized.GFP.rate.mean/GFP.rate/g' $MAINDIR/20-MSU-compare-output/rates.summary.csv
#sed -i .bak 's/normalized.GFP.rate.sd/GFP.rate.sd/g' $MAINDIR/20-MSU-compare-output/rates.summary.csv
##Now...Manually edit strains to be UTA or MSU
#$SCRIPTDIR/burden_summary.R -i $MAINDIR/20-MSU-compare-output/rates.summary.csv
#$SCRIPTDIR/igem2019_graph_normalized.R \
#  -i $MAINDIR/20-MSU-compare-output/msu.rates.all.csv \
#  -i $MAINDIR/20-MSU-compare-output/uta.rates.all.csv \
#  -o $MAINDIR/20-MSU-compare-output \
#  -m $MAINDIR/msu_plate_metadata.csv
