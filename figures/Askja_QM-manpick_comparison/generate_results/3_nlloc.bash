#!/bin/bash
# -------------------------------------------------------------------------------------
# Script to run NLLoc on the manually picked earthquakes being used in the manual-QM
# benchmark for the QuakeMigrate manuscript:

#     Winder, T., Bacon, C.A., Smith, J.D., Hudson, T.S., and White, R.S.
#     QuakeMigrate: a Python Package for Automatic Earthquake Detection and Location
#     Using Waveform Migration and Stacking. (to be submitted to Seismica).

# Author: Tom Winder, Oct 2024
# -------------------------------------------------------------------------------------


#_______________________ FUNCTION 1 --- CREATE CONTROL FILE ___________________________
#######################################################################################
create_control_file(){                                                                #
#######################################################################################

    controlFilename=$1
    phase=$2
    vel_model=$3
    time_grid=$4
    out_dir=$5
    obs_dir=$6
    loc_dir=$7
    loc_label=$8
    stationListFile=$9


######################################################################################
# Get generic lines of control file and write to file:                               #
#vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv

echo "# ----Generic Control File Statements---" > $controlFilename
echo "" >> $controlFilename
echo "CONTROL 1 54321" >> $controlFilename
echo "" >> $controlFilename
echo "TRANS LAMBERT WGS-84 65.1 -16.6 64.9 65.3 0.0" >> $controlFilename
echo "" >> $controlFilename
echo "# ----End of Generic Control File Statements---" >> $controlFilename
echo "" >> $controlFilename
echo "" >> $controlFilename
#_____________________________________________________________________________________
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^



######################################################################################
# Get/specify Vel2Grid control file statements:                                      #
#vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv

#vvvvvvvv

echo "# ----Vel2Grid Control File Statements---" >> $controlFilename
echo "VGOUT  $out_dir/model/$vel_model # Specified by user" >> $controlFilename
echo "" >> $controlFilename
#^^^^^^^


echo "VGTYPE $phase # Specified by user" >> $controlFilename


#vvvvvvvv

echo "" >> $controlFilename
echo "VGGRID 2 4000 900 0.0 0.0 -3.0 0.05 0.05 0.05 SLOW_LEN" >> $controlFilename
echo "" >> $controlFilename

# Greenfield et al. 2016 velocity model (1D average of 3D LET model)
# Greenfield, T., White, R.S. and Roecker, S., 2016. The magmatic plumbing system of the Askja central volcano, Iceland, as imaged by seismic tomography. Journal of Geophysical Research: Solid Earth, 121(10), pp.7211-7229.
echo "LAYER -2.0 3.18 0.05 1.85  0.020 0.0 0" >> $controlFilename
echo "LAYER -1.0 3.23 0.24 1.87  0.140 0.0 0" >> $controlFilename
echo "LAYER  0.0 3.47 0.32 2.01  0.160 0.0 0" >> $controlFilename
echo "LAYER  0.5 3.63 0.34 2.09  0.200 0.0 0" >> $controlFilename
echo "LAYER  1.0 3.80 0.32 2.19  0.180 0.0 0" >> $controlFilename
echo "LAYER  1.5 3.96 1.68 2.28  0.960 0.0 0" >> $controlFilename
echo "LAYER  2.0 4.80 0.54 2.76  0.300 0.0 0" >> $controlFilename
echo "LAYER  2.5 5.07 0.44 2.91  0.240 0.0 0" >> $controlFilename
echo "LAYER  3.0 5.29 1.14 3.03  0.640 0.0 0" >> $controlFilename
echo "LAYER  3.5 5.86 0.30 3.35  0.160 0.0 0" >> $controlFilename
echo "LAYER  4.0 6.01 0.34 3.43  0.200 0.0 0" >> $controlFilename
echo "LAYER  4.5 6.18 0.36 3.53  0.180 0.0 0" >> $controlFilename
echo "LAYER  5.0 6.36 0.22 3.62  0.120 0.0 0" >> $controlFilename
echo "LAYER  5.5 6.47 0.04 3.68  0.020 0.0 0" >> $controlFilename
echo "LAYER  6.0 6.49 0.00 3.69 -0.020 0.0 0" >> $controlFilename
echo "LAYER  6.5 6.49 0.04 3.68  0.020 0.0 0" >> $controlFilename
echo "LAYER  7.0 6.51 0.04 3.69  0.020 0.0 0" >> $controlFilename
echo "LAYER  7.5 6.53 0.08 3.70  0.020 0.0 0" >> $controlFilename
echo "LAYER  8.0 6.57 0.06 3.71  0.040 0.0 0" >> $controlFilename
echo "LAYER  8.5 6.60 0.06 3.73  0.020 0.0 0" >> $controlFilename
echo "LAYER  9.0 6.63 0.06 3.74  0.000 0.0 0" >> $controlFilename
echo "LAYER  9.5 6.66 0.02 3.74  0.000 0.0 0" >> $controlFilename
echo "LAYER 10.0 6.67 0.00 3.74  0.000 0.0 0" >> $controlFilename
echo "LAYER 10.5 6.67 0.00 3.74  0.003 0.0 0" >> $controlFilename
echo "LAYER 14.0 6.67 0.06 3.75  0.020 0.0 0" >> $controlFilename
echo "LAYER 14.5 6.70 0.04 3.76  0.040 0.0 0" >> $controlFilename
echo "LAYER 15.0 6.72 0.06 3.78  0.020 0.0 0" >> $controlFilename
echo "LAYER 15.5 6.75 0.02 3.79  0.020 0.0 0" >> $controlFilename
echo "LAYER 16.0 6.76 0.01 3.80  0.005 0.0 0" >> $controlFilename
echo "LAYER 18.0 6.78 0.06 3.81  0.040 0.0 0" >> $controlFilename
echo "LAYER 18.5 6.81 0.08 3.83  0.040 0.0 0" >> $controlFilename
echo "LAYER 19.0 6.85 0.06 3.85  0.040 0.0 0" >> $controlFilename
echo "LAYER 19.5 6.88 0.02 3.87  0.000 0.0 0" >> $controlFilename
echo "LAYER 20.0 6.89 0.00 3.87  0.000 0.0 0" >> $controlFilename
echo "LAYER 42.0 6.89 0.00 3.87  0.000 0.0 0" >> $controlFilename


echo "" >> $controlFilename
echo "# ----End of Vel2Grid Control File Statements---" >> $controlFilename
echo "" >> $controlFilename
echo "" >> $controlFilename
#^^^^^^^
#_____________________________________________________________________________________
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^



######################################################################################
# Get/specify Grid2Time control file statements:                                     #
#vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv


#vvvvvvvv
echo "# ----Grid2Time Control File Statements---" >> $controlFilename
echo "" >> $controlFilename
#^^^^^^^^


echo "GTFILES  $out_dir/model/$vel_model  $out_dir/time/$time_grid $phase # Specified by user" >> $controlFilename


#vvvvvvvv
echo "" >> $controlFilename
echo "GTMODE GRID2D ANGLES_YES" >> $controlFilename
echo "" >> $controlFilename
#^^^^^^^^


echo "#Stations:" >> $controlFilename
awk '{print $0}' ${stationListFile} >> $controlFilename

#vvvvvvvv
echo "" >> $controlFilename
echo "GT_PLFD  1.0e-3  0 # From Podvin and Lemcomte, 1991, finite difference scheme to calculate travel times" >> $controlFilename
echo "" >> $controlFilename
echo "# ----End of Grid2Time Control File Statements---" >> $controlFilename
#^^^^^^^^
#_____________________________________________________________________________________
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


######################################################################################
# Get/specify NLLoc control file statements:                                         #
#vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv

#vvvvvvvv
echo "# ----NLLoc Control File Statements---" >> $controlFilename
echo "" >> $controlFilename
echo "LOCSIG Tom Winder, University of Iceland" >> $controlFilename
echo "" >> $controlFilename
echo "LOCCOM replacecomment" >> $controlFilename
echo "" >> $controlFilename
#^^^^^^^^


echo "LOCFILES $obs_dir/*.nonlinloc  NLLOC_OBS $out_dir/time/$time_grid  $loc_dir/$loc_label # Specified by user" >> $controlFilename


#vvvvvvvv

echo "" >> $controlFilename
echo "LOCHYPOUT SAVE_NLLOC_ALL" >> $controlFilename
echo "LOCSEARCH  OCT 10 10 4 0.01 300000 5000 0 0 # Can be varied, check before major runs." >> $controlFilename
echo "LOCGRID  145 131 81  -35 -30 -3  0.5 0.5 0.5   PROB_DENSITY  SAVE" >> $controlFilename
echo "" >> $controlFilename
echo "LOCMETH EDT_OT_WT 9999.0 4 -1 -1 -1 6 -1.0 1" >> $controlFilename
echo "LOCGAU 0.1 0.0" >> $controlFilename
echo "LOCGAU2 0.03 0.00 1.5" >> $controlFilename
echo "LOCPHASEID  P   P p G PN PG" >> $controlFilename
echo "LOCPHASEID  S   S s G SN SG" >> $controlFilename
echo "LOCQUAL2ERR 0.01 0.02 0.05 0.1 99999.9" >> $controlFilename
echo "LOCPHSTAT 9999.0 -1 9999.0 1.0 1.0 9999.9 -9999.9 9999.9" >> $controlFilename
echo "LOCANGLES ANGLES_YES 5" >> $controlFilename
echo "" >> $controlFilename
echo "# ----End of NLLoc Control File Statements---" >> $controlFilename
#^^^^^^^^
#_____________________________________________________________________________________
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

######################################################################################
}                                                                                    #
######################################################################################




#________________________________ SPECIFY USER INPUTS _________________________________
#######################################################################################

# paths
in_dir=./inputs/NLLOC
out_dir=./outputs/NLLOC

# External station list file (in NLLoc GTSOURCE format):
stationListFile=$in_dir/askja_QM-manpick_stations_NLLOC.txt

# Specify velocity model name:
vel_model=layer.velocity_model_Greenfield2016

# Specify time grid name for Grid2Time run:
time_grid=layer.time_grid_LUT_Greenfield2016

# Specify Run overall label:
runOverallLabel=qm-man_matched

# input and output dirs
obs_dir=$in_dir/obs
loc_dir=$out_dir/loc

# Specify loc label
loc_label=loc.$runOverallLabel
mkdir -p ${loc_dir}

# Specify control file name:
controlFilename=$in_dir/control.in

# Create log file:
mkdir -p $out_dir/logs
run_time=`date +%Y%m%dT%H%M%S`
logfile=$out_dir/logs/nlloc_run_${run_time}.log
echo "Start: $run_time" > $logfile

# Count number of events
countMax=`ls $obs_dir/*nonlinloc | wc -l`
echo "TOTAL NUMBER OF EQS = $countMax"
echo "TOTAL NUMBER OF EQS = $countMax" >> $logfile


######################################################################################
######################################################################################



#____________________________ RUN Vel2Grid AND Grid2Time ______________________________
#######################################################################################
# Only have to do once for each velocity model! Slow.                                 #
#######################################################################################


# Run Vel2Grid for P and S:
mkdir -p $out_dir/model
create_control_file $controlFilename P $vel_model $time_grid $out_dir $obs_dir $loc_dir $loc_label $stationListFile
echo "Running Vel2Grid for P"
Vel2Grid $controlFilename  2>&1 | tee -a $logfile
create_control_file $controlFilename S $vel_model $time_grid $out_dir $obs_dir $loc_dir $loc_label $stationListFile
echo "Running Vel2Grid for S"
Vel2Grid $controlFilename 2>&1 | tee -a $logfile

# Run Grid2Time for P and S:
mkdir -p $out_dir/time
create_control_file $controlFilename P $vel_model $time_grid $out_dir $obs_dir $loc_dir $loc_label $stationListFile
echo "Running Grid2Time for P"
Grid2Time $controlFilename 2>&1 | tee -a $logfile
create_control_file $controlFilename S $vel_model $time_grid $out_dir $obs_dir $loc_dir $loc_label $stationListFile
echo "Running Grid2Time for S"
Grid2Time $controlFilename 2>&1 | tee -a $logfile

#######################################################################################
#######################################################################################




#___________________________ RUN EVENTS THROUGH NONLINLOC _____________________________
#######################################################################################
#######################################################################################

# Create control file
create_control_file $controlFilename P $vel_model $time_grid $out_dir $obs_dir $loc_dir $loc_label $stationListFile

# Run NonLinLoc
NLLoc $controlFilename 2>&1 | tee -a $logfile

# Combine individual hypocentre-phase summary files
cat ${loc_dir}/${loc_label}.2*.grid0.loc.hyp > ${loc_dir}/${loc_label}.summary.hyp

#
end_time=`date +%Y%m%dT%H%M%S`
echo "Finished: $end_time" >> $logfile

#######################################################################################
#######################################################################################


#END!
