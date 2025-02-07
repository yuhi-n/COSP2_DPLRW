# COSP2_DPLRW
Doppler radar signal simulator, implemented on COSP2.<br>
This is a beta version. The author is due for a code review by COSP2 development team.

Details of each file are as follows.<br>
The additional section is controlled by fortran preprocessor, named `OPT_DPLRW`.<br>
And other changes are indicated by the following comments: `modified by YN` or `added by YN`.

## 1. cosp.F90
Main program of COSP2 is modified to handle new output variables related to the Dopplet signal simulator.<br>
This calls `quickbeam.F90`.

## 2. quickbeam.F90
This calculates attenuation by gas and hydrometeor for stratiform and convective clouds separately at first.<br>
Statistics of Z<sub>e</sub> (dB mm<sup>6</sup>/m<sup>3</sup>), V<sub>d</sub> (m/s), and spectrum width (m/s) are produced.

## 3. quickbeam_optics.F90
This is modified to calculate droplet fall speeds multiplied by radar reflectivity and integrated over droplet size distribution.<br>
The `hydroclass init` routine is expanded to include parameters of droplet fall velocity, and is adjusted for control by NAMELIST.

## 4. cosp_config.F90
This contains parameters for min/max values and bin width for statistics of Z<sub>e</sub>, V<sub>d</sub>, and spectrum width.

## 4. cosp_utils.F90
The routine for precipitation flux conversion into mixing ratio is updated to refer to droplet fall speed.

## 5. cosp_miroc_drive.F90
This is an example COSP2 driver, which is for MIROC6 GCM.<br>
The vertcal air velocity and convective mass flux are new required variables.<br>
The variables and I/O routines for Doppler signals are added.<br>
This calls `cosp.F90` and subcolumn processes.
