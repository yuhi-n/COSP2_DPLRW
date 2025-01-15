# COSP2_DPLRW
Doppler radar signal simulator, implemented on COSP2.<br>
Details of each file are as follows.<br>
The additional section is controlled by fortran preprocessor, named `OPT_DPLRW`.
And other changes are indicated by the following comments: `modified by YN` or `added by YN`.

## 1. cosp.F90
Main program of COSP2.<br>
This calls `quickbeam.F90`.

## 2. quickbeam.F90
This calculates attenuation by gas and hydrometeor for stratiform and convective clouds separately at first.<br>
Statistics of Ze (dB mm $`{aa}^6`$/m $`{}^3`$)

$\sin x$

## 3. quickbeam_optics.F90
