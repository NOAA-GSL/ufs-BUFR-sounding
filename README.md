# ufs-BUFR-sounding
Creates BUFR sounding profiles from FV3 LAM output.

Matt Pyle at EMC has experience with this code and was able to help debug the compilation and run.
Ben Blake from EMC provide a path the to code on WCOSS.
Jeff Beck has experience with running this code on Hera/Jet.
Terra Ladwig followed the expertise of the above people to run the code on Jet for GSL RRFS/RTMA-RRFS runs during the Spring Experiment 2021.  Note that there are grid specific requirements, including the number of vertical levels being hard coded to 64.  

To compile this code, wrf libraries are required.  They are located on Jet at '/mnt/lfs4/BMC/nrtrr/wrflibs' .
Originally (code from WCOSS) these libraries were saved within the source code directory.

To run this code on GSL RRFS data, there are 3 fix files required.  These are available in '/fix'. Note that profdat file is grid specific and platform and/or model output differences may cause these files to need updating.  Matt Pyle provide the profdat file needed for the GSL 2021 RRFS grid.  
