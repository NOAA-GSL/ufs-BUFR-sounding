#! /bin/sh

set -x

source ./machine-setup.sh > /dev/null 2>&1

module purge

#module use ~rtrr/RRFS/sorc/EMC_post/modulefiles
module use /gpfs/dell2/emc/modeling/noscrub/Edward.Colon/NOAA_3drtma_bec/sorc/rtma_post.fd/modulefiles
module load post/v8.0.0-wcoss_dell_p3

module list

cd ./regional_bufr.fd

make clean
make FC="mpif90 -f90=ifort" NETCDF_INCLUDE="-I${NETCDF}/include"  \
  NETCDF_LDFLAGS="-L${NETCDF}/lib -lnetcdf -lnetcdff"


cd ../

