#! /bin/sh

set -x

source ./machine-setup.sh > /dev/null 2>&1

module purge

module use ~rtrr/RRFS/sorc/EMC_post/modulefiles
module load post/v8.0.0-jet

module list

cd ./regional_bufr.fd

make clean
make FC=mpiifort NETCDF_INCLUDE="-I${NETCDF}/include"  \
  NETCDF_LDFLAGS="-L${NETCDF}/lib -lnetcdf -lnetcdff"


cd ../

