#! /bin/sh

set -x

#. /apps/lmod/lmod/init/sh
source ./machine-setup.sh > /dev/null 2>&1

module purge
#module use -a ./regional_bufr.fd/modulefiles
#module load modulefile.nam_post0_wcoss_p3_d

module use ~rtrr/RRFS/sorc/EMC_post/modulefiles
module load post/v8.0.0-jet

module list

cd ./regional_bufr.fd

make clean
#make
make FC=mpiifort NETCDF_INCLUDE="-I${NETCDF}/include" WRF_PATH=/lfs4/BMC/wrfruc/beck \
  NETCDF_LDFLAGS="-L${NETCDF}/lib -lnetcdf -lnetcdff"


cd ../

