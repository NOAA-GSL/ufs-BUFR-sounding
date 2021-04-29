#! /bin/sh

set -x

source ./machine-setup.sh > /dev/null 2>&1

module purge >& /dev/null

#. /usrx/local/prod/lmod/lmod/init/sh
#module purge

#module use -a ./regional_stnmlist.fd/modulefiles
#module load v4.0.0_build

module use ~rtrr/RRFS/sorc/EMC_post/modulefiles
module load post/v8.0.0-jet
module use /lfs4/HFIP/hfv3gfs/gwv/l0530/lib/modulefiles
module load bufr

module list


cd ./regional_stnmlist.fd

make clean
#make
make FC=mpiifort

cd ../
