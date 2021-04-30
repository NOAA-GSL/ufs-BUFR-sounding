#! /bin/sh

set -x

source ./machine-setup.sh > /dev/null 2>&1

module purge >& /dev/null

module use ~rtrr/RRFS/sorc/EMC_post/modulefiles
module load post/v8.0.0-jet
module use /lfs4/HFIP/hfv3gfs/gwv/l0530/lib/modulefiles
module load bufr

module list


cd ./regional_stnmlist.fd

make clean
make FC=mpiifort

cd ../
