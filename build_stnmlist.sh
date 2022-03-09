#! /bin/sh

set -x

source ./machine-setup.sh > /dev/null 2>&1

module purge >& /dev/null

module use /gpfs/dell2/emc/modeling/noscrub/Edward.Colon/NOAA_3drtma_bec/sorc/rtma_post.fd/modulefiles
module load post/v8.0.0-wcoss_dell_p3
#module use /lfs4/HFIP/hfv3gfs/gwv/l0530/lib/modulefiles
#module load bufr
module load bufr/11.3.1 
module list


cd ./regional_stnmlist.fd

make clean
make FC="mpif90 -f90=ifort"

cd ../
