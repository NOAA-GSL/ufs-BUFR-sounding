#! /bin/sh

set -x

#. /usrx/local/prod/lmod/lmod/init/sh

source ./machine-setup.sh > /dev/null 2>&1

module purge >& /dev/null

#module use -a ./regional_sndp.fd/modulefiles
#module load v4.0.0_build

module use ~rtrr/RRFS/sorc/EMC_post/modulefiles
module load post/v8.0.0-jet
module use /lfs4/HFIP/hfv3gfs/gwv/l0530/lib/modulefiles
module load bufr

module list


cd ./regional_sndp.fd

#make clean
#make

make delete
make FC=mpiifort

cd ../




