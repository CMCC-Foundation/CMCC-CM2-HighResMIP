#!/bin/sh
#BSUB -J postrun          # Name of the job.
#BSUB -o postrun_%J.out   # Appends std output to file %J.out.
#BSUB -e postrun_%J.out   # Appends std error to  file %J.out.
#BSUB -q serial_6h             # queue
#BSUB -n 4                 # Number of CPUs
#BSUB -R "span[ptile=4]"   #
#BSUB -N 
#BSUB -W 01:30
#BSUB -L /bin/sh

#######################################################################

set -xv

#######################################################################

which module
stat=$?
set -e
if [ ${stat} -ne 0 ]; then
  source /usr/share/Modules/init/sh
fi
module purge
module purge   # Needed! Not an error!
#
module load SZIP/szip-2.1_int15
module load UDUNITS/udunits-1.12.11_int15
module load UDUNITS/udunits-2.2.19
module load GRIB_API/grib_api-1.13.1
module load HDF5/hdf5-1.8.15-patch1
module load NETCDF/netcdf-C_4.3.3.1-F_4.4.2_C++_4.2.1
module load ESMF/esmf-6.3.0rp1-mpiuni-64-O_int15
module load MAGICS/magics-2.24.7
module load EMOS/emos_000392
module load EMOS/libemos-4.0.6
module load CDO/cdo-1.7.0rc2
module load NCO/nco-4.4.9
module unload INTEL/intel_xe_2013
module unload INTEL/intel_xe_2013.5.192
module load INTEL/intel_xe_2015.3.187
module load IMPI/intel_mpi_5.0.3.048

#######################################################################

OUTDIR="DOUT_S_ROOT"
AUXDIR=`echo ${OUTDIR} | sed -e "s/\/archive\//\/postproc\//" -e "s/\/hist//"`
EXPID="CASE"
CONF="ORCA1"   # redefined later

CDO="cdo"

REBUILD="./rebuild_nemo"

#######################################################################

module list
printenv | sort

#######################################################################

NCPU=`echo ${LSB_HOSTS} | wc -w`

y1="YYYYMMDD1"
y2="YYYYMMDD2"
freqs="1h 6h 1d 5d 1m 1y"
fouts="grid_W grid_T ptrc_T grid_U grid_V diaptr"

nproc=0
for freq in ${freqs} ; do
  for v in ${fouts} ; do

    nproc=`ps -f -u ${USER} | grep "${REBUILD}" | grep ${y1} | wc -l`
    while [ ${nproc} -ge ${NCPU} ]; do
      sleep 10
      nproc=`ps -f -u ${USER} | grep "${REBUILD}" | grep ${y1} | wc -l`
    done

    if [ -f ${EXPID}_${freq}_${y1}_${y2}_${v}.nc -a -f ${EXPID}_${freq}_${y1}_${y2}_${v}_0000.nc ]; then
      rm -f ${EXPID}_${freq}_${y1}_${y2}_${v}.nc
    fi
    if [ -f ${EXPID}_${freq}_${y1}_${y2}_${v}_0000.nc ]; then
      npes=`ls ${EXPID}_${freq}_${y1}_${y2}_${v}_[0-9][0-9][0-9][0-9].nc | wc -w`
        ${REBUILD} ${EXPID}_${freq}_${y1}_${y2}_${v} ${npes} && \
      rm -f ${EXPID}_${freq}_${y1}_${y2}_${v}_[0-9][0-9][0-9][0-9].nc || exit 1 &
    fi
  done
done

wait

if [ -f ${EXPID}_1m_${y1}_${y2}_grid_T.nc ]; then
  fname="${EXPID}_1m_${y1}_${y2}_grid_T.nc"
  freq="1m"
elif [ -f ${EXPID}_1d_${y1}_${y2}_grid_T.nc ]; then
  fname="${EXPID}_1d_${y1}_${y2}_grid_T.nc"
  freq="1d"
elif [ -f ${EXPID}_1y_${y1}_${y2}_grid_T.nc ]; then
  fname="${EXPID}_1y_${y1}_${y2}_grid_T.nc"
  freq="1y"
else
  echo "ERROR: unknown resolution! nlon = ${nlon}"
  exit 1
fi
nlon=`${CDO} griddes ${fname} | grep -i xsize | cut -d'=' -f2 | tr -d '[:blank:]'`

case ${nlon} in
182)
  CONF="ORCA2"
  ;;
362)
  CONF="ORCA1"
  ;;
1442)
  CONF="ORCA025"
  ;;
*)
  echo "ERROR: unknown resolution! nlon = ${nlon}"
  exit 1
  ;;
esac


if [ ! -f ${CONF}_grid_T_griddes.txt -a -f ${EXPID}_${freq}_${y1}_${y2}_grid_T.nc ]; then
    ${CDO} -f nc griddes -selname,sosstsst \
    ${EXPID}_${freq}_${y1}_${y2}_grid_T.nc >${CONF}_grid_T_griddes.txt
fi
if [ ! -f ${CONF}_grid_U_griddes.txt -a -f ${EXPID}_1m_${y1}_${y2}_grid_U.nc ]; then
    ${CDO} -f nc griddes -selname,vozocrtx \
    ${EXPID}_1m_${y1}_${y2}_grid_U.nc >${CONF}_grid_U_griddes.txt
fi
if [ ! -f ${CONF}_grid_V_griddes.txt -a -f ${EXPID}_1m_${y1}_${y2}_grid_V.nc ]; then
    ${CDO} -f nc griddes -selname,vomecrty \
    ${EXPID}_1m_${y1}_${y2}_grid_V.nc >${CONF}_grid_V_griddes.txt
fi
if [ ! -f ${CONF}_grid_T_zaxisdes.txt -a -f ${EXPID}_1m_${y1}_${y2}_grid_T.nc ]; then
    ${CDO} -f nc zaxisdes -selname,votemper \
    ${EXPID}_1m_${y1}_${y2}_grid_T.nc | \
    sed -e "s/generic/depth_below_sea/" >${CONF}_grid_T_zaxisdes.txt
fi
if [ ! -f ${CONF}_grid_W_zaxisdes.txt -a -f ${EXPID}_1m_${y1}_${y2}_grid_W.nc ]; then
    ${CDO} -f nc zaxisdes -selname,vovecrtz \
    ${EXPID}_1m_${y1}_${y2}_grid_W.nc | \
    sed -e "s/generic/depth_below_sea/" >${CONF}_grid_W_zaxisdes.txt
fi

#######################################################################

meshmask="${EXPID}_mesh_mask.nc"

if [ ! -f ${meshmask} ]; then
  if [ -f mesh_mask_0000.nc ]; then
    npes=`ls mesh_mask_[0-9][0-9][0-9][0-9].nc | wc -w`
    ${REBUILD} mesh_mask ${npes} && \
    mv mesh_mask.nc ${meshmask} && \
    rm -f mesh_mask_[0-9][0-9][0-9][0-9].nc || exit 1
  fi
fi

for v in tmaskutil umaskutil vmaskutil ; do
  if [ ! -f ${EXPID}_${v}.nc -a -f ${meshmask} ]; then
    grid=`echo ${v} | cut -c1 | tr '[:lower:]' '[:upper:]'`
    ${CDO} -b F32 setgrid,${CONF}_grid_${grid}_griddes.txt -selname,${v} ${meshmask} ${EXPID}_${v}.nc
  fi
done

for v in tmask umask vmask ; do
  if [ ! -f ${EXPID}_${v}.nc -a -f ${meshmask} ]; then
    grid=`echo ${v} | cut -c1 | tr '[:lower:]' '[:upper:]'`
    ${CDO} -f nc setzaxis,${CONF}_grid_T_zaxisdes.txt -setgrid,${CONF}_grid_${grid}_griddes.txt \
      -mul -selname,${v} ${meshmask} ${EXPID}_${v}util.nc ${EXPID}_${v}.nc
  fi
done

for v in t u v ; do
  g=`echo ${v} | tr '[:lower:]' '[:upper:]'`
  if [ ! -f ${EXPID}_${v}area.nc -a -f ${meshmask} ]; then
    ${CDO} selname,e1${v} ${meshmask} ${EXPID}_e1${v}.nc
    ${CDO} setmissval,-1.e30 -ifthen ${EXPID}_${v}maskutil.nc ${EXPID}_e1${v}.nc ${EXPID}_e1${v}_masked.nc
    ${CDO} selname,e2${v} ${meshmask} ${EXPID}_e2${v}.nc
    ${CDO} setmissval,-1.e30 -ifthen ${EXPID}_${v}maskutil.nc ${EXPID}_e2${v}.nc ${EXPID}_e2${v}_masked.nc
    ${CDO} -f nc -setgrid,${CONF}_grid_${g}_griddes.txt -chname,e1${v},cell_area -mul \
      ${EXPID}_e1${v}.nc ${EXPID}_e2${v}.nc ${EXPID}_${v}area.nc
    ${CDO} -f nc setmissval,-1.e30 -ifthen ${EXPID}_${v}maskutil.nc \
      ${EXPID}_${v}area.nc ${EXPID}_${v}area_masked.nc
    ${CDO} fldsum ${EXPID}_${v}area_masked.nc ${EXPID}_${v}area_int.nc
  fi

  if [ ! -f ${EXPID}_${v}volume.nc -a -f ${meshmask} ]; then
    ${CDO} selname,e3${v} ${meshmask} ${EXPID}_e3${v}.nc
    ${CDO} setmissval,-1.e30 -ifthen ${EXPID}_${v}maskutil.nc ${EXPID}_e3${v}.nc ${EXPID}_e3${v}_masked.nc
    ${CDO} -f nc -setzaxis,${CONF}_grid_T_zaxisdes.txt -setgrid,${CONF}_grid_${g}_griddes.txt \
      -chname,e3${v},cell_volume -mul ${EXPID}_e3${v}.nc ${EXPID}_${v}area.nc ${EXPID}_${v}volume.nc
    ${CDO} -f nc setmissval,-1.e30 -ifthen ${EXPID}_${v}mask.nc \
      ${EXPID}_${v}volume.nc ${EXPID}_${v}volume_masked.nc
    ${CDO} vertsum -fldsum ${EXPID}_${v}volume_masked.nc ${EXPID}_${v}volume_int.nc
  fi
done

while [ -f ../${EXPID}.st_archive.lock ]; do
  sleep 10
done
#

if [ ! -d ${AUXDIR} ]; then
  mkdir -p ${AUXDIR}
fi

mv ${EXPID}_*${y1}_${y2}*[z-aA-Z].nc ${OUTDIR}/ || exit 1

for v in T_griddes U_griddes V_griddes T_zaxisdes W_zaxisdes ; do
  if [ ! -f ${AUXDIR}/${CONF}_grid_${v}.txt ]; then
    cp ${CONF}_grid_${v}.txt ${AUXDIR}/
  fi
done
#
if [ ! -f ${AUXDIR}/${meshmask} ]; then
  cp ${meshmask} ${AUXDIR}/
fi
#
for v in tmaskutil umaskutil vmaskutil tmask umask vmask ; do
  if [ ! -f ${AUXDIR}/${EXPID}_${v}.nc ]; then
    cp ${EXPID}_${v}.nc ${AUXDIR}/
  fi
done
#
for v in t u v ; do
  if [ ! -f ${AUXDIR}/${EXPID}_${v}area.nc -a -f ${EXPID}_${v}area.nc ]; then
    cp ${EXPID}_${v}area.nc ${AUXDIR}/
  fi
  if [ ! -f ${AUXDIR}/${EXPID}_${v}area_masked.nc -a -f ${EXPID}_${v}area_masked.nc ]; then
    cp ${EXPID}_${v}area_masked.nc ${AUXDIR}/
  fi
  if [ ! -f ${AUXDIR}/${EXPID}_${v}area_int.nc -a -f ${EXPID}_${v}area_int.nc ]; then
    cp ${EXPID}_${v}area_int.nc ${AUXDIR}/
  fi
  if [ ! -f ${AUXDIR}/${EXPID}_${v}volume.nc -a -f ${EXPID}_${v}volume.nc ]; then
    cp ${EXPID}_${v}volume.nc ${AUXDIR}/
  fi
  if [ ! -f ${AUXDIR}/${EXPID}_${v}volume_masked.nc -a -f ${EXPID}_${v}volume_masked.nc ]; then
    cp ${EXPID}_${v}volume_masked.nc ${AUXDIR}/
  fi
  if [ ! -f ${AUXDIR}/${EXPID}_${v}volume_int.nc -a -f ${EXPID}_${v}volume_int.nc ]; then
    cp ${EXPID}_${v}volume_int.nc ${AUXDIR}/
  fi
  for m in e1 e2 e3 ; do
    if [ ! -f ${AUXDIR}/${EXPID}_${m}${v}.nc -a -f ${EXPID}_${m}${v}.nc ]; then
      cp ${EXPID}_${m}${v}.nc ${AUXDIR}/
    fi
    if [ ! -f ${AUXDIR}/${EXPID}_${m}${v}_masked.nc -a -f ${EXPID}_${m}${v}_masked.nc ]; then
      cp ${EXPID}_${m}${v}_masked.nc ${AUXDIR}/
    fi
  done
done
#

