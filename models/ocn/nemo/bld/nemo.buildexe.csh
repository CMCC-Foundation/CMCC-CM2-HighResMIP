#! /bin/csh -fv

set objdir = $OBJROOT/ocn/obj; cd ${objdir}
set comp = 'unknown'
if (${COMP_INTERFACE} == 'MCT' ) set comp = mct
if (${COMP_INTERFACE} == 'ESMF') set comp = esmf

#######################################################################
#
# NEMO CPP keys
#
#######################################################################

# Resolution independent CPP keys
set nemodefs = "-Dkey_diaar5 -Dkey_diahth -Dkey_dynspg_flt -Dkey_ldfslp -Dkey_trabbl -Dkey_traldf_c2d -Dkey_zdfddm -Dkey_zdftke -Dkey_dtatem -Dkey_dtasal -Dkey_mpp_mpi -Dkey_iomput -Dkey_coupled -DCCSMCOUPLED"

# Resolution dependent CPP keys
set res_cpp =
switch (${OCN_GRID})
case "tn2v?":
  set res_cpp = "-Dkey_orca_r2 -Dkey_diaeiv -Dkey_dynldf_c3d -Dkey_traldf_eiv -Dkey_zdftmx"
  breaksw
case "tn1v?":
  set res_cpp = "-Dkey_orca_r1 -Dkey_diaeiv -Dkey_dynldf_c3d -Dkey_traldf_eiv -Dkey_zdftmx"
  breaksw
case "tn0.5v?":
  set res_cpp = "-Dkey_orca_r05"
  echo "ORCA05 configuration not yet implemented in CESM!"
  exit 2
  breaksw
case "tn0.25v?":
  set res_cpp = "-Dkey_orca_r025 -Dkey_dynldf_c2d"
  breaksw
default
  echo "Unknown resolution ${OCN_GRID}"
  exit 1
  breaksw
endsw

# TOP (passive tracers) module
set top_cpp
if ( "${NEMO_TOP_MODULES}" != "" ) then
  set top_cpp = "-Dkey_top"
  foreach m (`echo ${NEMO_TOP_MODULES} | tr ',' ' '`)
    switch (${m})
      case "age":
        set top_cpp = "${top_cpp} -Dkey_my_trc"   # Ideal age tracer active
	breaksw
      case "cfc":
        set top_cpp = "${top_cpp} -Dkey_cfc"      # CFCs tracers active
	breaksw
      default:
        echo "ERROR: unknown NEMO TOP module: ${m}!"
	exit 1
        breaksw
    endsw
endif

# NEMO CPP keys
# TODO: add XML variables to handle NEMO CPPs
set nemodefs = "${res_cpp} ${top_cpp} ${nemodefs}"

#######################################################################

# Build the NEMO library

\cat >! Filepath << EOF1
$CASEROOT/SourceMods/src.$COMP_OCN
$CODEROOT/ocn/nemo/drivers/cpl_$comp
$CODEROOT/ocn/nemo/drivers/cpl_share
$CODEROOT/ocn/nemo/NEMOGCM/EXTERNAL/IOIPSL/src
$CODEROOT/ocn/nemo/NEMOGCM/EXTERNAL/XMLIO_SERVER/src/XMLIO
$CODEROOT/ocn/nemo/NEMOGCM/EXTERNAL/XMLIO_SERVER/src/IOSERVER
$CODEROOT/ocn/nemo/NEMOGCM/EXTERNAL/XMLF90/src/cml
$CODEROOT/ocn/nemo/NEMOGCM/EXTERNAL/XMLF90/src/dom
$CODEROOT/ocn/nemo/NEMOGCM/EXTERNAL/XMLF90/src/sax
$CODEROOT/ocn/nemo/NEMOGCM/EXTERNAL/XMLF90/src/strings
$CODEROOT/ocn/nemo/NEMOGCM/EXTERNAL/XMLF90/src/wxml
$CODEROOT/ocn/nemo/NEMOGCM/EXTERNAL/XMLF90/src/xpath
$CODEROOT/ocn/nemo/NEMOGCM/NEMO/OPA_SRC
$CODEROOT/ocn/nemo/NEMOGCM/NEMO/OPA_SRC/DIA
$CODEROOT/ocn/nemo/NEMOGCM/NEMO/OPA_SRC/DOM
$CODEROOT/ocn/nemo/NEMOGCM/NEMO/OPA_SRC/ASM
$CODEROOT/ocn/nemo/NEMOGCM/NEMO/OPA_SRC/BDY
$CODEROOT/ocn/nemo/NEMOGCM/NEMO/OPA_SRC/C1D
$CODEROOT/ocn/nemo/NEMOGCM/NEMO/OPA_SRC/DTA
$CODEROOT/ocn/nemo/NEMOGCM/NEMO/OPA_SRC/DYN
$CODEROOT/ocn/nemo/NEMOGCM/NEMO/OPA_SRC/FLO
$CODEROOT/ocn/nemo/NEMOGCM/NEMO/OPA_SRC/IOM
$CODEROOT/ocn/nemo/NEMOGCM/NEMO/OPA_SRC/LBC
$CODEROOT/ocn/nemo/NEMOGCM/NEMO/OPA_SRC/LDF
$CODEROOT/ocn/nemo/NEMOGCM/NEMO/OPA_SRC/OBC
$CODEROOT/ocn/nemo/NEMOGCM/NEMO/OPA_SRC/OBS
$CODEROOT/ocn/nemo/NEMOGCM/NEMO/OPA_SRC/SBC
$CODEROOT/ocn/nemo/NEMOGCM/NEMO/OPA_SRC/SOL
$CODEROOT/ocn/nemo/NEMOGCM/NEMO/OPA_SRC/TRA
$CODEROOT/ocn/nemo/NEMOGCM/NEMO/OPA_SRC/TRD
$CODEROOT/ocn/nemo/NEMOGCM/NEMO/OPA_SRC/ZDF
EOF1

# Add TOP modules path if TOP is active
if ( "${NEMO_TOP_MODULES}" != "" ) then
\cat >> Filepath << EOF1
$CODEROOT/ocn/nemo/NEMOGCM/NEMO/TOP_SRC
$CODEROOT/ocn/nemo/NEMOGCM/NEMO/TOP_SRC/TRP
$CODEROOT/ocn/nemo/NEMOGCM/NEMO/TOP_SRC/CFC
$CODEROOT/ocn/nemo/NEMOGCM/NEMO/TOP_SRC/C14b
$CODEROOT/ocn/nemo/NEMOGCM/NEMO/TOP_SRC/SED
$CODEROOT/ocn/nemo/NEMOGCM/NEMO/TOP_SRC/LOBSTER
$CODEROOT/ocn/nemo/NEMOGCM/NEMO/TOP_SRC/PISCES
$CODEROOT/ocn/nemo/NEMOGCM/NEMO/TOP_SRC/MY_TRC
EOF1
endif


${GMAKE} complib -j ${GMAKE_J} MODEL=nemo COMPLIB=${LIBROOT}/libocn.a -f ${CASETOOLS}/Makefile USER_CPPDEFS="${nemodefs}" MACFILE=${CASEROOT}/Macros.${MACH} || exit 2

# Compile NEMO REBUILD tool

if !( -d ${EXEROOT}/${COMP_OCN}/REBUILD_NEMO ) then
  cp -r ${CODEROOT}/ocn/${COMP_OCN}/NEMOGCM/EXTERNAL/REBUILD_NEMO ${EXEROOT}/ocn/ || exit 2
endif
cd ${EXEROOT}/ocn/REBUILD_NEMO/src  || exit 2
if (! -f Macros ) then
  ${CCSM_MACHDIR}/configure -mpilib mpi-serial -mach ${MACH} || exit 2
endif
if (-e build.csh ) then
  ./build.csh default || exit 2
else
  echo "Unable to build NEMO rebuild tool"
  exit 2
endif

