#! /bin/csh -f

########################################################################
#
# NEMO nemo.buildnml.csh
#
# P.G. Fogli, CMCC - Jun 2012: CESM 1.1.1, initial version for calyspo
# P.G. Fogli, CMCC - Oct 2013: CESM 1.2.0, updated version, athena port
#
########################################################################

########################################################################
# 1. NEMO namelist & XIOS definitions
# TODO: In CESM namelists & I/O related files should be generated
# by the script build-namelist based on the specific compset/resolution
# while user specific changes should go in the file user_nl_nemo. TO BE DONE!
# For now NEMO namelist and XIOS xml files are copied from existing
# template files in ${CODEROOT}/ocn/${COMP_OCN}/bld to
# ${CASEBUILD}/nemoconf and ${RUNDIR} .
# The user should modify the files in ${CASEBUILD}/nemoconf after case
# build and before case run!
# (for run-time namelist changes see below).
########################################################################

# First make a local copy in ${CASEBUILD}/nemoconf which the user can modify

if !(-d ${CASEBUILD}/nemoconf) mkdir -p ${CASEBUILD}/nemoconf

cd ${CASEBUILD}/nemoconf || exit 1

# Copy NEMO template namelist (modified later at run-time, see below)
if !( -f namelist ) then
  set f = "${CODEROOT}/ocn/${COMP_OCN}/bld/namelist.${OCN_GRID}"
  if ( -f ${f} ) then
    cp ${f} namelist
  else
    echo "ERROR: file not found ${f}"
    exit 1
  endif
endif

# TOP passive tracers module
set NEMO_TOP_AGE = 0
set NEMO_TOP_CFC = 0
set NEMO_TOP_C14 = 0
if ( "${NEMO_TOP_MODULES}" != "" ) then
  foreach m (`echo ${NEMO_TOP_MODULES} | tr ',' ' '`)
    switch (${m})
      case "age":
        # Copy TOP module template namelist
        set NEMO_TOP_AGE = 1
	if !( -f namelist_top ) then
          set f = "${CODEROOT}/ocn/${COMP_OCN}/bld/namelist_top.${OCN_GRID}"
          if ( -f ${f} ) then
            cp ${f} namelist_top
          else
            echo "ERROR: file not found ${f}"
            exit 1
          endif
	endif
        breaksw
      case "cfc":
        # Copy CFCs module template namelist
        set NEMO_TOP_CFC = 1
	if !( -f namelist_cfc ) then
          set f = "${CODEROOT}/ocn/${COMP_OCN}/bld/namelist_cfc.${OCN_GRID}"
          if ( -f ${f} ) then
            cp ${f} namelist_cfc
          else
            echo "ERROR: file not found ${f}"
            exit 1
          endif
	endif
        breaksw
# C14 not implemented/tested yet
#      case "c14":
#        set NEMO_TOP_C14 = 1
#        if !( -f namelist_c14 ) then
#          set f = "${CODEROOT}/ocn/${COMP_OCN}/bld/namelist_c14.${OCN_GRID}"
#          if ( -f ${f} ) then
#            cp ${f} namelist_c14
#          else
#            echo "ERROR: file not found ${f}"
#            exit 1
#          endif
#        endif
#        breaksw
      default:
        echo "ERROR: unknown or unsupported NEMO TOP module: ${m}!"
	exit 1
        breaksw
    endsw
  end
endif

# Copy XIO server files
if !( -f xmlio_server.def ) then
  set f = "${CODEROOT}/ocn/${COMP_OCN}/bld/xmlio_server.def"
  if ( -f ${f} ) then
    cp ${CODEROOT}/ocn/${COMP_OCN}/bld/xmlio_server.def .
  else
    echo "ERROR: file not found ${f}"
    exit 1
  endif
endif
# NEMO XIO definition (XIO Server XML specification file)
# Edit ${CASEBUILD}/nemo.iodef.xml iodef.xml  as you need!!!
if !( -f iodef.xml ) then
  set f = "${CODEROOT}/ocn/${COMP_OCN}/bld/iodef.xml"
  if ( -f ${f} ) then
    # Mostly monthly output frequency
    cp ${CODEROOT}/ocn/${COMP_OCN}/bld/iodef.xml .
  else
    echo "ERROR: file not found ${f}"
    exit 1
  endif
endif

# Now copy namelist & XIO server files to the run dir ${RUNDIR}

cd ${RUNDIR} || exit 1

cp ${CASEBUILD}/nemoconf/namelist* . || exit 1
cp ${CASEBUILD}/nemoconf/xmlio_server.def . || exit 1
cp ${CASEBUILD}/nemoconf/iodef.xml . || exit 1

# NEMO postrun template script for output files rebuild
if !( -f ${CASEBUILD}/nemoconf/postrun.tpl ) then
  cp ${CODEROOT}/ocn/${COMP_OCN}/bld/postrun.tpl ${CASEBUILD}/nemoconf/ || exit 1
endif

########################################################################
# 2. IC, BC, grids/masks & restart files management
# TODO: I.C., B.C., restart files and grids/masks are copied/linked in
# the run dir by this script. As of version 3.3 NEMO doesn't support
# IC/BC and restart file in a directory different from the working dir
# so we can't use the default CESM nemo.input_data_list file method.
########################################################################

set NEMO_IN = "${DIN_LOC_ROOT}/ocn/${COMP_OCN}/${OCN_GRID}"

# I.C., B.C., grids coordinates & metrics
# Comment out unused files!
switch (${OCN_GRID})
case "tn2v?":
  set AHMCOEF       = "forcing/ahmcoef"
  set COORDINATES   = "grid/coordinates"
  set SUBBASINS     = "grid/subbasins"
  set GEOTHERMAL    = "forcing/geothermal_heating"
  set K1ROWDRG      = "forcing/K1rowdrg"
  set M2ROWDRG      = "forcing/M2rowdrg"
  set MASK_ITF      = "grid/mask_itf"
#  set DISTCOAST     = "dist.coast"
  set CHLOROPHYLL   = "forcing/chlorophyll"
# WOA 2009/PHC3.0
  set POTEMPERATURE = "ic/data_1m_potential_temperature_nomask"
  set SALINITY      = "ic/data_1m_salinity_nomask"
  set BATHYMETRY    = "grid/bathy_meter_closea_no_NAGL"
#  set SSSR          = ""
  breaksw
case "tn1v1":
  set AHMCOEF       = "forcing/ahmcoef"
  set COORDINATES   = "grid/coordinates_ukorca1"
  set SUBBASINS     = "grid/basinmask_050308_UKMO"
#  set GEOTHERMAL    = ""
  set K1ROWDRG      = "forcing/K1rowdrg_R1_modif"
  set M2ROWDRG      = "forcing/M2rowdrg_R1_modif"
  set MASK_ITF      = "grid/mask_itf_ORCA1_new"
#  set DISTCOAST     = ""
  set CHLOROPHYLL   = "forcing/chlorophyll-ORCA1_1m"
# Levitus 2001
#  set POTEMPERATURE = "ic/potemp_1m_z46_nomask"
#  set SALINITY      = "ic/salin_1m_z46_nomask"
# WOA 2009
#  set POTEMPERATURE = "ic/WOA2009_ORCA1_L46_1m_potential_temperature_nomask"
#  set SALINITY      = "ic/WOA2009_ORCA1_L46_1m_salinity_nomask"
# WOA 2009/PHC3.0
  set POTEMPERATURE = "ic/WOA2009_PHC3.0_ORCA1_L46_1m_potential_temperature_nomask"
  set SALINITY      = "ic/WOA2009_PHC3.0_ORCA1_L46_1m_salinity_nomask"
  set BATHYMETRY    = "grid/bathy_meter_121126_CMCC"
#  set SSSR          = "forcing/ORCA1_PHC2_1m_sss_nomask_patched"
#  set SSSR          = "forcing/ORCA1_PHC2_1m_sss_nomask"
  breaksw
case "tn0.5v?":
  echo "Error: Resolution ${OCN_GRID} not implemented yet!"
  exit 1
  breaksw
case "tn0.25v?":
#  set AHMCOEF       = "forcing/"
  set COORDINATES   = "grid/orca025_coordinates_280809"
#  set SUBBASINS     = "grid/"
#  set GEOTHERMAL    = "forcing/"
#  set K1ROWDRG      = "forcing/"
#  set M2ROWDRG      = "forcing/"
#  set MASK_ITF      = "grid/"
#  set DISTCOAST     = "dist.coast"
  set CHLOROPHYLL   = "forcing/chlorophyll"
# Levitus 98/PHC2.1
  set POTEMPERATURE = "ic/data_1m_potential_temperature_nomask"
  set SALINITY      = "ic/data_1m_salinity_nomask"
  set BATHYMETRY    = "grid/bathy_meter"
#  set SSSR          = ""
  breaksw
default:
  echo "Error: Resolution not supported: ${OCN_GRID}!"
  exit 1
  breaksw
endsw
#
# Resolution independent files
if ( ${NEMO_TOP_CFC} == 1 ) then
  set CFC = "../resolution_independent/forcing/cfc1112"
endif
#
# Copy/link files
if (-e data_1m_potential_temperature_nomask.nc == 0 && ${?POTEMPERATURE} == 1 && "${CONTINUE_RUN}" == "FALSE" ) then
  ln -s ${NEMO_IN}/${POTEMPERATURE}.nc    data_1m_potential_temperature_nomask.nc
endif
if (-e data_1m_salinity_nomask.nc == 0 && ${?SALINITY} == 1 ) then
  ln -s ${NEMO_IN}/${SALINITY}.nc    data_1m_salinity_nomask.nc
endif
if (-e chlorophyll.nc == 0 && ${?CHLOROPHYLL} == 1 ) then
  ln -s ${NEMO_IN}/${CHLOROPHYLL}.nc chlorophyll.nc
endif
if (-e ahmcoef == 0 && ${?AHMCOEF} == 1 ) then
  ln -s ${NEMO_IN}/${AHMCOEF}        ahmcoef
endif
if (-e bathy_meter.nc == 0 && ${?BATHYMETRY} == 1 ) then
  ln -s ${NEMO_IN}/${BATHYMETRY}.nc  bathy_meter.nc
endif
if (-e coordinates.nc == 0 && ${?COORDINATES} == 1 ) then
  ln -s ${NEMO_IN}/${COORDINATES}.nc coordinates.nc
endif
if (-e subbasins.nc == 0 && ${?SUBBASINS} == 1 ) then
  ln -s ${NEMO_IN}/${SUBBASINS}.nc   subbasins.nc
endif
if (-e geothermal_heating.nc == 0 && ${?GEOTHERMAL} == 1 ) then
  ln -s ${NEMO_IN}/${GEOTHERMAL}.nc  geothermal_heating.nc
endif
if (-e K1rowdrg.nc == 0 && ${?K1ROWDRG} == 1 ) then
  ln -s ${NEMO_IN}/${K1ROWDRG}.nc    K1rowdrg.nc
endif
if (-e M2rowdrg.nc == 0 && ${?M2ROWDRG} == 1 ) then
  ln -s ${NEMO_IN}/${M2ROWDRG}.nc    M2rowdrg.nc
endif
if (-e mask_itf.nc == 0 && ${?MASK_ITF} == 1 ) then
  ln -s ${NEMO_IN}/${MASK_ITF}.nc    mask_itf.nc
endif
if (-e distcoast.nc == 0 && ${?DISTCOAST} == 1 ) then
  ln -s ${NEMO_IN}/${DISTCOAST}.nc   distcoast.nc
endif
if (-e sss_1m.nc == 0 && ${?SSSR} == 1 ) then
  ln -s ${NEMO_IN}/${SSSR}.nc   sss_1m.nc
endif
if (-e cfc1112.atm == 0 && ${?CFC} == 1 ) then
  ln -s ${NEMO_IN}/${CFC}.atm   .
endif

#######################################################################
# last timestep and last date, assume initial run
# FIXME: should come from env_*.xml ?
set lastkt = 0
set ndastp = `echo ${RUN_STARTDATE} | sed -e "s/[^0-9]//g"`
#
set hyb_rest = 1   # 1 => hybrid run, starts from restart file
if ("${CONTINUE_RUN}" == "TRUE") then
  set rstart = ".true."
  set rstctl = 2
  set msh = 0        # msh>0 => write out a mesh mask file
else
  set rstart = ".false."
  set rstctl = 0
  set msh = 1        # msh>0 => write out a mesh mask file
endif
#
# Restart files
#
if ("${CONTINUE_RUN}" == "TRUE" || "${RUN_TYPE}" != "startup") then
  #
  if ("${CONTINUE_RUN}" == "FALSE" && "${RUN_TYPE}" != "startup") then
    # branch or hybrid run && CONTINUE_RUN=FALSE
    set EXPID = ${RUN_REFCASE}
  else
    # startup, branch or hybrid run & CONTINUE_RUN=TRUE
    set EXPID = ${CASE}
  endif
  #
  # Looks for restart files (1 file per PE case)
  set flist = `ls -1rt ${EXPID}_[0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9]_restart_[0-9][0-9][0-9][0-9].nc | tail -1`
  if (${?flist} == 1 && "x${flist}" != "x" ) then
    set flist = `ls -1rt ${EXPID}_[0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9]_restart_[0-9][0-9][0-9][0-9].nc`
    # TODO: check consistency with $NTASKS_OCN
    foreach f (${flist})
      set pe = `echo ${f} | sed -e "s/\.nc//" | awk -F '_' '{ print $NF ; }'`
#      if (-e restart_${pe}.nc) then
        rm -f restart_${pe}.nc
#      endif
      ln -s ${f} restart_${pe}.nc
    end
    set rest_file = "restart_0000.nc"
  else
    # Looks for restart file (single restart file case)
    set flist = `ls -1rt ${EXPID}_[0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9]_restart.nc | tail -1`
    if (${?flist} == 1 && "x${flist}" != "x") then
      rm -f restart.nc
      ln -s ${flist} restart.nc
      set rest_file = "restart.nc"
    else
      # No restart file found
      if ("${CONTINUE_RUN}" == "TRUE") then
        # No restart file & CONTINUE_RUN==TRUE ==> ERROR
        echo "ERROR: CONTINUE_RUN=${CONTINUE_RUN} (RUN_TYPE=${RUN_TYPE}) BUT NO RESTART FILE FOUND in ${RUNDIR}"
	exit 1
      else if ("${RUN_TYPE}" == "hybrid") then
        # No restart file & hybrid run ==> assume NEMO starts from initial conditions
        echo "WARNING: RUN_TYPE=${RUN_TYPE}, CONTINUE_RUN=${CONTINUE_RUN} AND NO RESTART FILE FOUND in ${RUNDIR}"
	echo " ******  Assume NEMO starts from initial conditions  ******"
	set hyb_rest = 0   # 0 => hybrid run, starts from initial conditions
      else
        # No restart file & branch run ==> ERROR
        echo "ERROR: RUN_TYPE=${RUN_TYPE}, CONTINUE_RUN=${CONTINUE_RUN} BUT NO RESTART FILE FOUND in ${RUNDIR}"
        exit 1
      endif
    endif
    #
  endif
  #
  # Restart files for TOP module
  if ( "${NEMO_TOP_MODULES}" != "" ) then
    # Looks for restart files (1 file per PE case)
    set flist = `ls -1rt ${EXPID}_[0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9]_restart_trc_[0-9][0-9][0-9][0-9].nc | tail -1`
    if (${?flist} == 1 && "x${flist}" != "x") then
      set flist = `ls -1rt ${EXPID}_[0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9]_restart_trc_[0-9][0-9][0-9][0-9].nc`
      # TODO: check consistency with $NTASKS_OCN
      foreach f (${flist})
        set pe = `echo ${f} | sed -e "s/\.nc//" | awk -F '_' '{ print $NF ; }'`
#        if (-e restart_trc_${pe}.nc) then
          rm -f restart_trc_${pe}.nc
#        endif
        ln -s ${f} restart_trc_${pe}.nc
      end
    else
      # Looks for restart file (single restart file case)
      set flist = `ls -1rt ${EXPID}_[0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9]_restart_trc.nc | tail -1`
      if (${?flist} == 1 && "x${flist}" != "x") then
        rm -f restart_trc.nc
        ln -s ${flist} restart_trc.nc
      else
        # No restart file found
        if ("${CONTINUE_RUN}" == "TRUE") then
          # No restart file & CONTINUE_RUN==TRUE ==> ERROR
          echo "ERROR: CONTINUE_RUN=${CONTINUE_RUN} (RUN_TYPE=${RUN_TYPE}) & NEMO_TOP_MODULES=${NEMO_TOP_MODULES} BUT NO TRACER RESTART FILE FOUND in ${RUNDIR}"
  	  exit 1
        else if ("${RUN_TYPE}" == "hybrid") then
          # No restart file & hybrid run ==> assume NEMO starts from initial conditions
          echo "WARNING: RUN_TYPE=${RUN_TYPE}, CONTINUE_RUN=${CONTINUE_RUN} & NEMO_TOP_MODULES=${NEMO_TOP_MODULES} BUT NO TRACER RESTART FILE FOUND in ${RUNDIR}"
  	  echo " ******  Assume NEMO starts from initial conditions  ******"
        else
          # No restart file & branch run ==> ERROR
          echo "ERROR: RUN_TYPE=${RUN_TYPE}, CONTINUE_RUN=${CONTINUE_RUN} & NEMO_TOP_MODULES=${NEMO_TOP_MODULES} BUT NO TRACER RESTART FILE FOUND in ${RUNDIR}"
          exit 1
        endif
      endif
      #
    endif
  endif
  # last timestep and last date
  # FIXME: should come from env_*.xml ?
  if (${?rest_file} == 1) then
    set lastkt = `ncdump -v kt ${rest_file} | grep "kt.*=" | cut -d'=' -f2 | sed -e "s/[^0-9]//g"`
    set ndastp = `ncdump -v ndastp ${rest_file} | grep "ndastp.*=" | cut -d'=' -f2 | sed -e "s/[^0-9]//g"`
    set rstart = ".true."
    set rstctl = 2
    unset rest_file
  endif
  # branch/hybrid run: RUN_REFDATE overrides RUN_STARTDATE & restart file date
  if ("${CONTINUE_RUN}" == "FALSE") then
    set ndastp = `echo ${RUN_REFDATE} | sed -e "s/[^0-9]//g"`
    set rstctl = 1
  endif
  #
endif

########################################################################
# 3. date management
########################################################################

# get last day of previous run
set pday = `expr ${ndastp} % 100`
set pmonth = `expr ${ndastp} / 100 % 100`
set pyear = `expr ${ndastp} / 100 / 100`
#
# Compute last day of the previous & current months
set lastday  = 31
set nlastday = 31
set leapyear = 0
if (${pmonth} == 2) then
   set lastday  = 28
   set nlastday = 31
   if ("${CALENDAR}" == "GREGORIAN") then
      # TODO: leap year calendar not implemented yet!
      set leapyear = 0
   endif
else if (${pmonth} == 4 || ${pmonth} == 6 || ${pmonth} == 9 || ${pmonth} == 11) then
  set lastday  = 30
  set nlastday = 31
else if (${pmonth} == 3 || ${pmonth} == 5 || ${pmonth} == 8 || ${pmonth} == 10) then
  set nlastday = 30
endif
@ lastday = ${lastday} + ${leapyear}
#
# Compute next model day (start date)
@ nday = ${pday}
if ( "${CONTINUE_RUN}" == "TRUE" ) then
  @ nday = ${nday} + 1
endif
set nmonth = ${pmonth}
set nyear = ${pyear}
if (${nday} > ${lastday}) then
   set nday = 1
   @ nmonth = ${pmonth} + 1
endif
if (${nmonth} > 12) then
   set nmonth = 1
   @ nyear = ${pyear} + 1
endif
# compute date0
set date0 = `printf "%04d%02d%02d" ${nyear} ${nmonth} ${nday}`

# Compute next model day + 1 (start date + 1)
@ nnday = ${nday} + 1
set nnmonth = ${nmonth}
set nnyear = ${nyear}
if (${nnday} > ${nlastday}) then
   set nnday = 1
   @ nnmonth = ${nmonth} + 1
endif
if (${nnmonth} > 12) then
   set nnmonth = 1
   @ nnyear = ${nyear} + 1
endif
# compute date1
set date1 = `printf "%04d%02d%02d" ${nnyear} ${nnmonth} ${nnday}`

#######################################################################
# 4. NEMO time step length, resolution dependent
#######################################################################

switch (${OCN_GRID})
case "tn2v?":
#  @ DTSEC = 5760   # Default ORCA2 time step length
  @ DTSEC = 5400   # ORCA2 time step length for coupling
  breaksw
case "tn1v?":
  @ DTSEC = 3600   # Default ORCA1 time step length
  breaksw
case "tn0.5v?":
#  @ DTSEC = 5760
  echo "Error: Resolution ${OCN_GRID} not implemented yet!"
  exit 1
  breaksw
case "tn0.25v?":
  @ DTSEC = 1080
#  echo "Error: Resolution ${OCN_GRID} not implemented yet!"
#  exit 1
  breaksw
default:
  echo "Error: Resolution not supported: ${OCN_GRID}!"
  exit 1
endsw

# number of time steps per day
@ ts_per_day = 86400 / ${DTSEC}
@ ts_mod = 86400 % ${DTSEC}
if (${ts_mod} != 0) then
  echo "Error: the time step length ${DTSEC} is not an integer divisor of the day length!"
  exit 1
endif

# Coupling time step related options
@ nnfsbc = 0
switch (${NCPL_BASE_PERIOD})
case "*day*":
  switch (${OCN_NCPL})
  case "1":
    # Daily OCN-CPL coupling
    @ nnfsbc = ${ts_per_day}
    breaksw
  case "4":
    if ("${OCN_GRID}" != "tn2v1" && "${OCN_GRID}" != "tn1v1") then
      echo "Error: wrong OCN_NCPL value (${OCN_NCPL}) for grid ${OCN_GRID} with time step ${DTSEC}!"
      exit 1
    endif
    # 6 hours OCN-CPL coupling
    @ nnfsbc = ${ts_per_day} / ${OCN_NCPL}
    breaksw
  case "6":
    if ("${OCN_GRID}" != "tn2v1" && "${OCN_GRID}" != "tn1v1") then
      echo "Error: wrong OCN_NCPL value (${OCN_NCPL}) for grid ${OCN_GRID} with time step ${DTSEC}!"
      exit 1
    endif
    # 4 hours OCN-CPL coupling
    @ nnfsbc = ${ts_per_day} / ${OCN_NCPL}
    breaksw
  case "8":
    if ("${OCN_GRID}" != "tn2v1" && "${OCN_GRID}" != "tn1v1" && "${OCN_GRID}" != "tn0.25v1") then
      echo "Error: wrong OCN_NCPL value (${OCN_NCPL}) for grid ${OCN_GRID} with time step ${DTSEC}!"
      exit 1
    endif
    # 3 hours OCN-CPL coupling
    @ nnfsbc = ${ts_per_day} / ${OCN_NCPL}
    breaksw
  case "12":
    if ("${OCN_GRID}" != "tn1v1") then
      echo "Error: wrong OCN_NCPL value (${OCN_NCPL}) for grid ${OCN_GRID} with time step ${DTSEC}!"
      exit 1
    endif
    # 2 hours OCN-CPL coupling
    @ nnfsbc = ${ts_per_day} / ${OCN_NCPL}
    breaksw
  case "16":
    if ("${OCN_GRID}" != "tn0.25v1") then
      echo "Error: wrong OCN_NCPL value (${OCN_NCPL}) for grid ${OCN_GRID} with time step ${DTSEC}!"
      exit 1
    endif
    # 1.5 hours OCN-CPL coupling
    @ nnfsbc = ${ts_per_day} / ${OCN_NCPL}
    breaksw
  case "24":
    if ("${OCN_GRID}" != "tn1v1") then
      echo "Error: wrong OCN_NCPL value (${OCN_NCPL}) for grid ${OCN_GRID} with time step ${DTSEC}!"
      exit 1
    endif
    # 1 hour OCN-CPL coupling
    @ nnfsbc = ${ts_per_day} / ${OCN_NCPL}
    breaksw
  default:
    echo "Error: OCN_NCPL ${OCN_NCPL} invalid or not implemented yet!"
    exit 1
  endsw
  breaksw
default:
  echo "Error: NCPL_BASE_PERIOD ${OCN_GRID} not implemented yet!"
  exit 1
endsw


########################################################################
# 5. Compute the run length
# TODO: leap year calendar not implemented yet!
# TODO: date/time computation requires a more advanced language than the shell
# (g)awk, perl or what else?
########################################################################

@ nit000 = ${lastkt} + 1
switch (${STOP_OPTION})
case "*day*":
  @ nitend = ${ts_per_day} * ${STOP_N}
  breaksw
case "*month*":
  @ mon = ${nmonth}
  @ days_run = 0
  @ n = 1
  while (${n} <= ${STOP_N})
    if (${mon} == 2) then
      @ days_run = ${days_run} + 28
    else if (${mon} == 4 || ${mon} == 6 || ${mon} == 9 || ${mon} == 11) then
      @ days_run = ${days_run} + 30
    else
      @ days_run = ${days_run} + 31
    endif
    @ mon = ${mon} + 1
    if (${mon} > 12) then
      @ mon = 1
    endif
    @ n = ${n} + 1
  end
  @ nitend = ${ts_per_day} * ${days_run}
  breaksw
case "*year*":
  @ nitend = ${ts_per_day} * ${STOP_N} * 365
  breaksw
endsw
@ nstock = ${nitend}
@ nitend = ${nitend} + ${lastkt}

# In an initial run if daily coupling frequency the ocean starts one day later than the other models
if ( "${CONTINUE_RUN}" == "FALSE" && ${nnfsbc} == ${ts_per_day}) then
  if ("${RUN_TYPE}" == "startup" || ${hyb_rest} == 0) then
    @ nitend = ${nitend} - ${nnfsbc}
    @ nstock = ${nstock} - ${nnfsbc}
    set date0 = ${date1}
  endif
endif

########################################################################
# 6. NEMO run-time namelist settings
########################################################################

# Frewsh water budget control
if ("${COMP_ATM}" == "datm") then
  set nnfwb = 1
else if ("${COMP_ATM}" == "cam") then
  set nnfwb = 0
endif

sed -e 's/^[[:blank:]]*cn_exp .*$/cn_exp = "'${CASE}'"/' \
    -e 's/^[[:blank:]]*nn_closea .*$/nn_closea = 1/' \
    -e 's/^[[:blank:]]*ln_zco .*$/ln_zco = .false./' \
    -e 's/^[[:blank:]]*ln_zps .*$/ln_zps = .true./' \
    -e 's/^[[:blank:]]*ln_sco .*$/ln_sco = .false./' \
    -e 's/^[[:blank:]]*nn_fwb .*$/nn_fwb = '"${nnfwb}"'/' \
    -e 's/^[[:blank:]]*nn_ice .*$/nn_ice = 0/' \
    -e 's/^[[:blank:]]*ln_rnf .*$/ln_rnf = .false./' \
    -e 's/^[[:blank:]]*ln_ssr .*$/ln_ssr = .false./' \
    -e 's/^[[:blank:]]*rn_shlat .*$/rn_shlat = 0.0/' \
    -e 's/^[[:blank:]]*nn_evdm .*$/nn_evdm = 1/' \
    -e 's/^[[:blank:]]*rn_avevd .*$/rn_avevd = 10.0/' \
    -e 's/^[[:blank:]]*nn_fsbc .*$/nn_fsbc = '"${nnfsbc}"'/' \
    -e 's/^[[:blank:]]*ln_cpl .*$/ln_cpl = .true./' \
    -e 's/^[[:blank:]]*ln_blk_core .*$/ln_blk_core = .false./' \
    -e 's/^[[:blank:]]*rn_rdt .*$/rn_rdt = '"${DTSEC}"'./' \
    -e 's/^[[:blank:]]*nn_it000 .*$/nn_it000 = '"${nit000}"'/' \
    -e 's/^[[:blank:]]*nn_itend .*$/nn_itend = '"${nitend}"'/' \
    -e 's/^[[:blank:]]*nn_date0 .*$/nn_date0 = '"${date0}"'/' \
    -e 's/^[[:blank:]]*nn_stock .*$/nn_stock = '"${nstock}"'/' \
    -e 's/^[[:blank:]]*nn_write .*$/nn_write = '"${nstock}"'/' \
    -e 's/^[[:blank:]]*nn_fwri .*$/nn_fwri = '"${ts_per_day}"'/' \
    -e 's/^[[:blank:]]*nn_rstctl .*$/nn_rstctl = '"${rstctl}"'/' \
    -e 's/^[[:blank:]]*ln_rstart .*$/ln_rstart = '"${rstart}"'/' \
    -e 's/^[[:blank:]]*nn_msh .*$/nn_msh = '"${msh}"'/' \
    namelist > namelist.tmp && \
    mv -f namelist.tmp namelist

########################################################################
# 7. TOP run-time namelist settings
########################################################################

if ( "${NEMO_TOP_MODULES}" != "" ) then
  sed -e 's/^[[:blank:]]*nn_dttrc .*$/nn_dttrc = 1/' \
      -e 's/^[[:blank:]]*nn_writetrc .*$/nn_writetrc = '"${ts_per_day}"'/' \
      -e 's/^[[:blank:]]*nn_writedia .*$/nn_writedia = '"${ts_per_day}"'/' \
      -e 's/^[[:blank:]]*ln_trcrad .*$/ln_trcrad = .true./' \
      -e 's/^[[:blank:]]*nn_rstctl .*$/nn_rstctl = '"${rstctl}"'/' \
      -e 's/^[[:blank:]]*ln_rsttr .*$/ln_rsttr = '"${rstart}"'/' \
      namelist_top > namelist_top.tmp && \
      mv -f namelist_top.tmp namelist_top

endif

########################################################################

