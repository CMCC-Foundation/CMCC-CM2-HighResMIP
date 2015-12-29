module ocn_comp_mct

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!BOP
! !MODULE: ocn_comp_mct
! !INTERFACE:

! !DESCRIPTION:
!  This is the main driver for the NEMO ocean model
!
! !REVISION HISTORY:
!  SVN:$Id:
!
! !USES:
   use par_kind
   use par_oce
   use dom_oce
   use oce
   use phycst
   use lib_mpp
   use lbclnk
   use sbccpl_cesm
   use in_out_manager
   use geo2ocean
   use step
   use restart
   use nemogcm,    ONLY: cform_aaa, nemo_cesm_init, nemo_closefile
   use lib_fortran
   use iom

   use qflxice

   use mct_mod
   use esmf
   use seq_flds_mod
   use seq_cdata_mod
   use seq_infodata_mod
   use seq_timemgr_mod
   use shr_file_mod 
   use shr_cal_mod 
   use shr_sys_mod
   use perf_mod

   use NEMO_CplIndices

! !PUBLIC MEMBER FUNCTIONS:
  implicit none
  public :: ocn_init_mct
  public :: ocn_run_mct
  public :: ocn_final_mct
  SAVE
  private                              ! By default make data private

!
! ! PUBLIC DATA:
!
! !REVISION HISTORY:
! Author: Mariana Vertenstein
!         Pier Giuseppe Fogli, CMCC, Italy - NEMO version
!
!EOP
! !PRIVATE MODULE FUNCTIONS:
  private :: ocn_export_mct
  private :: ocn_import_mct
  private :: ocn_SetGSMap_mct
  private :: ocn_domain_mct
  private :: ocn_sum_buffer
  private :: lice_form
!
! !PRIVATE MODULE VARIABLES

   logical, parameter ::           &
      lnohalo = .true.     ! logical flag: avoid to exchange the
                           ! extra halo with the driver

   ! The following variables are local copies of the respective NEMO
   ! variables (without the l prefix)
   ! They are identical only when lnohalo=.false.
   integer (i4)  ::     &
      ljpiglo, ljpjglo     ! local copy of the global domain dimensions 

   integer (i4)  ::     &
      lnldi,   lnlei,   &  ! local copy of the interior domain indices
      lnldj,   lnlej

   integer (i4)  ::     &  ! local copy of the global domain location
      lnimpp, lnjmpp       ! of the lower left corner of the current subdomain

   integer (i4)  ::     &
      lnlon, lnlat         ! local domain dimensions 

   !
   real (wp), save ::  &
      area

   real (wp), dimension(:,:,:), allocatable ::  &
      SBUFF_SUM           ! accumulated sum of send buffer quantities
                          ! for averaging before being sent

   integer (i4)  ::   &
      istp,           & ! time step counter (from nit000 to nitend)
      nn_ncpl           ! # of model time steps in 1 coupling time step

   real (wp) ::  &
      tlast_coupled

   integer (i4)  ::   &
      nsend, nrecv

   character(len=lc) :: &
      runtype,          &
      message

   type(seq_infodata_type), pointer :: &
      infodata

   logical ::     &
      ldiag_cpl = .false.

! !PRIVATE MODULE PARAMETERS

   real(wp), parameter ::     &
      c0 = 0.0_wp

!=======================================================================
#  include "vectopt_loop_substitute.h90"

contains

!***********************************************************************
!BOP
!
! !IROUTINE: ocn_init_mct
!
! !INTERFACE:
  subroutine ocn_init_mct( EClock, cdata_o, x2o_o, o2x_o, NLFilename )
!
! !DESCRIPTION:
! Initialize NEMO 
!
! !INPUT/OUTPUT PARAMETERS:

    type(ESMF_Clock)            , intent(in)    :: EClock
    type(seq_cdata)             , intent(inout) :: cdata_o
    type(mct_aVect)             , intent(inout) :: x2o_o, o2x_o
    character(len=*), optional  , intent(in)    :: NLFilename ! Namelist filename
!
! !REVISION HISTORY:
! Author: Mariana Vertenstein
!         Pier Giuseppe Fogli, CMCC, Italy - NEMO version
!EOP
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

    character(len=*), parameter  :: &
         SubName = "ocn_init_mct"

    integer(i4) ::  &
       OCNID,       &
       mpicom_o,    &
       lsize,       &
       start_ymd,   &
       start_tod,   &
       start_year,  &
       start_day,   &
       start_month, &
       start_hour,  &
       stop_ymd,    &
       stop_tod,    &
       cur_ymd,     &
       cur_tod,     &
       nstp,        &
       nemo_cpl_dt, &
       ocn_cpl_dt,  &
       shrlogunit,  &  ! old values
       shrloglev       ! old values

    type(mct_gsMap), pointer :: &
       gsMap_o

    type(mct_gGrid), pointer :: &
       dom_o

    integer (i4) :: &
       errorCode         ! error code

    integer(i4) :: info_debug

    integer (i4) :: iam, ierr, lbnum
    integer (i4) :: n
    character(len=32)  :: starttype          ! infodata start type
    character(len=32)  :: flux_epbal         ! infodata flux_epbal
    INTEGER :: hh,mm,ss

!-----------------------------------------------------------------------
!
!  set cdata pointers
!
!-----------------------------------------------------------------------

    call seq_cdata_setptrs(cdata_o, ID=OCNID, mpicom=mpicom_o, &
         gsMap=gsMap_o, dom=dom_o, infodata=infodata)

#if (defined _MEMTRACE)
    call MPI_comm_rank(mpicom_o,iam,ierr)
    if(iam == 0) then
        lbnum=1
        call memmon_dump_fort('memmon.out','ocn_init_mct:start::',lbnum) 
    endif
#endif

!-----------------------------------------------------------------------
!
!  initialize the model run 
!
!-----------------------------------------------------------------------

    call NEMO_CplIndicesSet()

    call seq_infodata_GetData( infodata, case_name=cn_exp )
   
    call seq_infodata_GetData( infodata, start_type=starttype)

    if (     trim(starttype) == trim(seq_infodata_start_type_start)) then
       runtype = "initial"
    else if (trim(starttype) == trim(seq_infodata_start_type_cont) ) then
       runtype = "continue"
    else if (trim(starttype) == trim(seq_infodata_start_type_brnch)) then
       runtype = "branch"
    else
       call shr_sys_abort('ocn_comp_mct: ocn_init_mct: unknown starttype '//runtype)
    end if


    call t_startf ('nemo_cesm_init')

    ! NEMO initialization
    ! NOTE: we can't use any NEMO variable/routine before the following call
    call nemo_cesm_init(mpicom_o)

    ! check that all process are still there...
    ! If some process have an error, they will never enter in step and other
    ! processes will wait until the end of the cpu time!
    if ( lk_mpp ) call mpp_max( nstop )
    if (nstop /= 0) then
       if (lwp) then   ! error print
          write(numout,cform_err)
          write(numout,*) nstop, ' error have been found' 
       end if
       if (lk_mpp) call mppsync       ! sync PEs
       call shr_sys_abort('ocn_comp_mct: ocn_init_mct: '//SubName//': nstop>0 !')
    endif

    if (lwp) write(numout,cform_aaa)   ! Flag AAAAAAA

    ! Initialize the model time step
    istp = nit000

!-----------------------------------------------------------------------
!
!   allocate space for received fields and initialize to zero
!
!-----------------------------------------------------------------------

    call sbc_cpl_cesm_init

!-----------------------------------------------------------------------
!
!   initialize sea ice formation/melting computation
!
!-----------------------------------------------------------------------

    call init_qflxice

    call t_stopf ('nemo_cesm_init')

    if (lwp) then
      write(numout,FMT='(A,I)') ' ocn_comp_mct:ocn_init_mct: NEMO ID ', OCNID
      write(numout,FMT='(A,I)') ' ocn_comp_mct:ocn_init_mct: NEMO MPICOM ', mpicom_o
    end if

!----------------------------------------------------------------------------
!
! reset shr logging to my log file
!
!----------------------------------------------------------------------------

    call shr_file_getLogUnit (shrlogunit)
    call shr_file_getLogLevel(shrloglev)
    call shr_file_setLogUnit (numout)
   
!-----------------------------------------------------------------------
!
!  check for consistency of nemo and sync clock initial time
!
!-----------------------------------------------------------------------

    if (runtype == 'initial') then

       ! check for consistency of nemo ln_rstart and runtype from seq_infodata
       if (ln_rstart) then
!          if (lk_mpp) call mppsync       ! sync PEs
!          call shr_sys_abort('ocn_comp_mct: ocn_init_mct: '//&
!            'CESM runtype=initial BUT NEMO ln_rstart=.T. !')
          if (lwp) then
             write(numout,FMT='(A,I)') 'ocn_comp_mct: ocn_init_mct: '//&
                 'CESM runtype=initial BUT NEMO ln_rstart=.T. !'
          end if
       endif

       call seq_timemgr_EClockGetData(EClock, &
            start_ymd=start_ymd, start_tod=start_tod, &
            stop_ymd=stop_ymd, stop_tod=stop_tod,     &
            curr_ymd=cur_ymd, curr_tod=cur_tod, StepNo=nstp)
       call shr_cal_date2ymd(start_ymd,start_year,start_month,start_day)
!       if (lwp) write(numout,*) '1 start_ymd=', start_ymd, ' start_tod=', &
!       start_tod, ' stop_ymd=', stop_ymd, ' stop_tod=', stop_tod,        &
!       ' curr_ymd=', cur_ymd, ' curr_tod=', cur_tod, ' StepNo=', nstp
!       hh = nsec_day/3600
!       ss = MOD(nsec_day,3600)
!       mm = ss/60
!       ss = MOD(ss,60)
!       if (lwp) WRITE(numout,FMT='(A,i8,A,i4.4,A,i2.2,A,i2.2,A,i2.2,A,i2.2,Ai2.2)') &
!       'ocn_init_mct: kt=', istp, ' Y/M/D=', nyear,'/',nmonth,'/',nday, ' H/M/S=', &
!       hh, '/', mm, '/', ss
!       if (lwp) write(numout,*) '1 ndastp=', ndastp, ' fjulday=', fjulday
!       if (lwp) call shr_sys_flush(numout)

       if (nyear /= start_year) then
          if (lwp) then
             write(numout,fmt='(A,I)') 'ocn_init_mct: nyear      ', nyear
             write(numout,fmt='(A,I)') 'ocn_init_mct: start_year ', start_year
          endif
          call shr_sys_abort('ocn_comp_mct: ocn_init_mct: nyear does not match start_year')
       end if
       if (nmonth /= start_month) then
          if (lwp) then
             write(numout,fmt='(A,I)') 'ocn_init_mct: nmonth      ', nmonth
             write(numout,fmt='(A,I)') 'ocn_init_mct: start_month ', start_month
          endif
          call shr_sys_abort('ocn_comp_mct: ocn_init_mct: nmonth does not match start_month')
       end if
       if (nday /= start_day) then
          if (lwp) then
             write(numout,fmt='(A,I)') 'ocn_init_mct: nday      ', nday
             write(numout,fmt='(A,I)') 'ocn_init_mct: start_day ', start_day
          endif
       end if
#ifndef _HIRES 
       if (nsec_day /= start_tod) then
          if (lwp) then
             write(numout,fmt='(A,I)') 'ocn_init_mct: nsec      ', nsec_day
             write(numout,fmt='(A,I)') 'ocn_init_mct: start_tod ', start_tod
          end if
       end if
#endif
    else
       ! check for consistency of nemo ln_rstart and runtype from seq_infodata
       if (.not. ln_rstart) then
          if (lk_mpp) call mppsync       ! sync PEs
          call shr_sys_abort('ocn_comp_mct: ocn_init_mct: '//&
            'CESM runtype='//TRIM(runtype)//' BUT NEMO ln_rstart=.F. !')
       endif

    end if

!-----------------------------------------------------------------------
!
!  initialize MCT attribute vectors and indices
!
!-----------------------------------------------------------------------

    call t_startf ('nemo_mct_init')

    call ocn_SetGSMap_mct( mpicom_o, OCNID, gsMap_o )
    lsize = mct_gsMap_lsize(gsMap_o, mpicom_o)

    ! Initialize mct ocn domain (needs ocn initialization info)
   
    call ocn_domain_mct( lsize, gsMap_o, dom_o )
   
    ! Inialize mct attribute vectors
   
    call mct_aVect_init(x2o_o, rList=seq_flds_x2o_fields, lsize=lsize)
    call mct_aVect_zero(x2o_o)
   
    call mct_aVect_init(o2x_o, rList=seq_flds_o2x_fields, lsize=lsize) 
    call mct_aVect_zero(o2x_o)
   
    nsend = mct_avect_nRattr(o2x_o)
    nrecv = mct_avect_nRattr(x2o_o)

!-----------------------------------------------------------------------
!
!   set debugging level
!
!-----------------------------------------------------------------------

    call seq_infodata_GetData( infodata, info_debug=info_debug)
    if (info_debug >= 2) then
       ldiag_cpl = .true. 
    else
       ldiag_cpl = .false.
    endif

!-----------------------------------------------------------------------
!
!   initialize necessary  coupling info
!
!-----------------------------------------------------------------------

    ! TODO: NEMO doesn't have any information about coupling frequency
    ! For now just check that the coupling frequency provided by cpl
    ! is an exact multiple of the ocean time step
    call seq_timemgr_EClockGetData(EClock, dtime=ocn_cpl_dt)
    if (MOD(ocn_cpl_dt, INT(rdt)) /= 0) then
       write(numout,*)' ocn_cpl_dt= ',ocn_cpl_dt, &
                      '        rdt= ',INT(rdt)   
       if (lk_mpp) call mppsync       ! sync PEs
       call shr_sys_abort('ocn_comp_mct: ocn_init_mct: ocn_cpl_dt must be an exact multiple of rdt')
    end if
    ! # of model times step in one coupling time step
    nn_ncpl = ocn_cpl_dt/INT(rdt)
    if (nn_ncpl /= nn_fsbc) then
       write(numout,*)' nn_ncpl= ',nn_ncpl, '  nn_fsbc= ', nn_fsbc
       if (lk_mpp) call mppsync       ! sync PEs
       call shr_sys_abort('ocn_comp_mct: ocn_init_mct: nn_ncpl_dt must be equal to nn_fsbc!')
    end if
    if (ldiag_cpl .AND. lwp) then
      write(numout,*) 'ocn_comp_mct: ocn_init_mct: coupling time step (sec)',  &
      ocn_cpl_dt
      write(numout,*) 'ocn_comp_mct: ocn_init_mct: coupling freq. (steps)  ',  &
      nn_ncpl
    end if

    ! FIXME: used by POP when OCN_COUPLING=partial . Needed here ???
    call seq_infodata_GetData( infodata, flux_epbal=flux_epbal)
    if (lwp) &
      write(numout,FMT='(A)') 'ocn_comp_mct: ocn_init_mct: flux_epbal '//trim(flux_epbal)
    IF (trim(flux_epbal) == 'ocn') THEN
      if (lwp) &
        write(numout,FMT='(A)') 'ocn_comp_mct: ocn_init_mct: send precip_fact to cpl'
      ! From cesm1.2 precip_fact=1.0 ; Previous releases was precip_fact=1.0e6 !!!
      call seq_infodata_PutData( infodata, precip_fact=1.0e6_wp)
    END IF

!-----------------------------------------------------------------------
!
!   allocate space for send buffer
!
!-----------------------------------------------------------------------

    if (.not. allocated(SBUFF_SUM)) allocate (SBUFF_SUM(jpi,jpj,nsend))

!-----------------------------------------------------------------------
!
!  send intial state to driver
!
!-----------------------------------------------------------------------

    tlast_coupled = c0

    call ocn_sum_buffer

    call ocn_export_mct(o2x_o, errorCode)  
    if (errorCode /= 0) then
       if (lk_mpp) call mppsync       ! sync PEs
       call shr_sys_abort('ocn_comp_mct: ocn_init_mct: in ocn_export_mct')
    endif

    call t_stopf ('nemo_mct_init')

    call seq_infodata_PutData( infodata, &
         ocn_nx = ljpiglo , ocn_ny = ljpjglo)
    call seq_infodata_PutData( infodata, &
         ocn_prognostic=.true., ocnrof_prognostic=.true.)

!----------------------------------------------------------------------------
!
! If coupling frequency is higher than daily advance NEMO 1 coupling time step
! Needed to synchronize NEMO, which assumes that the starting time is 00:00.
!
!----------------------------------------------------------------------------

    if (runtype == 'initial' .and. .not. ln_rstart) then
       ! Number of model time step per coupling time step
       ! NOTE: if daily coupling nstp=0, no need to synchronize NEMO
       nstp = MOD(ocn_cpl_dt, INT(rday))/INT(rdt)

       if (nstp /= 0 .and. nstp /= nn_fsbc) then
          write(numout,*) ' Coupling frequency different from SBC/ICE frequency!'
          if (lk_mpp) call mppsync       ! sync PEs
          call shr_sys_abort('ocn_comp_mct: '//SubName//': nstp != nn_fsbc !')
       end if

       if (lwp .and. nstp/=0) then
          write(numout,*) ' Coupling frequency: ', nstp, ' model time steps'
          write(numout,*) ' Advance NEMO 1 coupling time step for coupling synchronization'
       end if

       do n=1, nstp

          lice_form_ts = lice_form(istp)
          lice_cpl_ts  = (MOD(istp-nit000+1, nn_ncpl)==0)
          call increment_tlast_ice(istp)

          call stp(istp)

          ! check for errors
          if (lk_mpp) call mpp_max( nstop )
          if (nstop /= 0) then
             if (lwp) then   ! error print
                write(numout,cform_err)
                write(numout,*) nstop, ' error have been found' 
             end if
             if (lk_mpp) call mppsync       ! sync PEs
             call shr_sys_abort('ocn_comp_mct: '//SubName//': nstop>0 !')
          endif

          call ocn_sum_buffer
          if (MOD(istp-nit000+1, nn_ncpl)==0) then
             call ocn_export_mct(o2x_o, errorCode, lexport=.false.)
          endif

          istp = istp+1

       end do

    end if ! ln_rstart

!----------------------------------------------------------------------------
!
! Reset shr logging to original values
!
!----------------------------------------------------------------------------
    if (lwp) write(numout,FMT='(A,I)') ' END of ocn_init_mct, istp= ', istp
    if (lwp) write(numout,cform_aaa)   ! Flag AAAAAAA

    call shr_sys_flush(numout)
    call shr_file_setLogUnit (shrlogunit)
    call shr_file_setLogLevel(shrloglev)

#if (defined _MEMTRACE)
    if(iam  == 0) then
!        write(6,*) 'ocn_init_mct:end::'
        lbnum=1
        call memmon_dump_fort('memmon.out','ocn_init_mct:end::',lbnum) 
        call memmon_reset_addr()
    endif
#endif

!-----------------------------------------------------------------------
!EOC

 end subroutine ocn_init_mct

!***********************************************************************
!BOP
!
! !IROUTINE: ocn_run_mct
!
! !INTERFACE:
  subroutine ocn_run_mct( EClock, cdata_o, x2o_o, o2x_o)
!
! !DESCRIPTION:
! Run NEMO for a coupling interval
!
! !INPUT/OUTPUT PARAMETERS:
    type(ESMF_Clock)            , intent(in)    :: EClock
    type(seq_cdata)             , intent(inout) :: cdata_o
    type(mct_aVect)             , intent(inout) :: x2o_o
    type(mct_aVect)             , intent(inout) :: o2x_o

!
! !REVISION HISTORY:
! Author: Mariana Vertenstein
!         Pier Giuseppe Fogli, CMCC, Italy - NEMO version
!EOP
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

    integer(i4) :: & 
         errorCode           ! error flag

    integer(i4) :: & 
       shrlogunit,  &  ! old values
       shrloglev       ! old values

    character(len=*), parameter  :: &
         SubName = "ocn_run_mct"

    logical :: &
         rstwr          ! true => write restart at end of day

    integer :: lbnum

    integer :: nstp, reftod, refymd, ymd_sync, tod_sync, nitrst_old

    logical :: lsnd, lrcv

!-----------------------------------------------------------------------
nitrst_old = 0
#if (defined _MEMTRACE)
    if (nproc==0) then
       lbnum=1
       call memmon_dump_fort('memmon.out',SubName//':start::',lbnum) 
    endif
#endif

!-----------------------------------------------------------------------
!
! reset shr logging to my log file
!
!----------------------------------------------------------------------------

    errorCode = 0

    call shr_file_getLogUnit (shrlogunit)
    call shr_file_getLogLevel(shrloglev)
    call shr_sys_flush(shrlogunit)
    call shr_file_setLogUnit (numout)
    call shr_sys_flush(numout)

    call seq_cdata_setptrs(cdata_o, infodata=infodata)

    if (lwp .AND. ldiag_cpl) write(numout,FMT='(A)') ' BEGINNING of ocn_init_run'
!----------------------------------------------------------------------------
!
! restart flag (rstwr) will assume only an eod restart for now
!
!----------------------------------------------------------------------------

    rstwr = seq_timemgr_RestartAlarmIsOn(EClock)
!    if (rstwr) then
!       call override_time_flag(cpl_write_restart,value=.true.)
!       write(message,'(6a)') 'driver requests restart file at eod  ',  &
!            cyear,'/',cmonth,'/',cday
!       call document ('ocn_comp_mct(run):', message)
!      if (lwp) write(numout,*) 'driver requests restart file at kt= ', istp
!    endif

   call seq_timemgr_EClockGetData( EClock, curr_ymd=ymd_sync, &
          curr_tod=tod_sync, StepNo=nstp, ref_ymd=refymd,     &
          ref_tod=reftod)
!   if (lwp) write(numout,*) '2 cur_ymd=', ymd_sync, ' cur_tod=', tod_sync, &
!      ' StepNo=', nstp
!   if (lwp) write(numout,*) '2 ndastp=', ndastp, ' fjulday=', fjulday
!   if (lwp) call shr_sys_flush(numout)

!-----------------------------------------------------------------------
!
!  advance the model in time over coupling interval
!
!-----------------------------------------------------------------------

    advance: do 

       lsnd = .false.
       lrcv = .false.

       ! obtain import state from driver
       if (MOD(istp-nit000, nn_ncpl)==0) then
          call ocn_import_mct(x2o_o, errorCode)   

          if (errorCode /= 0) then
             if (lk_mpp) call mppsync       ! sync PEs
             call shr_sys_abort('ocn_comp_mct: '//SubName//': in ocn_import_mct')
          endif
          lrcv = .true.
       end if
       
       ! open restart file if requested from driver
       if (rstwr .and. nn_ncpl>1 .and. (MOD(istp-nit000+2, nn_ncpl)==0)) then
          if (lwp) write(numout,*) 'driver requests restart file at kt= ', istp+1
          if (istp+1 /= nitrst) THEN
             if (lwp) write(numout,*) 'nitrst=', nitrst, ' set nitrst to ', istp+1
             nitrst_old = nitrst  ! Save original value
             nitrst = istp+1
          endif
          if (lwp) write(numout,*) &
            'open restart file at kt= ', istp, ' ndastp=', ndastp
          if (lwp) call shr_sys_flush(numout)
       end if

       ! Compute sea ice formation
       lice_form_ts = lice_form(istp)
!       lice_form_ts = .false.
       if (lice_form_ts) then
          if (lwp .AND. istp-nit000+1<48) then
            write(numout,*)  &
            SubName, ': time step (istp) ', istp, ' lice_form_ts=', lice_form_ts
          endif
       end if

       ! Compute sea ice formation/melting flux at coupling time step
       lice_cpl_ts = (MOD(istp-nit000+1, nn_ncpl)==0)
!       lice_cpl_ts = .false.
       if (lice_cpl_ts) then
          if (lwp .AND. istp-nit000+1<48) then
            write(numout,*)  &
            SubName, ': time step (istp) ', istp, ' lice_cpl_ts=', lice_cpl_ts
          endif
       end if
       call shr_sys_flush(numout)

       call increment_tlast_ice(istp)

       ! advance NEMO 1 model time step
       call stp(istp)

       if (lrcv .and. ldiag_cpl .and. lwp)     &
          write(numout,*)  &
          SubName, ': ocn_import_mct called at time step (istp) ', istp, &
          ' ndastp=', ndastp

       if (ldiag_cpl .and. lwp)     &
          write(numout,*)  &
          SubName, ': stp() called at time step (istp) ', istp, ' ndastp=', ndastp

       ! check for errors
       if (lk_mpp) call mpp_max( nstop )
       if (nstop /= 0) then
          if (lwp) then   ! error print
             write(numout,cform_err)
             write(numout,*) nstop, ' error have been found' 
          end if
          if (lk_mpp) call mppsync       ! sync PEs
          call shr_sys_abort('ocn_comp_mct: '//SubName//': nstop>0 !')
       endif

       ! return export state to driver
       call ocn_sum_buffer
       if (MOD(istp-nit000+1, nn_ncpl)==0) then
          call ocn_export_mct(o2x_o, errorCode)
          if (errorCode /= 0) then
             if (lk_mpp) call mppsync       ! sync PEs
             call shr_sys_abort('ocn_comp_mct: '//SubName//': in ocn_export_mct')
          endif
          if (ldiag_cpl .and. lwp)     &
             write(numout,*)  &
             SubName, ': ocn_export_mct called at time step (istp) ', istp, &
             ' ndastp=', ndastp
          istp = istp+1
          exit advance
       end if

       istp = istp+1
       
    enddo advance

!--------------------------------------------------------------------
!
! check that internal clock is in sync with master clock
!
!--------------------------------------------------------------------

    if (rstwr) then
       if (lwp) write(numout,*) 'nitrst=', nitrst, ' restore nitrst to ', nitrst_old
       nitrst = nitrst_old  ! Restore original value
    end if

!   if (lwp) write(numout,*) '3 ndastp=', ndastp, ' fjulday=', fjulday
!   if (lwp) write(numout,*) '3 cur_ymd=', ymd_sync, ' cur_tod=', tod_sync, &
!      ' StepNo=', nstp
!   if (lwp) call shr_sys_flush(numout)


!    ymd = iyear*10000 + imonth*100 + iday
!    tod = ihour*seconds_in_hour + iminute*seconds_in_minute + isecond
!    if ( .not. seq_timemgr_EClockDateInSync( EClock, ymd, tod ) )then
!       call seq_timemgr_EClockGetData( EClock, curr_ymd=ymd_sync, &
!          curr_tod=tod_sync )
!       write(stdout,*)' pop2 ymd=',ymd     ,'  pop2 tod= ',tod
!       write(stdout,*)' sync ymd=',ymd_sync,'  sync tod= ',tod_sync
!       write(stdout,*)' Internal pop2 clock not in sync with Sync Clock'
!       call shr_sys_abort( SubName// &
!          ":: Internal pop2 clock not in sync with Sync Clock")
!    end if
   
   if (lwp .AND. ldiag_cpl) write(numout,FMT='(A)') ' END of ocn_init_run'
   call shr_sys_flush(numout)
!----------------------------------------------------------------------------
!
! Reset shr logging to original values
!
!----------------------------------------------------------------------------

   call shr_file_setLogUnit (shrlogunit)
   call shr_file_setLogLevel(shrloglev)
   call shr_sys_flush(shrlogunit)

#if (defined _MEMTRACE)
    if(nproc == 0) then
       lbnum=1
       call memmon_dump_fort('memmon.out',SubName//':end::',lbnum) 
       call memmon_reset_addr()
    endif
#endif
!-----------------------------------------------------------------------
!EOC

  end subroutine ocn_run_mct

!***********************************************************************
!BOP
!
! !IROUTINE: ocn_final_mct
!
! !INTERFACE:
  subroutine ocn_final_mct( EClock, cdata_o, x2o_o, o2x_o)
!
! !DESCRIPTION:
! Finalize NEMO
!
! !USES:
!
! !ARGUMENTS:
    type(ESMF_Clock)            , intent(in)    :: EClock
    type(seq_cdata)             , intent(inout) :: cdata_o
    type(mct_aVect)             , intent(inout) :: x2o_o
    type(mct_aVect)             , intent(inout) :: o2x_o
!
!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

!   if ( lk_diaobs ) call dia_obs_wri

   !                            !------------------------!
   !                            !==  finalize the run  ==!
   !                            !------------------------!
   if (lwp) write(numout,cform_aaa)   ! Flag AAAAAAA
   !
   if ( nstop /= 0 .and. lwp ) then   ! error print
      write(numout,cform_err)
      write(numout,*) nstop, 'error have been found'
   endif
   !
   call nemo_closefile

   ! Free allocated space
   call sbc_cpl_cesm_finalize
   if (allocated(SBUFF_SUM)) deallocate (SBUFF_SUM)

!-----------------------------------------------------------------------
  end subroutine ocn_final_mct

!***********************************************************************
!BOP
!IROUTINE: ocn_SetGSMap_mct
! !INTERFACE:

 subroutine ocn_SetGSMap_mct( mpicom_ocn, OCNID, gsMap_ocn )

! !DESCRIPTION:
!  This routine mct global seg maps for the NEMO decomposition
!
! !REVISION HISTORY:
!  same as module

! !INPUT/OUTPUT PARAMETERS:

    implicit none
    integer        , intent(in)    :: mpicom_ocn
    integer        , intent(in)    :: OCNID
    type(mct_gsMap), intent(inout) :: gsMap_ocn

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

    integer,allocatable :: &
      gindex(:)

    integer (i4) ::     &
      lnictls, lnictle, &
      lnjctls, lnjctle

    integer (i4) ::   &
      i, j, n,        &
      lsize, gsize,   &
      ier
    logical      ::   &
      mpp_ocn
 

!-----------------------------------------------------------------------
!  Build the NEMO grid numbering for MCT
!  NOTE:  Numbering scheme is: West to East and South to North starting
!  at the south pole.  Should be the same as what's used in SCRIP
!
!  *** NOTE *** : To be consistent with domain and mapping files the
!  overlapping grid points belonging to the East-West cyclic boundary
!  and to the North folding are removed (lnohalo=.true.).
!  This is equivalent to remove from the *global* domain the first and
!  the last colums and the last row.
!
!  2015-10-12 TL: add control to check if nemo is initialized with ocean
!                 domain decomposition only and bypass check on total
!                 number of grid points coherence. 
!-----------------------------------------------------------------------

    ! Currently only global cyclic east-west boundary condition is implemented
    if (jperio /= 4 .and. jperio /= 6) &
       call shr_sys_abort('ocn_comp_mct: ocn_SetGSMap_mct: unexpected lateral domain '// &
       'boundary condition (jperio /= [4, 6])!')

    ! Check if ocean domain decomposition only or not
    mpp_ocn = .false.
    if ( jpnij < jpni * jpnj ) then
       mpp_ocn = .true.
    endif

    ! Global domain size
    ljpiglo = jpiglo
    ljpjglo = jpjglo
    if (lnohalo) then                ! remove the extra halo
      ljpiglo = ljpiglo-2*jpreci     ! East-West halo
      ljpjglo = ljpjglo-jprecj       ! North halo
    end if
    gsize = ljpiglo*ljpjglo   ! Number of grid points for the global domain

    ! Local domain interior indeces
    lnldi = nldi
    lnldj = nldj
    lnlei = nlei
    lnlej = nlej

    ! Global domain location of the interior indeces
    lnictls = nimppt(narea) + nldit(narea) - 1
    lnictle = nimppt(narea) + nleit(narea) - 1
    lnjctls = njmppt(narea) + nldjt(narea) - 1
    lnjctle = njmppt(narea) + nlejt(narea) - 1

    if (lnohalo) then      ! remove the extra halo
      if (lnictls==1)      lnldi=lnldi+1  ! West  halo
      if (lnictle==jpiglo) lnlei=lnlei-1  ! East  halo
      if (lnjctle==jpjglo) lnlej=lnlej-1  ! North halo
    end if
!    write(numout,*) nproc,':  ', lnldi,':',lnlei, ' - ', lnldj,':',lnlej

    ! Global domain location of the lower left corner for
    ! the current subdomain
    lnimpp = nimpp
    lnjmpp = njmpp
    if (lnohalo) lnimpp = lnimpp-1 ! remove the extra halo
!    write(numout,*) nproc,':  ',nimpp,' ',njmpp,' ',lnimpp,' ',lnjmpp
!    call shr_sys_flush(numout)

    ! Local domain size
    lnlon = lnlei-lnldi+1
    lnlat = lnlej-lnldj+1
    lsize = lnlon*lnlat     ! Number of grid points for the local domain

    if ((lsize<=0).or.(lsize>gsize)) then
       write(message,FMT='(A,I)')     &
          'ocn_comp_mct: ocn_SetGSMap_mct: wrong local size: lsize=', lsize
       call shr_sys_abort(message)
    end if

    n = lsize
    call mpp_sum(n)
    if (gsize /= n .and. .NOT. mpp_ocn ) then
       write(numout,FMT='(A)') 'ocn_comp_mct: ocn_SetGSMap_mct: number of points in the global domain'
       write(numout,FMT='(2(A,I))') ' gsize=', gsize, ' mpp_sum(lsize)=', n
       if (lk_mpp) call mppsync       ! sync PEs
       call shr_sys_abort('ocn_comp_mct: ocn_SetGSMap_mct: mpp_sum(lsize) /= gsize !')
    elseif (gsize /= n .and. mpp_ocn .and. lwp ) then
       write(numout,FMT='(A)') 'ocn_comp_mct: ocn_SetGSMap_mct: NEMO is using ocean only domains (mpp_init2)'
       write(numout,FMT='(2(A,I))') ' gsize=', gsize, ' mpp_sum(lsize)=', n
       write(numout,FMT='(2(A,I))') ' jpnij=', jpnij, ' < jpni * jpnj=', jpni * jpnj
    end if
    if (lwp) &
      write(*,*) 'ocn_comp_mct: ocn_SetGSMap_mct: mpp_sum(lsize)=', n

    allocate(gindex(lsize),stat=ier)
    if (ier/=0) then
       if (lk_mpp) call mppsync       ! sync PEs
       call shr_sys_abort('ocn_comp_mct: ocn_SetGSMap_mct: failed allocation (gindex)!')
    end if
    gindex(:) = 0

    n = 0
    do j=lnldj,lnlej
       do i=lnldi,lnlei
          n=n+1
          gindex(n) = (lnjmpp+j-2)*ljpiglo + (lnimpp+i-1)
       enddo
    enddo

    if (any(gindex(:)<0)) then
       if (lk_mpp) call mppsync       ! sync PEs
       call shr_sys_abort('ocn_comp_mct: ocn_SetGSMap_mct: in ocn_SetGSMap_mct gindex<0!')
    end if

    if (any(gindex(:)>gsize)) then
       write(*,*) gindex(:)
       if (lk_mpp) call mppsync       ! sync PEs
       call shr_sys_abort('ocn_comp_mct: ocn_SetGSMap_mct: in ocn_SetGSMap_mct gindex>gsize!')
    end if

    call mct_gsMap_init( gsMap_ocn, gindex, mpicom_ocn, OCNID, lsize, gsize )

    deallocate(gindex)

!-----------------------------------------------------------------------
!EOC

  end subroutine ocn_SetGSMap_mct

!***********************************************************************
!BOP
! !IROUTINE: ocn_domain_mct
! !INTERFACE:

 subroutine ocn_domain_mct( lsize, gsMap_o, dom_o )

! !DESCRIPTION:
!  This routine mct global seg maps for the NEMO decomposition
!
! !REVISION HISTORY:
!  same as module
!
! !INPUT/OUTPUT PARAMETERS:

    implicit none
    integer        , intent(in)    :: lsize
    type(mct_gsMap), intent(in)    :: gsMap_o
    type(mct_ggrid), intent(inout) :: dom_o     

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

    integer, pointer :: &
      idata(:)

    real(wp), pointer :: &
      data(:)

    integer (i4) ::   &
      i,j, n, &
      ier

!-------------------------------------------------------------------
!
!  initialize mct domain type, lat/lon in degrees,
!  area in radians^2, mask is 1 (ocean), 0 (non-ocean)
!
!-------------------------------------------------------------------

    call mct_gGrid_init( GGrid=dom_o, CoordChars=trim(seq_flds_dom_coord), &
       OtherChars=trim(seq_flds_dom_other), lsize=lsize )
    call mct_aVect_zero(dom_o%data)
    allocate(data(lsize))

!-------------------------------------------------------------------
!
! Determine global gridpoint number attribute, GlobGridNum, which is set automatically by MCT
!
!-------------------------------------------------------------------

    ! NOTE: idata is allocated inside the next call
    call mct_gsMap_orderedPoints(gsMap_o, nproc, idata)
    call mct_gGrid_importIAttr(dom_o,'GlobGridNum',idata,lsize)

!-------------------------------------------------------------------
!
! Determine domain (numbering scheme is: West to East and South to North to South pole)
! Initialize attribute vector with special value
!
!-------------------------------------------------------------------

    data(:) = -9999.0_wp 
    call mct_gGrid_importRAttr(dom_o,"lat"  ,data,lsize) 
    call mct_gGrid_importRAttr(dom_o,"lon"  ,data,lsize) 
    call mct_gGrid_importRAttr(dom_o,"area" ,data,lsize) 
    call mct_gGrid_importRAttr(dom_o,"aream",data,lsize) 
    data(:) = c0
    call mct_gGrid_importRAttr(dom_o,"mask",data,lsize) 
    call mct_gGrid_importRAttr(dom_o,"frac",data,lsize) 

!-------------------------------------------------------------------
!
! Fill in correct values for domain components
!
!-------------------------------------------------------------------

    n=0
    do j=lnldj,lnlej
       do i=lnldi,lnlei
          n=n+1
          data(n) = glamt(i,j)
       enddo
    enddo
    ! longitude shift: [-180:180] --> [0:360]
    ! consistent with domain and mapping files
    where (data(:)<c0)
       data(:) = data(:) + 360.0_wp
    end where
    call mct_gGrid_importRattr(dom_o,"lon",data,lsize) 

    n=0
    do j=lnldj,lnlej
       do i=lnldi,lnlei
          n=n+1
          data(n) = gphit(i,j)
       enddo
    enddo
    call mct_gGrid_importRattr(dom_o,"lat",data,lsize) 

    n=0
    do j=lnldj,lnlej
       do i=lnldi,lnlei
          n=n+1
          data(n) = e1e2t(i,j)/(ra*ra)
       enddo
    enddo
    call mct_gGrid_importRattr(dom_o,"area",data,lsize) 

    ! NOTE: grid folding (halo) is masked out (tmask_i)
    ! consistent with domain and mapping files
    n=0
    do j=lnldj,lnlej
       do i=lnldi,lnlei
          n=n+1
          data(n) = tmask_i(i,j)
          if (data(n) > 1.0_wp) data(n) = 1.0_wp
       enddo
    enddo
    call mct_gGrid_importRattr(dom_o,"mask",data,lsize) 
    call mct_gGrid_importRattr(dom_o,"frac",data,lsize) 

    deallocate(data)
    deallocate(idata)

    area = glob_sum(e1e2t(:,:))
    if (nproc==0) write(numout,'(a,1es28.19)') 'Global ocean area (m**2)', area

!-----------------------------------------------------------------------
!EOC

  end subroutine ocn_domain_mct

!***********************************************************************
!BOP
! !IROUTINE: ocn_import_mct
! !INTERFACE:

 subroutine ocn_import_mct(x2o_o, errorCode)

! !DESCRIPTION:
!-----------------------------------------------------------------------
!  This routine receives message from cpl7 driver
!
!    The following fields are always received from the coupler:
! 
!    o  taux   -- zonal wind stress (taux)                 (W/m2   )
!    o  tauy   -- meridonal wind stress (tauy)             (W/m2   )
!    o  snow   -- water flux due to snow                   (kg/m2/s)
!    o  rain   -- water flux due to rain                   (kg/m2/s)
!    o  evap   -- evaporation flux                         (kg/m2/s)
!    o  meltw  -- snow melt flux                           (kg/m2/s)
!    o  salt   -- salt                                     (kg(salt)/m2/s)
!    o  swnet  -- net short-wave heat flux                 (W/m2   )
!    o  sen    -- sensible heat flux                       (W/m2   )
!    o  lwup   -- longwave radiation (up)                  (W/m2   )
!    o  lwdn   -- longwave radiation (down)                (W/m2   )
!    o  melth  -- heat flux from snow&ice melt             (W/m2   )
!    o  ifrac  -- ice fraction
!    o  roff   -- river runoff flux                        (kg/m2/s)
!    o  ioff   -- ice runoff flux                          (kg/m2/s)
! 
!    The following fields are sometimes received from the coupler,
!      depending on model options:
! 
!    o  pslv   -- sea-level pressure                       (Pa)
!    o  duu10n -- 10m wind speed squared                   (m^2/s^2)
!    o  co2prog-- bottom atm level prognostic co2
!    o  co2diag-- bottom atm level diagnostic co2
! 
!-----------------------------------------------------------------------
!
! !REVISION HISTORY:
!  same as module

! !INPUT/OUTPUT PARAMETERS:

    type(mct_aVect)   , intent(inout) :: x2o_o

! !OUTPUT PARAMETERS:

   integer (i4), intent(out) :: &
      errorCode              ! returned error code

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   character (len=lc) ::   &
      label
 
   integer (i4) ::  &
      i,j,k,n

   real (wp), dimension(jpi,jpj) ::   &
      WORK1               ! local work space

   real (wp) ::  &
      sgn,       &
      gsum
 
   real (wp) ::  &
      rtmp1, rtmp2

!-----------------------------------------------------------------------
!
!  zero out padded cells 
!
!-----------------------------------------------------------------------

   errorCode = 0

   ! lrecv should be reset to false in sbc_cpl_cesm_rcv
   if (lrecv) then
     errorCode = 1
     return
   end if

!-----------------------------------------------------------------------
!
!  unpack fields received from cpl
!
!-----------------------------------------------------------------------

   n = 0
   do j=lnldj,lnlej
      do i=lnldi,lnlei
         n = n + 1
         taux_x2o(i,j)   = x2o_o%rAttr(index_x2o_Foxx_taux,n)
         tauy_x2o(i,j)   = x2o_o%rAttr(index_x2o_Foxx_tauy,n)
         !
         ! emp: E-P
         ! NOTE: evap_x2o does not include evap/subl over sea ice (see meltw)
         evap_x2o(i,j)   = x2o_o%rAttr(index_x2o_Foxx_evap,n)
         rain_x2o(i,j)   = x2o_o%rAttr(index_x2o_Faxa_rain,n)
         snow_x2o(i,j)   = x2o_o%rAttr(index_x2o_Faxa_snow,n)
         ! runoff
         roff_x2o(i,j)   = x2o_o%rAttr(index_x2o_Forr_roff,n)
         ioff_x2o(i,j)   = x2o_o%rAttr(index_x2o_Forr_ioff,n)
         !
         ! meltw: Fresh water flux from sea ice
         ! NOTE: includes evap over snow/sea ice (condensation/sublimation)
         !       and snow melting over sea ice
         !       The sea ice freezing/melting is not taken into account in
         !       CICE when coupled to NEMO (see Tartinville et al., 2001)
         meltw_x2o(i,j)  = x2o_o%rAttr(index_x2o_Fioi_meltw,n)
         !
         ! salt: salt flux from sea ice
         ! Converted later to an equivalent fresh water flux
         ! and added to emp to form emps, as in NEMO-LIM
         salt_x2o(i,j)   = x2o_o%rAttr(index_x2o_Fioi_salt,n)
         !
         swnet_x2o(i,j)  = x2o_o%rAttr(index_x2o_Foxx_swnet,n)
         sen_x2o(i,j)    = x2o_o%rAttr(index_x2o_Foxx_sen,n)
         lat_x2o(i,j)    = x2o_o%rAttr(index_x2o_Foxx_lat,n)
         lwup_x2o(i,j)   = x2o_o%rAttr(index_x2o_Foxx_lwup,n)
         lwdn_x2o(i,j)   = x2o_o%rAttr(index_x2o_Faxa_lwdn,n)
         melth_x2o(i,j)  = x2o_o%rAttr(index_x2o_Fioi_melth,n)
         !
         ifrac_x2o(i,j)  = x2o_o%rAttr(index_x2o_Si_ifrac,n)
         pslv_x2o(i,j)   = x2o_o%rAttr(index_x2o_Sa_pslv,n)
         duu10n_x2o(i,j) = x2o_o%rAttr(index_x2o_So_duu10n,n)
         !
      enddo
   enddo

!-----------------------------------------------------------------------
!
!  incoming data quality control
!
!-----------------------------------------------------------------------
#ifdef CCSMCOUPLED
      if ( any(ioff_x2o < c0) ) then
        if (lk_mpp) call mppsync       ! sync PEs
        call shr_sys_abort ('ocn_comp_mct: ocn_import_mct: incoming ioff_x2o is negative')
      endif
#endif

!-----------------------------------------------------------------------
!
!  unpack atmospheric CO2
!
!-----------------------------------------------------------------------
   if (index_x2o_Sa_co2prog > 0) then
      n = 0
      do j=lnldj,lnlej
         do i=lnldi,lnlei
            n = n + 1
            co2_x2o(i,j) = x2o_o%rAttr(index_x2o_Sa_co2prog,n)
         enddo
      enddo
   endif

   if (index_x2o_Sa_co2diag > 0) then
      n = 0
      do j=lnldj,lnlej
         do i=lnldi,lnlei
            n = n + 1
            co2_x2o(i,j) = x2o_o%rAttr(index_x2o_Sa_co2diag,n)
         enddo
      enddo
   endif

   ! coupling time step flag
   lrecv = .TRUE.
 
!-----------------------------------------------------------------------
!
!  diagnostics
!
!-----------------------------------------------------------------------

   if (ldiag_cpl) then
     if (lwp) write(numout,*) 'nemo_recv_from_coupler'
 
     do k = 1,nrecv

         WORK1 = c0
         n = 0
         do j=lnldj,lnlej
            do i=lnldi,lnlei
               n = n + 1
               WORK1(i,j) = x2o_o%rAttr(k,n)
            enddo
         enddo

         sgn = 1._wp
         if (k==index_x2o_Foxx_taux .or. k==index_x2o_Foxx_tauy) then
            sgn = -1._wp
         end if
         call lbc_lnk(WORK1, 'T', sgn)

         gsum = glob_sum(WORK1(:,:)*e1e2t(:,:))

         rtmp1=minval(WORK1, mask=(tmask(:,:,1)>0._wp))
         call mpp_min(rtmp1)
         rtmp2=maxval(WORK1, mask=(tmask(:,:,1)>0._wp))
         call mpp_max(rtmp2)

         if (nproc==0) then
            call seq_flds_getField(label,k,seq_flds_x2o_fields)
            write(numout,1100)'ocn','recv', TRIM(label), gsum, gsum/area, &
              rtmp1, rtmp2
         endif
      enddo
      call shr_sys_flush(numout)
   endif


1100  format ('comm_diag ', a3, 1x, a4, 1x, a16, 1x, 4es28.19:)

!-----------------------------------------------------------------------
!EOC

 end subroutine ocn_import_mct
!***********************************************************************
!BOP
! !IROUTINE: ocn_export_mct
! !INTERFACE:

 subroutine ocn_export_mct(o2x_o, errorCode, lexport)

! !DESCRIPTION:
!  This routine calls the routines necessary to send NEMO fields to
!  the CCSM cpl7 driver
!
! !REVISION HISTORY:
!  same as module
!
! !INPUT/OUTPUT PARAMETERS:

   type(mct_aVect)   , intent(inout) :: o2x_o

   logical, intent(in), optional :: lexport

! !OUTPUT PARAMETERS:

   integer (i4), intent(out) :: &
      errorCode              ! returned error code

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (i4) :: n
           
   character (len=lc)    :: label
 
   integer (i4) ::  &
      i,j,k

   real (wp), dimension(jpi,jpj) ::   &
      WORK1, WORK2,      &! local work space
      WORK3, WORK4

   real (wp), dimension(jpi,jpj) ::   &
        WORKA               ! local work space with full block dimension

   real (wp) ::   &
      sgn,        &
      gsum

   real (wp) ::   &
      rtmp1, rtmp2

   logical :: l_export

   logical, save :: lfirst = .true.

!-----------------------------------------------------------------------

   ! NOTE:
   ! We do not use iom_setkt before iom_put as in sbcmod because we're
   ! already at the right time step for SBC output (istp+nn_fsbc-1).
   ! We just avoid to call iom_put the first time this routine is called
   ! (from ocn_init_mct) for initialization pourpuses.

   l_export = .true.
   if (present(lexport)) l_export = lexport

!-----------------------------------------------------------------------
!
!  initialize control buffer
!
!-----------------------------------------------------------------------

   errorCode = 0

!-----------------------------------------------------------------------
!
!     interpolate onto T-grid points and rotate on T grid
!
!-----------------------------------------------------------------------

   WORK1(:,:) = c0
   WORK2(:,:) = c0
   WORK3(:,:) = c0
   WORK4(:,:) = c0

   ! Apply LBC
   SBUFF_SUM(:,:,index_o2x_So_u) = SBUFF_SUM(:,:,index_o2x_So_u)*umask(:,:,1)
   SBUFF_SUM(:,:,index_o2x_So_v) = SBUFF_SUM(:,:,index_o2x_So_v)*vmask(:,:,1)
   call lbc_lnk(SBUFF_SUM(:,:,index_o2x_So_u), 'U', -1._wp)
   call lbc_lnk(SBUFF_SUM(:,:,index_o2x_So_v), 'V', -1._wp)
   !
   ! (U,V) -> T
   do j = 2, jpj
      do i = 2, jpi
!          WORK3(i,j) = 0.5_wp * ( SBUFF_SUM(i,j,index_o2x_So_u) + &
!             SBUFF_SUM(i-1,j,index_o2x_So_u) )
!          WORK4(i,j) = 0.5_wp * ( SBUFF_SUM(i,j,index_o2x_So_v) + &
!             SBUFF_SUM(i,j-1,index_o2x_So_v) )
         rtmp1 = umask(i,j,1)+umask(i-1,j,1)
         if (rtmp1 > c0) then
            WORK3(i,j) = ( SBUFF_SUM(i,j,index_o2x_So_u)*umask(i,j,1) + &
               SBUFF_SUM(i-1,j,index_o2x_So_u)*umask(i-1,j,1) ) / rtmp1
         end if
         rtmp2 = vmask(i,j,1)+vmask(i,j-1,1)
         if (rtmp2 > c0) then
            WORK4(i,j) = ( SBUFF_SUM(i,j,index_o2x_So_v)*vmask(i,j,1) + &
               SBUFF_SUM(i,j-1,index_o2x_So_v)*vmask(i,j-1,1) ) / rtmp2
         end if
      end do
   end do
   WORK3(:,:) = WORK3(:,:)*tmask(:,:,1)
   WORK4(:,:) = WORK4(:,:)*tmask(:,:,1)
   ! Apply LBC
   call lbc_lnk( WORK3, 'T', -1._wp)
   call lbc_lnk( WORK4, 'T', -1._wp)

   if (.not. lfirst) then
      call iom_put('So_u_o2x', WORK3/tlast_coupled)
      call iom_put('So_v_o2x', WORK4/tlast_coupled)
   endif

   call rot_rep( WORK3, WORK4, 'T', 'ij->e', WORK1 )
   call rot_rep( WORK3, WORK4, 'T', 'ij->n', WORK2 )

   WORK3(:,:) = WORK1(:,:)*tmask(:,:,1)/tlast_coupled
   WORK4(:,:) = WORK2(:,:)*tmask(:,:,1)/tlast_coupled

   if (l_export) then
   n = 0
   do j=lnldj,lnlej
      do i=lnldi,lnlei
         n = n + 1
         o2x_o%rAttr(index_o2x_So_u,n) = WORK3(i,j)
         o2x_o%rAttr(index_o2x_So_v,n) = WORK4(i,j)
      enddo
   enddo
   end if

!-----------------------------------------------------------------------
!
!     convert and pack surface temperature
!
!-----------------------------------------------------------------------

   WORK1(:,:) = (SBUFF_SUM(:,:,index_o2x_So_t)/tlast_coupled + rt0)*tmask(:,:,1)
   call lbc_lnk( WORK1, 'T', 1._wp)

   if (l_export) then
   n = 0
   do j=lnldj,lnlej
      do i=lnldi,lnlei
         n = n + 1
         o2x_o%rAttr(index_o2x_So_t,n) = WORK1(i,j)
      enddo
   enddo
   end if

   if (.not. lfirst) then
      call iom_put('So_t_o2x', WORK1)
   endif

!-----------------------------------------------------------------------
!
!     convert and pack salinity
!
!-----------------------------------------------------------------------

   WORK1(:,:) = SBUFF_SUM(:,:,index_o2x_So_s)*tmask(:,:,1)/tlast_coupled
   call lbc_lnk( WORK1, 'T', 1._wp)

   if (l_export) then
   n = 0
   do j=lnldj,lnlej
      do i=lnldi,lnlei
         n = n + 1
         o2x_o%rAttr(index_o2x_So_s,n) = WORK1(i,j)
      enddo
   enddo
   end if

   if (.not. lfirst) then
      call iom_put('So_s_o2x', WORK1)
   endif

!-----------------------------------------------------------------------
!
!     interpolate onto T-grid points, then rotate on T grid
!
!-----------------------------------------------------------------------

   WORK1(:,:) = c0
   WORK2(:,:) = c0
   WORK3(:,:) = c0
   WORK4(:,:) = c0

   ! Apply LBC (not sure it's necessary)
   SBUFF_SUM(:,:,index_o2x_So_dhdx) = SBUFF_SUM(:,:,index_o2x_So_dhdx)*umask(:,:,1)
   SBUFF_SUM(:,:,index_o2x_So_dhdy) = SBUFF_SUM(:,:,index_o2x_So_dhdy)*vmask(:,:,1)
   call lbc_lnk(SBUFF_SUM(:,:,index_o2x_So_dhdx), 'U', -1._wp)
   call lbc_lnk(SBUFF_SUM(:,:,index_o2x_So_dhdy), 'V', -1._wp)
   !
   ! (U,V) -> T
   do j = 2, jpj
      do i = 2, jpi
!          WORK3(i,j) = 0.5_wp * ( SBUFF_SUM(i,j,index_o2x_So_dhdx) + &
!             SBUFF_SUM(i-1,j,index_o2x_So_dhdx) )
!          WORK4(i,j) = 0.5_wp * ( SBUFF_SUM(i,j,index_o2x_So_dhdy) + &
!             SBUFF_SUM(i,j-1,index_o2x_So_dhdy) )
         rtmp1 = umask(i,j,1)+umask(i-1,j,1)
         if (rtmp1>c0) then
            WORK3(i,j) = ( SBUFF_SUM(i,j,index_o2x_So_dhdx)*umask(i,j,1) + &
               SBUFF_SUM(i-1,j,index_o2x_So_dhdx)*umask(i-1,j,1) ) / rtmp1
         end if
         rtmp2 = vmask(i,j,1)+vmask(i,j-1,1)
         if (rtmp2>c0) then
            WORK4(i,j) = ( SBUFF_SUM(i,j,index_o2x_So_dhdy)*vmask(i,j,1) + &
               SBUFF_SUM(i,j-1,index_o2x_So_dhdy)*vmask(i,j-1,1) ) / rtmp2
         end if
      end do
   end do
   WORK3(:,:) = WORK3(:,:)*tmask(:,:,1)
   WORK4(:,:) = WORK4(:,:)*tmask(:,:,1)
   ! Apply LBC
   call lbc_lnk( WORK3, 'T', -1._wp)
   call lbc_lnk( WORK4, 'T', -1._wp)

   if (.not. lfirst) then
      call iom_put("So_dhdx_o2x", WORK3/tlast_coupled)
      call iom_put("So_dhdy_o2x", WORK4/tlast_coupled)
   endif

   call rot_rep( WORK3, WORK4, 'T', 'ij->e', WORK1 )
   call rot_rep( WORK3, WORK4, 'T', 'ij->n', WORK2 )

   WORK3(:,:) = WORK1(:,:)*tmask(:,:,1)/tlast_coupled
   WORK4(:,:) = WORK2(:,:)*tmask(:,:,1)/tlast_coupled

   if (l_export) then
   n = 0
   do j=lnldj,lnlej
      do i=lnldi,lnlei
         n = n + 1
         o2x_o%rAttr(index_o2x_So_dhdx,n) = WORK3(i,j)
         o2x_o%rAttr(index_o2x_So_dhdy,n) = WORK4(i,j)
      enddo
   enddo
   end if

!-----------------------------------------------------------------------
!
!     pack heat flux due to freezing/melting (W/m^2)
!     QFLUX computation and units conversion occurs in qflxice.F90
!
!-----------------------------------------------------------------------

   call lbc_lnk( QFLUX(:,:), 'T', 1._wp )

   if (l_export) then
   n = 0
   do j=lnldj,lnlej
      do i=lnldi,lnlei
         n = n + 1
         o2x_o%rAttr(index_o2x_Fioo_q,n) = QFLUX(i,j)
!         o2x_o%rAttr(index_o2x_Fioo_q,n) = c0
      enddo
   enddo
   end if

   if (.not. lfirst) then
      call iom_put('So_qflux_o2x', MAX(c0, QFLUX))
   endif

   tlast_ice  = c0
   AQICE(:,:) = c0
!   SALT_FREEZE(:,:) = c0
!   QICE(:,:)  = c0

!-----------------------------------------------------------------------
!
!     pack co2 flux, if requested (kg CO2/m^2/s)
!     units conversion occurs where co2 flux is computed
!
!-----------------------------------------------------------------------

   if (index_o2x_Faoo_fco2_ocn > 0 .AND. l_export) then
      n = 0
      do j=lnldj,lnlej
         do i=lnldi,lnlei
            n = n + 1
            o2x_o%rAttr(index_o2x_Faoo_fco2_ocn,n) = &
               SBUFF_SUM(i,j,index_o2x_Faoo_fco2_ocn)/tlast_coupled
         enddo
      enddo
   endif

!   if (.not. lfirst) CALL iom_setkt( istp )  ! iom_put outside of sbc is called at every time step

!-----------------------------------------------------------------------
!
!     diagnostics
!
!-----------------------------------------------------------------------

   if (ldiag_cpl .AND. l_export) then
      if (lwp) write(numout,*)'nemo_send_to_coupler'

      do k = 1,nsend
        n = 0
        WORKA(:,:) = c0
        do j=lnldj,lnlej
           do i=lnldi,lnlei
              n = n + 1
              WORKA(i,j) = o2x_o%rAttr(k,n)
           enddo
        enddo

        sgn = 1._wp
        if (k==index_o2x_So_u    .or. k==index_o2x_So_u     .or. &
            k==index_o2x_So_dhdx .or. k==index_o2x_So_dhdy) then
           sgn = -1._wp
        else if (k==index_o2x_Fioo_q) then
           WORKA(:,:)=MAX(c0, WORKA(:,:))
        end if
        call lbc_lnk(WORKA(:,:), 'T', sgn)

        gsum = glob_sum(WORKA(:,:)*e1e2t(:,:))

        rtmp1=minval(WORKA, mask=(tmask(:,:,1)>0._wp))
        call mpp_min(rtmp1)
        rtmp2=maxval(WORKA, mask=(tmask(:,:,1)>0._wp))
        call mpp_max(rtmp2)

        if (nproc==0) then
           call seq_flds_getField(label,k,seq_flds_o2x_fields)
           write(numout,1100)'ocn','send', TRIM(label), gsum, gsum/area, &
             rtmp1, rtmp2
        endif
      enddo ! k
      if (nproc==0) call shr_sys_flush(numout)
   endif
    
   tlast_coupled = c0

   lfirst = .false.

1100 format ('comm_diag ', a3, 1x, a4, 1x, a16, 1x, 4es28.19:)

!-----------------------------------------------------------------------
!EOC

 end subroutine ocn_export_mct

!***********************************************************************

!BOP
! !IROUTINE: ocn_sum_buffer
! !INTERFACE:

 subroutine ocn_sum_buffer

! !DESCRIPTION:
!  This routine accumulates sums for averaging fields to
!  be sent to the coupler
!
! !REVISION HISTORY:
!  same as module
! 
!EOP
!BOC

#ifdef CCSMCOUPLED
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   real (wp), dimension(jpi,jpj) ::  &
      WORK                ! local work arrays

   real (wp), dimension(jpi,jpj) ::  &
      ssgu, ssgv          ! sea surface gradient

   real (wp) ::   &
      delt                ! time interval since last step

   integer (i4) :: &
      sflux_co2_nf_ind = 0! named field index of fco2

   logical :: &
      first = .true.      ! only true for first call

   real (wp) ::   &
      ztmp                ! time interval since last step

   integer :: &
      ji, jj

!-----------------------------------------------------------------------
!
!  zero buffer if this is the first time after a coupling interval
!
!-----------------------------------------------------------------------

   if (tlast_coupled == c0) SBUFF_SUM = c0

!-----------------------------------------------------------------------
!
!  update time since last coupling
!
!-----------------------------------------------------------------------

   delt = rdt
   tlast_coupled = tlast_coupled + delt

!-----------------------------------------------------------------------
!  TODO: fco2 field handling in NEMO!
!
!  allow for fco2 field to not be registered on first call
!     because init_forcing is called before init_passive_tracers
!  use weight from previous timestep because flux used here is that
!     computed during the previous timestep
!
!-----------------------------------------------------------------------

!   if (index_o2x_Faoo_fco2_ocn > 0) then
!      if (sflux_co2_nf_ind == 0) then
!         call named_field_get_index('SFLUX_CO2', sflux_co2_nf_ind, &
!                                    exit_on_err=.not. first)
!      endif
!   endif

    ! T -> (U, V)
    ssgu(:,:) = c0
    ssgv(:,:) = c0
    DO jj = 1, jpjm1              ! Sea surface gradient (now)
       DO ji = 1, jpim1
          ssgu(ji,jj) = ( sshn(ji+1,jj) - sshn(ji,jj) ) / e1u(ji,jj)
          ssgv(ji,jj) = ( sshn(ji,jj+1) - sshn(ji,jj) ) / e2v(ji,jj)
       END DO 
    END DO 

!-----------------------------------------------------------------------
!
!  accumulate sums of U,V,T,S and GRADP
!  accumulate sum of co2 flux, if requested
!     implicitly use zero flux if fco2 field not registered yet
!  ice formation flux is handled separately in ice routine
!
!-----------------------------------------------------------------------

   SBUFF_SUM(:,:,index_o2x_So_u) =    &
      SBUFF_SUM(:,:,index_o2x_So_u) +  delt*un(:,:,1)

   SBUFF_SUM(:,:,index_o2x_So_v) =    &
      SBUFF_SUM(:,:,index_o2x_So_v) +  delt*vn(:,:,1)

   SBUFF_SUM(:,:,index_o2x_So_t) =    &
      SBUFF_SUM(:,:,index_o2x_So_t) + delt*tsn(:,:,1,jp_tem)

   SBUFF_SUM(:,:,index_o2x_So_s) =    &
      SBUFF_SUM(:,:,index_o2x_So_s) + delt*tsn(:,:,1,jp_sal)

   SBUFF_SUM(:,:,index_o2x_So_dhdx) =    &
      SBUFF_SUM(:,:,index_o2x_So_dhdx) + delt*ssgu(:,:)

   SBUFF_SUM(:,:,index_o2x_So_dhdy) =    &
      SBUFF_SUM(:,:,index_o2x_So_dhdy) + delt*ssgv(:,:)

!   if (index_o2x_Faoo_fco2_ocn > 0 .and. sflux_co2_nf_ind > 0) then
!      call named_field_get(sflux_co2_nf_ind, iblock, WORK(:,:))
!      SBUFF_SUM(:,:,index_o2x_Faoo_fco2_ocn) = &
!         SBUFF_SUM(:,:,index_o2x_Faoo_fco2_ocn) + delt*WORK(:,:)
!   endif

   first = .false.

#endif

!-----------------------------------------------------------------------
!EOC

 end subroutine ocn_sum_buffer
 
 function lice_form(kt)

   integer, intent(in) :: kt

   logical :: lice_form

   integer :: n

   lice_form = .false.

   do n=1,min(nn_nits,nn_ncpl)
      lice_form = lice_form .OR. (MOD(kt-nit000+n, nn_ncpl)==0)
      if (lice_form) exit
   end do

 end function lice_form

end module ocn_comp_mct

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
