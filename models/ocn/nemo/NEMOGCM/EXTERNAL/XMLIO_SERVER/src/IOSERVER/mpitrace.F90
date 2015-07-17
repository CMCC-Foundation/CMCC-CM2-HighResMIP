module mpitrace
  INTEGER,parameter :: VTcheck_buffer=91
  INTEGER,parameter :: VTprocess_event=92
  INTEGER,parameter :: VTbuffer_full=93
  
  INTEGER,parameter :: nb_inst=7
  INTEGER :: MPE_begin(nb_inst)
  INTEGER :: MPE_end(nb_inst)
  
contains

  subroutine InitVampir
    implicit none

#ifdef USE_VT
    include 'VT.inc'
    integer :: ierr
    
    call VTSYMDEF(VTcheck_buffer,"check_buffer","check_buffer",ierr)
    call VTSYMDEF(VTprocess_event,"process_event","process_event",ierr)
    call VTSYMDEF(VTbuffer_full,"buffer_full","buffer_full",ierr)
    PRINT *,'definition des symbole VT'
#endif

#ifdef USE_MPE
    include 'mpe_logf.h' 
    integer :: ierr,i
    
    DO i=1,nb_inst
      ierr = MPE_Log_get_state_eventIDs( MPE_begin(i), MPE_end(i) )
    ENDDO
    
    ierr = MPE_Describe_state( MPE_begin(VTcaldyn), MPE_end(VTcaldyn),"caldyn", "yellow" )
    ierr = MPE_Describe_state( MPE_begin(VTintegre), MPE_end(VTintegre),"integre", "blue" )
    ierr = MPE_Describe_state( MPE_begin(VTadvection), MPE_end(VTadvection),"advection", "green" )
    ierr = MPE_Describe_state( MPE_begin(VTdissipation), MPE_end(VTdissipation),"dissipation", "ivory" )
    ierr = MPE_Describe_state( MPE_begin(VThallo), MPE_end(VThallo),"hallo", "orange" )
    ierr = MPE_Describe_state( MPE_begin(VTphysiq), MPE_end(VTphysiq),"physiq", "purple" )
    ierr = MPE_Describe_state( MPE_begin(VTinca), MPE_end(VTinca),"inca", "LightBlue" )
#endif     
  end subroutine InitVampir

  subroutine VTb(number)
    implicit none
    INTEGER :: number
#ifdef USE_VT    
    include 'VT.inc'
    integer :: ierr
    
    call VTBEGIN(number,ierr)
#endif 
#ifdef USE_MPE
    include 'mpe_logf.h' 
    integer :: ierr,i
    ierr = MPE_Log_event( MPE_begin(number), 0, '' )
#endif

  end subroutine VTb

  subroutine VTe(number)
    implicit none
    INTEGER :: Number
#ifdef USE_VT    
    include 'VT.inc'
    integer :: ierr
   
    call VTEND(number,ierr)
#endif    

#ifdef USE_MPE
    include 'mpe_logf.h' 
    integer :: ierr,i
    ierr = MPE_Log_event( MPE_end(number), 0, '' )
#endif

  end subroutine VTe
  
end module mpitrace
