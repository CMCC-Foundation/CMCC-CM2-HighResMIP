MODULE wrk_nemo
   !!======================================================================
   !!                       ***  MODULE  wrk_nemo  ***
   !! NEMO work space:  define and allocate work-space arrays used in 
   !! all components of NEMO
   !!======================================================================
   !! History :  4.0  !  2011-01  (A Porter)  Original code
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   wrk_alloc         : define in memory the work space arrays
   !!   wrk_in_use, iwrk_in_use, wrk_in_use_xz : check the availability of a workspace 
   !!   wrk_not_released, iwrk_not_released, wrk_not_released_xz : release the workspace
   !!   print_in_use_list : print out the table holding which workspace arrays are currently marked as in use
   !!   get_next_arg      : get the next argument
   !!   wrk_stop          : act as local alternative to ctl_stop
   !!----------------------------------------------------------------------
   USE par_oce        ! ocean parameters

   IMPLICIT NONE
   PRIVATE

   PUBLIC   wrk_alloc   ! function called in nemogcm module (nemo_init routine)
   PUBLIC   wrk_in_use, iwrk_in_use, wrk_in_use_xz                     ! function called almost everywhere
   PUBLIC   wrk_not_released, iwrk_not_released, wrk_not_released_xz   ! function called almost everywhere

   INTEGER, PARAMETER :: num_1d_wrkspaces  = 27   ! No. of 1D workspace arrays ( MAX(jpi*jpj,jpi*jpk,jpj*jpk) )
   INTEGER, PARAMETER :: num_2d_wrkspaces  = 37   ! No. of 2D workspace arrays (jpi,jpj)
   INTEGER, PARAMETER :: num_3d_wrkspaces  = 15   ! No. of 3D workspace arrays (jpi,jpj,jpk)
   INTEGER, PARAMETER :: num_4d_wrkspaces  = 4    ! No. of 4D workspace arrays (jpi,jpj,jpk,jpts)

   INTEGER, PARAMETER :: num_xz_wrkspaces  = 4   ! No. of 2D, xz workspace arrays (jpi,jpk)

   INTEGER, PARAMETER :: num_1d_lwrkspaces = 0   ! No. of 1D logical workspace arrays
   INTEGER, PARAMETER :: num_2d_lwrkspaces = 3   ! No. of 2D logical workspace arrays
   INTEGER, PARAMETER :: num_3d_lwrkspaces = 1   ! No. of 3D logical workspace arrays
   INTEGER, PARAMETER :: num_4d_lwrkspaces = 0   ! No. of 4D logical workspace arrays

   INTEGER, PARAMETER :: num_1d_iwrkspaces = 0   ! No. of 1D integer workspace arrays
   INTEGER, PARAMETER :: num_2d_iwrkspaces = 1   ! No. of 2D integer workspace arrays
   INTEGER, PARAMETER :: num_3d_iwrkspaces = 0   ! No. of 3D integer workspace arrays
   INTEGER, PARAMETER :: num_4d_iwrkspaces = 0   ! No. of 4D integer workspace arrays
   ! Maximum no. of workspaces of any one dimensionality that can be
   ! requested - MAX(num_1d_wrkspaces, num_2d_wrkspaces, num_3d_wrkspaces, num_4d_wrkspaces) 
   INTEGER :: max_num_wrkspaces = 37

   ! If adding more arrays here, remember to increment the appropriate 
   ! num_Xd_wrkspaces parameter above and to allocate them in wrk_alloc()

   !                                                               !!**  1D, REAL(wp) workspaces  **
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:)      , TARGET, PUBLIC ::   wrk_1d_1 , wrk_1d_2 , wrk_1d_3 , wrk_1d_4 , wrk_1d_5
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:)      , TARGET, PUBLIC ::   wrk_1d_6 , wrk_1d_7 , wrk_1d_8 , wrk_1d_9 , wrk_1d_10
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:)      , TARGET, PUBLIC ::   wrk_1d_11, wrk_1d_12, wrk_1d_13, wrk_1d_14, wrk_1d_15
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:)      , TARGET, PUBLIC ::   wrk_1d_16, wrk_1d_17, wrk_1d_18, wrk_1d_19, wrk_1d_20
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:)      , TARGET, PUBLIC ::   wrk_1d_21, wrk_1d_22, wrk_1d_23, wrk_1d_24, wrk_1d_25
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:)      , TARGET, PUBLIC ::   wrk_1d_26, wrk_1d_27

   !                                                               !!**  2D, x-y, REAL(wp) workspaces  **
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:)    , TARGET, PUBLIC ::   wrk_2d_1 , wrk_2d_2 , wrk_2d_3 , wrk_2d_4 , wrk_2d_5
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:)    , TARGET, PUBLIC ::   wrk_2d_6 , wrk_2d_7 , wrk_2d_8 , wrk_2d_9 , wrk_2d_10
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:)    , TARGET, PUBLIC ::   wrk_2d_11, wrk_2d_12, wrk_2d_13, wrk_2d_14, wrk_2d_15
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:)    , TARGET, PUBLIC ::   wrk_2d_16, wrk_2d_17, wrk_2d_18, wrk_2d_19, wrk_2d_20
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:)    , TARGET, PUBLIC ::   wrk_2d_21, wrk_2d_22, wrk_2d_23, wrk_2d_24, wrk_2d_25
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:)    , TARGET, PUBLIC ::   wrk_2d_26, wrk_2d_27, wrk_2d_28, wrk_2d_29, wrk_2d_30
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:)    , TARGET, PUBLIC ::   wrk_2d_31, wrk_2d_32, wrk_2d_33, wrk_2d_34, wrk_2d_35
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:)    , TARGET, PUBLIC ::   wrk_2d_36, wrk_2d_37

   !                                                               !!**  2D, x-z, REAL(wp) workspaces  **
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:)            , PUBLIC ::   wrk_xz_1, wrk_xz_2, wrk_xz_3, wrk_xz_4 
   
   !                                                               !!**  3D, x-y-z, REAL(wp) workspaces  **
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:)  , TARGET, PUBLIC ::   wrk_3d_1 , wrk_3d_2 , wrk_3d_3 , wrk_3d_4 , wrk_3d_5
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:)  , TARGET, PUBLIC ::   wrk_3d_6 , wrk_3d_7 , wrk_3d_8 , wrk_3d_9 , wrk_3d_10
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:)  , TARGET, PUBLIC ::   wrk_3d_11, wrk_3d_12, wrk_3d_13, wrk_3d_14, wrk_3d_15

   !                                                               !!**  4D, x-y-z-tra, REAL(wp) workspaces  **
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:,:), TARGET, PUBLIC ::   wrk_4d_1, wrk_4d_2, wrk_4d_3, wrk_4d_4 
   
   !                                                               !!** 2D integer workspace  **
   INTEGER , ALLOCATABLE, SAVE, DIMENSION(:,:)            , PUBLIC ::   iwrk_2d_1

   LOGICAL, DIMENSION(num_1d_wrkspaces)  ::   in_use_1d     !: Flags to track which 1D workspace arrays are in use  
   LOGICAL, DIMENSION(num_2d_wrkspaces)  ::   in_use_2d     !: Flags to track which 2D workspace arrays are in use
   LOGICAL, DIMENSION(num_3d_wrkspaces)  ::   in_use_3d     !: Flags to track which 3D workspace arrays are in use
   LOGICAL, DIMENSION(num_4d_wrkspaces)  ::   in_use_4d     !: Flags to track which 4D workspace arrays are in use
   LOGICAL, DIMENSION(num_xz_wrkspaces)  ::   in_use_xz     !: Flags to track which 2D, xz workspace arrays are in use
   LOGICAL, DIMENSION(num_2d_lwrkspaces) ::   in_use_2dll   !: Flags to track which 2D, logical workspace arrays are in use
   LOGICAL, DIMENSION(num_3d_lwrkspaces) ::   in_use_3dll   !: Flags to track which 3D, logical workspace arrays are in use
   LOGICAL, DIMENSION(num_2d_iwrkspaces) ::   in_use_2di    !: Flags to track which 2D, integer workspace arrays are in use

   ! Labels for specifying workspace type in call to print_in_use_list()
   INTEGER, PARAMETER ::   INTEGER_TYPE = 0
   INTEGER, PARAMETER ::   LOGICAL_TYPE = 1
   INTEGER, PARAMETER ::   REAL_TYPE    = 2

   INTEGER :: kumout  ! Local copy of numout unit number for error/warning messages
   LOGICAL :: llwp    ! Local copy of lwp - whether we are master PE or not

   CHARACTER(LEN=*), PARAMETER ::   cform_err2 = "(/,' ===>>> : E R R O R',     /,'         ===========',/)"       !:
   CHARACTER(LEN=*), PARAMETER ::   cform_war2 = "(/,' ===>>> : W A R N I N G', /,'         ===============',/)"   !:

   !!----------------------------------------------------------------------
   !! NEMO/OPA 4.0 , NEMO Consortium (2011)
   !! $Id:$
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

  FUNCTION wrk_alloc(iunit, lwp_arg)
      !!----------------------------------------------------------------------
      !!                   ***  FUNCTION wrk_alloc  ***
      !!
      !! ** Purpose :   Define in memory once for all the NEMO 2D, 3D and 4d 
      !!                work space arrays
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   iunit         ! Unit no. to use for error/warning messages in this module
      LOGICAL, INTENT(in) ::   lwp_arg       ! Value of lwp
      !
      INTEGER ::   wrk_alloc   ! Return value
      INTEGER ::   extent_1d   ! Extent to allocate for 1D arrays
      INTEGER ::   ierror(6)   ! local integer
      !!----------------------------------------------------------------------
      !
      ! Save the unit number to use for err/warning messages
      kumout = iunit
      ! Save whether we are master PE or not (for output messages)
      llwp = lwp_arg
      !
      ! Extent to use for 1D work arrays - find the maximum product of 
      ! jpi*jpj, jpi*jpk and jpj*jpk and use that
      IF    ( jpi < jpj .AND. jpi < jpk ) THEN   ;   extent_1d = jpj*jpk
      ELSEIF( jpj < jpi .AND. jpj < jpk ) THEN   ;   extent_1d = jpi*jpk
      ELSE                                       ;   extent_1d = jpi*jpj
      ENDIF
      !
      ! Initialise the 'in use' flags for each work-space array
      in_use_1d  (:) = .FALSE.
      in_use_2d  (:) = .FALSE.
      in_use_3d  (:) = .FALSE.
      in_use_4d  (:) = .FALSE.
      in_use_xz  (:) = .FALSE.
      in_use_2dll(:) = .FALSE.
      in_use_3dll(:) = .FALSE.
      in_use_2di (:) = .FALSE.
      !
      ierror(:) = 0
      !
      ALLOCATE( wrk_1d_1 (extent_1d) , wrk_1d_2 (extent_1d) , wrk_1d_3 (extent_1d) , wrk_1d_4 (extent_1d) ,     &
         &      wrk_1d_5 (extent_1d) , wrk_1d_6 (extent_1d) , wrk_1d_7 (extent_1d) , wrk_1d_8 (extent_1d) ,     &
         &      wrk_1d_9 (extent_1d) , wrk_1d_10(extent_1d)                                               ,     &
         &      wrk_1d_11(extent_1d) , wrk_1d_12(extent_1d) , wrk_1d_13(extent_1d) , wrk_1d_14(extent_1d) ,     &
         &      wrk_1d_15(extent_1d) , wrk_1d_16(extent_1d) , wrk_1d_17(extent_1d) , wrk_1d_18(extent_1d) ,     &
         &      wrk_1d_19(extent_1d) , wrk_1d_20(extent_1d)                                               ,     &
         &      wrk_1d_21(extent_1d) , wrk_1d_22(extent_1d) , wrk_1d_23(extent_1d) , wrk_1d_24(extent_1d) ,     &
         &      wrk_1d_25(extent_1d) , wrk_1d_26(extent_1d) , wrk_1d_27(extent_1d)                        , STAT=ierror(1) )
         !
      ALLOCATE( wrk_2d_1 (jpi,jpj) , wrk_2d_2 (jpi,jpj) , wrk_2d_3 (jpi,jpj) , wrk_2d_4 (jpi,jpj) ,     & 
         &      wrk_2d_5 (jpi,jpj) , wrk_2d_6 (jpi,jpj) , wrk_2d_7 (jpi,jpj) , wrk_2d_8 (jpi,jpj) ,     &
         &      wrk_2d_9 (jpi,jpj) , wrk_2d_10(jpi,jpj)                                           ,     &
         &      wrk_2d_11(jpi,jpj) , wrk_2d_12(jpi,jpj) , wrk_2d_13(jpi,jpj) , wrk_2d_14(jpi,jpj) ,     &
         &      wrk_2d_15(jpi,jpj) , wrk_2d_16(jpi,jpj) , wrk_2d_17(jpi,jpj) , wrk_2d_18(jpi,jpj) ,     &
         &      wrk_2d_19(jpi,jpj) , wrk_2d_20(jpi,jpj)                                           ,     &
         &      wrk_2d_21(jpi,jpj) , wrk_2d_22(jpi,jpj) , wrk_2d_23(jpi,jpj) , wrk_2d_24(jpi,jpj) ,     &
         &      wrk_2d_25(jpi,jpj) , wrk_2d_26(jpi,jpj) , wrk_2d_27(jpi,jpj) , wrk_2d_28(jpi,jpj) ,     &
         &      wrk_2d_29(jpi,jpj) , wrk_2d_30(jpi,jpj)                                           ,     &
         &      wrk_2d_31(jpi,jpj) , wrk_2d_32(jpi,jpj) , wrk_2d_33(jpi,jpj) , wrk_2d_34(jpi,jpj) ,     &
         &      wrk_2d_35(jpi,jpj) , wrk_2d_36(jpi,jpj) , wrk_2d_37(jpi,jpj)                      , STAT=ierror(2) )
         !
      ALLOCATE( wrk_3d_1 (jpi,jpj,jpk) , wrk_3d_2 (jpi,jpj,jpk) , wrk_3d_3 (jpi,jpj,jpk) , wrk_3d_4 (jpi,jpj,jpk) ,     &
         &      wrk_3d_5 (jpi,jpj,jpk) , wrk_3d_6 (jpi,jpj,jpk) , wrk_3d_7 (jpi,jpj,jpk) , wrk_3d_8 (jpi,jpj,jpk) ,     &
         &      wrk_3d_9 (jpi,jpj,jpk) , wrk_3d_10(jpi,jpj,jpk)                                                   ,     & 
         &      wrk_3d_11(jpi,jpj,jpk) , wrk_3d_12(jpi,jpj,jpk) , wrk_3d_13(jpi,jpj,jpk) , wrk_3d_14(jpi,jpj,jpk) ,     & 
         &      wrk_3d_15(jpi,jpj,jpk)                                                                            , STAT=ierror(3) )
         !
      ALLOCATE( wrk_4d_1(jpi,jpj,jpk,jpts) , wrk_4d_2(jpi,jpj,jpk,jpts),     &
         &      wrk_4d_3(jpi,jpj,jpk,jpts) , wrk_4d_4(jpi,jpj,jpk,jpts), STAT=ierror(4) )
         !
      ALLOCATE( wrk_xz_1(jpi,jpk) , wrk_xz_2(jpi,jpk) , wrk_xz_3(jpi,jpk) , wrk_xz_4(jpi,jpk) , STAT=ierror(5) )
         !
      ALLOCATE( iwrk_2d_1(jpi,jpj)      , STAT=ierror(6) )
      !
      wrk_alloc = MAXVAL( ierror )
      !
      ! Calling routine, nemo_alloc(), checks for errors and takes 
      ! appropriate action - we just print a warning message
      IF( wrk_alloc /= 0 ) THEN
         WRITE(kumout,cform_war2)
         WRITE(kumout,*) 'wrk_alloc: allocation of workspace arrays failed'
      ENDIF
      !
   END FUNCTION wrk_alloc


   FUNCTION wrk_in_use( kdim,    index1,  index2,  index3,  index4,    &
      &                 index5,  index6,  index7,  index8,  index9,    &
      &                 index10, index11, index12, index13, index14,   &
      &                 index15, index16, index17, index18, index19,   &
      &                 index20, index21, index22, index23, index24,   &
      &                 index25, index26, index27)
      !!----------------------------------------------------------------------
      !!                   ***  FUNCTION wrk_in_use  ***
      !!
      !! ** Purpose :   Request a set of KIND(wp) workspaces to use. Returns 
      !!                .TRUE. if any of those requested are already in use, 
      !!                .FALSE. otherwise. 
      !!
      !! ** Method  :   Sets internal flags to signal that requested workspaces are in use.
      !!                key_no_workspace_check defined ==> always return FALSE
      !!----------------------------------------------------------------------
      INTEGER          , INTENT(in) ::   kdim        ! Dimensionality of requested workspace(s)
      INTEGER          , INTENT(in) ::   index1      ! Index of first requested workspace
      INTEGER, OPTIONAL, INTENT(in) ::             index2 , index3 , index4 , index5 , index6 , index7 , index8 , index9, index10
      INTEGER, OPTIONAL, INTENT(in) ::   index11, index12, index13, index14, index15, index16, index17, index18, index19, index20
      INTEGER, OPTIONAL, INTENT(in) ::   index21, index22, index23, index24, index25, index26, index27
      !
      LOGICAL ::   wrk_in_use      ! Return value
      INTEGER ::   iarg, iptr   ! local integer
      !!----------------------------------------------------------------------
      !
      wrk_in_use = .FALSE.
      !
#if ! defined   key_no_workspace_check   ||   ! defined   key_agrif
      ! NB: check not available with AGRIF
      !
      iptr    = index1
      iarg    = 1
      !
      DO WHILE( .NOT. wrk_in_use .AND. iarg <= max_num_wrkspaces )
         !
         IF( kdim == 1 ) THEN
            IF( iptr > num_1d_wrkspaces ) THEN
               CALL wrk_stop('wrk_in_use - more 1D workspace arrays requested than defined in wrk_nemo module')
               wrk_in_use = .TRUE.
               EXIT
            ELSEIF( in_use_1d(iptr) ) THEN
               wrk_in_use = .TRUE.
               CALL print_in_use_list(1, REAL_TYPE, in_use_1d)
            ENDIF
            in_use_1d(iptr) = .TRUE.
            !
         ELSEIF( kdim == 2 ) THEN
            IF( iptr > num_2d_wrkspaces ) THEN
               CALL wrk_stop('wrk_in_use - more 2D workspace arrays requested than defined in wrk_nemo module')
               wrk_in_use = .TRUE.
               EXIT
            ELSEIF( in_use_2d(iptr) ) THEN
               wrk_in_use = .TRUE.
               CALL print_in_use_list(2, REAL_TYPE, in_use_2d)
            ENDIF
            in_use_2d(iptr) = .TRUE.
            !
         ELSEIF( kdim == 3 ) THEN
            IF( iptr > num_3d_wrkspaces ) THEN
               CALL wrk_stop( 'wrk_in_use - more 3D workspace arrays requested than defined in wrk_nemo module' )
               wrk_in_use = .TRUE.
               EXIT
            ELSEIF( in_use_3d(iptr) ) THEN
               wrk_in_use = .TRUE.
               CALL print_in_use_list(3, REAL_TYPE, in_use_3d)
            ENDIF
            in_use_3d(iptr) = .TRUE.
            !
         ELSEIF( kdim == 4 ) THEN
            IF(iptr > num_4d_wrkspaces)THEN
               CALL wrk_stop( 'wrk_in_use - more 4D workspace arrays requested than defined in wrk_nemo module' )
               wrk_in_use = .TRUE.
               EXIT
            ELSEIF( in_use_4d(iptr) ) THEN
               wrk_in_use = .TRUE.
               CALL print_in_use_list( 4, REAL_TYPE, in_use_4d )
            ENDIF
            !
            in_use_4d(iptr) = .TRUE.
            !
         ELSE 
            IF(llwp) WRITE(kumout,*) 'wrk_in_use: unsupported value of kdim = ',kdim
            CALL wrk_stop( 'wrk_in_use: unrecognised value for number of dimensions' )
         ENDIF

         CALL get_next_arg( iarg  ,  iptr  ,  index2,  index3,  index4,    &
            &               index5,  index6,  index7,  index8,  index9,    &
            &               index10, index11, index12, index13, index14,   &
            &               index15, index16, index17, index18, index19,   &
            &               index20, index21, index22, index23, index24,   &
            &               index25, index26, index27)

         IF( iarg == -1 ) THEN      ! We've checked all of the arguments and are done
            EXIT
         ELSEIF( iarg == -99 ) THEN
            CALL wrk_stop( 'wrk_in_use : caught unexpected argument count - BUG' )
            EXIT
         ENDIF
         !
      END DO ! end of DO WHILE()
#endif
      !
   END FUNCTION wrk_in_use


   FUNCTION iwrk_in_use( kdim, index1, index2, index3, index4,   &
      &                        index5, index6, index7 )
      !!----------------------------------------------------------------------
      !!                   ***  FUNCTION iwrk_in_use  ***
      !!
      !! ** Purpose :   Request a set of INTEGER workspaces to use. Returns 
      !!                .TRUE. if any of those requested are already in use, 
      !!                .FALSE. otherwise. 
      !!
      !! ** Method  :   Sets internal flags to signal that requested workspaces
      !!                are in use.
      !!----------------------------------------------------------------------
      INTEGER          , INTENT(in) ::   kdim        ! Dimensionality of requested workspace(s)
      INTEGER          , INTENT(in) ::   index1      ! Index of first requested workspace
      INTEGER, OPTIONAL, INTENT(in) ::   index2, index3, index4, index5, index6, index7
      !
      LOGICAL ::   iwrk_in_use    ! Return value
      INTEGER ::   iarg, iptr
      !!----------------------------------------------------------------------
      !
      iwrk_in_use = .FALSE.
      !
#if ! defined   key_no_workspace_check   ||   ! defined   key_agrif
      ! NB: check not available with AGRIF
      !
      iptr     = index1
      iarg     = 1
      !
      DO WHILE( .NOT.iwrk_in_use .AND. iarg <= max_num_wrkspaces )
         !
         IF( kdim == 2 ) THEN
            IF( iptr > num_2d_wrkspaces ) THEN
               CALL wrk_stop( 'wrk_in_use - more 2D workspace arrays requested than defined in wrk_nemo module' )
               iwrk_in_use = .TRUE.
            ELSEIF( in_use_2di(iptr) ) THEN
               iwrk_in_use = .TRUE.
               CALL print_in_use_list( 2, INTEGER_TYPE, in_use_2di )
            ENDIF
            in_use_2di(iptr) = .TRUE.
            !
         ELSE
            IF(llwp) WRITE(kumout,*) 'iwrk_in_use: unsupported value of kdim = ',kdim
            CALL wrk_stop('iwrk_in_use: unsupported value for number of dimensions')
         ENDIF
         !
         SELECT CASE (iarg)         ! Move on to next optional argument
         CASE ( 1 )
            IF( .NOT. PRESENT(index2) ) THEN   ;   EXIT
            ELSE                               ;   iarg = 2   ;   iptr = index2
            ENDIF
         CASE ( 2 )
            IF( .NOT. PRESENT(index3) ) THEN   ;   EXIT
            ELSE                               ;   iarg = 3   ;   iptr = index3
            ENDIF
         CASE ( 3 )
            IF( .NOT. PRESENT(index4) ) THEN   ;   EXIT
            ELSE                               ;   iarg = 4   ;   iptr = index4
            ENDIF
         CASE ( 4 )
            IF( .NOT. PRESENT(index5) ) THEN   ;   EXIT
            ELSE                               ;   iarg = 5   ;   iptr = index5
            ENDIF
         CASE ( 5 )
            IF( .NOT. PRESENT(index6) ) THEN   ;   EXIT
            ELSE                               ;   iarg = 6   ;   iptr = index6
            ENDIF
         CASE ( 6 )
            IF( .NOT. PRESENT(index7) ) THEN   ;   EXIT
            ELSE                               ;   iarg = 7   ;   iptr = index7
            ENDIF
         CASE ( 7 )
            EXIT
         CASE DEFAULT
            CALL wrk_stop( 'iwrk_in_use : caught unexpected argument count - BUG' )
            EXIT
         END SELECT
         !
      END DO ! end of DO WHILE()
#endif
      !
   END FUNCTION iwrk_in_use


   FUNCTION wrk_in_use_xz( index1, index2, index3, index4,   &
      &                    index5, index6, index7, index8, index9 )
      !!----------------------------------------------------------------------
      !!                   ***  FUNCTION wrk_in_use_xz  ***
      !!
      !! ** Purpose :   Request a set of 2D, xz (jpi,jpk) workspaces to use. 
      !!                Returns .TRUE. if any of those requested are already in
      !!                use, .FALSE. otherwise. 
      !!
      !! ** Method  :   Sets internal flags to signal that requested workspaces
      !!                are in use.
      !!----------------------------------------------------------------------
      INTEGER          , INTENT(in) ::   index1      ! Index of first requested workspace
      INTEGER, OPTIONAL, INTENT(in) ::   index2, index3, index4, index5
      INTEGER, OPTIONAL, INTENT(in) ::   index6, index7, index8, index9
      !
      LOGICAL ::   wrk_in_use_xz   ! Return value
      INTEGER ::   iarg, iptr      ! local integer
      !!----------------------------------------------------------------------
      !
      wrk_in_use_xz = .FALSE.
      !
#if ! defined   key_no_workspace_check   ||   ! defined   key_agrif
      ! NB: check not available with AGRIF
      !
      iptr = index1
      iarg = 1
      !
      DO WHILE( .NOT. wrk_in_use_xz .AND. iarg <= max_num_wrkspaces )
         !
         IF( iptr > num_xz_wrkspaces ) THEN
            CALL wrk_stop('wrk_in_use_xz - more 2D xz workspace arrays requested than defined in wrk_nemo module')
            wrk_in_use_xz = .TRUE.
            EXIT
         ELSE IF( in_use_xz(iptr) ) THEN
            wrk_in_use_xz = .TRUE.
            CALL print_in_use_list(2, REAL_TYPE, in_use_xz) !ARPDBG - bug
         ENDIF
         !
         in_use_xz(iptr) = .TRUE.
         !
         CALL get_next_arg( iarg  , iptr  , index2, index3, index4,   &
            &               index5, index6, index7, index8, index9 )
         !
         IF( iarg == -1 ) THEN      ! We've checked all of the arguments and are done
            EXIT
         ELSEIF( iarg == -99 ) THEN
            CALL wrk_stop( 'wrk_in_use_xz : caught unexpected argument count - BUG' )   ;   EXIT
         ENDIF
         !
      END DO ! while( (.NOT. wrk_in_use_xz) .AND. iarg <= max_num_wrkspaces)
#endif
      !
   END FUNCTION wrk_in_use_xz


   FUNCTION wrk_not_released( kdim,    index1,  index2,  index3,  index4,  &
      &                       index5,  index6,  index7,  index8,  index9,  &
      &                       index10, index11, index12, index13, index14, &
      &                       index15, index16, index17, index18, index19, &
      &                       index20, index21, index22, index23, index24, &
      &                       index25, index26, index27)
      !!----------------------------------------------------------------------
      !!                 ***  FUNCTION wrk_not_released  ***
      !!
      !! ** Purpose :   Flag that the specified workspace arrays are no-longer
      !!                in use.
      !!----------------------------------------------------------------------
      INTEGER          , INTENT(in) ::   kdim             ! Dimensionality of workspace(s)
      INTEGER          , INTENT(in) ::   index1           ! Index of 1st workspace to release
      INTEGER, OPTIONAL, INTENT(in) ::            index2 , index3 , index4 , index5 , index6 , index7 , index8 , index9 , index10
      INTEGER, OPTIONAL, INTENT(in) ::   index11, index12, index13, index14, index15, index16, index17, index18, index19, index20
      INTEGER, OPTIONAL, INTENT(in) ::   index21, index22, index23, index24, index25, index26, index27
      !
      LOGICAL ::   wrk_not_released   ! Return value
      INTEGER ::   iarg, iptr
      !!----------------------------------------------------------------------
      !
      wrk_not_released = .FALSE.
      !
#if ! defined   key_no_workspace_check   ||   ! defined   key_agrif
      ! NB: check not available with AGRIF
      !
      iptr = index1
      iarg = 1
      !
      DO WHILE( iarg <= max_num_wrkspaces )
         !
         IF( kdim == 1 ) THEN
            IF( iptr > num_1d_wrkspaces ) THEN
               CALL wrk_stop( 'wrk_not_released : attempt to release a non-existent 1D workspace array' )
               wrk_not_released = .TRUE.
            ELSE
               in_use_1d(iptr) = .FALSE.
            ENDIF
            !
         ELSE IF(kdim == 2) THEN
            IF( iptr > num_2d_wrkspaces ) THEN
               CALL wrk_stop( 'wrk_not_released : attempt to release a non-existent 2D workspace array' )
               wrk_not_released = .TRUE.
            ENDIF
            in_use_2d(iptr) = .FALSE.
            !
         ELSEIF( kdim == 3 ) THEN
            IF( iptr > num_3d_wrkspaces ) THEN
               CALL wrk_stop('wrk_not_released : attempt to release a non-existent 3D workspace array')
               wrk_not_released = .TRUE.
            ENDIF
            in_use_3d(iptr) = .FALSE.
            !
          ELSEIF( kdim == 4 ) THEN
            IF( iptr > num_4d_wrkspaces ) THEN
               CALL wrk_stop('wrk_not_released : attempt to release a non-existent 4D workspace array')
               wrk_not_released = .TRUE.
            ENDIF
            in_use_4d(iptr) = .FALSE.
            !
         ELSE 
            IF(llwp) WRITE(kumout,*) 'wrk_not_released: unsupported value of kdim = ',kdim
            CALL wrk_stop('wrk_not_released: unrecognised value for number of dimensions')
         ENDIF
         !
         ! Move on to next optional argument
         CALL get_next_arg( iarg  ,  iptr  ,  index2,  index3,  index4,   &
            &               index5,  index6,  index7,  index8,  index9,   &
            &               index10, index11, index12, index13,           &
            &               index14, index15, index16, index17,           &
            &               index18, index19, index20, index21,           &
            &               index22, index23, index24, index25,           &
            &               index26, index27 )

         IF( iarg == -1 ) THEN      ! We've checked all of the arguments and are done
            EXIT
         ELSEIF( iarg == -99 ) THEN
             CALL wrk_stop('wrk_not_released - caught unexpected argument count - BUG')   ;   EXIT
         ENDIF
         !
      END DO ! end of DO WHILE()
#endif
      !
   END FUNCTION wrk_not_released


   FUNCTION iwrk_not_released( kdim, index1, index2, index3, index4,   &
      &                              index5, index6, index7 )
      !!----------------------------------------------------------------------
      !!                 ***  FUNCTION iwrk_not_released  ***
      !!
      !! ** Purpose :   Flag that the specified INTEGER workspace arrays are
      !!                no-longer in use.
      !!----------------------------------------------------------------------
      INTEGER          , INTENT(in) ::   kdim             ! Dimensionality of workspace(s)
      INTEGER          , INTENT(in) ::   index1           ! Index of 1st workspace to release
      INTEGER, OPTIONAL, INTENT(in) ::   index2, index3, index4, index5, index6, index7
      !
      LOGICAL :: iwrk_not_released   ! Return value
      INTEGER :: iarg, iptr          ! local integer
      !!----------------------------------------------------------------------
      !
      iwrk_not_released = .FALSE.
      !
#if ! defined   key_no_workspace_check   ||   ! defined   key_agrif
      ! NB: check not available with AGRIF
      !
      iptr = index1
      iarg = 1
      !
      DO WHILE(iarg <= max_num_wrkspaces)
         !
         IF( kdim == 2 ) THEN
            IF( iptr > num_2d_iwrkspaces ) THEN
               CALL wrk_stop('iwrk_not_released : attempt to release a non-existant 2D workspace array')
               iwrk_not_released = .TRUE.
            ENDIF
            in_use_2di(iptr) = .FALSE.
         ELSE 
            IF(llwp) WRITE(kumout,*) 'iwrk_not_released: unsupported value of kdim = ',kdim
            CALL wrk_stop('iwrk_not_released: unsupported value for number of dimensions')
         ENDIF
         !
         ! Move on to next optional argument
         SELECT CASE (iarg)
         CASE ( 1 )
            IF( .NOT. PRESENT(index2) ) THEN   ;   EXIT
            ELSE                               ;   iarg = 2   ;   iptr = index2
            ENDIF
         CASE ( 2 )
            IF( .NOT. PRESENT(index3) ) THEN   ;   EXIT
            ELSE                               ;   iarg = 3   ;   iptr = index3
            ENDIF
         CASE ( 3 )
            IF( .NOT. PRESENT(index4) ) THEN   ;   EXIT
            ELSE                               ;   iarg = 4   ;   iptr = index4
            ENDIF
         CASE ( 4 )
            IF( .NOT. PRESENT(index5) ) THEN   ;   EXIT
            ELSE                               ;   iarg = 5   ;   iptr = index5
            ENDIF
         CASE ( 5 )
            IF( .NOT. PRESENT(index6) ) THEN   ;   EXIT
            ELSE                               ;   iarg = 6   ;   iptr = index6
            ENDIF
         CASE ( 6 )
            IF( .NOT. PRESENT(index7) ) THEN   ;   EXIT
            ELSE                               ;   iarg = 7   ;   iptr = index7
            ENDIF
         CASE ( 7 )
            EXIT
         CASE DEFAULT
            CALL wrk_stop( 'iwrk_not_released : caught unexpected argument count - BUG' )
            EXIT
         END SELECT
         !
      END DO ! end of DO WHILE()
#endif
      !
   END FUNCTION iwrk_not_released


   FUNCTION wrk_not_released_xz( index1, index2, index3, index4, index5,   &
      &                          index6, index7, index8, index9 )
      !!----------------------------------------------------------------------
      !!                 ***  FUNCTION wrk_not_released_xz  ***
      !!
      !!----------------------------------------------------------------------
      INTEGER          , INTENT(in) ::   index1   ! Index of 1st workspace to release
      INTEGER, OPTIONAL, INTENT(in) ::   index2, index3, index4, index5, index6, index7, index8, index9
      !
      LOGICAL ::   wrk_not_released_xz   ! Return value
      INTEGER ::   iarg, iptr            ! local integer
      !!----------------------------------------------------------------------
      !
      wrk_not_released_xz = .FALSE.
      !
#if ! defined   key_no_workspace_check   ||   ! defined   key_agrif
      ! NB: check not available with AGRIF
      !
      iptr           = index1
      iarg           = 1
      !
      DO WHILE( iarg <= max_num_wrkspaces )
         !
         IF( iptr > num_xz_wrkspaces ) THEN
            CALL wrk_stop('wrk_not_released_xz : attempt to release a non-existant 2D xz workspace array')
            wrk_not_released_xz = .TRUE.
            EXIT
         ENDIF
         in_use_xz(iptr) = .FALSE.
         !
         ! Move on to next optional argument
         CALL get_next_arg( iarg, iptr, index2, index3, index4,   &
            &                           index5, index6, index7, index8, index9)
         !
         IF(  iarg == -1 ) THEN     ! We've checked all of the arguments and are done
            EXIT
         ELSEIF( iarg == -99 ) THEN
            CALL wrk_stop('wrk_not_released_xz : caught unexpected argument count - BUG')
            EXIT
         ENDIF
         !
      END DO ! while (iarg <= max_num_wrkspaces)
#endif
      !
   END FUNCTION wrk_not_released_xz


   SUBROUTINE print_in_use_list( kdim, itype, in_use_list )
      !!----------------------------------------------------------------------
      !!                 *** ROUTINE print_in_use_list ***
      !!
      !! ** Purpose:   to print out the table holding which workspace arrays
      !!             are currently marked as in use.
      !!----------------------------------------------------------------------
      INTEGER,               INTENT(in) :: kdim
      INTEGER,               INTENT(in) :: itype
      LOGICAL, DIMENSION(:), INTENT(in) :: in_use_list
      !
      INTEGER          ::   ji, icount
      CHARACTER(LEN=7) ::   type_string
      !!----------------------------------------------------------------------
      !
      IF(.NOT. llwp)   RETURN
      !
      SELECT CASE ( kdim )
      !
      CASE (1)
         SELECT CASE (itype)
         CASE (INTEGER_TYPE)   ;   icount = num_1d_iwrkspaces
         CASE (LOGICAL_TYPE)   ;   icount = num_1d_lwrkspaces
         CASE (REAL_TYPE   )   ;   icount = num_1d_wrkspaces
         END SELECT
         !
      CASE (2)
         SELECT CASE (itype)
         CASE (INTEGER_TYPE)   ;   icount = num_2d_iwrkspaces
         CASE (LOGICAL_TYPE)   ;   icount = num_2d_lwrkspaces
         CASE (REAL_TYPE   )   ;   icount = num_2d_wrkspaces
         END SELECT
         !
      CASE (3)
         SELECT CASE (itype)
         CASE (INTEGER_TYPE)   ;   icount = num_3d_iwrkspaces
         CASE (LOGICAL_TYPE)   ;   icount = num_3d_lwrkspaces
         CASE (REAL_TYPE   )   ;   icount = num_3d_wrkspaces
         END SELECT
         !
      CASE (4)
         SELECT CASE (itype)
         CASE (INTEGER_TYPE)   ;   icount = num_4d_iwrkspaces
         CASE (LOGICAL_TYPE)   ;   icount = num_4d_lwrkspaces
         CASE (REAL_TYPE   )   ;   icount = num_4d_wrkspaces
         END SELECT
         !
      CASE DEFAULT   ;   RETURN
      !
      END SELECT
      !
      ! Set character string with type of workspace
      SELECT CASE (itype)
      CASE (INTEGER_TYPE)   ;   type_string = "INTEGER" 
      CASE (LOGICAL_TYPE)   ;   type_string = "LOGICAL"
      CASE (REAL_TYPE   )   ;   type_string = "REAL" 
      END SELECT
      !
      WRITE(kumout,*)
      WRITE(kumout,"('------------------------------------------')")
      WRITE(kumout,"('Table of ',I1,'D ',(A),' workspaces currently in use:')") kdim, TRIM(type_string)
      WRITE(kumout,"('Workspace   In use')")
      DO ji = 1, icount, 1
         WRITE(kumout,"(4x,I2,8x,L1)") ji, in_use_list(ji)
      END DO
      WRITE(kumout,"('------------------------------------------')")
      WRITE(kumout,*)
      !
   END SUBROUTINE print_in_use_list


   SUBROUTINE get_next_arg( iargidx, iargval, index2,  index3,  index4,  &
      &                     index5 , index6,  index7,  index8,  index9,  &
      &                     index10, index11, index12, index13, index14, &
      &                     index15, index16, index17, index18, index19, &
      &                     index20, index21, index22, index23, index24, &
      &                     index25, index26, index27 )
      !!----------------------------------------------------------------------
      INTEGER          , INTENT(inout) ::   iargidx   ! Index of current arg
      INTEGER          , INTENT(inout) ::   iargval   ! Value of current arg
      INTEGER, OPTIONAL, INTENT(in   ) ::            index2 , index3 , index4 , index5 , index6 , index7 , index8 , index9 , index10
      INTEGER, OPTIONAL, INTENT(in   ) ::   index11, index12, index13, index14, index15, index16, index17, index18, index19, index20
      INTEGER, OPTIONAL, INTENT(in   ) ::   index21, index22, index23, index24, index25, index26, index27
      !!----------------------------------------------------------------------
      !
      SELECT CASE (iargidx)       ! Move on to next optional argument
      CASE ( 1 )
         IF( .NOT. PRESENT(index2 ) ) THEN   ;   iargidx = -1
         ELSE                                ;   iargidx =  2   ;   iargval = index2
         ENDIF
      CASE ( 2 )
         IF( .NOT. PRESENT(index3 ) ) THEN   ;   iargidx = -1
         ELSE                                ;   iargidx =  3   ;   iargval = index3
         ENDIF
      CASE ( 3 )
         IF( .NOT. PRESENT(index4 ) ) THEN   ;   iargidx = -1
         ELSE                                ;   iargidx =  4   ;   iargval = index4
         ENDIF
      CASE ( 4 )
         IF( .NOT. PRESENT(index5 ) ) THEN   ;   iargidx = -1
         ELSE                                ;   iargidx =  5   ;   iargval = index5
         ENDIF
      CASE ( 5 )
         IF( .NOT. PRESENT(index6 ) ) THEN   ;   iargidx = -1
         ELSE                                ;   iargidx =  6   ;   iargval = index6
         ENDIF
      CASE ( 6 )
         IF( .NOT. PRESENT(index7 ) ) THEN   ;   iargidx = -1
         ELSE                                ;   iargidx =  7   ;   iargval = index7
         ENDIF
      CASE ( 7 )
         IF( .NOT. PRESENT(index8 ) ) THEN   ;   iargidx = -1
         ELSE                                ;   iargidx =  8   ;   iargval = index8
         ENDIF
      CASE ( 8 )
         IF( .NOT. PRESENT(index9 ) ) THEN   ;   iargidx = -1
         ELSE                                ;   iargidx =  9   ;   iargval = index9
         ENDIF
      CASE ( 9 )
         IF( .NOT. PRESENT(index10) ) THEN   ;   iargidx = -1
         ELSE                                ;   iargidx = 10   ;   iargval = index10
         ENDIF
      CASE ( 10 )
         IF( .NOT. PRESENT(index11) ) THEN   ;   iargidx = -1
         ELSE                                ;   iargidx = 11   ;   iargval = index11
         ENDIF
      CASE ( 11 )
         IF( .NOT. PRESENT(index12) ) THEN   ;   iargidx = -1
         ELSE                                ;   iargidx = 12   ;   iargval = index12
         ENDIF
      CASE ( 12 )
         IF( .NOT. PRESENT(index13) ) THEN   ;   iargidx = -1
         ELSE                                ;   iargidx =  13   ;   iargval = index13
         ENDIF
      CASE ( 13 )
         IF( .NOT. PRESENT(index14) ) THEN   ;   iargidx = -1
         ELSE                                ;   iargidx = 14   ;   iargval = index14
         ENDIF
      CASE ( 14 )
         IF( .NOT. PRESENT(index15) ) THEN   ;   iargidx = -1
         ELSE                                ;   iargidx = 15   ;   iargval = index15
         ENDIF
      CASE ( 15 )
         IF( .NOT. PRESENT(index16) ) THEN   ;   iargidx = -1
         ELSE                                ;   iargidx = 16   ;   iargval = index16
         ENDIF
      CASE ( 16 )
         IF( .NOT. PRESENT(index17) ) THEN   ;   iargidx = -1
         ELSE                                ;   iargidx = 17   ;   iargval = index17
         ENDIF
      CASE ( 17 )
         IF( .NOT. PRESENT(index18) ) THEN   ;   iargidx = -1
         ELSE                                ;   iargidx = 18   ;   iargval = index18
         ENDIF
      CASE ( 18 )
         IF( .NOT. PRESENT(index19) ) THEN   ;   iargidx = -1
         ELSE                                ;   iargidx = 19   ;   iargval = index19
         ENDIF
      CASE ( 19 )
         IF( .NOT. PRESENT(index20) ) THEN   ;   iargidx = -1
         ELSE                                ;   iargidx = 20   ;   iargval = index20
         ENDIF
      CASE ( 20 )
         IF( .NOT. PRESENT(index21) ) THEN   ;   iargidx = -1
         ELSE                                ;   iargidx = 21   ;   iargval = index21
         ENDIF
      CASE ( 21 )
         IF( .NOT. PRESENT(index22) ) THEN   ;   iargidx = -1
         ELSE                                ;   iargidx = 22   ;   iargval = index22
         ENDIF
      CASE ( 22 )
         IF( .NOT. PRESENT(index23) ) THEN   ;   iargidx = -1
         ELSE                                ;   iargidx = 23   ;   iargval = index23
         ENDIF
      CASE ( 23 )
         IF( .NOT. PRESENT(index24) ) THEN   ;   iargidx = -1
         ELSE                                ;   iargidx = 24   ;   iargval = index24
         ENDIF
      CASE ( 24 )
         IF( .NOT. PRESENT(index25) ) THEN   ;   iargidx = -1
         ELSE                                ;   iargidx = 25   ;   iargval = index25
         ENDIF
      CASE ( 25 )
         IF( .NOT. PRESENT(index26) ) THEN   ;   iargidx = -1
         ELSE                                ;   iargidx = 26   ;   iargval = index26
         ENDIF
      CASE ( 26 )
         IF( .NOT. PRESENT(index27) ) THEN   ;   iargidx = -1
         ELSE                                ;   iargidx = 27   ;   iargval = index27
         ENDIF
      CASE ( 27 )
         iargidx = -1
      CASE DEFAULT
         ! BUG - iargidx shouldn't take any other values!
         ! Flag error for calling routine
         iargidx = -99
      END SELECT
      !
   END SUBROUTINE get_next_arg


   SUBROUTINE wrk_stop(cmsg)
      !!----------------------------------------------------------------------
      !!               ***  ROUTINE wrk_stop  ***
      !! ** Purpose :   to act as local alternative to ctl_stop. 
      !!                Avoids dependency on in_out_manager module.
      !!----------------------------------------------------------------------
      CHARACTER(LEN=*), INTENT(in) :: cmsg
      !!----------------------------------------------------------------------
      !
      WRITE(kumout, cform_err2)
      WRITE(kumout,*) TRIM(cmsg)
      ! ARPDBG - would like to call mppstop here to force a stop but that
      ! introduces a dependency on lib_mpp. Could call mpi_abort() directly
      ! but that's fairly brutal. Better to rely on calling routine to
      ! deal with the error passed back from the wrk_X routine?
      !CALL mppstop
      !
   END SUBROUTINE wrk_stop

   !!=====================================================================
END MODULE wrk_nemo
