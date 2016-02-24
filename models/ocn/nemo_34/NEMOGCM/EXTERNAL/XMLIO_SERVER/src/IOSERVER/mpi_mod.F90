MODULE mpi_mod

#if ! defined key_mpp_mpi
    LOGICAL, PARAMETER :: using_mpi=.FALSE.
    INCLUDE 'ios_mpif.h'
#else
    LOGICAL, PARAMETER :: using_mpi=.TRUE.
    INCLUDE 'mpif.h'
#endif

END MODULE mpi_mod
