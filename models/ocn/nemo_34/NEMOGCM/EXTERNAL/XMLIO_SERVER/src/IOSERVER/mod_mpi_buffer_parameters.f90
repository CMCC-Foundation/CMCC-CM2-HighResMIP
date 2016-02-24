MODULE mod_mpi_buffer_parameters
  USE mod_ioserver_namelist, ONLY : global_mpi_buffer_size
  INTEGER,PARAMETER :: max_request = 10000
  REAL, PARAMETER   :: frac_alert  = 0.9
  INTEGER           :: tag_iocomm  = 4000
  INTEGER           :: tag_terminate = 4001
END MODULE mod_mpi_buffer_parameters
