MODULE mod_time_parameters

  INTEGER,POINTER  :: initial_timestep
  REAL,POINTER     :: initial_date
  REAL,POINTER     :: timestep_value 
  INTEGER,POINTER  :: timestep_number

CONTAINS 

  SUBROUTINE time_parameters__swap_context(saved_initial_timestep,saved_initial_date, &
                                           saved_timestep_value,saved_timestep_number)
  IMPLICIT NONE
    INTEGER,POINTER  :: saved_initial_timestep
    REAL,POINTER     :: saved_initial_date
    REAL,POINTER     :: saved_timestep_value 
    INTEGER,POINTER  :: saved_timestep_number
   
    initial_timestep => saved_initial_timestep
    initial_date     => saved_initial_date
    timestep_value   => saved_timestep_value 
    timestep_number  => saved_timestep_number    
    
  END SUBROUTINE time_parameters__swap_context  
    
                                           
END MODULE mod_time_parameters
