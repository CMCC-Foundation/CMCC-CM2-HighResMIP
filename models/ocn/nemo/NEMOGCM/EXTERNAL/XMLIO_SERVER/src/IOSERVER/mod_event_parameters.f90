MODULE mod_event_parameters
  INTEGER, PARAMETER :: event_id_parse_xml_file        = 100
  INTEGER, PARAMETER :: event_id_set_local_domain      = 101
  INTEGER, PARAMETER :: event_id_set_vert_axis         = 102
  INTEGER, PARAMETER :: event_id_set_time_parameters   = 103
  INTEGER, PARAMETER :: event_id_set_timestep          = 104
  INTEGER, PARAMETER :: event_id_write_field2d         = 105
  INTEGER, PARAMETER :: event_id_write_field3d         = 106
  INTEGER, PARAMETER :: event_id_close_io_definition   = 107
  INTEGER, PARAMETER :: event_id_set_grid_dimension    = 108
  INTEGER, PARAMETER :: event_id_set_grid_domain       = 109
  INTEGER, PARAMETER :: event_id_set_grid_type_nemo    = 110
  INTEGER, PARAMETER :: event_id_enable_field          = 111
  INTEGER, PARAMETER :: event_id_disable_field         = 112
  INTEGER, PARAMETER :: event_id_swap_context          = 113
  INTEGER, PARAMETER :: event_id_set_grid_type_lmdz    = 114
  INTEGER, PARAMETER :: event_id_write_field1d         = 115
  INTEGER, PARAMETER :: event_id_set_attribut          = 116
  INTEGER, PARAMETER :: event_id_set_calendar          = 117
  INTEGER, PARAMETER :: event_id_stop_ioserver         = 999
END MODULE mod_event_parameters
