MODULE mod_file_attribut
  USE mod_object
  
  INTEGER,PARAMETER :: file__begin = file_object + 10
    
  INTEGER,PARAMETER :: file__name = file_object + 11
  INTEGER,PARAMETER :: file__description = file_object + 12
  INTEGER,PARAMETER :: file__output_freq = file_object + 13
  INTEGER,PARAMETER :: file__output_level = file_object + 14
  INTEGER,PARAMETER :: file__enabled = file_object + 15
  INTEGER,PARAMETER :: file__name_suffix = file_object + 16

  INTEGER,PARAMETER :: file__end = file_object + 17

END MODULE mod_file_attribut
  
