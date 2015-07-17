MODULE mod_field_attribut
  USE mod_object

  INTEGER,PARAMETER :: field__begin = field_object + 10
  
  INTEGER,PARAMETER :: field__name = field_object + 11
  INTEGER,PARAMETER :: field__description = field_object + 12
  INTEGER,PARAMETER :: field__unit = field_object + 13
  INTEGER,PARAMETER :: field__operation = field_object + 14
  INTEGER,PARAMETER :: field__freq_op = field_object + 15
  INTEGER,PARAMETER :: field__axis_ref = field_object + 16
  INTEGER,PARAMETER :: field__grid_ref = field_object + 17
  INTEGER,PARAMETER :: field__zoom_ref = field_object + 18
  INTEGER,PARAMETER :: field__level = field_object + 19
  INTEGER,PARAMETER :: field__prec = field_object + 20
  INTEGER,PARAMETER :: field__field_ref = field_object + 21
  INTEGER,PARAMETER :: field__enabled = field_object + 22

  INTEGER,PARAMETER :: field__end = field_object + 23
  
END MODULE mod_field_attribut
