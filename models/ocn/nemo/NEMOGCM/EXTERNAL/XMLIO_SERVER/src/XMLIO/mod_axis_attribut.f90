MODULE mod_axis_attribut
  USE mod_object
  
  INTEGER,PARAMETER :: axis__begin = axis_object + 10
  
  INTEGER,PARAMETER :: axis__name = axis_object + 11
  INTEGER,PARAMETER :: axis__description = axis_object + 12
  INTEGER,PARAMETER :: axis__unit = axis_object + 13
  INTEGER,PARAMETER :: axis__size = axis_object + 14
  INTEGER,PARAMETER :: axis__positive = axis_object + 15
  INTEGER,PARAMETER :: axis__values = axis_object + 16

  INTEGER,PARAMETER :: axis__end = axis_object + 17
  
END MODULE mod_axis_attribut
