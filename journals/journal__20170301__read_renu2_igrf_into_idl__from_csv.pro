;;2017/03/01
PRO JOURNAL__20170301__READ_RENU2_IGRF_INTO_IDL__FROM_CSV

  COMPILE_OPT IDL2

  datdir   = '/SPENCEdata/Research/database/RENU2/IGRF/'
  fpref    = 'RENU2_IGRF'
  infil    = fpref + '.csv'
  infil2   = fpref + '2.csv'
  outfil   = fpref + '.sav'

  ;; OPENW,lun,datdir+outfil,/GET_LUN
  
  IGRF     = READ_CSV(datdir+infil,HEADER=header)
  IGRF2    = READ_CSV(datdir+infil2,HEADER=header2)

  this     = CREATE_STRUCT(header[0],igrf.(0))
  FOR k=1,N_ELEMENTS(tag_names(igrf))-1 DO BEGIN
     this  = CREATE_STRUCT(header[k],igrf.(k),this)
  ENDFOR

  this2    = CREATE_STRUCT(header2[0],igrf2.(0)[0])
  FOR k=1,N_ELEMENTS(tag_names(igrf2))-1 DO BEGIN
     this2 = CREATE_STRUCT(header2[k],igrf2.(k)[0],this2)
  ENDFOR

  IGRF     = CREATE_STRUCT(this2,this)
  
  PRINT,"Saving to " + datdir + outfil
  SAVE,IGRF,FILENAME=datDir+outFil

  STOP

END
