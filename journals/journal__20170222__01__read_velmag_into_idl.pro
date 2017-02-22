;;2017/02/22
;;.csv file is generated by journal__20170217__read_velmag_into_python__output_csv.py
PRO JOURNAL__20170222__01__READ_VELMAG_INTO_IDL

  COMPILE_OPT IDL2

  fDir   = '/SPENCEdata/Research/database/RENU2/'
  fName  = 'velMAG.csv'
  fTmplt = 'velMAG_csv_tmplt.idl'

  outF   = 'velMAG.idl'

  ;; tmplt = ASCII_TEMPLATE(fDir+fName)
  ;; SAVE,tmplt,FILENAME=fDir+fTmplt
  RESTORE,fDir+fTmplt

  velMag = READ_ASCII(fDir+fName,TEMPLATE=tmplt)
  STR_ELEMENT,velMag,'vTot',SQRT(velMag.v1*velMag.v1+velMag.v2*velMag.v2+velMag.v3*velMag.v3),/ADD_REPLACE

  PRINT,"Saving velMag to " + outF + ' ...'
  ;; SAVE,velMag,FILENAME=fDir+outF
  STOP

END
