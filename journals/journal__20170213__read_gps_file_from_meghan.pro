;;2017/02/13
;;Is that all it took? Just call the routine?
PRO JOURNAL__20170213__READ_GPS_FILE_FROM_MEGHAN

  COMPILE_OPT IDL2

  saveFile = '/SPENCEdata/Research/database/RENU2/RENU2_GPS.sav'

  READ_RENU2_GPS_CSV,SAVEFILE=saveFile

END
