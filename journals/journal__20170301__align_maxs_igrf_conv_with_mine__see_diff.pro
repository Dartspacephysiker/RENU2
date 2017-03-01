;;2017/03/01
PRO JOURNAL__20170301__ALIGN_MAXS_IGRF_CONV_WITH_MINE__SEE_DIFF

  COMPILE_OPT IDL2

  datdir = '/SPENCEdata/Research/database/RENU2/'
  maxfile = 'IGRF/RENU2_IGRF.sav'
  mycoordfile = 'RENU2_coordinates.sav'
  myotherfile = 'RENU2_GPS.sav'

  RESTORE,datdir+maxfile
  RESTORE,datdir+mycoordfile
  RESTORE,datdir+myotherfile

  STOP

END
