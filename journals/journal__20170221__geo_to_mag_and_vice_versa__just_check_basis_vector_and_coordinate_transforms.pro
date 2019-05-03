;2017/02/21
PRO JOURNAL__20170221__GEO_TO_MAG_AND_VICE_VERSA__JUST_CHECK_BASIS_VECTOR_AND_COORDINATE_TRANSFORMS

  COMPILE_OPT IDL2

  orig_routineName = 'JOURNAL__20170221__GEO_TO_MAG_AND_VICE_VERSA__JUST_CHECK_BASIS_VECTOR_AND_COORDINATE_TRANSFORMS'
  R_E              = 6371.2D    ;Earth radius in km, from IGRFLIB_V2.pro

  outDir           = '/SPENCEdata/Research/database/RENU2/'

  inFile           = 'RENU2_GPS.sav'
  timeStrFile      = 'RENU2_GPS--timeStr.sav'
  outFile          = 'RENU2_coordinates.sav'


  RESTORE,outDir+inFile

  ;; GEOPACK_CONV_COORD
  ;; Description: Convert between a variety of commonly used coordinate systems.
  ;; Calling Sequence: geopack_conv_coord(_08), s1, s2, s3, d1, d2, d3.
  ;; Inputs: s1, s2, s3: Coordinates in system of origin.
  ;; Outputs: d1, d2, d3: Coordinates in target system.
  ;; Keywords: FROM_GEO: Specify source in geographic coordinates. 
  ;;  FROM_MAG: Specify source in geomagnetic coordinates.
  ;;  FROM_GEI: Specify source in geocentric equatorial inertial coordinates.
  ;;  FROM_SM: Specify source in solar magnetic coordinates.
  ;;  FROM_GSM: Specify source in geocentric solar magnetospheric
  ;;  coordinates.
  ;;  FROM_GSE: Specify source in geocentric solar ecliptic coordinates.
  ;;  TO_GEO: Specify destination in geopgraphic coordinates.
  ;;  TO_MAG: Specify destination in geomagnetic coordinates.
  ;;  TO_GEI: Specify destination in geocentric equatorial inertial coordinates.
  ;;  TO_SM: Specify destination in solar magnetic coordinates.
  ;;  TO_GSM: Specify destination in geocentric solar magnetospheric
  ;;  coordinates. 

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;;Times in CDF epoch time
  t1JulDay      = JULDAY(01,01,1997,00,00,00)
  t2JulDay      = JULDAY(12,31,2016,00,00,00)
  julArr        = TIMEGEN_EASIER(t1JulDay,t2JulDay,NMONTHS_PER_CALC=6)
  timeStr       = TIME_TO_STR(JULDAY_TO_UTC(julArr))

  time_epoch    = UTC_TO_CDF_EPOCH(JULDAY_TO_UTC(julArr))

  ;;PRINT,TIME_TO_STR(JULDAY_TO_UTC(julArr))
  ;;1997-01-01/00:00:00
  ;;1997-07-01/00:00:00
  ;;1998-01-01/00:00:00
  ;;...
  ;; 2016-01-01/00:00:00
  ;; 2016-07-01/00:00:00

  CONVERT_TIME_STRING_TO_YMDHMS_ARRAYS, $
     timeStr, $
     OUT_YEARARR=yearArr, $
     OUT_DOYARR=DOYArr, $
     OUT_MONTHARR=monthArr, $
     OUT_DAYARR=dayArr, $
     OUT_HOURARR=hourArr, $
     OUT_MINARR=minArr, $
     OUT_SECARR=secArr

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;;feed it to GEOPACK
  nTot           = N_ELEMENTS(timeStr)

  TiltArr        = !NULL

  ident          = IDENTITY(3)

  GEOin_arr      = MAKE_ARRAY(3,3,nTot,/FLOAT)

  GEO2MAG_coord  = MAKE_ARRAY(3,3,nTot,/FLOAT)
  GEO2GEI_coord  = MAKE_ARRAY(3,3,nTot,/FLOAT)

  MAG2GEO_coord  = MAKE_ARRAY(3,3,nTot,/FLOAT)
  MAG2GEI_coord  = MAKE_ARRAY(3,3,nTot,/FLOAT)

  GEO2MAG_vec    = MAKE_ARRAY(3,3,nTot,/FLOAT)
  GEO2GEI_vec    = MAKE_ARRAY(3,3,nTot,/FLOAT)

  MAG2GEO_vec    = MAKE_ARRAY(3,3,nTot,/FLOAT)
  MAG2GEI_vec    = MAKE_ARRAY(3,3,nTot,/FLOAT)

  GEO2MAG_sphC   = MAKE_ARRAY(3,3,nTot,/FLOAT)
  GEO2GEI_sphC   = MAKE_ARRAY(3,3,nTot,/FLOAT)

  MAG2GEO_sphC   = MAKE_ARRAY(3,3,nTot,/FLOAT)
  MAG2GEI_sphC   = MAKE_ARRAY(3,3,nTot,/FLOAT)

  GEO2MAG_sphV   = MAKE_ARRAY(3,3,nTot,/FLOAT)
  GEO2GEI_sphV   = MAKE_ARRAY(3,3,nTot,/FLOAT)

  MAG2GEO_sphV   = MAKE_ARRAY(3,3,nTot,/FLOAT)
  MAG2GEI_sphV   = MAKE_ARRAY(3,3,nTot,/FLOAT)

  ;; GEOSph_arr     = MAKE_ARRAY(3,3,nTot,/FLOAT)
  ;; GEOSph_arr     = MAKE_ARRAY(3,3,nTot,/FLOAT)
  ;; MAGSph_arr     = MAKE_ARRAY(3,3,nTot,/FLOAT)

  FOR i=0,nTot-1 DO BEGIN
     GEOPACK_RECALC_08,YearArr[i],MonthArr[i],DayArr[i],HourArr[i],MinArr[i],SecArr[i],/DATE 

     ;;do that dance

     ;;GEO to GEI
     GEOPACK_CONV_COORD_08,ident[0,0],ident[0,1],ident[0,2], $
                           geo2gei_x0,geo2gei_y0,geo2gei_z0, $
                           /FROM_GEO,/TO_GEI,EPOCH=time_epoch[i]

     GEOPACK_CONV_COORD_08,ident[1,0],ident[1,1],ident[1,2], $
                           geo2gei_x1,geo2gei_y1,geo2gei_z1, $
                           /FROM_GEO,/TO_GEI,EPOCH=time_epoch[i]

     GEOPACK_CONV_COORD_08,ident[2,0],ident[2,1],ident[2,2], $
                           geo2gei_x2,geo2gei_y2,geo2gei_z2, $
                           /FROM_GEO,/TO_GEI,EPOCH=time_epoch[i]

     ;;GEO to MAG
     GEOPACK_CONV_COORD_08,ident[0,0],ident[0,1],ident[0,2], $
                           mag_x0,mag_y0,mag_z0, $
                           /FROM_GEO,/TO_MAG,EPOCH=time_epoch[i]

     GEOPACK_CONV_COORD_08,ident[1,0],ident[1,1],ident[1,2], $
                           mag_x1,mag_y1,mag_z1, $
                           /FROM_GEO,/TO_MAG,EPOCH=time_epoch[i]

     GEOPACK_CONV_COORD_08,ident[2,0],ident[2,1],ident[2,2], $
                           mag_x2,mag_y2,mag_z2, $
                           /FROM_GEO,/TO_MAG,EPOCH=time_epoch[i]

     ;;MAG to GEI
     GEOPACK_CONV_COORD_08,ident[0,0],ident[0,1],ident[0,2], $
                           mag2gei_x0,mag2gei_y0,mag2gei_z0, $
                           /FROM_MAG,/TO_GEI,EPOCH=time_epoch[i]

     GEOPACK_CONV_COORD_08,ident[1,0],ident[1,1],ident[1,2], $
                           mag2gei_x1,mag2gei_y1,mag2gei_z1, $
                           /FROM_MAG,/TO_GEI,EPOCH=time_epoch[i]

     GEOPACK_CONV_COORD_08,ident[2,0],ident[2,1],ident[2,2], $
                           mag2gei_x2,mag2gei_y2,mag2gei_z2, $
                           /FROM_MAG,/TO_GEI,EPOCH=time_epoch[i]

     ;;MAG to GEO
     GEOPACK_CONV_COORD_08,ident[0,0],ident[0,1],ident[0,2], $
                           geo_x0,geo_y0,geo_z0, $
                           /FROM_MAG,/TO_GEO,EPOCH=time_epoch[i]

     GEOPACK_CONV_COORD_08,ident[1,0],ident[1,1],ident[1,2], $
                           geo_x1,geo_y1,geo_z1, $
                           /FROM_MAG,/TO_GEO,EPOCH=time_epoch[i]

     GEOPACK_CONV_COORD_08,ident[2,0],ident[2,1],ident[2,2], $
                           geo_x2,geo_y2,geo_z2, $
                           /FROM_MAG,/TO_GEO,EPOCH=time_epoch[i]

     ;;Store coord conversions
     GEO2MAG_coord[*,*,i]  = [[mag_x0,mag_y0,mag_z0], $
                              [mag_x1,mag_y1,mag_z1], $
                              [mag_x2,mag_y2,mag_z2]]
     GEO2GEI_coord[*,*,i]  = [[geo2gei_x0,geo2gei_y0,geo2gei_z0], $
                              [geo2gei_x1,geo2gei_y1,geo2gei_z1], $
                              [geo2gei_x2,geo2gei_y2,geo2gei_z2]]
     
     MAG2GEO_coord[*,*,i]  = [[geo_x0,geo_y0,geo_z0], $
                              [geo_x1,geo_y1,geo_z1], $
                              [geo_x2,geo_y2,geo_z2]]
     MAG2GEI_coord[*,*,i]  = [[mag2gei_x0,mag2gei_y0,mag2gei_z0], $
                              [mag2gei_x1,mag2gei_y1,mag2gei_z1], $
                              [mag2gei_x2,mag2gei_y2,mag2gei_z2]]

     GEO2MAG_vec[*,*,i]    = INVERT(GEO2MAG_coord[*,*,i])
     GEO2GEI_vec[*,*,i]    = INVERT(GEO2GEI_coord[*,*,i])
     MAG2GEO_vec[*,*,i]    = INVERT(MAG2GEO_coord[*,*,i])
     MAG2GEI_vec[*,*,i]    = INVERT(MAG2GEI_coord[*,*,i])
     

     GEOPACK_BCARSP_08,geo2gei_x0,geo2gei_y0,geo2gei_z0, $
                       GEO2GEI_vec[0,0,i],GEO2GEI_vec[0,1,i],GEO2GEI_vec[0,2,i], $
                       geo2geiVSph_r0,geo2geiVSph_theta0,geo2geiVSph_phi0
     GEOPACK_BCARSP_08,geo2gei_x1,geo2gei_y1,geo2gei_z1, $
                       GEO2GEI_vec[1,0,i],GEO2GEI_vec[1,1,i],GEO2GEI_vec[1,2,i], $
                       geo2geiVSph_r1,geo2geiVSph_theta1,geo2geiVSph_phi1
     GEOPACK_BCARSP_08,geo2gei_x2,geo2gei_y2,geo2gei_z2, $
                       GEO2GEI_vec[2,0,i],GEO2GEI_vec[2,1,i],GEO2GEI_vec[2,2,i], $
                       geo2geiVSph_r2,geo2geiVSph_theta2,geo2geiVSph_phi2

     GEOPACK_BCARSP_08,mag_x0,mag_y0,mag_z0, $
                       GEO2MAG_vec[0,0,i],GEO2MAG_vec[0,1,i],GEO2MAG_vec[0,2,i], $
                       magVSph_r0,magVSph_theta0,magVSph_phi0
     GEOPACK_BCARSP_08,mag_x1,mag_y1,mag_z1, $
                       GEO2MAG_vec[1,0,i],GEO2MAG_vec[1,1,i],GEO2MAG_vec[1,2,i], $
                       magVSph_r1,magVSph_theta1,magVSph_phi1
     GEOPACK_BCARSP_08,mag_x2,mag_y2,mag_z2, $
                       GEO2MAG_vec[2,0,i],GEO2MAG_vec[2,1,i],GEO2MAG_vec[2,2,i], $
                       magVSph_r2,magVSph_theta2,magVSph_phi2

     GEOPACK_BCARSP_08,geo_x0,geo_y0,geo_z0, $
                       MAG2GEO_vec[0,0,i],MAG2GEO_vec[0,1,i],MAG2GEO_vec[0,2,i], $
                       geoVSph_r0,geoVSph_theta0,geoVSph_phi0
     GEOPACK_BCARSP_08,geo_x1,geo_y1,geo_z1, $
                       MAG2GEO_vec[1,0,i],MAG2GEO_vec[1,1,i],MAG2GEO_vec[1,2,i], $
                       geoVSph_r1,geoVSph_theta1,geoVSph_phi1
     GEOPACK_BCARSP_08,geo_x2,geo_y2,geo_z2, $
                       MAG2GEO_vec[2,0,i],MAG2GEO_vec[2,1,i],MAG2GEO_vec[2,2,i], $
                       geoVSph_r2,geoVSph_theta2,geoVSph_phi2

     GEOPACK_BCARSP_08,mag2gei_x0,mag2gei_y0,mag2gei_z0, $
                       MAG2GEI_vec[0,0,i],MAG2GEI_vec[0,1,i],MAG2GEI_vec[0,2,i], $
                       mag2geiVSph_r0,mag2geiVSph_theta0,mag2geiVSph_phi0
     GEOPACK_BCARSP_08,mag2gei_x1,mag2gei_y1,mag2gei_z1, $
                       MAG2GEI_vec[1,0,i],MAG2GEI_vec[1,1,i],MAG2GEI_vec[1,2,i], $
                       mag2geiVSph_r1,mag2geiVSph_theta1,mag2geiVSph_phi1
     GEOPACK_BCARSP_08,mag2gei_x2,mag2gei_y2,mag2gei_z2, $
                       MAG2GEI_vec[2,0,i],MAG2GEI_vec[2,1,i],MAG2GEI_vec[2,2,i], $
                       mag2geiVSph_r2,mag2geiVSph_theta2,mag2geiVSph_phi2

     GEO2MAG_sphV[*,*,i]   = [[magVSph_r0,magVSph_theta0,magVSph_phi0], $
                              [magVSph_r1,magVSph_theta1,magVSph_phi1], $
                              [magVSph_r2,magVSph_theta2,magVSph_phi2]]
     GEO2GEI_sphV[*,*,i]   = [[geo2geiVSph_r0,geo2geiVSph_theta0,geo2geiVSph_phi0], $
                              [geo2geiVSph_r1,geo2geiVSph_theta1,geo2geiVSph_phi1], $
                              [geo2geiVSph_r2,geo2geiVSph_theta2,geo2geiVSph_phi2]]

     MAG2GEO_sphV[*,*,i]   = [[geoVSph_r0,geoVSph_theta0,geoVSph_phi0], $
                              [geoVSph_r1,geoVSph_theta1,geoVSph_phi1], $
                              [geoVSph_r2,geoVSph_theta2,geoVSph_phi2]]
     MAG2GEI_sphV[*,*,i]   = [[mag2geiVSph_r0,mag2geiVSph_theta0,mag2geiVSph_phi0], $
                              [mag2geiVSph_r1,mag2geiVSph_theta1,mag2geiVSph_phi1], $
                              [mag2geiVSph_r2,mag2geiVSph_theta2,mag2geiVSph_phi2]]
     

     ;;And spherical everything     
     GEOPACK_SPHCAR_08,mag_x0,mag_y0,mag_z0,mag_r0,mag_theta0,mag_phi0,/TO_SPHERE,/DEGREE
     GEOPACK_SPHCAR_08,mag_x1,mag_y1,mag_z1,mag_r1,mag_theta1,mag_phi1,/TO_SPHERE,/DEGREE
     GEOPACK_SPHCAR_08,mag_x2,mag_y2,mag_z2,mag_r2,mag_theta2,mag_phi2,/TO_SPHERE,/DEGREE

     GEOPACK_SPHCAR_08,geo2gei_x0,geo2gei_y0,geo2gei_z0,geo2gei_r0,geo2gei_theta0,geo2gei_phi0,/TO_SPHERE,/DEGREE
     GEOPACK_SPHCAR_08,geo2gei_x1,geo2gei_y1,geo2gei_z1,geo2gei_r1,geo2gei_theta1,geo2gei_phi1,/TO_SPHERE,/DEGREE
     GEOPACK_SPHCAR_08,geo2gei_x2,geo2gei_y2,geo2gei_z2,geo2gei_r2,geo2gei_theta2,geo2gei_phi2,/TO_SPHERE,/DEGREE

     GEOPACK_SPHCAR_08,geo_x0,geo_y0,geo_z0,geo_r0,geo_theta0,geo_phi0,/TO_SPHERE,/DEGREE
     GEOPACK_SPHCAR_08,geo_x1,geo_y1,geo_z1,geo_r1,geo_theta1,geo_phi1,/TO_SPHERE,/DEGREE
     GEOPACK_SPHCAR_08,geo_x2,geo_y2,geo_z2,geo_r2,geo_theta2,geo_phi2,/TO_SPHERE,/DEGREE

     GEOPACK_SPHCAR_08,mag2gei_x0,mag2gei_y0,mag2gei_z0,mag2gei_r0,mag2gei_theta0,mag2gei_phi0,/TO_SPHERE,/DEGREE
     GEOPACK_SPHCAR_08,mag2gei_x1,mag2gei_y1,mag2gei_z1,mag2gei_r1,mag2gei_theta1,mag2gei_phi1,/TO_SPHERE,/DEGREE
     GEOPACK_SPHCAR_08,mag2gei_x2,mag2gei_y2,mag2gei_z2,mag2gei_r2,mag2gei_theta2,mag2gei_phi2,/TO_SPHERE,/DEGREE

     ;;Update spherical
     GEO2MAG_sphC[*,*,i]  = [[mag_r0,mag_theta0,mag_phi0], $
                             [mag_r1,mag_theta1,mag_phi1], $
                             [mag_r2,mag_theta2,mag_phi2]]
     GEO2GEI_sphC[*,*,i]  = [[geo2gei_r0,geo2gei_theta0,geo2gei_phi0], $
                             [geo2gei_r1,geo2gei_theta1,geo2gei_phi1], $
                             [geo2gei_r2,geo2gei_theta2,geo2gei_phi2]]
     
     MAG2GEO_sphC[*,*,i]  = [[geo_r0,geo_theta0,geo_phi0], $
                             [geo_r1,geo_theta1,geo_phi1], $
                             [geo_r2,geo_theta2,geo_phi2]]
     MAG2GEI_sphC[*,*,i]  = [[mag2gei_r0,mag2gei_theta0,mag2gei_phi0], $
                             [mag2gei_r1,mag2gei_theta1,mag2gei_phi1], $
                             [mag2gei_r2,mag2gei_theta2,mag2gei_phi2]]

     PRINT,i
  ENDFOR

  STOP

END
