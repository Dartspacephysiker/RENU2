;2017/02/21
PRO JOURNAL__20170221__GEO_TO_MAG_AND_VICE_VERSA__JUST_CHECK_BASIS_VECTOR_AND_COORDINATE_TRANSFORMS

  COMPILE_OPT IDL2

  use_geopack_08   = 1

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
  ;; Keywords: FROM_GEO: Specify source in geopgraphic coordinates. 
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
  time_epoch    = UTC_TO_CDF_EPOCH(renu2.time.utc)

  IF FILE_TEST(outDir+timeStrFile) THEN BEGIN
     PRINT,"Restoring timeString file " + outDir+timeStrFile + ' ...'
     RESTORE,outDir+timeStrFile
  ENDIF ELSE BEGIN
     PRINT,"Making " + STRCOMPRESS(N_ELEMENTS(renu2.time.utc),/REMOVE_ALL) + " time strings ..."
     timeStr       = TIME_TO_STR(renu2.time.utc,/MS)

     PRINT,"Saving timeStrs ..."
     SAVE,timeStr,FILENAME=outDir+timeStrFile
  ENDELSE
  
  t1JulDay      = JULDAY(01,01,1997,00,00,00)
  t2JulDay      = JULDAY(12,31,2016,00,00,00)
  julArr        = TIMEGEN_EASIER(t1JulDay,t2JulDay,NMONTHS_PER_CALC=6)
  timeStr       = TIME_TO_STR(JULDAY_TO_UTC(julArr))
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

  basis_vecs     = 

  GEI_arr        = MAKE_ARRAY(3,nTot,/FLOAT)
  GEO_arr        = MAKE_ARRAY(3,nTot,/FLOAT)
  MAG_arr        = MAKE_ARRAY(3,nTot,/FLOAT)

  GEOSph_arr     = MAKE_ARRAY(3,nTot,/FLOAT)
  GEOvelSph_arr  = MAKE_ARRAY(3,nTot,/FLOAT)
  MAGSph_arr     = MAKE_ARRAY(3,nTot,/FLOAT)

  CASE 1 OF
     KEYWORD_SET(use_geopack_08): BEGIN
        FOR i=0,nTot-1 DO BEGIN
           GEOPACK_RECALC_08,YearArr[i],MonthArr[i],DayArr[i],HourArr[i],MinArr[i],SecArr[i],/DATE 

           ;;do that dance

           ;;To GEI
           GEOPACK_CONV_COORD_08,renu2.ecef.position.x[i],renu2.ecef.position.y[i],renu2.ecef.position.z[i], $
                              gei_x,gei_y,gei_z, $
                              /FROM_GEO,/TO_GEI,EPOCH=time_epoch[i]
           ;;To MAG
           GEOPACK_CONV_COORD_08,renu2.ecef.position.x[i],renu2.ecef.position.y[i],renu2.ecef.position.z[i], $
                              mag_x,mag_y,mag_z, $
                              /FROM_GEO,/TO_MAG,EPOCH=time_epoch[i]

           ;;To GSM for IGRF calc
           GEOPACK_CONV_COORD_08,renu2.ecef.position.x[i],renu2.ecef.position.y[i],renu2.ecef.position.z[i], $
                                 gsm_x,gsm_y,gsm_z, $
                              /FROM_GEO,/TO_GSW,EPOCH=time_epoch[i]


           ;;And spherical everything
           GEOPACK_SPHCAR_08,gei_x,gei_y,gei_z,gei_r,gei_theta,gei_phi,/TO_SPHERE,/DEGREE
           GEOPACK_SPHCAR_08,renu2.ecef.position.x[i],renu2.ecef.position.y[i],renu2.ecef.position.z[i], $
                          geo_r,geo_theta,geo_phi,/TO_SPHERE,/DEGREE
           GEOPACK_SPHCAR_08,renu2.ecef.velocity.x[i],renu2.ecef.velocity.y[i],renu2.ecef.velocity.z[i], $
                          geoVel_r,geoVel_theta,geoVel_phi,/TO_SPHERE,/DEGREE
           GEOPACK_SPHCAR_08,mag_x,mag_y,mag_z,mag_r,mag_theta,mag_phi,/TO_SPHERE,/DEGREE

           ;;Get IGRF
           gsm_xyz_R_E = [gsm_x,gsm_y,gsm_z]/1000.D/R_E
           ;; GEOPACK_RECALC_08,YearArr[i],MonthArr[i],DayArr[i],HourArr[i],MinArr[i],SecArr[i],/DATE 
           GEOPACK_IGRF_GEO_08,geo_r/1000./R_E,geo_theta,geo_phi,br,btheta,bphi,EPOCH=time_epoch[i],/DEGREE
           GEOPACK_BSPCAR_08  ,geo_theta,geo_phi,br,btheta,bphi,bx,by,bz,/DEGREE ; ,EPOCH=time_epoch[i]

           ;;alternate shot at IGRF
           GEOPACK_IGRF_GSW_08,gsm_xyz_R_E[0],gsm_xyz_R_E[1],gsm_xyz_R_E[2],bx_gsm,by_gsm,bz_gsm,EPOCH=time_epoch[i] ;,/DEGREE

           ;;Dipole, please?
           GEOPACK_DIP_08,gsm_xyz_R_E[1],gsm_xyz_R_E[1],gsm_xyz_R_E[2],bx_dip,by_dip,bz_dip,EPOCH=time_epoch[i] ;,/DEGREE


           ;;Update spherical
           GEISph_arr[*,i]    = [gei_theta,gei_phi,gei_r] 
           GEOSph_arr[*,i]    = [geo_theta,geo_phi,geo_r] ;Redundant, yes, but a check
           GEOvelSph_arr[*,i] = [geoVel_theta,geoVel_phi,geoVel_r]
           IGRFSph_arr[*,i]   = [btheta,bphi,br] 
           MAGSph_arr[*,i]    = [mag_theta,mag_phi,mag_r] 

           ;;Update not-spherical
           ;; TiltArr    = [TiltArr,tempTilt]
           GEI_arr[*,i]  = [gei_x,gei_y,gei_z]
           GSM_arr[*,i]  = [gsm_x,gsm_y,gsm_z]
           ;; GEO_arr[*,i]  = [geo_x,geo_y,geo_z]
           IGRF_arr[*,i] = [bx,by,bz]
           IGRFGSM_arr[*,i] = [bx_gsm,by_gsm,bz_gsm]
           IGRFDIP_arr[*,i] = [bx_dip,by_dip,bz_dip]
           MAG_arr[*,i]  = [mag_x,mag_y,mag_z]

           IF (i MOD 100) EQ 0 THEN PRINT,i
        ENDFOR
     END
     ELSE: BEGIN
        ;; FOR i=0,nTot-1 DO BEGIN
        ;;    GEOPACK_RECALC,YearArr[i],MonthArr[i],DayArr[i],HourArr[i],MinArr[i],SecArr[i],/DATE 

        ;;    ;;do that dance

        ;;    ;;To GEI
        ;;    GEOPACK_CONV_COORD,renu2.ecef.position.x[i],renu2.ecef.position.y[i],renu2.ecef.position.z[i], $
        ;;                       gei_x,gei_y,gei_z, $
        ;;                       /FROM_GEO,/TO_GEI,EPOCH=time_epoch[i]
        ;;    ;;To MAG
        ;;    GEOPACK_CONV_COORD,renu2.ecef.position.x[i],renu2.ecef.position.y[i],renu2.ecef.position.z[i], $
        ;;                       mag_x,mag_y,mag_z, $
        ;;                       /FROM_GEO,/TO_MAG,EPOCH=time_epoch[i]


        ;;    ;;And spherical everything
        ;;    GEOPACK_SPHCAR,gei_x,gei_y,gei_z,gei_r,gei_theta,gei_phi,/TO_SPHERE,/DEGREE
        ;;    GEOPACK_SPHCAR,renu2.ecef.position.x[i],renu2.ecef.position.y[i],renu2.ecef.position.z[i], $
        ;;                   geo_r,geo_theta,geo_phi,/TO_SPHERE,/DEGREE
        ;;    GEOPACK_SPHCAR,renu2.ecef.velocity.x[i],renu2.ecef.velocity.y[i],renu2.ecef.velocity.z[i], $
        ;;                   geoVel_r,geoVel_theta,geoVel_phi,/TO_SPHERE,/DEGREE
        ;;    GEOPACK_SPHCAR,mag_x,mag_y,mag_z,mag_r,mag_theta,mag_phi,/TO_SPHERE,/DEGREE

        ;;    ;;Get IGRF
        ;;    GEOPACK_IGRF_GEO,geo_r/R_E/1000.D,geo_theta,geo_phi,br,btheta,bphi,/DEGREE,EPOCH=time_epoch[i]
        ;;    GEOPACK_BSPCAR  ,geo_theta,geo_phi,br,btheta,bphi,bx,by,bz,/DEGREE ; ,EPOCH=time_epoch[i]

        ;;    ;;Update spherical
        ;;    GEISph_arr[*,i]    = [gei_theta,gei_phi,gei_r] 
        ;;    GEOSph_arr[*,i]    = [geo_theta,geo_phi,geo_r] ;Redundant, yes, but a check
        ;;    GEOvelSph_arr[*,i] = [geoVel_theta,geoVel_phi,geoVel_r]
        ;;    IGRFSph_arr[*,i]   = [btheta,bphi,br] 
        ;;    MAGSph_arr[*,i]    = [mag_theta,mag_phi,mag_r] 

        ;;    ;;Update not-spherical
        ;;    ;; TiltArr    = [TiltArr,tempTilt]
        ;;    GEI_arr[*,i]  = [gei_x,gei_y,gei_z]
        ;;    ;; GEO_arr[*,i]  = [geo_x,geo_y,geo_z]
        ;;    IGRF_arr[*,i] = [bx,by,bz]
        ;;    MAG_arr[*,i]  = [mag_x,mag_y,mag_z]

        ;;    IF (i MOD 100) EQ 0 THEN PRINT,i
        ;; ENDFOR
     END
  ENDCASE

  STOP

END
