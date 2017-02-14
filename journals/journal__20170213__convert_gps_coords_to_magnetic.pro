;;2017/02/13
PRO JOURNAL__20170213__CONVERT_GPS_COORDS_TO_MAGNETIC

  COMPILE_OPT IDL2

  use_geopack_08   = 1

  orig_routineName = 'JOURNAL__20170213__CONVERT_GPS_COORDS_TO_MAGNETIC'
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
  
  ;; GEOPACK_RECALC_08
  ;; GEOPACK_CONV_COORD_08,cuspLoc_MAG[0],cuspLoc_MAG[1],cuspLoc_MAG[2],clgeo_x,clgeo_y,clgeo_z,/FROM_MAG,/TO_GEO,EPOCH=time_epoch

  YearArr       = FIX(STRMID(timeStr,0,4))
  MonthArr      = FIX(STRMID(timeStr,5,2))
  DayArr        = FIX(STRMID(timeStr,8,2))
  HourArr       = FIX(STRMID(timeStr,11,2))
  MinArr        = FIX(STRMID(timeStr,14,2))
  SecArr        = FLOAT(STRMID(timeStr,17,6))

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;;feed it to GEOPACK
  nTot           = N_ELEMENTS(renu2.time.utc)

  TiltArr        = !NULL

  GEI_arr        = MAKE_ARRAY(3,nTot,/FLOAT)
  GSM_arr        = MAKE_ARRAY(3,nTot,/FLOAT)
  ;; GEO_arr        = MAKE_ARRAY(3,nTot,/FLOAT)
  IGRF_arr       = MAKE_ARRAY(3,nTot,/FLOAT)
  IGRFGSM_arr    = MAKE_ARRAY(3,nTot,/FLOAT)
  IGRFDIP_arr    = MAKE_ARRAY(3,nTot,/FLOAT)
  MAG_arr        = MAKE_ARRAY(3,nTot,/FLOAT)

  GEISph_arr     = MAKE_ARRAY(3,nTot,/FLOAT)
  GEOSph_arr     = MAKE_ARRAY(3,nTot,/FLOAT)
  GEOvelSph_arr  = MAKE_ARRAY(3,nTot,/FLOAT)
  IGRFSph_arr    = MAKE_ARRAY(3,nTot,/FLOAT)
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
        FOR i=0,nTot-1 DO BEGIN
           GEOPACK_RECALC,YearArr[i],MonthArr[i],DayArr[i],HourArr[i],MinArr[i],SecArr[i],/DATE 

           ;;do that dance

           ;;To GEI
           GEOPACK_CONV_COORD,renu2.ecef.position.x[i],renu2.ecef.position.y[i],renu2.ecef.position.z[i], $
                              gei_x,gei_y,gei_z, $
                              /FROM_GEO,/TO_GEI,EPOCH=time_epoch[i]
           ;;To MAG
           GEOPACK_CONV_COORD,renu2.ecef.position.x[i],renu2.ecef.position.y[i],renu2.ecef.position.z[i], $
                              mag_x,mag_y,mag_z, $
                              /FROM_GEO,/TO_MAG,EPOCH=time_epoch[i]


           ;;And spherical everything
           GEOPACK_SPHCAR,gei_x,gei_y,gei_z,gei_r,gei_theta,gei_phi,/TO_SPHERE,/DEGREE
           GEOPACK_SPHCAR,renu2.ecef.position.x[i],renu2.ecef.position.y[i],renu2.ecef.position.z[i], $
                          geo_r,geo_theta,geo_phi,/TO_SPHERE,/DEGREE
           GEOPACK_SPHCAR,renu2.ecef.velocity.x[i],renu2.ecef.velocity.y[i],renu2.ecef.velocity.z[i], $
                          geoVel_r,geoVel_theta,geoVel_phi,/TO_SPHERE,/DEGREE
           GEOPACK_SPHCAR,mag_x,mag_y,mag_z,mag_r,mag_theta,mag_phi,/TO_SPHERE,/DEGREE

           ;;Get IGRF
           GEOPACK_IGRF_GEO,geo_r/R_E/1000.D,geo_theta,geo_phi,br,btheta,bphi,/DEGREE,EPOCH=time_epoch[i]
           GEOPACK_BSPCAR  ,geo_theta,geo_phi,br,btheta,bphi,bx,by,bz,/DEGREE ; ,EPOCH=time_epoch[i]

           ;;Update spherical
           GEISph_arr[*,i]    = [gei_theta,gei_phi,gei_r] 
           GEOSph_arr[*,i]    = [geo_theta,geo_phi,geo_r] ;Redundant, yes, but a check
           GEOvelSph_arr[*,i] = [geoVel_theta,geoVel_phi,geoVel_r]
           IGRFSph_arr[*,i]   = [btheta,bphi,br] 
           MAGSph_arr[*,i]    = [mag_theta,mag_phi,mag_r] 

           ;;Update not-spherical
           ;; TiltArr    = [TiltArr,tempTilt]
           GEI_arr[*,i]  = [gei_x,gei_y,gei_z]
           ;; GEO_arr[*,i]  = [geo_x,geo_y,geo_z]
           IGRF_arr[*,i] = [bx,by,bz]
           MAG_arr[*,i]  = [mag_x,mag_y,mag_z]

           IF (i MOD 100) EQ 0 THEN PRINT,i
        ENDFOR
     END
  ENDCASE


  ;;Lat, long, height
  GEISph_arr2   = [ $
                  [90.-REFORM(GEISph_arr[0,*])], $
                  [REFORM(GEISph_arr[1,*])], $
                  [REFORM(GEISph_arr[2,*])-R_E] $ ;Convert to latitude from colatitude here
                  ]

  GEOSph_arr2   = [ $
                  [90.-REFORM(GEOSph_arr[0,*])], $
                  [REFORM(GEOSph_arr[1,*])], $
                  [REFORM(GEOSph_arr[2,*])-R_E] $ ;Convert to latitude from colatitude here
                  ]   

  MAGSph_arr2   = [ $
                  [90.-REFORM(MAGSph_arr[0,*])], $
                  [REFORM(MAGSph_arr[1,*])], $
                  [REFORM(MAGSph_arr[2,*])-R_E] $ ;Convert to latitude from colatitude here
                  ]   
  
  IGRF    = {sph : {bTheta : REFORM(IGRFSph_arr[0,*]), $
                    bPhi   : REFORM(IGRFSph_arr[1,*]), $
                    br     : REFORM(IGRFSph_arr[2,*])}, $
             car : {x      : REFORM(IGRF_arr[0,*]), $
                    y      : REFORM(IGRF_arr[1,*]), $
                    z      : REFORM(IGRF_arr[2,*])}, $
             cargsm : {x   : REFORM(IGRFGSM_arr[0,*]), $
                       y   : REFORM(IGRFGSM_arr[1,*]), $
                       z   : REFORM(IGRFGSM_arr[2,*])}, $
             cardip : {x   : REFORM(IGRFDIP_arr[0,*]), $
                       y   : REFORM(IGRFDIP_arr[1,*]), $
                       z   : REFORM(IGRFDIP_arr[2,*])}}
  
  GEI     = {ALT:GEISph_arr2[*,2], $
             LON:GEISph_arr2[*,1], $
             LAT:GEISph_arr2[*,0], $
             sph : {theta : REFORM(GEISph_arr[0,*]), $
                    phi   : REFORM(GEISph_arr[1,*]), $
                    r     : REFORM(GEISph_arr[2,*])}, $
             car : {x     : REFORM(GEI_arr[0,*]), $
                    y     : REFORM(GEI_arr[1,*]), $
                    z     : REFORM(GEI_arr[2,*])}}

  GEO     = {ALT:GEOSph_arr2[*,2], $
             LON:GEOSph_arr2[*,1], $
             LAT:GEOSph_arr2[*,0], $
             sph : {theta : REFORM(GEOSph_arr[0,*]), $
                    phi   : REFORM(GEOSph_arr[1,*]), $
                    r     : REFORM(GEOSph_arr[2,*])}, $
             ;; car : {x     : REFORM(GEO_arr[0,*]), $
             ;;        y     : REFORM(GEO_arr[1,*]), $
             ;;        z     : REFORM(GEO_arr[2,*])}}
             car : {x     : renu2.ecef.position.x, $
                    y     : renu2.ecef.position.y, $
                    z     : renu2.ecef.position.z}}

  GEOvel  = {sph : {theta : REFORM(GEOvelSph_arr[0,*]), $
                    phi   : REFORM(GEOvelSph_arr[1,*]), $
                    r     : REFORM(GEOvelSph_arr[2,*])}, $
             car : {x     : renu2.ecef.velocity.x, $
                    y     : renu2.ecef.velocity.y, $
                    z     : renu2.ecef.velocity.z}}



  MAG     = {ALT:MAGSph_arr2[*,2], $
             LON:MAGSph_arr2[*,1], $
             LAT:MAGSph_arr2[*,0], $
             sph : {theta : REFORM(MAGSph_arr[0,*]), $
                    phi   : REFORM(MAGSph_arr[1,*]), $
                    r     : REFORM(MAGSph_arr[2,*])}, $
             car : {x     : REFORM(MAG_arr[0,*]), $
                    y     : REFORM(MAG_arr[1,*]), $
                    z     : REFORM(MAG_arr[2,*])}}


  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;;make struct
  ;; renu2Coords = {GEO     : GEO_arr, $
  ;;                MAG     : MAG_arr, $
  ;;                GEI     : GEI_arr, $
  ;;                IGRF    : IGRF_arr, $
  ;;                CREATED : GET_TODAY_STRING(/DO_YYYYMMDD_FMT), $
  ;;                ORIGINATING_ROUTINE: orig_routineName}
  renu2Coords = {GEI     : GEI, $
                 GEO     : GEO, $
                 GEOvel  : GEOvel, $
                 MAG     : MAG, $
                 IGRF    : IGRF, $
                 CREATED : GET_TODAY_STRING(/DO_YYYYMMDD_FMT), $
                 ORIGINATING_ROUTINE: orig_routineName}

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;;Save it
  PRINT,'Saving ' + outDir + outFile + '...'
  save,renu2coords,FILENAME=outDir+outFile

  PRINT,"Did it!"

  geocarmagnitude  = SQRT(renu2.ecef.position.x*renu2.ecef.position.x+ $
                          renu2.ecef.position.y*renu2.ecef.position.y+ $
                          renu2.ecef.position.z*renu2.ecef.position.z)
  geocarmagnit2    = SQRT(geo.car.x*geo.car.x+geo.car.y*geo.car.y+geo.car.z*geo.car.z)
  ;; geosphmagnitude  = SQRT(geo.sph.r*geo.sph.r+geo.sph.theta*geo.sph.theta+geo.sph.phi*geo.sph.phi)
  geosphmagnitude  = geo.sph.r
  diffgeo          = geosphmagnitude - $
                     geocarmagnit2

  velcarMagnitude  = SQRT(geovel.car.x*geovel.car.x+geovel.car.y*geovel.car.y+geovel.car.z*geovel.car.z)
  ;; velsphMagnitude  = SQRT(geovel.sph.r*geovel.sph.r+geovel.sph.theta*geovel.sph.theta+geovel.sph.phi*geovel.sph.phi)
  velsphMagnitude  = geovel.sph.r ;;pff, idiot
  diffVel          = velsphMagnitude- $
                     velcarmagnitude

  bsphmagnitude    = SQRT(igrf.sph.br*igrf.sph.br+igrf.sph.btheta*igrf.sph.btheta+igrf.sph.bphi*igrf.sph.bphi)
  bcarmagnitude    = SQRT(igrf.car.x*igrf.car.x+igrf.car.y*igrf.car.y+igrf.car.z*igrf.car.z)
  bcargsmmagnitude = SQRT(igrf.cargsm.x*igrf.cargsm.x+igrf.cargsm.y*igrf.cargsm.y+igrf.cargsm.z*igrf.cargsm.z)
  diffB            = bsphmagnitude-bcarmagnitude

  magPercentDiff   = (bsphmagnitude[1:-1]-bsphmagnitude[0:-2])/bsphmagnitude[1:-1]*100.
  magPercentDiffGSM= (bcargsmmagnitude[1:-1]-bcargsmmagnitude[0:-2])/bcargsmmagnitude[1:-1]*100.

  STOP

  plots            = PLOT(renu2.time.flight,geosph_arr2[*,2],XTITLE='Time (s)',YTITLE='Alt (km)')
  CGHISTOPLOT,90.-geosph_arr[0,*]
  CGHISTOPLOT,geosph_arr[1,*]
  CGHISTOPLOT,90.-magsph_arr[0,*]
  CGHISTOPLOT,renu2coords.mag.sph.theta
  CGHISTOPLOT,renu2coords.mag.lon
  CGHISTOPLOT,renu2coords.mag.lat
  CGHISTOPLOT,renu2coords.geo.lat
  CGHISTOPLOT,renu2.lat
  CGHISTOPLOT,renu2coords.geo.lon
  CGHISTOPLOT,renu2.long
  CGHISTOPLOT,renu2.ecef.position.x
  CGHISTOPLOT,renu2coords.geo.car.x
  diffposcarx = renu2coords.geo.car.x-renu2.ecef.position.x
  saywhat     = where(ABS(renu2coords.geo.car.x-renu2.ecef.position.x) GT 0.1)
  CGHISTOPLOT,geovel.sph.r
  CGHISTOPLOT,geovel.sph.theta
  CGHISTOPLOT,geovel.sph.r

END


