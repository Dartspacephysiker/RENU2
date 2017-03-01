;2017/02/13
;
;2017/03/01 NOTE: velMagSph and velMag (in GEO coordinates) match to the fourth decimal place (so tenths of a meter?).
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
  nTot             = N_ELEMENTS(renu2.time.utc)

  TiltArr          = !NULL

  GEI_arr          = MAKE_ARRAY(3,nTot,/FLOAT)
  GSM_arr          = MAKE_ARRAY(3,nTot,/FLOAT)
  ;; GEO_arr       = MAKE_ARRAY(3,nTot,/FLOAT)
  IGRF_GEO_arr     = MAKE_ARRAY(3,nTot,/FLOAT)
  IGRF_GSM_arr     = MAKE_ARRAY(3,nTot,/FLOAT)
  IGRF_GSM_DIP_arr      = MAKE_ARRAY(3,nTot,/FLOAT)
  MAG_arr          = MAKE_ARRAY(3,nTot,/FLOAT)

  GEISph_arr       = MAKE_ARRAY(3,nTot,/FLOAT)
  GEOSph_arr       = MAKE_ARRAY(3,nTot,/FLOAT)
  vel_GEOSphArr    = MAKE_ARRAY(3,nTot,/DOUBLE)
  vel_MAGArr       = MAKE_ARRAY(3,nTot,/DOUBLE)
  vel_MAGSphArr    = MAKE_ARRAY(3,nTot,/DOUBLE)
  IGRF_GEO_sphArr  = MAKE_ARRAY(3,nTot,/FLOAT)
  IGRF_GSM_sphArr  = MAKE_ARRAY(3,nTot,/FLOAT)
  MAGSph_arr       = MAKE_ARRAY(3,nTot,/DOUBLE)

  ;;for vector transforms
  ident            = IDENTITY(3,/DOUBLE) 
  GEO2MAG_coord    = MAKE_ARRAY(3,3,nTot,/DOUBLE)
  GEO2MAG_sphC     = MAKE_ARRAY(3,3,nTot,/DOUBLE)
  GEO2MAG_vec      = MAKE_ARRAY(3,3,nTot,/DOUBLE)
  GEO2MAG_sphV     = MAKE_ARRAY(3,3,nTot,/DOUBLE)

  CASE 1 OF
     KEYWORD_SET(use_geopack_08): BEGIN
        FOR i=0,nTot-1 DO BEGIN
           GEOPACK_RECALC_08,YearArr[i],MonthArr[i],DayArr[i],HourArr[i],MinArr[i],SecArr[i],/DATE 

           ;;do that dance
           tmpGEOVel = [renu2.ecef.velocity.x[i],renu2.ecef.velocity.y[i],renu2.ecef.velocity.z[i]]

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
           ;; GEOPACK_SPHCAR_08,renu2.ecef.velocity.x[i],renu2.ecef.velocity.y[i],renu2.ecef.velocity.z[i], $
           ;; GEOPACK_SPHCAR_08,tmpGEOVel[0],tmpGEOVel[1],tmpGEOVel[2], $
           ;;                geoVel_r,geoVel_theta,geoVel_phi,/TO_SPHERE,/DEGREE
           GEOPACK_BCARSP_08,renu2.ecef.position.x[i],renu2.ecef.position.y[i],renu2.ecef.position.z[i], $
                             renu2.ecef.velocity.x[i],renu2.ecef.velocity.y[i],renu2.ecef.velocity.z[i], $
                             geoVel_r,geoVel_theta,geoVel_phi
           GEOPACK_SPHCAR_08,mag_x,mag_y,mag_z,mag_r,mag_theta,mag_phi,/TO_SPHERE,/DEGREE

           ;;Get IGRF
           gsm_xyz_R_E = [gsm_x,gsm_y,gsm_z]/1000.D/R_E ;div by 1000 to get to km, then by R_E (which is in units of km) 
           ;; GEOPACK_RECALC_08,YearArr[i],MonthArr[i],DayArr[i],HourArr[i],MinArr[i],SecArr[i],/DATE 
           GEOPACK_IGRF_GEO_08,geo_r/1000./R_E,geo_theta,geo_phi,br,btheta,bphi,EPOCH=time_epoch[i],/DEGREE
           GEOPACK_BSPCAR_08  ,geo_theta,geo_phi,br,btheta,bphi,bx,by,bz,/DEGREE ; ,EPOCH=time_epoch[i]

           ;;alternate shot at IGRF
           GEOPACK_IGRF_GSW_08,gsm_xyz_R_E[0],gsm_xyz_R_E[1],gsm_xyz_R_E[2],bx_gsm,by_gsm,bz_gsm,EPOCH=time_epoch[i] ;,/DEGREE
           GEOPACK_BCARSP_08  ,gsm_xyz_R_E[0],gsm_xyz_R_E[1],gsm_xyz_R_E[2],bx_gsm,by_gsm,bz_gsm,bmagr,bmagtheta,bmagphi

           ;;Dipole, please?
           GEOPACK_DIP_08,gsm_xyz_R_E[1],gsm_xyz_R_E[1],gsm_xyz_R_E[2],bx_gsm_dip,by_gsm_dip,bz_gsm_dip,EPOCH=time_epoch[i] ;,/DEGREE

           ;;Update spherical
           GEISph_arr[*,i]        = [gei_theta,gei_phi,gei_r] 
           GEOSph_arr[*,i]        = [geo_theta,geo_phi,geo_r] ;Redundant, yes, but a check
           vel_GEOSphArr[*,i]     = [geoVel_r,geoVel_theta,geoVel_phi] ;;NOTE!!!!!
           IGRF_GEO_sphArr[*,i]   = [btheta   ,bphi   ,br   ] 
           IGRF_GSM_sphArr[*,i]   = [bmagtheta,bmagphi,bmagr] 
           MAGSph_arr[*,i]        = [mag_theta,mag_phi,mag_r] 

           ;;Update not-spherical
           ;; TiltArr    = [TiltArr,tempTilt]
           GEI_arr[*,i]  = [gei_x,gei_y,gei_z]
           GSM_arr[*,i]  = [gsm_x,gsm_y,gsm_z]
           ;; GEO_arr[*,i]  = [geo_x,geo_y,geo_z]
           IGRF_GEO_arr[*,i] = [bx,by,bz]
           IGRF_GSM_arr[*,i] = [bx_gsm,by_gsm,bz_gsm]
           IGRF_GSM_DIP_arr[*,i] = [bx_gsm_dip,by_gsm_dip,bz_gsm_dip]
           MAG_arr[*,i]  = [mag_x,mag_y,mag_z]

           ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
           ;;Now transform matrices
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


           GEO2MAG_coord[*,*,i]  = [[mag_x0,mag_y0,mag_z0], $
                                    [mag_x1,mag_y1,mag_z1], $
                                    [mag_x2,mag_y2,mag_z2]]
           GEO2MAG_vec[*,*,i]    = INVERT(GEO2MAG_coord[*,*,i])


           ;;Convert Cartesian vector transform matrices to spherical vector transform matrices

           GEOPACK_BCARSP_08,mag_x0,mag_y0,mag_z0, $
           ;; GEOPACK_BCARSP_08,ident[0,0],ident[0,1],ident[0,2], $
                             ;; GEO2MAG_vec[0,0,i],GEO2MAG_vec[0,1,i],GEO2MAG_vec[0,2,i], $
                             ident[0,0],ident[0,1],ident[0,2], $
                             magVSph_r0,magVSph_theta0,magVSph_phi0
           GEOPACK_BCARSP_08,mag_x1,mag_y1,mag_z1, $
           ;; GEOPACK_BCARSP_08,ident[1,0],ident[1,1],ident[1,2], $
                             ;; GEO2MAG_vec[1,0,i],GEO2MAG_vec[1,1,i],GEO2MAG_vec[1,2,i], $
                             ident[1,0],ident[1,1],ident[1,2], $
                             magVSph_r1,magVSph_theta1,magVSph_phi1
           GEOPACK_BCARSP_08,mag_x2,mag_y2,mag_z2, $
           ;; GEOPACK_BCARSP_08,ident[2,0],ident[2,1],ident[2,2], $
                             ;; GEO2MAG_vec[2,0,i],GEO2MAG_vec[2,1,i],GEO2MAG_vec[2,2,i], $
                             ident[2,0],ident[2,1],ident[2,2], $
                             magVSph_r2,magVSph_theta2,magVSph_phi2

           GEO2MAG_sphV[*,*,i]   = [[magVSph_r0,magVSph_theta0,magVSph_phi0], $
                                    [magVSph_r1,magVSph_theta1,magVSph_phi1], $
                                    [magVSph_r2,magVSph_theta2,magVSph_phi2]]

           ;;Now update spherical coordinate transform matrices
           GEOPACK_SPHCAR_08,mag_x0,mag_y0,mag_z0,mag_r0,mag_theta0,mag_phi0,/TO_SPHERE,/DEGREE
           GEOPACK_SPHCAR_08,mag_x1,mag_y1,mag_z1,mag_r1,mag_theta1,mag_phi1,/TO_SPHERE,/DEGREE
           GEOPACK_SPHCAR_08,mag_x2,mag_y2,mag_z2,mag_r2,mag_theta2,mag_phi2,/TO_SPHERE,/DEGREE

           GEO2MAG_sphC[*,*,i]  = [[mag_r0,mag_theta0,mag_phi0], $
                                   [mag_r1,mag_theta1,mag_phi1], $
                                   [mag_r2,mag_theta2,mag_phi2]]


           ;;get Cartesian and spherical velocity in MAG coordinates
           vel_MAGArr[*,i]       = GEO2MAG_vec[*,*,i] # tmpGEOVel
           vel_MAGSphArr[*,i]    = GEO2MAG_sphV[*,*,i] # vel_GEOSphArr[*,i]

           ;;... And, a TEST. These should come out the same:
           GEOPACK_BSPCAR_08  ,mag_theta,mag_phi, $
                               vel_MAGSphArr[0,i],vel_MAGSphArr[1,i],vel_MAGSphArr[2,i], $
                               testVelx_MAG,testVely_MAG,testVelz_MAG,/DEGREE ; ,EPOCH=time_epoch[i]
           GEOPACK_BCARSP_08  ,gsm_xyz_R_E[0],gsm_xyz_R_E[1],gsm_xyz_R_E[2], $
                               vel_MAGArr[0,i],vel_MAGArr[1,i],vel_MAGArr[2,i], $
                               testVelr_MAG,testVeltheta_MAG,testVelphi_MAG


           velMagCar_GEO    = SQRT(renu2.ecef.velocity.x[i]*renu2.ecef.velocity.x[i] + $
                                   renu2.ecef.velocity.y[i]*renu2.ecef.velocity.y[i] + $
                                   renu2.ecef.velocity.z[i]*renu2.ecef.velocity.z[i])
           velMagSph_GEO    = SQRT(geoVel_r*geoVel_r + $
                                   geoVel_theta*geoVel_theta + $
                                   geoVel_phi*geoVel_phi)

           ;; IF ABS(velMagCar_GEO-velMagSph_GEO) GT 0.0001 THEN STOP ;Precision to the fourth decimal place. This will make it quit.
           IF ABS(velMagCar_GEO-velMagSph_GEO) GT 0.001 THEN STOP
           
           conv_velMagSph_MAG    = SQRT(vel_MAGSphArr[0,i]*vel_MAGSphArr[0,i] + $
                                        vel_MAGSphArr[1,i]*vel_MAGSphArr[1,i] + $
                                        vel_MAGSphArr[2,i]*vel_MAGSphArr[2,i])
           
           conv_velMagCar_MAG    = SQRT(vel_MAGArr[0,i]*vel_MAGArr[0,i] + $
                                        vel_MAGArr[1,i]*vel_MAGArr[1,i] + $
                                        vel_MAGArr[2,i]*vel_MAGArr[2,i])

           IF ABS(velMagCar_GEO-conv_velMagCar_MAG) GT 0.001 THEN STOP
           IF ABS(velMagSph_GEO-conv_velMagSph_MAG) GT 0.001 THEN STOP
           
           convConv_velMagSph_MAG = SQRT(testVelr_MAG*testVelr_MAG + $
                                         testVeltheta_MAG*testVeltheta_MAG + $
                                         testVelphi_MAG*testVelphi_MAG)
           
           convConv_velMagCar_MAG = SQRT(testVelx_MAG*testVelx_MAG + $
                                         testVely_MAG*testVely_MAG + $
                                         testVelz_MAG*testVelz_MAG)
           
           velMagnitudeDiff      = conv_velMagSph_MAG-conv_velMagCar_MAG
           velSphMagnitudeDiff   = conv_velMagSph_MAG-convConv_velMagSph_MAG
           velCarMagnitudeDiff   = conv_velMagCar_MAG-convConv_velMagCar_MAG

           ;; IF (ABS(velCarMagnitudeDiff) GT 1) OR (ABS(velSphMagnitudeDiff) GT 1) OR (ABS(velMagnitudeDiff) GT 1) THEN STOP

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
        ;;    vel_GEOSphArr[*,i] = [geoVel_theta,geoVel_phi,geoVel_r]
        ;;    IGRF_GEO_sphArr[*,i]   = [btheta,bphi,br] 
        ;;    MAGSph_arr[*,i]    = [mag_theta,mag_phi,mag_r] 

        ;;    ;;Update not-spherical
        ;;    ;; TiltArr    = [TiltArr,tempTilt]
        ;;    GEI_arr[*,i]  = [gei_x,gei_y,gei_z]
        ;;    ;; GEO_arr[*,i]  = [geo_x,geo_y,geo_z]
        ;;    IGRF_GEO_arr[*,i] = [bx,by,bz]
        ;;    MAG_arr[*,i]  = [mag_x,mag_y,mag_z]

        ;;    IF (i MOD 100) EQ 0 THEN PRINT,i
        ;; ENDFOR
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
  
  DIPOLE  = {gsm : {car : {x   : REFORM(IGRF_GSM_DIP_arr[0,*]), $
                           y   : REFORM(IGRF_GSM_DIP_arr[1,*]), $
                           z   : REFORM(IGRF_GSM_DIP_arr[2,*])}}}

  IGRF    = {geo : {sph : {bTheta : REFORM(IGRF_GEO_sphArr[0,*]), $
                           bPhi   : REFORM(IGRF_GEO_sphArr[1,*]), $
                           br     : REFORM(IGRF_GEO_sphArr[2,*])}, $
                    car : {x      : REFORM(IGRF_GEO_arr[0,*]), $
                           y      : REFORM(IGRF_GEO_arr[1,*]), $
                           z      : REFORM(IGRF_GEO_arr[2,*])}}, $
             gsm : {sph : {bTheta : REFORM(IGRF_GSM_sphArr[0,*]), $
                           bPhi   : REFORM(IGRF_GSM_sphArr[1,*]), $
                           br     : REFORM(IGRF_GSM_sphArr[2,*])}, $
                    car : {x      : REFORM(IGRF_GSM_arr[0,*]), $
                           y      : REFORM(IGRF_GSM_arr[1,*]), $
                           z      : REFORM(IGRF_GSM_arr[2,*])}}}
  
  GEI     = {ALT:GEISph_arr2[*,2], $
             LON:GEISph_arr2[*,1], $
             LAT:GEISph_arr2[*,0], $
             sph : {theta : REFORM(GEISph_arr[0,*]), $
                    phi   : REFORM(GEISph_arr[1,*]), $
                    r     : REFORM(GEISph_arr[2,*])}, $
             car : {x     : REFORM(GEI_arr[0,*]), $
                    y     : REFORM(GEI_arr[1,*]), $
                    z     : REFORM(GEI_arr[2,*])}}


  ;; vel     = {GEO : {sph : {theta : REFORM(vel_GEOSphArr[0,*]), $
  ;;                          phi   : REFORM(vel_GEOSphArr[1,*]), $
  ;;                          r     : REFORM(vel_GEOSphArr[2,*])}, $
  vel     = {GEO : {sph : {theta : REFORM(vel_GEOSphArr[1,*]), $
                           phi   : REFORM(vel_GEOSphArr[2,*]), $
                           r     : REFORM(vel_GEOSphArr[0,*])}, $
                    car : {x     : renu2.ecef.velocity.x, $
                           y     : renu2.ecef.velocity.y, $
                           z     : renu2.ecef.velocity.z}}, $
             MAG : {sph : {theta : REFORM(vel_MAGSphArr[0,*]), $
                           phi   : REFORM(vel_MAGSphArr[1,*]), $
                           r     : REFORM(vel_MAGSphArr[2,*])}, $
                    car : {x     : REFORM(vel_MAGArr[0,*]), $
                           y     : REFORM(vel_MAGArr[1,*]), $
                           z     : REFORM(vel_MAGArr[2,*])}}}


  pos      = {GEO : {ALT:GEOSph_arr2[*,2], $
                     LON:GEOSph_arr2[*,1], $
                     LAT:GEOSph_arr2[*,0], $
                     sph : {theta : REFORM(GEOSph_arr[0,*]), $
                            phi   : REFORM(GEOSph_arr[1,*]), $
                            r     : REFORM(GEOSph_arr[2,*])}, $
                     car : {x     : renu2.ecef.position.x, $
                            y     : renu2.ecef.position.y, $
                            z     : renu2.ecef.position.z}}, $
              MAG : {ALT:MAGSph_arr2[*,2], $
                     LON:MAGSph_arr2[*,1], $
                     LAT:MAGSph_arr2[*,0], $
                     sph : {theta : REFORM(MAGSph_arr[0,*]), $
                            phi   : REFORM(MAGSph_arr[1,*]), $
                            r     : REFORM(MAGSph_arr[2,*])}, $
                     car : {x     : REFORM(MAG_arr[0,*]), $
                            y     : REFORM(MAG_arr[1,*]), $
                            z     : REFORM(MAG_arr[2,*])}}}


  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;;make struct
  ;; renu2Coords = {GEO     : GEO_arr, $
  ;;                MAG     : MAG_arr, $
  ;;                GEI     : GEI_arr, $
  ;;                IGRF    : IGRF_GEO_arr, $
  ;;                CREATED : GET_TODAY_STRING(/DO_YYYYMMDD_FMT), $
  ;;                ORIGINATING_ROUTINE: orig_routineName}
  coords = {pos     : pos, $
                 vel     : vel, $
                 ;; MAG     : MAG, $
                 IGRF    : IGRF, $
                 DIPOLE  : DIPOLE, $
                 CREATED : GET_TODAY_STRING(/DO_YYYYMMDD_FMT), $
                 ORIGINATING_ROUTINE: orig_routineName}

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;;Save it
  PRINT,'Saving ' + outDir + outFile + '...'
  save,coords,FILENAME=outDir+outFile

  PRINT,"Did it!"

  geocarmagnitude   = SQRT(renu2.ecef.position.x*renu2.ecef.position.x+ $
                           renu2.ecef.position.y*renu2.ecef.position.y+ $
                           renu2.ecef.position.z*renu2.ecef.position.z)
  geocarmagnit2     = SQRT(pos.geo.car.x*pos.geo.car.x+pos.geo.car.y*pos.geo.car.y+pos.geo.car.z*pos.geo.car.z)
  ;; geosphmagnitude  = SQRT(geo.sph.r*geo.sph.r+geo.sph.theta*geo.sph.theta+geo.sph.phi*geo.sph.phi)
  geosphmagnitude   = pos.geo.sph.r
  diffgeo           = geosphmagnitude - $
                      geocarmagnit2

  velcarMagnitude   = SQRT(vel.geo.car.x*vel.geo.car.x+vel.geo.car.y*vel.geo.car.y+vel.geo.car.z*vel.geo.car.z)
  velsphMagnitude  = SQRT(vel.geo.sph.r*vel.geo.sph.r+vel.geo.sph.theta*vel.geo.sph.theta+vel.geo.sph.phi*vel.geo.sph.phi)
  ;; velsphMagnitude   = vel.geo.sph.r ;;pff, idiot
  diffVel           = velsphMagnitude- $
                      velcarmagnitude

  bgeosphmagnitude  = SQRT(igrf.geo.sph.br*igrf.geo.sph.br+igrf.geo.sph.btheta*igrf.geo.sph.btheta+igrf.geo.sph.bphi*igrf.geo.sph.bphi)
  bgeocarmagnitude  = SQRT(igrf.geo.car.x*igrf.geo.car.x+igrf.geo.car.y*igrf.geo.car.y+igrf.geo.car.z*igrf.geo.car.z)
  bgsmcarmagnitude  = SQRT(igrf.gsm.car.x*igrf.gsm.car.x+igrf.gsm.car.y*igrf.gsm.car.y+igrf.gsm.car.z*igrf.gsm.car.z)
  diffB             = bgeosphmagnitude-bgeocarmagnitude

  magPercentDiffGEO = (bgeosphmagnitude[1:-1]-bgeosphmagnitude[0:-2])/bgeosphmagnitude[1:-1]*100.
  magPercentDiffGSM = (bgsmcarmagnitude[1:-1]-bgsmcarmagnitude[0:-2])/bgsmcarmagnitude[1:-1]*100.

  STOP

  plots             = PLOT(renu2.time.flight,geosph_arr2[*,2],XTITLE='Time (s)',YTITLE='Alt (km)')
  CGHISTOPLOT,90.-geosph_arr[0,*]
  CGHISTOPLOT,geosph_arr[1,*]
  CGHISTOPLOT,90.-magsph_arr[0,*]
  CGHISTOPLOT,coords.mag.sph.theta
  CGHISTOPLOT,coords.mag.lon
  CGHISTOPLOT,coords.mag.lat
  CGHISTOPLOT,coords.geo.lat
  CGHISTOPLOT,renu2.lat
  CGHISTOPLOT,coords.geo.lon
  CGHISTOPLOT,renu2.long
  CGHISTOPLOT,renu2.ecef.position.x
  CGHISTOPLOT,coords.geo.car.x
  diffposcarx = coords.geo.car.x-renu2.ecef.position.x
  saywhat     = where(ABS(coords.geo.car.x-renu2.ecef.position.x) GT 0.1)
  CGHISTOPLOT,vel.geo.sph.r
  CGHISTOPLOT,vel.geo.sph.theta
  CGHISTOPLOT,vel.geo.sph.r

END


