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

  pos_GEI_arr      = MAKE_ARRAY(3,nTot,/FLOAT)
  pos_GSM_arr      = MAKE_ARRAY(3,nTot,/FLOAT)
  ;; GEO_arr       = MAKE_ARRAY(3,nTot,/FLOAT)

  pos_MAG_arr      = MAKE_ARRAY(3,nTot,/FLOAT)
  pos_GEISph_arr   = MAKE_ARRAY(3,nTot,/FLOAT)
  pos_GEOSph_arr   = MAKE_ARRAY(3,nTot,/FLOAT)
  pos_GSMSph_arr   = MAKE_ARRAY(3,nTot,/DOUBLE)
  pos_NEDArr       = MAKE_ARRAY(3,nTot,/DOUBLE)

  vel_GEOSphArr    = MAKE_ARRAY(3,nTot,/DOUBLE)
  vel_MAGArr       = MAKE_ARRAY(3,nTot,/DOUBLE)
  vel_MAGSphArr    = MAKE_ARRAY(3,nTot,/DOUBLE)
  vel_NEDArr       = MAKE_ARRAY(3,nTot,/DOUBLE)

  IGRF_GEO_arr     = MAKE_ARRAY(3,nTot,/FLOAT)
  IGRF_GEO_sphArr  = MAKE_ARRAY(3,nTot,/FLOAT)
  IGRF_GSM_arr     = MAKE_ARRAY(3,nTot,/FLOAT)
  IGRF_GSM_sphArr  = MAKE_ARRAY(3,nTot,/FLOAT)
  IGRF_GSM_DIP_arr = MAKE_ARRAY(3,nTot,/FLOAT)

  pos_MAGSph_arr   = MAKE_ARRAY(3,nTot,/DOUBLE)

  ;;for vector transforms
  ident            = IDENTITY(3,/DOUBLE) 

  GEO2GEI_coord    = MAKE_ARRAY(3,3,nTot,/DOUBLE)
  GEO2GEI_vec      = MAKE_ARRAY(3,3,nTot,/DOUBLE)

  GEO2MAG_coord    = MAKE_ARRAY(3,3,nTot,/DOUBLE)
  GEO2MAG_sphC     = MAKE_ARRAY(3,3,nTot,/DOUBLE)
  GEO2MAG_vec      = MAKE_ARRAY(3,3,nTot,/DOUBLE)
  GEO2MAG_sphV     = MAKE_ARRAY(3,3,nTot,/DOUBLE)

  ;;North-East-Down
  GEO2NED_coord    = MAKE_ARRAY(3,3,nTot,/DOUBLE)
  GEO2NED_vec      = MAKE_ARRAY(3,3,nTot,/DOUBLE)

  CASE 1 OF
     KEYWORD_SET(use_geopack_08): BEGIN
        FOR i=0,nTot-1 DO BEGIN
           GEOPACK_RECALC_08,YearArr[i],MonthArr[i],DayArr[i],HourArr[i],MinArr[i],SecArr[i],/DATE 

           ;;do that dance
           tmpPos_GEO              = [renu2.ecef.position.x[i],renu2.ecef.position.y[i],renu2.ecef.position.z[i]]
           tmpVel_GEO              = [renu2.ecef.velocity.x[i],renu2.ecef.velocity.y[i],renu2.ecef.velocity.z[i]]

           ;;To GEI
           GEOPACK_CONV_COORD_08,tmpPos_GEO[0],tmpPos_GEO[1],tmpPos_GEO[2], $
                                 pos_gei_x,pos_gei_y,pos_gei_z, $
                                 /FROM_GEO,/TO_GEI,EPOCH=time_epoch[i]
           ;;To MAG
           GEOPACK_CONV_COORD_08,tmpPos_GEO[0],tmpPos_GEO[1],tmpPos_GEO[2], $
                                 pos_mag_x,pos_mag_y,pos_mag_z, $
                                 /FROM_GEO,/TO_MAG,EPOCH=time_epoch[i]

           ;;To GSM for IGRF calc
           GEOPACK_CONV_COORD_08,tmpPos_GEO[0],tmpPos_GEO[1],tmpPos_GEO[2], $
                                 pos_gsm_x,pos_gsm_y,pos_gsm_z, $
                                 /FROM_GEO,/TO_GSW,EPOCH=time_epoch[i]


           ;;And spherical everything
           GEOPACK_SPHCAR_08,pos_gei_x,pos_gei_y,pos_gei_z, $
                             pos_gei_r,pos_gei_theta,pos_gei_phi,/TO_SPHERE,/DEGREE
           GEOPACK_SPHCAR_08,tmpPos_GEO[0],tmpPos_GEO[1],tmpPos_GEO[2], $
                             pos_geo_r,pos_geo_theta,pos_geo_phi,/TO_SPHERE,/DEGREE
           GEOPACK_SPHCAR_08,pos_gsm_x,pos_gsm_y,pos_gsm_z, $
                             pos_gsm_r,pos_gsm_theta,pos_gsm_phi,/TO_SPHERE,/DEGREE
           GEOPACK_SPHCAR_08,pos_mag_x,pos_mag_y,pos_mag_z, $
                             pos_mag_r,pos_mag_theta,pos_mag_phi,/TO_SPHERE,/DEGREE

           ;;velocity vector in spherical GEO coords
           GEOPACK_BCARSP_08,tmpPos_GEO[0],tmpPos_GEO[1],tmpPos_GEO[2], $
                             tmpVel_GEO[0],tmpVel_GEO[1],tmpVel_GEO[2], $
                             vel_GEO_r,vel_GEO_theta,vel_GEO_phi

           ;;Get IGRF
           pos_gsm_xyz_R_E         = [pos_gsm_x,pos_gsm_y,pos_gsm_z]/1000.D/R_E ;div by 1000 to get to km, then by R_E (which is in units of km) 
           ;; GEOPACK_RECALC_08,YearArr[i],MonthArr[i],DayArr[i],HourArr[i],MinArr[i],SecArr[i],/DATE 
           GEOPACK_IGRF_GEO_08,pos_geo_r/1000./R_E,pos_geo_theta,pos_geo_phi,br,btheta,bphi,EPOCH=time_epoch[i],/DEGREE
           GEOPACK_BSPCAR_08  ,pos_geo_theta,pos_geo_phi,br,btheta,bphi,bx,by,bz,/DEGREE ; ,EPOCH=time_epoch[i]

           ;;alternate shot at IGRF
           GEOPACK_IGRF_GSW_08,pos_gsm_xyz_R_E[0],pos_gsm_xyz_R_E[1],pos_gsm_xyz_R_E[2],bx_gsm,by_gsm,bz_gsm,EPOCH=time_epoch[i] ;,/DEGREE
           GEOPACK_BCARSP_08  ,pos_gsm_xyz_R_E[0],pos_gsm_xyz_R_E[1],pos_gsm_xyz_R_E[2],bx_gsm,by_gsm,bz_gsm,bmagr,bmagtheta,bmagphi

           ;;Dipole, please?
           GEOPACK_DIP_08,pos_gsm_xyz_R_E[1],pos_gsm_xyz_R_E[1],pos_gsm_xyz_R_E[2],bx_gsm_dip,by_gsm_dip,bz_gsm_dip,EPOCH=time_epoch[i] ;,/DEGREE

           ;;Update not-spherical
           ;; TiltArr              = [TiltArr,tempTilt]
           pos_GEI_arr[*,i]        = [pos_gei_x,pos_gei_y,pos_gei_z]
           pos_GSM_arr[*,i]        = [pos_gsm_x,pos_gsm_y,pos_gsm_z]
           pos_MAG_arr[*,i]        = [pos_mag_x,pos_mag_y,pos_mag_z]

           ;;Update spherical
           pos_GEISph_arr[*,i]     = [pos_gei_theta,pos_gei_phi,pos_gei_r] 
           pos_GEOSph_arr[*,i]     = [pos_geo_theta,pos_geo_phi,pos_geo_r] ;Redundant, yes, but a check
           pos_GSMSph_arr[*,i]     = [pos_gsm_theta,pos_gsm_phi,pos_gsm_r] ;Redundant, yes, but a check
           pos_MAGSph_arr[*,i]     = [pos_mag_theta,pos_mag_phi,pos_mag_r] 

           vel_GEOSphArr[*,i]      = [vel_GEO_r,vel_GEO_theta,vel_GEO_phi]    ;;NOTE!!!!!

           IGRF_GEO_sphArr[*,i]    = [btheta   ,bphi   ,br   ] 
           IGRF_GSM_sphArr[*,i]    = [bmagtheta,bmagphi,bmagr] 

           ;;Update IGRFs
           IGRF_GSM_DIP_arr[*,i]   = [bx_gsm_dip,by_gsm_dip,bz_gsm_dip]
           IGRF_GEO_arr[*,i]       = [bx,by,bz]
           IGRF_GSM_arr[*,i]       = [bx_gsm,by_gsm,bz_gsm]

           ;;Where is magnetic center of earth?
           GEOPACK_CONV_COORD_08,0D,0D,0D, $
                                 pos_magCenter_x_GEI,pos_magCenter_y_GEI,pos_magCenter_z_GEI, $
                                 /FROM_MAG,/TO_GEI,EPOCH=time_epoch[i]
           pos_magCenter_GEI       = [TEMPORARY(pos_magCenter_x_GEI),TEMPORARY(pos_magCenter_y_GEI),TEMPORARY(pos_magCenter_z_GEI)]

           ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
           ;;Now transform matrices

           ;;GEO to GEI
           GEOPACK_CONV_COORD_08,ident[0,0],ident[0,1],ident[0,2], $
                                 gei_x0,gei_y0,gei_z0, $
                                 /FROM_GEO,/TO_GEI,EPOCH=time_epoch[i]

           GEOPACK_CONV_COORD_08,ident[1,0],ident[1,1],ident[1,2], $
                                 gei_x1,gei_y1,gei_z1, $
                                 /FROM_GEO,/TO_GEI,EPOCH=time_epoch[i]

           GEOPACK_CONV_COORD_08,ident[2,0],ident[2,1],ident[2,2], $
                                 gei_x2,gei_y2,gei_z2, $
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


           GEO2GEI_coord[*,*,i]    = [[gei_x0,gei_y0,gei_z0], $
                                      [gei_x1,gei_y1,gei_z1], $
                                      [gei_x2,gei_y2,gei_z2]]
           GEO2GEI_vec[*,*,i]      = INVERT(GEO2GEI_coord[*,*,i])

           GEO2MAG_coord[*,*,i]    = [[mag_x0,mag_y0,mag_z0], $
                                      [mag_x1,mag_y1,mag_z1], $
                                      [mag_x2,mag_y2,mag_z2]]
           GEO2MAG_vec[*,*,i]      = INVERT(GEO2MAG_coord[*,*,i])


           ;;North-East-Down
           GEOPACK_CONV_COORD_08,ident[2,0],ident[2,1],ident[2,2], $
                                 mag_x2,mag_y2,mag_z2, $
                                 /FROM_GEO,/TO_MAG,EPOCH=time_epoch[i]

           tmpLambda               = renu2.long[i]*!DTOR
           tmpPhi                  = renu2.lat[i]*!DTOR


           ;; GEOPACK_GEODGEO_08,0.D,tmpPhi,refAlt_GEO,refTheta_GEO,(KEYWORD_SET(GEOD_TO_GEO) ? 1 : -1) ;see GEOPACK_2008 documentation
           GEOPACK_GEODGEO_08,0.D,tmpPhi,refAlt_GEO,refTheta_GEO,/TO_GEOCENTRIC ;/TO_GEODETIC (see GEOPACK_2008 documentation, or type GEOPACK_HELP)
           GEOPACK_SPHCAR_08,refAlt_GEO,refTheta_GEO,tmpLambda,refAlt_x,refAlt_y,refAlt_z,/TO_RECT

           tmpGEODRef              = [refAlt_x,refAlt_y,refAlt_z]
           tmpNEDPos               = [tmpPos_GEO[0]-tmpGEODRef[0],tmpPos_GEO[1]-tmpGEODRef[1],tmpPos_GEO[2]-tmpGEODRef[2]]

           GEO2NED_coord[*,*,i]    = [[-SIN(tmpPhi)*COS(tmpLambda),-SIN(tmpPhi)*SIN(tmpLambda),COS(tmpPhi)], $
                                      [-SIN(tmpLambda)            ,COS(tmpLambda)             ,0.D         ], $
                                      [-COS(tmpPhi)*COS(tmpLambda),-COS(tmpPhi)*SIN(tmpLambda),-SIN(tmpPhi)]]
           GEO2NED_vec[*,*,i]      = INVERT(GEO2NED_coord[*,*,i])

           pos_NEDArr[*,i]         = GEO2NED_coord[*,*,i] # tmpNEDPos
           vel_NEDArr[*,i]         = GEO2NED_vec[*,*,i]   # tmpVel_GEO

           ;;Convert Cartesian vector transform matrices to spherical vector transform matrices
           ;; GEOPACK_BCARSP_08,mag_x0,mag_y0,mag_z0, $
           GEOPACK_BCARSP_08,pos_mag_x,pos_mag_y,pos_mag_z, $
                             ;; GEOPACK_BCARSP_08,ident[0,0],ident[0,1],ident[0,2], $
                             ;; GEO2MAG_vec[0,0,i],GEO2MAG_vec[0,1,i],GEO2MAG_vec[0,2,i], $
                             ident[0,0],ident[0,1],ident[0,2], $
                             magVSph_r0,magVSph_theta0,magVSph_phi0
           ;; GEOPACK_BCARSP_08,mag_x1,mag_y1,mag_z1, $
           GEOPACK_BCARSP_08,pos_mag_x,pos_mag_y,pos_mag_z, $
                             ;; GEOPACK_BCARSP_08,ident[1,0],ident[1,1],ident[1,2], $
                             ;; GEO2MAG_vec[1,0,i],GEO2MAG_vec[1,1,i],GEO2MAG_vec[1,2,i], $
                             ident[1,0],ident[1,1],ident[1,2], $
                             magVSph_r1,magVSph_theta1,magVSph_phi1
           ;; GEOPACK_BCARSP_08,mag_x2,mag_y2,mag_z2, $
           GEOPACK_BCARSP_08,pos_mag_x,pos_mag_y,pos_mag_z, $
                             ;; GEOPACK_BCARSP_08,ident[2,0],ident[2,1],ident[2,2], $
                             ;; GEO2MAG_vec[2,0,i],GEO2MAG_vec[2,1,i],GEO2MAG_vec[2,2,i], $
                             ident[2,0],ident[2,1],ident[2,2], $
                             magVSph_r2,magVSph_theta2,magVSph_phi2

           GEO2MAG_sphV[*,*,i]     = [[magVSph_r0,magVSph_theta0,magVSph_phi0], $
                                      [magVSph_r1,magVSph_theta1,magVSph_phi1], $
                                      [magVSph_r2,magVSph_theta2,magVSph_phi2]]
           
           

           ;;Now update spherical coordinate transform matrices
           GEOPACK_SPHCAR_08,mag_x0,mag_y0,mag_z0,mag_r0,mag_theta0,mag_phi0,/TO_SPHERE,/DEGREE
           GEOPACK_SPHCAR_08,mag_x1,mag_y1,mag_z1,mag_r1,mag_theta1,mag_phi1,/TO_SPHERE,/DEGREE
           GEOPACK_SPHCAR_08,mag_x2,mag_y2,mag_z2,mag_r2,mag_theta2,mag_phi2,/TO_SPHERE,/DEGREE

           GEO2MAG_sphC[*,*,i]     = [[mag_r0,mag_theta0,mag_phi0], $
                                      [mag_r1,mag_theta1,mag_phi1], $
                                      [mag_r2,mag_theta2,mag_phi2]]


           ;;get Cartesian and spherical velocity in MAG coordinates
           vel_MAGArr[*,i]         = GEO2MAG_vec[*,*,i] # tmpVel_GEO
           ;; vel_MAGSphArr[*,i]   = GEO2MAG_sphV[*,*,i] # vel_GEOSphArr[*,i]
           ;; vel_MAGSphArr[*,i]   = GEO2MAG_sphV[*,*,i] # vel_MAGArr[*,i]

           
           GEOPACK_BCARSP_08,pos_mag_x,pos_mag_y,pos_mag_z, $
                             vel_MAGArr[0,i],vel_MAGArr[1,i],vel_MAGArr[2,i], $
                             magVSph_r,magVSph_theta,magVSph_phi
           vel_MAGSphArr[*,i]      = [magVSph_r,magVSph_theta,magVSph_phi]

           ;;... And, a TEST. These should come out the same:
           GEOPACK_BSPCAR_08  ,pos_mag_theta,pos_mag_phi, $
                               vel_MAGSphArr[0,i],vel_MAGSphArr[1,i],vel_MAGSphArr[2,i], $
                               testVelx_MAG,testVely_MAG,testVelz_MAG,/DEGREE ; ,EPOCH=time_epoch[i]
           GEOPACK_BCARSP_08  ,pos_gsm_xyz_R_E[0],pos_gsm_xyz_R_E[1],pos_gsm_xyz_R_E[2], $
                               vel_MAGArr[0,i],vel_MAGArr[1,i],vel_MAGArr[2,i], $
                               testVelr_MAG,testVeltheta_MAG,testVelphi_MAG


           velMagCar_GEO           = SQRT(tmpVel_GEO[0]*tmpVel_GEO[0] + $
                                          tmpVel_GEO[1]*tmpVel_GEO[1] + $
                                          tmpVel_GEO[2]*tmpVel_GEO[2])
           velMagSph_GEO           = SQRT(vel_GEO_r*vel_GEO_r + $
                                          vel_GEO_theta*vel_GEO_theta + $
                                          vel_GEO_phi*vel_GEO_phi)

           ;; IF ABS(velMagCar_GEO-velMagSph_GEO) GT 0.0001 THEN STOP ;Precision to the fourth decimal place. This will make it quit.
           IF ABS(velMagCar_GEO-velMagSph_GEO) GT 0.001 THEN STOP
           
           conv_velMagSph_MAG      = SQRT(vel_MAGSphArr[0,i]*vel_MAGSphArr[0,i] + $
                                          vel_MAGSphArr[1,i]*vel_MAGSphArr[1,i] + $
                                          vel_MAGSphArr[2,i]*vel_MAGSphArr[2,i])
           
           conv_velMagCar_MAG      = SQRT(vel_MAGArr[0,i]*vel_MAGArr[0,i] + $
                                          vel_MAGArr[1,i]*vel_MAGArr[1,i] + $
                                          vel_MAGArr[2,i]*vel_MAGArr[2,i])

           IF ABS(velMagCar_GEO-conv_velMagCar_MAG) GT 0.001 THEN STOP
           IF ABS(velMagSph_GEO-conv_velMagSph_MAG) GT 0.001 THEN STOP
           
           convConv_velMagSph_MAG  = SQRT(testVelr_MAG*testVelr_MAG + $
                                          testVeltheta_MAG*testVeltheta_MAG + $
                                          testVelphi_MAG*testVelphi_MAG)
           
           convConv_velMagCar_MAG  = SQRT(testVelx_MAG*testVelx_MAG + $
                                          testVely_MAG*testVely_MAG + $
                                          testVelz_MAG*testVelz_MAG)
           
           velMagnitudeDiff        = conv_velMagSph_MAG-conv_velMagCar_MAG
           velSphMagnitudeDiff     = conv_velMagSph_MAG-convConv_velMagSph_MAG
           velCarMagnitudeDiff     = conv_velMagCar_MAG-convConv_velMagCar_MAG

           IF (ABS(velCarMagnitudeDiff) GT 1) OR (ABS(velSphMagnitudeDiff) GT 1) OR (ABS(velMagnitudeDiff) GT 1) THEN STOP

           IF (i MOD 100) EQ 0 THEN PRINT,i
        ENDFOR
     END
     ELSE: BEGIN
        ;; FOR i=0,nTot-1 DO BEGIN
        ;;    GEOPACK_RECALC,YearArr[i],MonthArr[i],DayArr[i],HourArr[i],MinArr[i],SecArr[i],/DATE 

        ;;    ;;do that dance

        ;;    ;;To GEI
        ;;    GEOPACK_CONV_COORD,tmpPos_GEO[0],tmpPos_GEO[1],tmpPos_GEO[2], $
        ;;                       pos_gei_x,pos_gei_y,pos_gei_z, $
        ;;                       /FROM_GEO,/TO_GEI,EPOCH=time_epoch[i]
        ;;    ;;To MAG
        ;;    GEOPACK_CONV_COORD,tmpPos_GEO[0],tmpPos_GEO[1],tmpPos_GEO[2], $
        ;;                       pos_mag_x,pos_mag_y,pos_mag_z, $
        ;;                       /FROM_GEO,/TO_MAG,EPOCH=time_epoch[i]


        ;;    ;;And spherical everything
        ;;    GEOPACK_SPHCAR,pos_gei_x,pos_gei_y,pos_gei_z,pos_gei_r,pos_gei_theta,pos_gei_phi,/TO_SPHERE,/DEGREE
        ;;    GEOPACK_SPHCAR,tmpPos_GEO[0],tmpPos_GEO[1],tmpPos_GEO[2], $
        ;;                   pos_geo_r,pos_geo_theta,pos_geo_phi,/TO_SPHERE,/DEGREE
        ;;    GEOPACK_SPHCAR,tmpVel_GEO[0],renu2.ecef.velocity.y[i],tmpVel_GEO[2], $
        ;;                   vel_GEO_r,vel_GEO_theta,vel_GEO_phi,/TO_SPHERE,/DEGREE
        ;;    GEOPACK_SPHCAR,pos_mag_x,pos_mag_y,pos_mag_z,pos_mag_r,pos_mag_theta,pos_mag_phi,/TO_SPHERE,/DEGREE

        ;;    ;;Get IGRF
        ;;    GEOPACK_IGRF_GEO,pos_geo_r/R_E/1000.D,pos_geo_theta,pos_geo_phi,br,btheta,bphi,/DEGREE,EPOCH=time_epoch[i]
        ;;    GEOPACK_BSPCAR  ,pos_geo_theta,pos_geo_phi,br,btheta,bphi,bx,by,bz,/DEGREE ; ,EPOCH=time_epoch[i]

        ;;    ;;Update spherical
        ;;    pos_GEISph_arr[*,i]    = [pos_gei_theta,pos_gei_phi,pos_gei_r] 
        ;;    pos_GEOSph_arr[*,i]    = [pos_geo_theta,pos_geo_phi,pos_geo_r] ;Redundant, yes, but a check
        ;;    vel_GEOSphArr[*,i]     = [vel_GEO_theta,vel_GEO_phi,vel_GEO_r]
        ;;    IGRF_GEO_sphArr[*,i]   = [btheta,bphi,br] 
        ;;    pos_MAGSph_arr[*,i]    = [pos_mag_theta,pos_mag_phi,pos_mag_r] 

        ;;    ;;Update not-spherical
        ;;    ;; TiltArr    = [TiltArr,tempTilt]
        ;;    pos_GEI_arr[*,i]  = [pos_gei_x,pos_gei_y,pos_gei_z]
        ;;    ;; GEO_arr[*,i]  = [geo_x,geo_y,geo_z]
        ;;    IGRF_GEO_arr[*,i] = [bx,by,bz]
        ;;    pos_MAG_arr[*,i]  = [pos_mag_x,pos_mag_y,pos_mag_z]

        ;;    IF (i MOD 100) EQ 0 THEN PRINT,i
        ;; ENDFOR
     END
  ENDCASE


  ;;Lat, long, height
  GEISph_arr2   = [ $
                  [90.-REFORM(pos_GEISph_arr[0,*])], $
                  [REFORM(pos_GEISph_arr[1,*])], $
                  [REFORM(pos_GEISph_arr[2,*])-R_E] $ ;Convert to latitude from colatitude here
                  ]

  GEOSph_arr2   = [ $
                  [90.-REFORM(pos_GEOSph_arr[0,*])], $
                  [REFORM(pos_GEOSph_arr[1,*])], $
                  [REFORM(pos_GEOSph_arr[2,*])-R_E] $ ;Convert to latitude from colatitude here
                  ]   

  MAGSph_arr2   = [ $
                  [90.-REFORM(pos_MAGSph_arr[0,*])], $
                  [REFORM(pos_MAGSph_arr[1,*])], $
                  [REFORM(pos_MAGSph_arr[2,*])-R_E] $ ;Convert to latitude from colatitude here
                  ]   
  
  DIPOLE  = {gsm : {car : {x   : REFORM(IGRF_GSM_DIP_arr[0,*]), $
                           y   : REFORM(IGRF_GSM_DIP_arr[1,*]), $
                           z   : REFORM(IGRF_GSM_DIP_arr[2,*])}}}

  IGRF    = {GEO : {car : {x      : REFORM(IGRF_GEO_arr[0,*]), $
                           y      : REFORM(IGRF_GEO_arr[1,*]), $
                           z      : REFORM(IGRF_GEO_arr[2,*])}, $
                    sph : {bTheta : REFORM(IGRF_GEO_sphArr[0,*]), $
                           bPhi   : REFORM(IGRF_GEO_sphArr[1,*]), $
                           br     : REFORM(IGRF_GEO_sphArr[2,*])}}, $
             GSM : {car : {x      : REFORM(IGRF_GSM_arr[0,*]), $
                           y      : REFORM(IGRF_GSM_arr[1,*]), $
                           z      : REFORM(IGRF_GSM_arr[2,*])}, $
                    sph : {bTheta : REFORM(IGRF_GSM_sphArr[0,*]), $
                           bPhi   : REFORM(IGRF_GSM_sphArr[1,*]), $
                           br     : REFORM(IGRF_GSM_sphArr[2,*])}}}
  
  ;;Position struct
  pos      = {GEI : {car : {x     : REFORM(pos_GEI_arr[0,*]), $
                            y     : REFORM(pos_GEI_arr[1,*]), $
                            z     : REFORM(pos_GEI_arr[2,*])}, $
                     sph : {theta : REFORM(pos_GEISph_arr[0,*]), $
                            phi   : REFORM(pos_GEISph_arr[1,*]), $
                            r     : REFORM(pos_GEISph_arr[2,*])}}, $
              GEO : {ALT : renu2.alt, $
                     LON : renu2.long, $
                     LAT : renu2.lat, $
                     car : {x     : renu2.ecef.position.x, $
                            y     : renu2.ecef.position.y, $
                            z     : renu2.ecef.position.z}, $
                     sph : {theta : REFORM(pos_GEOSph_arr[0,*]), $
                            phi   : REFORM(pos_GEOSph_arr[1,*]), $
                            r     : REFORM(pos_GEOSph_arr[2,*])}}, $
              GSM : {car : {x     : REFORM(pos_GSM_arr[0,*]), $
                            y     : REFORM(pos_GSM_arr[1,*]), $
                            z     : REFORM(pos_GSM_arr[2,*])}}, $
              MAG : {ALT : MAGSph_arr2[*,2], $
                     LON : MAGSph_arr2[*,1], $
                     LAT : MAGSph_arr2[*,0], $
                     car : {x     : REFORM(pos_MAG_arr[0,*]), $
                            y     : REFORM(pos_MAG_arr[1,*]), $
                            z     : REFORM(pos_MAG_arr[2,*])}, $
                     sph : {theta : REFORM(pos_MAGSph_arr[0,*]), $
                            phi   : REFORM(pos_MAGSph_arr[1,*]), $
                            r     : REFORM(pos_MAGSph_arr[2,*])}}}

  ;;Velocity struct
  vel     = {GEO : {car : {x     : renu2.ecef.velocity.x, $
                           y     : renu2.ecef.velocity.y, $
                           z     : renu2.ecef.velocity.z}, $
                    sph : {theta : REFORM(vel_GEOSphArr[1,*]), $
                           phi   : REFORM(vel_GEOSphArr[2,*]), $
                           r     : REFORM(vel_GEOSphArr[0,*])}}, $
             MAG : {car : {x     : REFORM(vel_MAGArr[0,*]), $
                           y     : REFORM(vel_MAGArr[1,*]), $
                           z     : REFORM(vel_MAGArr[2,*])}, $
                    sph : {theta : REFORM(vel_MAGSphArr[0,*]), $
                           phi   : REFORM(vel_MAGSphArr[1,*]), $
                           r     : REFORM(vel_MAGSphArr[2,*])}}}


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

  plots             = PLOT(renu2.time.flight,GEOSph_arr2[*,2],XTITLE='Time (s)',YTITLE='Alt (km)')
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


