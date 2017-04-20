;2017/04/20
FUNCTION DOTP,v1,v2
  RETURN,(TRANSPOSE(v1) # v2)[0]
END
FUNCTION VECNORM,vec
  RETURN,(SQRT(TRANSPOSE(vec) # vec))[0]
END
FUNCTION VNORMALIZE,vec
  ;; RETURN,[vec[0],vec[1],vec[2]]/SQRT(vec[0]*vec[0]+vec[1]*vec[1]+vec[2]*vec[2])
  RETURN,([vec[0],vec[1],vec[2]]/VECNORM(vec))
END
FUNCTION CROSSP_NORMED,v1,v2
  tmp = CROSSP(v1,v2)
  RETURN,VNORMALIZE(tmp)
END

PRO JOURNAL__20170420__ROTVEC_MAG_TO_GEO_FOR_ROB

  COMPILE_OPT IDL2,STRICTARRSUBS

  ;; long = 267.0D
  ;; lat  = 70.0D

  ;;For practice 
  long_MAG = 0.D
  lat_MAG  = 0.D


                                ;rad, colat
  pos_sph_MAG = [ 1, ((90.D - lat_MAG) * !DTOR) , long_MAG * !DTOR ]

  timeStr    = '2017-03-02/13:00:00'
  date       = TIME_TO_STR(timeStr,/MS)
  
  time_epoch = UTC_TO_CDF_EPOCH(date)

  CONVERT_TIME_STRING_TO_YMDHMS_ARRAYS,timeStr, $
                                       OUT_YEARARR=yearArr, $
                                       OUT_DOYARR=DOYArr, $
                                       OUT_MONTHARR=monthArr, $
                                       OUT_DAYARR=dayArr, $
                                       OUT_HOURARR=hourArr, $
                                       OUT_MINARR=minArr, $
                                       OUT_SECARR=secArr

  nTot             = N_ELEMENTS(timeStr)

  TiltArr          = !NULL

  ;; pos_GEI_arr      = MAKE_ARRAY(3,nTot,/FLOAT)
  ;; pos_GSM_arr      = MAKE_ARRAY(3,nTot,/FLOAT)
  ;; pos_GEISph_arr   = MAKE_ARRAY(3,nTot,/FLOAT)
  pos_GEOSph_arr   = MAKE_ARRAY(3,nTot,/FLOAT)
  pos_GEO_arr      = MAKE_ARRAY(3,nTot,/FLOAT)
  ;; pos_GSMSph_arr   = MAKE_ARRAY(3,nTot,/DOUBLE)
  pos_MAG_arr      = MAKE_ARRAY(3,nTot,/FLOAT)
  pos_MAGSph_arr   = MAKE_ARRAY(3,nTot,/DOUBLE)
  ;; pos_NED_arr      = MAKE_ARRAY(3,nTot,/DOUBLE)
  
  pos_MAGSph_arr[*,0] = pos_sph_MAG

  ;; vel_FAC_arr      = MAKE_ARRAY(3,nTot,/DOUBLE)
  ;; vel_FACV_arr     = MAKE_ARRAY(3,nTot,/DOUBLE)
  ;; vel_GEI_arr      = MAKE_ARRAY(3,nTot,/DOUBLE)
  ;; vel_GEOSph_arr   = MAKE_ARRAY(3,nTot,/DOUBLE)
  ;; vel_MAG_arr      = MAKE_ARRAY(3,nTot,/DOUBLE)
  ;; vel_MAGSph_arr   = MAKE_ARRAY(3,nTot,/DOUBLE)
  ;; vel_NED_arr      = MAKE_ARRAY(3,nTot,/DOUBLE)
  ;; vel_VDH_arr      = MAKE_ARRAY(3,nTot,/DOUBLE)

  ;; IGRF_FAC_arr     = MAKE_ARRAY(3,nTot,/FLOAT)
  ;; IGRF_GEI_arr     = MAKE_ARRAY(3,nTot,/FLOAT)
  ;; IGRF_GEO_arr     = MAKE_ARRAY(3,nTot,/FLOAT)
  ;; IGRF_GEO_sphArr  = MAKE_ARRAY(3,nTot,/FLOAT)
  ;; IGRF_GSM_arr     = MAKE_ARRAY(3,nTot,/FLOAT)
  ;; IGRF_GSM_sphArr  = MAKE_ARRAY(3,nTot,/FLOAT)
  ;; IGRF_GSM_DIP_arr = MAKE_ARRAY(3,nTot,/FLOAT)
  ;; IGRF_VDH_arr     = MAKE_ARRAY(3,nTot,/FLOAT)

  ;;for vector transforms
  ident            = IDENTITY(3,/DOUBLE)

  ;; GEO2FAC_vec      = MAKE_ARRAY(3,3,nTot,/DOUBLE)
  ;; GEO2FACV_vec     = MAKE_ARRAY(3,3,nTot,/DOUBLE)

  ;; GEO2GEI_coord    = MAKE_ARRAY(3,3,nTot,/DOUBLE)
  ;; GEO2GEI_vec      = MAKE_ARRAY(3,3,nTot,/DOUBLE)

  ;; GEO2MAG_coord    = MAKE_ARRAY(3,3,nTot,/DOUBLE)
  ;; GEO2MAG_vec      = MAKE_ARRAY(3,3,nTot,/DOUBLE)

  MAG2GEO_coord    = MAKE_ARRAY(3,3,nTot,/DOUBLE)
  MAG2GEO_vec      = MAKE_ARRAY(3,3,nTot,/DOUBLE)

  ;;North-East-Down
  ;; GEO2NED_coord    = MAKE_ARRAY(3,3,nTot,/DOUBLE)
  ;; GEO2NED_vec      = MAKE_ARRAY(3,3,nTot,/DOUBLE)

  ;; GEI2VDH_vec      = MAKE_ARRAY(3,3,nTot,/DOUBLE)

  FOR i=0,nTot-1 DO BEGIN

     GEOPACK_RECALC_08,YearArr[i],MonthArr[i],DayArr[i],HourArr[i],MinArr[i],SecArr[i],/DATE

     ;;do that dance
     tmpPos_MAGSph              = pos_MAGSph_arr[*,i]

     ;; tmpPos_GEO              = DOUBLE([renu2.ecef.position.x[i],renu2.ecef.position.y[i],renu2.ecef.position.z[i]])
     ;; tmpVel_GEO              = DOUBLE([renu2.ecef.velocity.x[i],renu2.ecef.velocity.y[i],renu2.ecef.velocity.z[i]])

     ;;Convert position in MAG spherical to position in MAG Cartesian 
     GEOPACK_SPHCAR_08,tmpPos_MagSph[0],tmpPos_MagSph[1],tmpPos_MagSph[2], $
                       tmpPos_MAG_x,tmpPos_MAG_y,tmpPos_MAG_z, $
                       /TO_RECT

     tmpPos_MAG                  = [tmpPos_MAG_x,tmpPos_MAG_y,tmpPos_MAG_z]

     ;;Now get (LOCAL) MAGSph basis vectors in MAG Cartesian basis
     tmp_MAGSph2MAG              = ident * 0.D
     FOR k=0,2 DO BEGIN
        ;;k indexes 
        GEOPACK_BSPCAR_08,tmpPos_MagSph[1],tmpPos_MagSph[2], $
                          ident[k,0],ident[k,1],ident[k,2], $
                          tmpX,tmpY,tmpZ

        tmp_MAGSph2MAG[k,*] = [TEMPORARY(tmpX),TEMPORARY(tmpY),TEMPORARY(tmpZ)]

     ENDFOR

     ;;rHat, thetaHat, and phiHat in MAG Cartesian coordinates
     ;;Note, -thetaHat is our NORTH_MAG vector
     rHat_MAG = REFORM(tmp_MAGSph2MAG[0,*])
     tHat_MAG = REFORM(tmp_MAGSph2MAG[1,*])
     pHat_MAG = REFORM(tmp_MAGSph2MAG[2,*])

     ;;So this is the (LOCAL) NORTH_MAG vector
     MAGNorthHat_MAG = -tHat_MAG

     ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
     ;;Now let's convert our position in MAG_Cartesian to GEO_Cartesian
     GEOPACK_CONV_COORD_08,tmpPos_MAG[0],tmpPos_MAG[1],tmpPos_MAG[2], $
                           tmpPos_GEO_x,tmpPos_GEO_y,tmpPos_GEO_z, $
                           /FROM_MAG,/TO_GEO,EPOCH=time_epoch[i]

     tmpPos_GEO   = [TEMPORARY(tmpPos_GEO_x),TEMPORARY(tmpPos_GEO_y),TEMPORARY(tmpPos_GEO_z)]

     ;;Use GEO_Cartesian position to get GEO_GEODETIC, from which we'll

     ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
     ;;MAG to GEO (coords) and MAG to GEO (vec)
     GEOPACK_CONV_COORD_08,ident[0,0],ident[0,1],ident[0,2], $
                           geo_x0,geo_y0,geo_z0, $
                           /FROM_MAG,/TO_GEO,EPOCH=time_epoch[i]

     GEOPACK_CONV_COORD_08,ident[1,0],ident[1,1],ident[1,2], $
                           geo_x1,geo_y1,geo_z1, $
                           /FROM_MAG,/TO_GEO,EPOCH=time_epoch[i]

     GEOPACK_CONV_COORD_08,ident[2,0],ident[2,1],ident[2,2], $
                           geo_x2,geo_y2,geo_z2, $
                           /FROM_MAG,/TO_GEO,EPOCH=time_epoch[i]


     MAG2GEO_coord[*,*,i] = [[geo_x0,geo_y0,geo_z0], $
                             [geo_x1,geo_y1,geo_z1], $
                             [geo_x2,geo_y2,geo_z2]]

     MAG2GEO_vec[*,*,i]   = INVERT(MAG2GEO_coord[*,*,i])


     ;;This and the former tmpPos_GEO should give the same answer
     tmpPos_GEO2          = MAG2GEO_coord[*,*,i] # tmpPos_MAG

     testDiff = tmpPos_GEO2 - tmpPos_GEO
     PRINT,testDiff
     IF (WHERE(ABS(testDiff) GT 0.005))[0] NE -1 THEN STOP

     ;;Now get MAGNorthHat_MAG in GEO coords
     MAGNorthHat_GEO      = MAG2GEO_vec[*,*,i] # MAGNorthHat_MAG

     ;; PRINT,DOTP(

     ;;Position in GEO to GEO_sph
     GEOPACK_SPHCAR_08,tmpPos_GEO[0],tmpPos_GEO[1],tmpPos_GEO[2], $
                       tmpPos_GEO_r,tmpPos_GEO_theta,tmpPos_GEO_phi, $
                       /TO_SPHERE

     tmpPos_GEOSph = [tmpPos_GEO_r,tmpPos_GEO_theta,tmpPos_GEO_phi]

     ;;Now from GEOCENTRIC to GEODETIC
     ;;Wait, no. Just use spherical
     ;; GEOPACK_GEODGEO_08,tmpPos_GEOSph[0],tmpPos_GEOSph[1], $
     ;;                    refAlt_GEO,refTheta_GEO,/TO_GEOCENTRIC ;/TO_GEODETIC (see GEOPACK_2008 documentation, or type GEOPACK_HELP)
     
     ;;Get local GEONorthHat in GEO coordinates
     ;; GEONorthHat_
     ;; GEOPACK_BSPCAR_08,tmpPos_GEOSph[1],tmpPos_GEOSph[2], $
                       

     ;;Now get (LOCAL) MAGSph basis vectors in MAG Cartesian basis
     tmp_GEOSph2GEO              = ident * 0.D
     FOR k=0,2 DO BEGIN
        ;;k indexes 
        GEOPACK_BSPCAR_08,tmpPos_GEOSph[1],tmpPos_GEOSph[2], $
                          ident[k,0],ident[k,1],ident[k,2], $
                          tmpX,tmpY,tmpZ

        tmp_GEOSph2GEO[k,*] = [TEMPORARY(tmpX),TEMPORARY(tmpY),TEMPORARY(tmpZ)]

     ENDFOR

     ;;rHat, thetaHat, and phiHat in GEO Cartesian coordinates
     ;;Note, -thetaHat is our NORTH_GEO vector
     rHat_GEO = REFORM(tmp_GEOSph2GEO[0,*])
     tHat_GEO = REFORM(tmp_GEOSph2GEO[1,*])
     pHat_GEO = REFORM(tmp_GEOSph2GEO[2,*])

     ;;So this is the (LOCAL) NORTH_GEO vector
     GEONorthHat_GEO = -tHat_GEO

     PRINT,DOTP(GEONorthHat_GEO,MAGNorthHat_GEO)

     PRINT,"Norm- notNorm GEONorthHat : ",VNORMALIZE(GEONorthHat_GEO)-GEONorthHat_GEO
     PRINT,"Theta:::::: ",ACOS(DOTP(GEONorthHat_GEO,MAGNorthHat_GEO))*!RADEG
     STOP

     GEOPACK_CONV_COORD_08,tmpPos_MAG[0],tmpPos_MAG[1],tmpPos_MAG[2], $
                           pos_x_MAG,pos_y_MAG,pos_z_MAG, $
                           /FROM_MAG,/TO_GEO,EPOCH=time_epoch[i]


     ;;And spherical everything
     ;; GEOPACK_SPHCAR_08,pos_x_GEI,pos_y_GEI,pos_z_GEI, $
     ;;                   pos_r_GEI,pos_theta_GEI,pos_phi_GEI,/TO_SPHERE,/DEGREE
     ;; GEOPACK_SPHCAR_08,tmpPos_GEO[0],tmpPos_GEO[1],tmpPos_GEO[2], $
     ;;                   pos_r_GEO,pos_theta_GEO,pos_phi_GEO,/TO_SPHERE,/DEGREE
     ;; GEOPACK_SPHCAR_08,pos_x_GSM,pos_y_GSM,pos_z_GSM, $
     ;;                   pos_r_GSM,pos_theta_GSM,pos_phi_GSM,/TO_SPHERE,/DEGREE
     ;; GEOPACK_SPHCAR_08,pos_x_MAG,pos_y_MAG,pos_z_MAG, $
     ;;                   pos_r_MAG,pos_theta_MAG,pos_phi_MAG,/TO_SPHERE,/DEGREE

     ;;velocity vector in spherical GEO coords
     GEOPACK_BCARSP_08,tmpPos_GEO[0],tmpPos_GEO[1],tmpPos_GEO[2], $
                       tmpVel_GEO[0],tmpVel_GEO[1],tmpVel_GEO[2], $
                       vel_GEO_r,vel_GEO_theta,vel_GEO_phi

     ;;Get IGRF
     pos_gsm_xyz_R_E         = [pos_x_GSM,pos_y_GSM,pos_z_GSM]/1000.D/R_E ;div by 1000 to get to km, then by R_E (which is in units of km)
     ;; GEOPACK_RECALC_08,YearArr[i],MonthArr[i],DayArr[i],HourArr[i],MinArr[i],SecArr[i],/DATE
     GEOPACK_IGRF_GEO_08,pos_r_GEO/1000./R_E,pos_theta_GEO,pos_phi_GEO, $
                         br_GEO,btheta_GEO,bphi_GEO, $
                         EPOCH=time_epoch[i],/DEGREE
     GEOPACK_BSPCAR_08  ,pos_theta_GEO,pos_phi_GEO, $
                         br_GEO,btheta_GEO,bphi_GEO, $
                         bx_GEO,by_GEO,bz_GEO,/DEGREE ; ,EPOCH=time_epoch[i]

     ;;alternate shot at IGRF
     GEOPACK_IGRF_GSW_08,pos_gsm_xyz_R_E[0],pos_gsm_xyz_R_E[1],pos_gsm_xyz_R_E[2], $
                         bx_GSM,by_GSM,bz_GSM, $
                         EPOCH=time_epoch[i]
     GEOPACK_BCARSP_08  ,pos_gsm_xyz_R_E[0],pos_gsm_xyz_R_E[1],pos_gsm_xyz_R_E[2], $
                         bx_GSM,by_GSM,bz_GSM, $
                         br_GSM,btheta_GSM,bphi_GSM

     ;;Dipole, please?
     GEOPACK_DIP_08,pos_gsm_xyz_R_E[1],pos_gsm_xyz_R_E[1],pos_gsm_xyz_R_E[2],bx_GSM_dip,by_GSM_dip,bz_GSM_dip,EPOCH=time_epoch[i] ;,/DEGREE

     ;;Update not-spherical
     ;; TiltArr              = [TiltArr,tempTilt]
     pos_GEI_arr[*,i]        = [pos_x_GEI,pos_y_GEI,pos_z_GEI]
     pos_GSM_arr[*,i]        = [pos_x_GSM,pos_y_GSM,pos_z_GSM]
     pos_MAG_arr[*,i]        = [pos_x_MAG,pos_y_MAG,pos_z_MAG]

     ;;Update spherical
     pos_GEISph_arr[*,i]     = [pos_theta_GEI,pos_phi_GEI,pos_r_GEI]
     pos_GEOSph_arr[*,i]     = [pos_theta_GEO,pos_phi_GEO,pos_r_GEO] ;Redundant, yes, but a check
     pos_GSMSph_arr[*,i]     = [pos_theta_GSM,pos_phi_GSM,pos_r_GSM]
     pos_MAGSph_arr[*,i]     = [pos_theta_MAG,pos_phi_MAG,pos_r_MAG]

     vel_GEOSph_arr[*,i]     = [vel_GEO_r,vel_GEO_theta,vel_GEO_phi]

     IGRF_GEO_sphArr[*,i]    = [btheta_GEO,bphi_GEO,br_GEO]
     IGRF_GSM_sphArr[*,i]    = [btheta_GSM,bphi_GSM,br_GSM]

     ;;Update IGRFs
     IGRF_GSM_DIP_arr[*,i]   = [bx_GSM_dip,by_GSM_dip,bz_GSM_dip]
     IGRF_GEO_arr[*,i]       = [bx_GEO,by_GEO,bz_GEO]
     IGRF_GSM_arr[*,i]       = [bx_GSM,by_GSM,bz_GSM]

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

     ;;GEO to NED (North-East-Down)
     ;; tmpLambda               = renu2.long[i]*!DTOR
     ;; tmpPhi                  = renu2.lat[i]*!DTOR

     ;;**Get colatitude (refTheta) and altitude in GEO coordinates for geodetic reference altitude (sea level)
     GEOPACK_GEODGEO_08,0.D,tmpPhi, $
                        refAlt_GEO,refTheta_GEO,/TO_GEOCENTRIC ;/TO_GEODETIC (see GEOPACK_2008 documentation, or type GEOPACK_HELP)
     GEOPACK_SPHCAR_08,refAlt_GEO,refTheta_GEO,tmpLambda, $
                       refAlt_x,refAlt_y,refAlt_z,/TO_RECT

     tmpGEODRef              = [refAlt_x,refAlt_y,refAlt_z]
     tmpNEDRefPos_GEO        = [tmpPos_GEO[0]-tmpGEODRef[0],tmpPos_GEO[1]-tmpGEODRef[1],tmpPos_GEO[2]-tmpGEODRef[2]]

     ;;GEI to VDH (Vertical [vehicle?]-Down-Horizontal)

     ;;**First, Get vector pointing east in GEO coordinates (later, we'll transform to GEI coordinates with GEO2GEI_arr)
     GEOPACK_BSPCAR_08,refTheta_GEO,tmpLambda, $
                       0.D,0.D,1.D, $
                       eEast_x_GEO,eEast_y_GEO,eEast_z_GEO
     eEast_Car_GEO           = [TEMPORARY(eEast_x_GEO),TEMPORARY(eEast_y_GEO),TEMPORARY(eEast_z_GEO)]

     ;;**Also, vector along magnetic pole (northward) in GEI coordinates, which is H_hat (if normalized) in VDH system
     GEOPACK_CONV_COORD_08,0.D,0.D,1.D, $
                           magPole_x0_GEI,magPole_y0_GEI,magPole_z0_GEI, $
                           /FROM_MAG,/TO_GEI,EPOCH=time_epoch[i]
     magPole_GEI             = [TEMPORARY(magPole_x0_GEI),TEMPORARY(magPole_y0_GEI),TEMPORARY(magPole_z0_GEI)]

     ;;**Here are some other ingredients for VDH transform
     hHat_GEI                = magPole_GEI
     ;; rPos_norm_GEI           = SQRT(pos_x_GEI*pos_x_GEI+pos_y_GEI*pos_y_GEI+pos_z_GEI*pos_z_GEI)
     ;; rPosHat_GEI             = [pos_x_GEI,pos_y_GEI,pos_z_GEI]/rPos_norm_GEI
     rPosHat_Car_GEI         = VNORMALIZE([pos_x_GEI,pos_y_GEI,pos_z_GEI])
     dHat_GEI                = CROSSP_NORMED(hHat_GEI,rPosHat_Car_GEI)
     vHat_GEI                = CROSSP_NORMED(dHat_GEI,hHat_GEI)

           ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
     ;;Field-aligned coordinate systems

     ;;FAC system 1 (ripped from UCLA_MAG_DESPIN)
     ;;   "Field-aligned coordinates defined as:   "
     ;;   "z-along B, y-east (BxR), x-nominally out"

     bHat_GEO                = VNORMALIZE([bx_GEO,by_GEO,bz_GEO]) ;along-b unit vector (or zHat)
     rPosHat_Car_GEO         = VNORMALIZE(tmpPos_GEO)
     eHat_GEO                = CROSSP_NORMED(bHat_GEO,rPosHat_Car_GEO) ;eastward unit vector (or yHat)
     oHat_GEO                = CROSSP_NORMED(eHat_GEO,bHat_GEO)        ;"nominally out" unit vector (or xHat)

     ;;FAC system 2 (ripped from UCLA_MAG_DESPIN)
     ;;   "Field-aligned velocity-based coordinates defined as: "
     ;;   "z-along B, y-cross track (BxV), x-along track ((BxV)xB)."

     vHat_GEO                = VNORMALIZE(tmpVel_GEO)
     cHat_GEO                = CROSSP_NORMED(bHat_GEO,vHat_GEO) ;cross-track unit vector (or yHat)
     aHat_GEO                = CROSSP_NORMED(cHat_GEO,bHat_GEO) ;along-track unit vector (or xHat)

           ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
     ;;Update arrays of transformation matrices
     GEO2GEI_coord[*,*,i]    = [[gei_x0,gei_y0,gei_z0], $
                                [gei_x1,gei_y1,gei_z1], $
                                [gei_x2,gei_y2,gei_z2]]
     GEO2MAG_coord[*,*,i]    = [[mag_x0,mag_y0,mag_z0], $
                                [mag_x1,mag_y1,mag_z1], $
                                [mag_x2,mag_y2,mag_z2]]
     GEO2NED_coord[*,*,i]    = [[-SIN(tmpPhi)*COS(tmpLambda),-SIN(tmpPhi)*SIN(tmpLambda),COS(tmpPhi)], $
                                [-SIN(tmpLambda)            ,COS(tmpLambda)             ,0.D         ], $
                                [-COS(tmpPhi)*COS(tmpLambda),-COS(tmpPhi)*SIN(tmpLambda),-SIN(tmpPhi)]]

     GEO2FACV_vec[*,*,i]     = INVERT([[aHat_GEO],[cHat_GEO],[bHat_GEO]])
     GEO2FAC_vec[*,*,i]      = INVERT([[oHat_GEO],[eHat_GEO],[bHat_GEO]])
     GEO2GEI_vec[*,*,i]      = INVERT(GEO2GEI_coord[*,*,i])
     GEO2MAG_vec[*,*,i]      = INVERT(GEO2MAG_coord[*,*,i])
     GEO2NED_vec[*,*,i]      = INVERT(GEO2NED_coord[*,*,i])

     ;;get velocity and magnetic field (as a check) in FAC and FACV coordinates
     vel_FAC_arr[*,i]        = GEO2FAC_vec[*,*,i]    # tmpVel_GEO
     vel_FACV_arr[*,i]       = GEO2FACV_vec[*,*,i]   # tmpVel_GEO

     IGRF_FAC_arr[*,i]       = GEO2FAC_vec[*,*,i]    # [bx_GEO,by_GEO,bz_GEO]

     ;;get velocity and magnetic field in GEI coordinates
     vel_GEI_arr[*,i]        = GEO2GEI_vec[*,*,i]    # tmpVel_GEO
     IGRF_GEI_arr[*,i]       = GEO2GEI_vec[*,*,i]    # [bx_GEO,by_GEO,bz_GEO]

     ;;get velocity in MAG coordinates (Cartesian and spherical)
     vel_MAG_arr[*,i]        = GEO2MAG_vec[*,*,i]    # tmpVel_GEO
     GEOPACK_BCARSP_08,pos_x_MAG,pos_y_MAG,pos_z_MAG, $
                       vel_MAG_arr[0,i],vel_MAG_arr[1,i],vel_MAG_arr[2,i], $
                       magVSph_r,magVSph_theta,magVSph_phi
     vel_MAGSph_arr[*,i]     = [magVSph_r,magVSph_theta,magVSph_phi]


     pos_NED_arr[*,i]        = GEO2NED_coord[*,*,i]  # TEMPORARY(tmpNEDRefPos_GEO)
     vel_NED_arr[*,i]        = GEO2NED_vec[*,*,i]    # tmpVel_GEO

     ;;get velocity and magnetic field in VDH

     eEast_Car_GEI           = GEO2GEI_vec[*,*,i]    # eEast_Car_GEO
     ;; IF (TRANSPOSE(dHat_GEI) # eEast_Car_GEI) LT 0 THEN BEGIN
     IF DOTP(dHat_GEI,eEast_Car_GEI) LT 0 THEN BEGIN
        dHat_GEI            *= (-1D)
        vHat_GEI            *= (-1D)
     ENDIF

     GEI2VDH_vec[*,*,i]      = [[vHat_GEI],[dHat_GEI],[hHat_GEI]]

     vel_VDH_arr[*,i]        = GEI2VDH_vec[*,*,i]    # vel_GEI_arr[*,i]
     IGRF_VDH_arr[*,i]       = GEI2VDH_vec[*,*,i]    # IGRF_GEI_arr[*,i]

     ;;Convert Cartesian vector transform matrices to spherical vector transform matrices
     ;; GEOPACK_BCARSP_08,mag_x0,mag_y0,mag_z0, $
     ;; GEOPACK_BCARSP_08,pos_x_MAG,pos_y_MAG,pos_z_MAG, $
     ;;                   ;; GEOPACK_BCARSP_08,ident[0,0],ident[0,1],ident[0,2], $
     ;;                   ;; GEO2MAG_vec[0,0,i],GEO2MAG_vec[0,1,i],GEO2MAG_vec[0,2,i], $
     ;;                   ident[0,0],ident[0,1],ident[0,2], $
     ;;                   magVSph_r0,magVSph_theta0,magVSph_phi0
     ;; ;; GEOPACK_BCARSP_08,mag_x1,mag_y1,mag_z1, $
     ;; GEOPACK_BCARSP_08,pos_x_MAG,pos_y_MAG,pos_z_MAG, $
     ;;                   ;; GEOPACK_BCARSP_08,ident[1,0],ident[1,1],ident[1,2], $
     ;;                   ;; GEO2MAG_vec[1,0,i],GEO2MAG_vec[1,1,i],GEO2MAG_vec[1,2,i], $
     ;;                   ident[1,0],ident[1,1],ident[1,2], $
     ;;                   magVSph_r1,magVSph_theta1,magVSph_phi1
     ;; ;; GEOPACK_BCARSP_08,mag_x2,mag_y2,mag_z2, $
     ;; GEOPACK_BCARSP_08,pos_x_MAG,pos_y_MAG,pos_z_MAG, $
     ;;                   ;; GEOPACK_BCARSP_08,ident[2,0],ident[2,1],ident[2,2], $
     ;;                   ;; GEO2MAG_vec[2,0,i],GEO2MAG_vec[2,1,i],GEO2MAG_vec[2,2,i], $
     ;;                   ident[2,0],ident[2,1],ident[2,2], $
     ;;                   magVSph_r2,magVSph_theta2,magVSph_phi2

     ;; GEO2MAG_sphV[*,*,i]     = [[magVSph_r0,magVSph_theta0,magVSph_phi0], $
     ;;                            [magVSph_r1,magVSph_theta1,magVSph_phi1], $
     ;;                            [magVSph_r2,magVSph_theta2,magVSph_phi2]]



     ;;Now update spherical coordinate transform matrices
     ;; GEOPACK_SPHCAR_08,mag_x0,mag_y0,mag_z0,mag_r0,mag_theta0,mag_phi0,/TO_SPHERE,/DEGREE
     ;; GEOPACK_SPHCAR_08,mag_x1,mag_y1,mag_z1,mag_r1,mag_theta1,mag_phi1,/TO_SPHERE,/DEGREE
     ;; GEOPACK_SPHCAR_08,mag_x2,mag_y2,mag_z2,mag_r2,mag_theta2,mag_phi2,/TO_SPHERE,/DEGREE

     ;; GEO2MAG_sphC[*,*,i]     = [[mag_r0,mag_theta0,mag_phi0], $
     ;;                            [mag_r1,mag_theta1,mag_phi1], $
     ;;                            [mag_r2,mag_theta2,mag_phi2]]

     ;;... And, TESTS. These should come out the same:
     GEOPACK_BSPCAR_08  ,pos_theta_MAG,pos_phi_MAG, $
                         vel_MAGSph_arr[0,i],vel_MAGSph_arr[1,i],vel_MAGSph_arr[2,i], $
                         testVelx_MAG,testVely_MAG,testVelz_MAG,/DEGREE ; ,EPOCH=time_epoch[i]
     GEOPACK_BCARSP_08  ,pos_gsm_xyz_R_E[0],pos_gsm_xyz_R_E[1],pos_gsm_xyz_R_E[2], $
                         vel_MAG_arr[0,i],vel_MAG_arr[1,i],vel_MAG_arr[2,i], $
                         testVelr_MAG,testVeltheta_MAG,testVelphi_MAG


     velMagCar_GEO           = SQRT(tmpVel_GEO[0]*tmpVel_GEO[0] + $
                                    tmpVel_GEO[1]*tmpVel_GEO[1] + $
                                    tmpVel_GEO[2]*tmpVel_GEO[2])
     velMagSph_GEO           = SQRT(vel_GEO_r*vel_GEO_r + $
                                    vel_GEO_theta*vel_GEO_theta + $
                                    vel_GEO_phi*vel_GEO_phi)

     ;; IF ABS(velMagCar_GEO-velMagSph_GEO) GT 0.0001 THEN STOP ;Precision to the fourth decimal place. This will make it quit.

     conv_velMagSph_MAG      = SQRT(vel_MAGSph_arr[0,i]*vel_MAGSph_arr[0,i] + $
                                    vel_MAGSph_arr[1,i]*vel_MAGSph_arr[1,i] + $
                                    vel_MAGSph_arr[2,i]*vel_MAGSph_arr[2,i])

     conv_velMagCar_MAG      = SQRT(vel_MAG_arr[0,i]*vel_MAG_arr[0,i] + $
                                    vel_MAG_arr[1,i]*vel_MAG_arr[1,i] + $
                                    vel_MAG_arr[2,i]*vel_MAG_arr[2,i])

     convConv_velMagSph_MAG  = SQRT(testVelr_MAG*testVelr_MAG + $
                                    testVeltheta_MAG*testVeltheta_MAG + $
                                    testVelphi_MAG*testVelphi_MAG)

     convConv_velMagCar_MAG  = SQRT(testVelx_MAG*testVelx_MAG + $
                                    testVely_MAG*testVely_MAG + $
                                    testVelz_MAG*testVelz_MAG)

     velMagnitudeDiff        = conv_velMagSph_MAG-conv_velMagCar_MAG
     velSphMagnitudeDiff     = conv_velMagSph_MAG-convConv_velMagSph_MAG
     velCarMagnitudeDiff     = conv_velMagCar_MAG-convConv_velMagCar_MAG

     IF ABS(velMagCar_GEO-velMagSph_GEO) GT 0.001 THEN STOP
     IF ABS(velMagCar_GEO-conv_velMagCar_MAG) GT 0.001 THEN STOP
     IF ABS(velMagSph_GEO-conv_velMagSph_MAG) GT 0.001 THEN STOP

     IF (ABS(velCarMagnitudeDiff) GT 1) OR (ABS(velSphMagnitudeDiff) GT 1) OR (ABS(velMagnitudeDiff) GT 1) THEN STOP

     ;;More tests! Determinants.
     maxDiff = 0.001
     IF ABS(DETERM(GEO2FAC_vec  [*,*,i]) - 1.) GT maxDiff THEN STOP
     IF ABS(DETERM(GEO2FACV_vec [*,*,i]) - 1.) GT maxDiff THEN STOP

     IF ABS(DETERM(GEO2GEI_coord[*,*,i]) - 1.) GT maxDiff THEN STOP
     IF ABS(DETERM(GEO2GEI_vec  [*,*,i]) - 1.) GT maxDiff THEN STOP

     IF ABS(DETERM(GEO2MAG_coord[*,*,i]) - 1.) GT maxDiff THEN STOP
     IF ABS(DETERM(GEO2MAG_vec  [*,*,i]) - 1.) GT maxDiff THEN STOP

     IF ABS(DETERM(GEO2NED_coord[*,*,i]) - 1.) GT maxDiff THEN STOP
     IF ABS(DETERM(GEO2NED_vec  [*,*,i]) - 1.) GT maxDiff THEN STOP

     IF ABS(DETERM(GEI2VDH_vec  [*,*,i]) - 1.) GT maxDiff THEN STOP

     IF (i MOD 100) EQ 0 THEN PRINT,i
  ENDFOR


END
