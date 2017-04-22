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
;Rob's question: What is the angle between MAG north and GEO north at a particular location?
PRO JOURNAL__20170420__ROTVEC_MAG_TO_GEO_FOR_ROB,PRACTICE=practice

  COMPILE_OPT IDL2,STRICTARRSUBS

  long_MAG = 267.0D
  lat_MAG  = 70.0D

  ;;For practice 
  IF KEYWORD_SET(practice) THEN BEGIN
     long_MAG = 0.D
     lat_MAG  = 0.D
  ENDIF

  ;;rad, colat
  pos_MAGSph = [ 1, ((90.D - lat_MAG) * !DTOR) , long_MAG * !DTOR ]

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
  
  pos_MAGSph_arr[*,0] = pos_MAGSph

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
     ;;Wait, no. Just use spherical. Unless you want GEODETIC
     GEOPACK_GEODGEO_08,tmpPos_GEOSph[0],tmpPos_GEOSph[1], $
                        refAlt_GEO,refLat_GEO,/TO_GEODETIC
     
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

     ;;More tests! Determinants.
     maxDiff = 0.001
     IF ABS(DETERM(tmp_GEOSph2GEO  [*,*,i]) - 1.) GT maxDiff THEN STOP

     IF ABS(DETERM(MAG2GEO_coord[*,*,i]) - 1.) GT maxDiff THEN STOP
     IF ABS(DETERM(MAG2GEO_vec  [*,*,i]) - 1.) GT maxDiff THEN STOP

     IF (i MOD 100) EQ 0 THEN PRINT,i

     PRINT,FORMAT='(A0,T20,": ",3(F6.2,:,","))',"Lat,long,alt (MAG)",90.-tmpPos_MAGSph[1]*!RADEG,tmpPos_MAGSph[2]*!RADEG,tmpPos_MAGSph[0]
     PRINT,FORMAT='(A0,T20,": ",3(F6.2,:,","))',"Lat,long,alt (GEO)",refLat_GEO*!RADEG,tmpPos_GEOSph[2]*!RADEG,refAlt_GEO

     PRINT,"Norm- notNorm GEONorthHat : ",VNORMALIZE(GEONorthHat_GEO)-GEONorthHat_GEO
     PRINT,"Dot prod::: ",DOTP(GEONorthHat_GEO,MAGNorthHat_GEO)
     PRINT,"Theta:::::: ",ACOS(DOTP(GEONorthHat_GEO,MAGNorthHat_GEO))*!RADEG
     STOP

     GEOPACK_CONV_COORD_08,tmpPos_MAG[0],tmpPos_MAG[1],tmpPos_MAG[2], $
                           pos_x_MAG,pos_y_MAG,pos_z_MAG, $
                           /FROM_MAG,/TO_GEO,EPOCH=time_epoch[i]


  ENDFOR


END
