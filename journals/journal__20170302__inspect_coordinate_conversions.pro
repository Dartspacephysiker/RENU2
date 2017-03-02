;2017/03/02
;See what we really got from JOURNAL__20170213__CONVERT_GPS_COORDS_TO_MAGNETIC
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

PRO JOURNAL__20170302__INSPECT_COORDINATE_CONVERSIONS,NORMALIZE=normalize, $
   PLOT_MAGNITUDE=plot_magnitude, $
   BUFFER=buffer

  COMPILE_OPT IDL2

  outDir           = '/SPENCEdata/Research/database/RENU2/'

  inFile           = 'RENU2_GPS.sav'
  timeStrFile      = 'RENU2_GPS--timeStr.sav'
  outFile          = 'RENU2_coordinates.sav'
  
  RESTORE,outDir+inFile
  RESTORE,outDir+outFile


  ;; xComp             = coords.igrf.vdh.car.v
  ;; yComp             = coords.igrf.vdh.car.d
  ;; zComp             = coords.igrf.vdh.car.h
  ;; title                = "B-field"
  ;; names                = ['V','D','H']

  ;; xComp             = coords.igrf.fac.car.o
  ;; yComp             = coords.igrf.fac.car.e
  ;; zComp             = coords.igrf.fac.car.b
  ;; title                = "B-field"
  ;; names                = ['Out','East','Along B','Magnitude']

  ;; xComp                = coords.vel.fac.car.o/1000.D
  ;; yComp                = coords.vel.fac.car.e/1000.D
  ;; zComp                = coords.vel.fac.car.b/1000.D
  ;; title                = "Velocity (km/s)"
  ;; names                = ['Out','East','Along B','Speed']

  ;; xComp                = coords.vel.ned.car.x/1000.D
  ;; yComp                = coords.vel.ned.car.y/1000.D
  ;; zComp                = coords.vel.ned.car.z/1000.D
  ;; title                = "Velocity (km/s)"
  ;; names                = ['North','East','Down','Speed']
  ;; legPos               = [600,-0.5]
  ;; legCoord_data        = 1
  ;; legCoord_normal      = 0

  xComp                = coords.vel.vdh.car.v/1000.D
  yComp                = coords.vel.vdh.car.d/1000.D
  zComp                = coords.vel.vdh.car.h/1000.D
  title                = "Velocity (km/s)"
  names                = ['Vehicle','Dipole','Horizon','Speed']
  legPos               = [50,-1.5]
  legCoord_data        = 1
  legCoord_normal      = 0

  ;; xComp                = coords.vel.facv.car.a/1000.D
  ;; yComp                = coords.vel.facv.car.c/1000.D
  ;; zComp                = coords.vel.facv.car.b/1000.D
  ;; title                = "Velocity (km/s)"
  ;; names                = ['Along track','Cross track','Along B','Speed']
  ;; legPos               = [600,-0.5]
  ;; legCoord_data        = 1
  ;; legCoord_normal      = 0

  lineStyle            = ['-','--',':']
  thick                = [2.0,2.0,2.0]
  IF KEYWORD_SET(normalize) THEN BEGIN
     IGRFVDHNorm       = SQRT(xComp*xComp + $
                              yComp*yComp + $
                              zComp*zComp)

     xComp            /= IGRFVDHNorm
     yComp            /= IGRFVDHNorm
     zComp            /= IGRFVDHNorm
  ENDIF ELSE BEGIN
     IGRFVDHNorm       = 1.D
  ENDELSE

  ;; yLims = MINMAX([xComp,yComp,zComp]/[IGRFVDHNorm,IGRFVDHNorm,IGRFVDHNorm])
  xQuant = renu2.time.flight
  xTitle = 'Flight time (s)'

  nPlots = 3 + KEYWORD_SET(plot_magnitude)
  plotArr = MAKE_ARRAY(nPlots,/OBJ)

  window = WINDOW(DIMENSIONS=KEYWORD_SET(buffer) ? !NULL : [1000,800])

  plotArr[0] = PLOT(xQuant,xComp, $
                    NAME=names[0], $
                    COLOR='brown', $
                    LINESTYLE=lineStyle[0], $
                    THICK=thick[0], $
                    ;; SYMBOL='*', $
                    YRANGE=yLims, $
                    YTITLE=title, $
                    XTITLE=xTitle, $
                    /CURRENT)

  plotArr[1] = PLOT(xQuant,yComp, $
                    NAME=names[1], $
                    COLOR='Red', $
                    LINESTYLE=lineStyle[1], $
                    THICK=thick[1], $
                    ;; SYMBOL='*', $
                    /OVERPLOT)

  plotArr[2] = PLOT(xQuant,zComp, $
                    NAME=names[2], $
                    COLOR='Blue', $
                    LINESTYLE=lineStyle[2], $
                    THICK=thick[2], $
                    ;; SYMBOL='*', $
                    /OVERPLOT)

  IF KEYWORD_SET(plot_magnitude) THEN BEGIN
     magnitude  = SQRT(xComp*xComp+yComp*yComp+zComp*zComp)
     plotArr[3] = PLOT(xQuant,magnitude, $
                       NAME=names[3], $
                       /OVERPLOT)
  ENDIF
  leg   = LEGEND(TARGET=plotArr,POSITION=legPos,DATA=legCoord_data,NORMAL=legCoord_normal)

  STOP


END
