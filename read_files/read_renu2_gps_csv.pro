;;02/13/17
;;Why is the CSV file lying? It gives GPS weeks as 851, but RENU2 launched Dec 13, 2015 at 0734 UTC, which corresponds to GPS week 1875.
;;saveFile = '/SPENCEdata/Research/database/RENU2/RENU2_GPS.sav'
PRO READ_RENU2_GPS_CSV,ASCII_tmplt, $
                       SAVEFILE=saveFile

  COMPILE_OPT IDL2

  renu2dir    = '/SPENCEdata/Research/database/RENU2/'
  renu2fil    = 'RENU2_GPS.csv'
  
  IF ~FILE_TEST(renu2dir+renu2fil) THEN STOP

  navner = ['Flight_time'   , 'Delta_Time'               , 'ECEF_X_POS'       , 'ECEF_Y_POS'             , 'ECEF_Z_POS'         , 'ECEF_X_VEL'         , $
            'ECEF_Y_VEL'    , 'ECEF_Z_VEL'               , 'Lat'              , 'Long'                   , 'Alt'                , 'rV_Sample_Number'   , $
            'SVs_Used'      , 'Position_Valid'           , 'Velocity_Valid'   , 'PDOP'                   , 'Position_Coord_Ref' , 'Pos_Comp_Mode'      , $
            'Vel_Comp_Mode' , 'RMS_Pos_Error'            , 'RMS_Vel_Error'    , 'JNS_Datum_Number'       , 'SNR_133'            , 'SNR_135'            , $
            'SNR_138'       , 'Number_of_rM_SVD_Records' , 'rM_Sample_Number' , 'GPS_Time__mSec_of_week' , 'GPS_Week'           , 'Leap_Second_Status' , $
            'Time_Scale_ID' , 'Clock_Offset_Sec'         , 'Clock_Offset_ID'  , 'Clock_Offset__Delta_Time']

  ;; ASCII_tmplt = ASCII_TEMPLATE(renu2dir+renu2fil)

  ;; STOP
  fieldTypes      = [4, 4, 4, 4, 4, 4, $
                     4, 4, 4, 4, 4, 3, $
                     3, 3, 3, 4, 3, 3, $
                     3, 4, 4, 3, 3, 3, $
                     3, 3, 3, 3, 3, 3, $
                     3, 4, 3, 3]

  fieldLocations  = [   6,  23,  28,  43,  58,  75, $
                        88, 101, 111, 123, 137, 161, $
                        175, 191, 206, 211, 234, 252, $
                        270, 282, 298, 324, 334, 344, $
                        354, 383, 402, 427, 444, 468, $
                        483, 492, 520, 533]

  fieldGroups     = [ 0,  1,  2,  3,  4,  5, $
                      6,  7,  8,  9, 10, 11, $
                      12, 13, 14, 15, 16, 17, $
                      18, 19, 20, 21, 22, 23, $
                      24, 25, 26, 27, 28, 29, $
                      30, 31, 32, 33]

  ascii_tmplt      = {VERSION        : 1.0            , $
                      DATASTART      : 1              , $
                      DELIMITER      : BYTE(32)       , $
                      MISSINGVALUE   : !VALUES.F_NaN  , $
                      COMMENTSYMBOL  : ''             , $
                      FIELDCOUNT     : 34             , $
                      FIELDTYPES     : fieldTypes     , $
                      FIELDNAMES     : navner         , $
                      FIELDLOCATIONS : fieldLocations , $
                      FIELDGROUPS    : fieldGroups}

  gps              = READ_ASCII(renu2dir+renu2fil,TEMPLATE=ascii_tmplt)

  pos_ecef         = {x          : gps.ecef_x_pos     , $
                      y          : gps.ecef_y_pos     , $
                      z          : gps.ecef_z_pos     , $
                      valid      : gps.position_valid , $
                      comp_mode  : gps.pos_comp_mode  , $
                      rms_error  : gps.rms_pos_error  , $
                      coord_ref  : gps.position_coord_ref}

  vel_ecef         = {x          : gps.ecef_x_vel     , $
                      y          : gps.ecef_y_vel     , $
                      z          : gps.ecef_z_vel     , $
                      valid      : gps.velocity_valid , $
                      comp_mode  : gps.vel_comp_mode  , $
                      rms_error  : gps.rms_vel_error}

  snr              = {_133 : gps.snr_133, $
                      _135 : gps.snr_135, $
                      _138 : gps.snr_138}

  clock_offset     = {sec        : gps.clock_offset_sec        , $
                      id         : gps.clock_offset_id         , $
                      delta_time : gps.clock_offset__delta_time}

  time             = {utc          : GPS_TO_UTC(1875,gps.gps_time__msec_of_week), $
                      flight       : gps.flight_time                            , $
                      flight_dt    : gps.delta_time                             , $
                      gps          : {msec_of_week    : gps.gps_time__mSec_of_week, $
                                     week            : gps.gps_week              , $
                                     leap_sec_status : gps.leap_second_status    , $
                                     time_scale_id   : gps.time_scale_id}      , $
                      clock_offset : TEMPORARY(clock_offset)}

  ephem            = {lat        : gps.lat, $
                      long       : gps.long, $
                      alt        : gps.alt}

  STR_ELEMENT,gps,'ECEF_X_POS',/DELETE
  STR_ELEMENT,gps,'ECEF_Y_POS',/DELETE
  STR_ELEMENT,gps,'ECEF_Z_POS',/DELETE
  STR_ELEMENT,gps,'POSITION_VALID',/DELETE
  STR_ELEMENT,gps,'POS_COMP_MODE',/DELETE
  STR_ELEMENT,gps,'RMS_POS_ERROR',/DELETE
  STR_ELEMENT,gps,'POSITION_COORD_REF',/DELETE
  STR_ELEMENT,gps,'ECEF_X_VEL',/DELETE
  STR_ELEMENT,gps,'ECEF_Y_VEL',/DELETE
  STR_ELEMENT,gps,'ECEF_Z_VEL',/DELETE
  STR_ELEMENT,gps,'VELOCITY_VALID',/DELETE
  STR_ELEMENT,gps,'VEL_COMP_MODE',/DELETE
  STR_ELEMENT,gps,'RMS_VEL_ERROR',/DELETE
  STR_ELEMENT,gps,'SNR_133',/DELETE
  STR_ELEMENT,gps,'SNR_135',/DELETE
  STR_ELEMENT,gps,'SNR_138',/DELETE
  STR_ELEMENT,gps,'CLOCK_OFFSET_SEC',/DELETE
  STR_ELEMENT,gps,'CLOCK_OFFSET_ID',/DELETE
  STR_ELEMENT,gps,'CLOCK_OFFSET__DELTA_TIME',/DELETE
  STR_ELEMENT,gps,'FLIGHT_TIME',/DELETE
  STR_ELEMENT,gps,'DELTA_TIME',/DELETE
  STR_ELEMENT,gps,'GPS_TIME__MSEC_OF_WEEK',/DELETE
  STR_ELEMENT,gps,'GPS_WEEK',/DELETE
  STR_ELEMENT,gps,'TIME_SCALE_ID',/DELETE
  STR_ELEMENT,gps,'LEAP_SECOND_STATUS',/DELETE
  STR_ELEMENT,gps,'LAT',/DELETE
  STR_ELEMENT,gps,'LONG',/DELETE
  STR_ELEMENT,gps,'ALT',/DELETE

  ancillary  = TEMPORARY(gps)

  renu2            = {time       : TEMPORARY(time), $
                      lat        : ephem.lat, $
                      long       : ephem.long, $
                      alt        : ephem.alt , $
                      ecef       : {position : TEMPORARY(pos_ecef), $
                                    velocity : TEMPORARY(vel_ecef)}, $
                      snr        : TEMPORARY(snr), $
                      ancillary  : TEMPORARY(ancillary)}

  IF KEYWORD_SET(saveFile) THEN BEGIN
     IF FILE_TEST(saveFile) THEN BEGIN
        PRINT,"File exists! You sure?"
        STOP
     ENDIF

     PRINT,'Saving to ' + saveFile + ' ...'
     SAVE,renu2,FILENAME=saveFile
  ENDIF

END
