;;2017/03/03
PRO JOURNAL__20170303__OUTPUT_COORDINATE_CONVERSIONS_TO_CSV

  COMPILE_OPT IDL2

  outDir           = '/SPENCEdata/Research/database/RENU2/'

  inFile           = 'RENU2_GPS.sav'
  timeStrFile      = 'RENU2_GPS--timeStr.sav'
  outFile          = 'RENU2_coordinates.sav'
  
  RESTORE,outDir+inFile
  RESTORE,outDir+outFile

  allTags       = TAG_NAMES_R(coords)
  nTags         = N_ELEMENTS(allTags)
  splitTagsList = STRSPLIT(allTags,'.',/EXTRACT)

  level = 0
  done = 0B
  levelSelection = !NULL
  WHILE ~done DO BEGIN

     ;;Get all unique tag names on this level 
     theseTags = !NULL
     ;; nnTags    = 0
     FOR k=0,nTags-1 DO BEGIN
        IF N_ELEMENTS(splitTagsList[k]) GT level THEN BEGIN
           theseTags = [theseTags,splitTagsList[k,level]]
           ;; nnTags++
        ENDIF
     ENDFOR

     ;;Get uniquers
     theseTags = theseTags[UNIQ(theseTags,SORT(theseTags))]
     nnTags    = N_ELEMENTS(theseTags)

     PRINT,"Which would you like? (Or type 'quit')"
     
     FOR kk=0,nnTags-1 DO PRINT,FORMAT='(I3,") ",A0)',kk,theseTags[kk]

     cont = 0B
     tries = 0B
     read = ''
     WHILE ~cont DO BEGIN

        READ,read

        STOP
        IF STRMATCH(read,'quit',/FOLD_CASE) THEN BEGIN

           cont = 1B
           broken = 1B

        ENDIF ELSE BEGIN

           bro = WHERE(STRMATCH(theseTags,read,/FOLD_CASE),nMatch)

           CASE nMatch OF
              0: BEGIN

                 ;;Maybe it's an int?
                 CATCH,error_status

                 tryNum = FIX(read)

                 isNum  = 1B
                 isValid = 1B
                 IF error_status EQ 0 THEN BEGIN

                    IF FIX(read) LT nnTags THEN BEGIN
                       
                       levelSelection = [levelSelection,theseTags[read]]
                       
                       nextLevel      = WHERE(STRMATCH(allTags,'*'+theseTags[read]+'*',/FOLD_CASE),nTags)
                       IF nTags GT 0 THEN BEGIN

                          allTags       = allTags[nextLevel]
                          splitTagsList = STRSPLIT(allTags,'.',/EXTRACT)

                          level++
                          cont        = 1B
                          
                       ENDIF ELSE BEGIN
                          PRINT,"DEATH"
                          STOP
                       ENDELSE

                    ENDIF ELSE BEGIN
                       isValid = 0B
                    ENDELSE

                 ENDIF ELSE BEGIN

                    isNum = 0B
                    PRINT,"errNum: ",error_status
                    PRINT,!ERROR_STATE.MSG
                    CATCH,/CANCEL


                 ENDELSE


                 IF isNum THEN BEGIN

                    STOP

                 ENDIF ELSE BEGIN
                    CASE 1 OF
                       (tries LT 5): BEGIN
                          PRINT,"No! Try again."
                       END
                       ELSE: BEGIN
                          PRINT,"You hosed it."
                          broken = 1B
                          cont   = 1B
                       ENDELSE

                    ENDCASE
                    tries++
                    
                 ENDELSE

              END
              1: BEGIN
                 PRINT,"Match: ",theseTags[bro]
                 STOP
              END
              ELSE: BEGIN
                 PRINT,"manyMatch: ",theseTags[bro]
                 STOP
              END
           ENDCASE

        ENDELSE

     ENDWHILE

     IF KEYWORD_SET(broken) THEN break

  ENDWHILE

  STOP

END
