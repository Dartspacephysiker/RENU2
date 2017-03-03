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
  levelSelNum  = !NULL
  levelSelName = !NULL
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
     uniqTag_i = UNIQ(theseTags,SORT(theseTags))
     uniqTags  = theseTags[uniqTag_i]
     nnTags    = N_ELEMENTS(uniqTags)

     CASE nnTags OF
        1: BEGIN

           read  = ''
           cont  = 0B
           tries = 0B
           WHILE ~cont DO BEGIN

              levelSelName = [levelSelName,uniqTags[0]]
              levelSelNum  = [levelSelNum,uniqTag_i[0]]

              ;;Nous sommes arrivés
              PRINT,levelSelName
              PRINT,levelSelNum
              PRINT,"Look OK? (y/n/q)"
              READ,read

              CASE 1 OF
                 STRMATCH(read,'y*',/FOLD_CASE) OR STRMATCH(read,'1',/FOLD_CASE): BEGIN

                    STOP
                    
                    allTags         = TAG_NAMES_R(coords)
                    splitTagsList   = STRSPLIT(allTags,'.',/EXTRACT)
                    matchTags       = WHERE(STRMATCH(allTags, $
                                                     '*' + STRJOIN(levelSelName,'.') + '*', $
                                                     /FOLD_CASE), $
                                            nMatch)

                    mems            = !NULL
                    FOR k=0,nMatch-1 DO mems = [mems,splitTagsList[matchTags[k],-1]]
                    PRINT,FORMAT='("There are ",I0," members here: ",' + $
                          STRCOMPRESS(nMatch,/REMOVE_ALL) + '(A0,:,"  "))',nMatch,mems

                 END
                 STRMATCH(read,'n*',/FOLD_CASE) OR STRMATCH(read,'0',/FOLD_CASE): BEGIN
                    PRINT,"Trying again ..."

                    allTags         = TAG_NAMES_R(coords)
                    nTags           = N_ELEMENTS(allTags)
                    splitTagsList   = STRSPLIT(allTags,'.',/EXTRACT)

                    level           = 0
                    levelSelName    = !NULL
                    levelSelNum     = !NULL
                    cont            = 1B

                 END
                 STRMATCH(read,'q*',/FOLD_CASE) : BEGIN
                    cont            = 1B
                    doQuit          = 1B

                 END
                 ELSE: BEGIN

                    IF tries GT maxTries THEN BEGIN
                       PRINT,"Can't understand a word comin' out of his mouth ..."
                       cont         = 1B
                       doQuit       = 1B
                    ENDIF ELSE BEGIN
                       PRINT,"Huh?"
                       tries++
                    ENDELSE

                 END
              ENDCASE

           ENDWHILE

        END
        ELSE: BEGIN

           cont  = 0B
           tries = 0B
           read  = ''
           WHILE ~cont DO BEGIN

              PRINT,"Which would you like? (Or type 'quit')"
              
              FOR kk=0,nnTags-1 DO PRINT,FORMAT='(I3,") ",A0)',kk,uniqTags[kk]
              READ,read

              ;; STOP
              CASE 1 OF
                 (read EQ ''): BEGIN

                 END
                 STRMATCH(read,'quit',/FOLD_CASE): BEGIN

                    PRINT,'Quitting ...'

                    cont = 1B
                    doQuit = 1B

                 END
                 ELSE: BEGIN

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
                                
                                levelSelName = [levelSelName,uniqTags[read]]
                                levelSelNum  = [levelSelNum,uniqTag_i[read]]
                                
                                nextLevel      = WHERE(STRMATCH(allTags,'*'+uniqTags[read]+'*',/FOLD_CASE),nTags)
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

                             ;;What to do? Keep it movin'
                             
                          ENDIF ELSE BEGIN
                             CASE 1 OF
                                (tries LT 5): BEGIN
                                   PRINT,"No! Try again."
                                END
                                ELSE: BEGIN
                                   PRINT,"You hosed it."
                                   doQuit = 1B
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
              ENDCASE

           ENDWHILE

        END
     ENDCASE

     IF KEYWORD_SET(doQuit) THEN BREAK

  ENDWHILE

  STOP

END
