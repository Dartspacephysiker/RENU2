;;2017/03/03
PRO PRINT_CSV,outMe,mem,levelSelName,fileName,CSV_dir,doQuit

  PRINT,"Outputting " + fileName + ' ...'

  ;;Get formats for output
  defaultOutNames = TAG_NAMES(outMe)
  nOut            = N_ELEMENTS(defaultOutNames)

  confirmed       = 0B

  WHILE ~confirmed DO BEGIN

     defaultFormatStr = !NULL

     FOR k=0,nOut-1 DO BEGIN
        intTypes = [2,3,12,13,14,15]
        floatTypes  = [4]
        doubleTypes = [5]
        stringTypes = [7]

        ;;First, format string
        tmpType  = SIZE(outMe.(k),/TYPE)
        CASE 1 OF
           ((WHERE(tmpType EQ intTypes))[0] NE -1): BEGIN
              defaultFormatStr = [defaultFormatStr,'I0']
           END
           ((WHERE(tmpType EQ floatTypes))[0] NE -1): BEGIN
              defaultFormatStr = [defaultFormatStr,'F0.10']
           END
           ((WHERE(tmpType EQ doubleTypes))[0] NE -1): BEGIN
              defaultFormatStr = [defaultFormatStr,'F0.15']
           END
           ((WHERE(tmpType EQ stringTypes))[0] NE -1): BEGIN
              defaultFormatStr = [defaultFormatStr,'A0']
           END
           ELSE: BEGIN
              PRINT,"Just got " + STRCOMPRESS(tmpType,/REMOVE_ALL) + " for " + defaultOutNames[k] + ". What is this type? "
              PRINT,"Kidding. Type how you want this to appear (no need for quotation marks):"

              cont = 0B
              tries = 0
              doQuit = 0B
              WHILE ~cont DO BEGIN

                 IF tries GT 5 THEN BEGIN
                    doQuit = 1B
                    BREAK
                 ENDIF

                 PRINT,(tries EQ 0 ? "Got it." : "K, try again.") + " Like this? '" + read + "' (y/n/q)"
                 read2 = ''
                 READ,read2

                 CASE 1 OF
                    STRMATCH(read2,'y*',/FOLD_CASE) OR STRMATCH(read,'1',/FOLD_CASE): BEGIN
                       defaultFormatStr = [defaultFormatStr,read]
                       cont = 1B
                    END
                    STRMATCH(read2,'n*',/FOLD_CASE) OR STRMATCH(read,'1',/FOLD_CASE): BEGIN
                       tries++
                    END
                    STRMATCH(read2,'q*',/FOLD_CASE): BEGIN
                       doQuit = 1B
                    END
                 ENDCASE
              ENDWHILE
           END
        ENDCASE

        IF KEYWORD_SET(doQuit) THEN BREAK

     ENDFOR

     ;;Now output names
     defaultHeader   = STRJOIN(defaultOutNames,",")
     confirmedHeader = 0B
     PRINT,"Default header for .csv file: " + defaultHeader

     cont = 0B
     tries = 0
     doQuit = 0B
     WHILE ~cont DO BEGIN

        IF tries GT 5 THEN BEGIN
           doQuit = 1B
           BREAK
        ENDIF

        PRINT,(tries EQ 0 ? "You likee?" : "Now really, tell me.") + " (y/n/q)"
        read2 = ''
        READ,read2

        CASE 1 OF
           STRMATCH(read2,'y*',/FOLD_CASE) OR STRMATCH(read2,'1',/FOLD_CASE): BEGIN
              header          = defaultHeader
              confirmedHeader = 1B
              cont = 1B
           END
           STRMATCH(read2,'n*',/FOLD_CASE) OR STRMATCH(read2,'1',/FOLD_CASE): BEGIN
              enterTheVortex = 1B
           END
           STRMATCH(read2,'q*',/FOLD_CASE): BEGIN
              doQuit = 1B
           END
        ENDCASE

        IF KEYWORD_SET(doQuit) THEN BREAK

     ENDWHILE

     IF KEYWORD_SET(enterTheVortex) THEN BEGIN

        PRINT,"K, enter a comma-delimited header string (you need " + STRCOMPRESS(nOut,/REMOVE_ALL) + " entries)"

        cont = 0B
        tries = 0
        doQuit = 0B
        WHILE ~cont DO BEGIN

           IF tries GT 5 THEN BEGIN
              doQuit = 1B
              BREAK
           ENDIF

           read = ''
           READ,read

           IF N_ELEMENTS(STRSPLIT(read,',')) EQ nOut THEN BEGIN
              PRINT,"K, got it"
              cont            = 1
              header          = read
              confirmedHeader = 1B
           ENDIF ELSE BEGIN
              PRINT,"Nope, that was " + STRCOMPRESS(N_ELEMENTS(STRSPLIT(read,',')),/REMOVE_ALL) + " elements."
              tries++
           ENDELSE

           IF KEYWORD_SET(doQuit) THEN BREAK
           
        ENDWHILE

     ENDIF

     formatString = '(' + STRJOIN(defaultFormatStr,',",",') + ')'
     PRINT,"FormatString will be " + formatString
     PRINT,"Header will be " + header

     PRINT,(tries EQ 0 ? "Can I pull the trigger?" : "Now really, tell me.") + " (y/n/q)"
     read2 = ''
     READ,read2

     CASE 1 OF
        STRMATCH(read2,'y*',/FOLD_CASE) OR STRMATCH(read2,'1',/FOLD_CASE): BEGIN
           cont      = 1B
           confirmed = 1
        END
        STRMATCH(read2,'n*',/FOLD_CASE) OR STRMATCH(read2,'1',/FOLD_CASE): BEGIN
           tries++
        END
        STRMATCH(read2,'q*',/FOLD_CASE): BEGIN
           doQuit = 1B
        END
     ENDCASE

     IF KEYWORD_SET(doQuit) THEN BREAK
           
  ENDWHILE

  IF KEYWORD_SET(doQuit) THEN RETURN
  
  nHere = N_ELEMENTS(outMe.(0))
  PRINT,"Printing " + STRCOMPRESS(nHere,/REMOVE_ALL) + " elements to file ..."

  structString = ''
  structString = STRING(FORMAT='(20("(outMe.(",I0,"))[k]",:,","))',INDGEN(nOut))

  OPENW,outLun,CSV_dir+fileName,/GET_LUN
  
  PRINTF,outLun,header
  execStr = 'FOR k=0,nHere-1 DO PRINTF,outLun,FORMAT=formatString,' + structString
  IF ~EXECUTE(execStr) THEN STOP

  CLOSE,outLun
  FREE_LUN,outLun

  PRINT,"Outputted " + CSV_dir+fileName

  STOP
END
PRO READ_TAG_SELECTION,theseTags,uniqTags,uniqTag_i,nnTags,level,levelSelName,levelSelNum,nTags,allTags,splitTagsList

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
PRO JOURNAL__20170303__OUTPUT_COORDINATE_CONVERSIONS_TO_CSV, $
   FILEPREF=filePref, $
   CSV_DIR=CSV_dir

  COMPILE_OPT IDL2

  outDir              = '/SPENCEdata/Research/database/RENU2/'

  inFile              = 'RENU2_GPS.sav'
  timeStrFile         = 'RENU2_GPS--timeStr.sav'
  outFile             = 'RENU2_coordinates.sav'
  
  RESTORE,outDir+inFile
  RESTORE,outDir+outFile

  IF N_ELEMENTS(filePref) EQ 0 THEN BEGIN
     filePref         = 'RENU2_'
  ENDIF

  IF N_ELEMENTS(CSV_dir) EQ 0 THEN BEGIN
     CSV_dir          = outDir + 'csv/'
  ENDIF

  allTags             = TAG_NAMES_R(coords)
  nTags               = N_ELEMENTS(allTags)
  splitTagsList       = STRSPLIT(allTags,'.',/EXTRACT)

  level               = 0
  levelArr            = !NULL
  done                = 0B
  levelSelNum         = !NULL
  levelSelName        = !NULL
  maxTries            = 5
  WHILE ~done DO BEGIN

     ;;Get all unique tag names on this level 
     theseTags        = !NULL
     depths           = !NULL
     ;; nnTags        = 0
     FOR k=0,nTags-1 DO BEGIN
        IF N_ELEMENTS(splitTagsList[k]) GT level THEN BEGIN
           theseTags  = [theseTags,splitTagsList[k,level]]
           depths     = [depths,N_ELEMENTS(splitTagsList[k])]
        ENDIF
     ENDFOR

     ;;Get uniquers
     uniqTag_i  = UNIQ(theseTags,SORT(theseTags))
     uniqTags   = theseTags[uniqTag_i]
     uniqDepths = depths[uniqTag_i]
     nnTags     = N_ELEMENTS(uniqTags)

     PRINT,uniqDepths
     CASE 1 OF
        (nnTags EQ 1) OR ((MAX(uniqDepths)-level) LE 2): BEGIN

           ;;One last read, if need be
           IF nnTags GT 1 THEN BEGIN
              READ_TAG_SELECTION, $
                 theseTags,uniqTags,uniqTag_i,nnTags,level,levelSelName,levelSelNum,nTags,allTags,splitTagsList
           ENDIF ELSE BEGIN
              levelSelName = [levelSelName,uniqTags[0]]
              levelSelNum  = [levelSelNum,uniqTag_i[0]]
           ENDELSE
              
           read  = ''
           cont  = 0B
           tries = 0B
           WHILE ~cont DO BEGIN

              ;; levelSelName = [levelSelName,uniqTags[0]]
              ;; levelSelNum  = [levelSelNum,uniqTag_i[0]]

              ;;Nous sommes arrivés
              PRINT,levelSelName
              PRINT,levelSelNum
              PRINT,"Look OK? (y/n/q)"
              READ,read

              CASE 1 OF
                 STRMATCH(read,'y*',/FOLD_CASE) OR STRMATCH(read,'1',/FOLD_CASE): BEGIN

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


                    ;;Get struct indices
                    levelTags           = TAG_NAMES(coords)
                    FOR k=0,N_ELEMENTS(levelSelName)-1 DO BEGIN
                       newLevel = WHERE(STRMATCH(levelTags,'*'+levelSelName[k]+'*',/FOLD_CASE),nNewLevel)

                       IF N_ELEMENTS(nNewLevel) GT 0 THEN BEGIN
                          newLevel = WHERE(STRMATCH(levelTags,levelSelName[k],/FOLD_CASE),nNewLevel)

                          IF nNewLevel EQ 0 THEN BEGIN
                             PRINT,"Can't figure this one out"
                             STOP
                          ENDIF
                       ENDIF
                       
                       CASE nNewLevel OF
                          0: BEGIN
                             PRINT,"FLUKE"
                             STOP
                          END
                          1: BEGIN
                             levelArr = [levelArr,newLevel[0]]
                          END
                          ELSE: BEGIN
                             PRINT,"Qu'est-ce que tu veux dire?"
                             STOP
                          END
                       ENDCASE

                       ;;Get tags on this level
                       execStr   = STRING(FORMAT='("levelTags = TAG_NAMES(coords",10(:,".(",I0,")"))',levelArr)+')'
                       IF ~EXECUTE(execStr) THEN STOP

                    ENDFOR

                    execStr = 'outMe = coords'
                    FOR k=0,N_ELEMENTS(levelArr)-1 DO BEGIN
                       execStr += STRING(FORMAT='(".(levelArr[",I0,"])")',k)
                    ENDFOR

                    IF ~EXECUTE(execStr) THEN STOP

                    fileName = filePref + STRJOIN(levelSelName,'_') + '.csv'
                    PRINT_CSV,outMe,mem,levelSelName,fileName,CSV_dir,doQuit

                    STOP
                    
                    PRINT,"Do another?"

                    READ,read
                    CASE 1 OF
                       STRMATCH(read,'y*',/FOLD_CASE) OR STRMATCH(read,'1',/FOLD_CASE): BEGIN

                          allTags         = TAG_NAMES_R(coords)
                          nTags           = N_ELEMENTS(allTags)
                          splitTagsList   = STRSPLIT(allTags,'.',/EXTRACT)

                          level           = 0
                          levelSelName    = !NULL
                          levelSelNum     = !NULL
                          cont            = 1B

                       END
                       STRMATCH(read,'n*',/FOLD_CASE) OR STRMATCH(read,'0',/FOLD_CASE): BEGIN
                          
                          cont            = 1B
                          doQuit          = 1B

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

           READ_TAG_SELECTION,theseTags,uniqTags,uniqTag_i,nnTags,level,levelSelName,levelSelNum,nTags,allTags,splitTagsList

        END
     ENDCASE

     IF KEYWORD_SET(doQuit) THEN BREAK

  ENDWHILE

END
