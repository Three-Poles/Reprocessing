      SUBROUTINE UPDAT(JDATE,KHR,MDATE) 
C 
C        FEBRUARY 1993   GLAHN   TDL   MOS-2000
C        JULY     1997   GLAHN   SKIP CALCULATIONS WHEN KHR = 0.
C
C        PURPOSE 
C            TO ADD KHR HOURS TO DATE/TIME IN JDATE, AND TO
C            PROVIDE THE UPDATED DATE/TIME IN MDATE. 
C            JDATE AND MDATE ARE OF FORMAT YYYYMMDDHH.
C            ADAPTED FROM HP CHGDT. 
C 
C        DATA SET USE 
C            NONE. 
C 
C        VARIABLES 
C 
C            INPUT 
C               JDATE = BASIC DATE OF FORMAT YYYYMMDDHH.
C                 KHR = HOURS TO ADD TO BASIC DATE, CAN BE 
C                       POSITIVE OR NEGATIVE. 
C 
C            OUTPUT 
C               MDATE = UPDATED DATE IN SAME FORMAT AS JDATE
C 
C            INTERNAL
C                  JY = YEAR (E.G., 1993). 
C                  JM = MONTH (E.G., 12). 
C                  JD = DAY (E.G., 01). 
C                  JH = HOUR (E.G., 08). 
C            MTEST(J) = LAST DAY OF EACH MONTH (J=1,12). 
C 
C        NONSYSTEM SUBROUTINES CALLED 
C            NONE. 
C 
      DIMENSION MTEST(12)
      DATA MTEST/31,28,31,30,31,30,31,31,30,31,30,31/ 
C 
      IF(KHR.EQ.0)THEN
         MDATE=JDATE
         GO TO 260
      ENDIF
C
C        BREAK DATE INTO COMPONENT PARTS. 
C 
      JY=JDATE/1000000
      JM=JDATE/10000-JY*100
      JD=JDATE/100-JY*10000-JM*100
      JH=MOD(JDATE,100)+KHR
C 
C        HANDLE LEAP YEAR, WHICH INCLUDES YEAR 2000
C        BUT NOT OTHER CENTENIAL YEARS. 
C 
 110  MTEST(2)=28 
      MFEB=MOD(JY,4) 
      IF(MFEB.EQ.0)MTEST(2)=29 
C 
C        CHANGE HOURS IF NECESSARY. 
C 
 120  IF(JH-23)140,250,130 
 130  JH=JH-24 
      JD=JD+1 
      IF(JD.LE.MTEST(JM))GO TO 120 
      JD=1 
      JM=JM+1 
      IF(JM.LE.12)GO TO 120 
      JY=JY+1 
      JM=1 
      GO TO 110 
C 
 140  IF(JH.GE.0)GO TO 250 
      JH=JH+24 
      JD=JD-1 
      IF(JD.GT.0)GO TO 120 
      JM=JM-1 
      IF(JM.LE.0)GO TO 170 
      JD=MTEST(JM) 
C        CHECK FOR LEAP YEAR. 
      GO TO 110 
C 
 170  JM=12 
      JD=31 
      JY=JY-1 
      GO TO 120 
C 
 250  MDATE=JY*1000000+JM*10000+JD*100+JH 
 260  RETURN 
      END 
