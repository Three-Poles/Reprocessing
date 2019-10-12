!  f2py -c -m verify_bss_raw_qmapped_dressed verify_bss_raw_qmapped_dressed.f90 read_prob_forecasts_raw_qmapped_dressed.f90 read_prob_forecasts_raw_qmapped_dressed_mme.f90 ran3.f sort.f check.f90 -I/usr/local/gfortran/include -L/usr/local/gfortran/lib  -L/opt/local/lib -lnetcdff

SUBROUTINE verify_bss_raw_qmapped_dressed(cleade, cmodel, cftype, rthresh, &
    date_list_anal, apcp_anal_t, validdays, nxa, nya, ndates, &
    bss_raw, bss_qmapped, bss_dressed, &
    bs_raw_daily, bs_qmapped_daily, bs_dressed_daily, bs_climo_daily)

! purpose:  generate reliability information, frequency of usage, and Brier
!  skill scores for raw, quantile-mapped, and dressed data for input model

PARAMETER (nthreshes = 7)

CHARACTER*3, INTENT(IN) :: cleade
CHARACTER*5, INTENT(IN) :: cmodel
CHARACTER*(*), INTENT(IN) :: cftype
INTEGER, INTENT(IN) :: nxa, nya, ndates
REAL, INTENT(IN) :: rthresh ! precip threshold amount
REAL, DIMENSION(nxa,nya,ndates), INTENT(IN) :: apcp_anal_t
CHARACTER*10, DIMENSION(ndates), INTENT(IN) :: date_list_anal
INTEGER, INTENT(INOUT), DIMENSION(ndates) :: validdays

REAL, INTENT(OUT) :: bss_raw, bss_qmapped, bss_dressed 
REAL*8, DIMENSION(ndates), INTENT(OUT) :: bs_raw_daily
REAL*8, DIMENSION(ndates), INTENT(OUT)  :: bs_qmapped_daily
REAL*8, DIMENSION(ndates), INTENT(OUT)  :: bs_dressed_daily
REAL*8, DIMENSION(ndates), INTENT(OUT)  :: bs_climo_daily

! f2py intent(in) cleade, cmodel, cftype, rthresh, date_list_anal, apcp_anal_t
! f2py intent(in) nxa, nya, ndates  
! f2py intent(in,out) validdays
! f2py intent(out) bss_raw, bss_raw_daily
! f2py intent(out) bss_qmapped, bss_qmapped_daily
! f2py intent(out) bss_dressed, bss_dressed_daily

! f2py depend(nxa,nya,ndates) apcp_anal_t
! f2py depend(ndates) date_list_anal, validdays
! f2py depend(ndates) bs_raw_daily, bs_qmapped_daily, 
! f2py depend(ndates) bs_dressed_daily, bs_climo_daily

! --- now local variables

INTEGER*2, DIMENSION(nxa,nya) :: conusmask

REAL, DIMENSION(nxa,nya) :: climo_prob
REAL, DIMENSION(nxa,nya) :: rlonsa
REAL, DIMENSION(nxa,nya) :: rlatsa
REAL, DIMENSION(nxa,nya) :: prob_forecast
REAL, DIMENSION(nxa,nya) :: prob_forecast_raw
REAL, DIMENSION(nxa,nya) :: prob_forecast_qmapped
REAL, DIMENSION(nxa,nya) :: prob_forecast_dressed
REAL, DIMENSION(nthreshes) :: pthreshes

REAL*8 :: bs
REAL*8 :: bs_climo

REAL*8 :: bs_raw
REAL*8 :: bs_qmapped
REAL*8 :: bs_dressed

REAL*8, DIMENSION(ndates) :: bs_daily

CHARACTER*120 infile
CHARACTER*10 cyyyymmddhh
CHARACTER*3 clead_use


bs_climo = 0.

bs_raw = 0.
bs_qmapped = 0.
bs_dressed = 0.

bs_raw_daily(:) = 0.0
bs_qmapped_daily(:) = 0.0
bs_dressed_daily(:) = 0.0
bs_climo_daily(:) = 0.0


pid180 = 3.1415926/180.
clead_use = cleade
IF (clead_use .eq. '012') clead_use = '12'
IF (clead_use .eq. '024') clead_use = '24'
IF (clead_use .eq. '036') clead_use = '36'
IF (clead_use .eq. '048') clead_use = '48'
IF (clead_use .eq. '060') clead_use = '60'
IF (clead_use .eq. '072') clead_use = '72'
IF (clead_use .eq. '084') clead_use = '84'
IF (clead_use .eq. '096') clead_use = '96'

PRINT *,'max(apcp_anal_t) = ', maxval(apcp_anal_t)

! ---- loop thru days and verify

DO idate = 1, ndates
!DO idate = 12, 12

    cyyyymmddhh = date_list_anal(idate) 

    ! ---- read in the forecast for this particular initial date/time and forecast lead and resolution

    IF (TRIM(cmodel) .eq. 'MME' .or. TRIM(cmodel) .eq. 'mme') THEN
        CALL read_prob_forecasts_raw_qmapped_dressed_mme(nxa, nya, cyyyymmddhh, cftype, &
            clead_use, rthresh, nthreshes, pthreshes, conusmask, rlonsa, rlatsa, &
            prob_forecast_raw, prob_forecast_qmapped, prob_forecast_dressed, climo_prob)
    ELSE
        CALL read_prob_forecasts_raw_qmapped_dressed(nxa, nya, cyyyymmddhh, cmodel, cftype, &
            clead_use, rthresh, nthreshes, pthreshes, conusmask, rlonsa, rlatsa, &
            prob_forecast_raw, prob_forecast_qmapped, prob_forecast_dressed, climo_prob)
    ENDIF

    pmax_raw = MAXVAL(prob_forecast_raw)
    pmax_qmapped = MAXVAL(prob_forecast_qmapped)
    pmax_dressed = MAXVAL(prob_forecast_dressed)
    
    ! ---- only verify if there is good data for all three models.
    
    contab = 0.
    bs = 0.
    bs_daily(idate) = 0.
    bs_climo_daily(idate) = 0.
    
    IF (pmax_raw .lt. 0.0 .or. pmax_qmapped .lt. 0.0 .or. pmax_dressed .lt. 0.0) THEN
        PRINT *,'*********   Some missing data for ',cyyyymmddhh
        PRINT *,'            pmax_raw, _qmapped, _dressed = ', &
            pmax_raw, pmax_qmapped, pmax_dressed
    ENDIF
    
    IF (pmax_raw .ge. 0.0 .and. pmax_qmapped .ge. 0.0 .and. pmax_dressed .ge. 0.0 &
    .and. validdays(idate) .eq. 1) THEN

	    DO itype = 1, 3
		    contab = 0.
		    bs = 0.
            bs_daily(idate) = 0.
            
		    IF (itype .eq. 1) THEN ! MME
			    prob_forecast(:,:) = prob_forecast_raw
		    ELSE IF (itype .eq. 2) THEN
                prob_forecast(:,:) = prob_forecast_qmapped
		    ELSE 
                prob_forecast(:,:) = prob_forecast_dressed
            ENDIF
	   
	   	    ! ---- now let's verify probability forecast

		    DO i = 1, nxa
          	    DO j = 1, nya
             	    IF (conusmask(i,j) .gt. 0 .and. prob_forecast(i,j) .GE. 0.0 .and. &
             	    prob_forecast(i,j) .le. 1.0 .and. apcp_anal_t(i,j,idate) .GE. 0.0) THEN
            	 	    cfac = cos(rlatsa(i,j)*pid180)  ! cos of latitude to acct for unequal grid box size
            	 	    pclimo = climo_prob(i,j)
            	 	    p      = prob_forecast(i,j)
            	 	    ipcat  = nint(p*20)
            	 	    v      = apcp_anal_t(i,j,idate)
				 
            	 	    IF (v .GE. rthresh) THEN
                    	    bs = bs + cfac*(1.-p)**2
                            bs_daily(idate) = bs_daily(idate) + cfac*(1.-p)**2
                    	    IF (itype .eq. 1) THEN
                                bs_climo = bs_climo + cfac*(1.-pclimo)**2
                                bs_climo_daily(idate) = bs_climo_daily(idate) + cfac*(1.-pclimo)**2
                            ENDIF
            	 	    ELSE
                    	    bs = bs + cfac * p**2
                            bs_daily(idate) = bs_daily(idate) + cfac * p**2
                    	    IF (itype .eq. 1) THEN
                                bs_climo = bs_climo + cfac*pclimo**2
                                bs_climo_daily(idate) = bs_climo_daily(idate) + cfac * pclimo**2
                            ENDIF
                 	    ENDIF
             	    ENDIF ! conusmask
      	  	    END DO   ! j = 1, nya
   	   	    END DO      ! i = 1, nxa
	   
	   	    IF (itype .eq. 1) THEN
		   	    bs_raw = bs_raw + bs
                bs_raw_daily(idate) = bs_daily(idate)
	   	    ELSE IF (itype .eq. 2) THEN
		   	    bs_qmapped = bs_qmapped + bs
                bs_qmapped_daily(idate) = bs_daily(idate)
	        ELSE IF (itype .eq. 3) THEN
		   	    bs_dressed = bs_dressed + bs
                bs_dressed_daily(idate) = bs_daily(idate)
	        ENDIF 
	   
   	    END DO  ! itype
        
    ELSE

        
        validdays(idate) = 0
        bs_climo_daily(idate) = 0.
        bs_raw_daily(idate) = 0.
        bs_qmapped_daily(idate) = 0.
        bs_dressed_daily(idate) = 0.
    ENDIF !(iex)
    !PRINT *,'  bs climo raw qmapped dressed = ', bs_climo_daily(idate), &
    !    bs_raw_daily(idate), bs_qmapped_daily(idate), bs_dressed_daily(idate)   
  
END DO ! idate

DO itype = 1, 3
	
	IF (itype .eq. 1) THEN
		PRINT *,'Raw'
		bs = bs_raw
	ELSE IF (itype .eq. 2) THEN
		PRINT *, 'Quantile-mapped'
		bs = bs_qmapped
	ELSE IF (itype .eq. 3) THEN
		print *, 'Dressed'
		bs = bs_dressed
    ENDIF
	
	! ---- with tallied contingency tables, now set mean reliability and frequency of use for mean

	bss = 1. - bs / bs_climo
	!PRINT *,'  bss = ',bss

	! ---- copy to output arrays
		
	IF (itype .eq. 1) THEN
		bss_raw = bss
	ELSE IF (itype .eq. 2) THEN
		bss_qmapped = bss
	ELSE IF (itype .eq. 3) THEN
		bss_dressed = bss
	ENDIF
	
END DO  ! itype

RETURN
END SUBROUTINE verify_bss_raw_qmapped_dressed

