!  f2py -c -m verify_relia_bss_qmapped verify_relia_bss_qmapped.f90 read_prob_forecasts_raw_qmapped_dressed.f90 ran3.f sort.f check.f90 -I/usr/local/gfortran/include -L/usr/local/gfortran/lib  -L/opt/local/lib -lnetcdff

SUBROUTINE verify_relia_bss_qmapped(cleade, cmodel, nclasses, rthresh, &
    date_list_anal, apcp_anal_t, nxa, nya, ndates, &
    relia_dressed, relia_dressed_05, relia_dressed_95, frequse_dressed, bss_dressed, &
    relia_raw, relia_raw_05, relia_raw_95, frequse_raw, bss_raw, &
    relia_qmapped, relia_qmapped_05, relia_qmapped_95, frequse_qmapped, bss_qmapped)

! purpose:  generate reliability information, frequency of usage, and Brier
!  skill scores for post-processed NCEP, CMC, ECMWF, and MME forecasts

PARAMETER (nresa = 1000)
PARAMETER (nthreshes = 7)
CHARACTER*3, INTENT(IN) :: cleade
CHARACTER*(*), INTENT(IN) :: cmodel
INTEGER, INTENT(IN) :: nclasses, nxa, nya, ndates
REAL, INTENT(IN) :: rthresh ! precip threshold amount
REAL, DIMENSION(nxa,nya,ndates), INTENT(IN) :: apcp_anal_t
CHARACTER*10, DIMENSION(ndates), INTENT(IN) :: date_list_anal

REAL, DIMENSION(nclasses), INTENT(OUT) :: relia_dressed, relia_dressed_05, relia_dressed_95
REAL, DIMENSION(nclasses), INTENT(OUT) :: relia_raw, relia_raw_05, relia_raw_95
REAL, DIMENSION(nclasses), INTENT(OUT) :: relia_qmapped, relia_qmapped_05, relia_qmapped_95	
REAL, DIMENSION(nclasses), INTENT(OUT) :: frequse_dressed, frequse_raw, frequse_qmapped
    
REAL, INTENT(OUT) :: bss_dressed, bss_raw, bss_qmapped
REAL*4, DIMENSION(ndates) :: bss_dressed_daily
REAL*4, DIMENSION(ndates) :: bss_raw_daily
REAL*4, DIMENSION(ndates) :: bss_qmapped_daily



! f2py intent(in) nxa, nya, ndates, cleade, cmodel, nclasses 
! f2py intent(in) rthresh, apcp_anal_t, date_list_anal
! f2py depend(ndates) date_list_anal
! f2py depend(nxa,nya,ndates) apcp_anal_t
! f2py intent(out) relia_dressed, relia_dressed_05, relia_dressed_95
! f2py intent(out) relia_raw, relia_raw_05, relia_raw_95
! f2py intent(out) relia_qmapped, relia_qmapped_05, relia_qmapped_95
! f2py intent(out) frequse_dressed, frequse_raw, frequse_qmapped
! f2py intent(out) bss_dressed, bss_raw, bss_qmapped
! f2py depend(nclasses) relia_raw, relia_raw_05, relia_raw_95
! f2py depend(nclasses) relia_dressed, relia_dressed_05, relia_dressed_95
! f2py depend(nclasses) relia_qmapped, relia_qmapped_05, relia_qmapped_95
! f2py depend(nclasses) frequse_dressed, frequse_raw, frequse_qmapped
! f2py depend(nthreshes) bss_dressed, bss_raw, bss_qmapped

! --- now local variables

INTEGER*2, DIMENSION(nxa,nya) :: conusmask

REAL, DIMENSION(nxa,nya) :: climo_prob
REAL, DIMENSION(nxa,nya) :: rlonsa
REAL, DIMENSION(nxa,nya) :: rlatsa
REAL, DIMENSION(nxa,nya) :: prob_forecast_dressed
REAL, DIMENSION(nxa,nya) :: prob_forecast_raw
REAL, DIMENSION(nxa,nya) :: prob_forecast_qmapped

REAL, DIMENSION(nthreshes) :: pthreshes

REAL*8, DIMENSION(0:nclasses-1,2) :: contab
REAL*8, DIMENSION(0:nclasses-1,2) :: contab_sum

REAL*8, DIMENSION(0:nclasses-1,2) :: contab_raw
REAL*8, DIMENSION(0:nclasses-1,2) :: contab_qmapped
REAL*8, DIMENSION(0:nclasses-1,2) :: contab_dressed

REAL*8, DIMENSION(ndates, 0:nclasses-1,2) :: contab_daily
REAL*8, DIMENSION(ndates, 0:nclasses-1,2) :: contab_raw_daily
REAL*8, DIMENSION(ndates, 0:nclasses-1,2) :: contab_dressed_daily
REAL*8, DIMENSION(ndates, 0:nclasses-1,2) :: contab_qmapped_daily

REAL*8, DIMENSION(0:nclasses-1) :: relia, relia_05, relia_95
REAL*8, DIMENSION(nresa, 0:nclasses-1) :: relia_resa
REAL*8, DIMENSION(0:nclasses-1) :: frequse

REAL*8 :: bs
REAL*8 :: bs_climo

REAL*8 :: bs_raw
REAL*8 :: bs_qmapped
REAL*8 :: bs_dressed

REAL*8, DIMENSION(ndates) :: bs_daily
REAL*8, DIMENSION(ndates) :: bs_climo_daily
REAL*8, DIMENSION(ndates) :: bs_raw_daily
REAL*8, DIMENSION(ndates) :: bs_qmapped_daily
REAL*8, DIMENSION(ndates) :: bs_dressed_daily

REAL, DIMENSION(nresa) :: rsamps

CHARACTER*120 infile
CHARACTER*10 cyyyymmddhh
CHARACTER*3 clead_use

contab_raw = 0.
contab_qmapped = 0.
contab_dressed = 0.

contab_raw_daily = 0.
contab_dressed_daily = 0.
contab_qmapped_daily = 0.

bs_climo = 0.

bs_raw = 0.
bs_qmapped = 0.
bs_dressed = 0.

bs_raw_daily(:) = 0.0
bs_dressed_daily(:) = 0.0
bs_qmapped_daily(:) = 0.0

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

    cyyyymmddhh = date_list_anal(idate) 

    ! ---- read in the forecast for this particular initial date/time and forecast lead and resolution

    CALL read_prob_forecasts_raw_qmapped_dressed(nxa, nya, cyyyymmddhh, cmodel, &
        clead_use, rthresh, nthreshes, pthreshes, conusmask, rlonsa, rlatsa, &
        prob_forecast_raw, prob_forecast_qmapped, prob_forecast_dressed, climo_prob)

    pmax_raw = MAXVAL(prob_forecast_raw)
    pmax_qmapped = MAXVAL(prob_forecast_qmapped)
    pmax_dressed = MAXVAL(prob_forecast_dressed)
    
    ! ---- only verify if there is good data for all three models.
    
    IF (pmax_raw .ge. 0.0 .and. pmax_qmapped .ge. 0.0 .and. pmax_dressed .ge. 0.0) THEN

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
               	    	    contab(ipcat,2) = contab(ipcat,2) + cfac
                    	    bs = bs + cfac*(1.-p)**2
                            bs_daily(idate) = bs_daily(idate) + cfac*(1.-p)**2
                    	    IF (itype .eq. 1) THEN
                                bs_climo = bs_climo + cfac*(1.-pclimo)**2
                                bs_climo_daily(idate) = bs_climo_daily(idate) + cfac*(1.-pclimo)**2
                            ENDIF
            	 	    ELSE
                    	    contab(ipcat,1) = contab(ipcat,1) + cfac
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
		   	    contab_raw = contab_raw + contab
			    contab_raw_daily(idate,:,:) = contab(:,:)
		   	    bs_raw = bs_raw + bs
                bs_raw_daily(idate) = bs_daily(idate)
	   	    ELSE IF (itype .eq. 2) THEN
		   	    contab_qmapped = contab_qmapped + contab
			    contab_qmapped_daily(idate,:,:) = contab(:,:)
		   	    bs_qmapped = bs_qmapped + bs
                bs_qmapped_daily(idate) = bs_daily(idate)
	        ELSE IF (itype .eq. 3) THEN
		   	    contab_dressed = contab_dressed + contab
			    contab_dressed_daily(idate,:,:) = contab(:,:)
		   	    bs_dressed = bs_dressed + bs
                bs_dressed_daily(idate) = bs_daily(idate)
	        ENDIF 
	   
   	    END DO  ! itype
        
    ELSE
        bs_climo_daily(idate) = 0.
        bs_raw_daily(idate) = 0.
        bs_qmapped_daily(idate) = 0.
        bs_dressed_daily(idate) = 0.
    ENDIF !(iex)
  
END DO ! idate

DO itype = 1, 3
	
	IF (itype .eq. 1) THEN
		PRINT *,'Raw'
		contab = contab_raw
		contab_daily = contab_raw_daily
		bs = bs_NCEP
	ELSE IF (itype .eq. 2) THEN
		PRINT *, 'Quantile-mapped'
		contab = contab_qmapped
		contab_daily = contab_qmapped_daily
		bs = bs_qmapped
	ELSE IF (itype .eq. 3) THEN
		print *, 'Dressed'
		contab = contab_dressed
		contab_daily = contab_dressed_daily
		bs = bs_dressed
    ENDIF
	
	! ---- with tallied contingency tables, now set mean reliability and frequency of use for mean

	!print *,'bs, bs_climo = ',bs, bs_climo
	ctot = SUM(contab)
	bss = 1. - bs / bs_climo
	relia(:) = -99.9999
	PRINT *,'  bss = ',bss
	PRINT *,'  p   reliability   freq of usage'
	DO icat = 0,20
		frequse(icat) = (contab(icat,2)  + contab(icat,1)) / ctot
		IF ((contab(icat,1) + contab(icat,2)) .gt. 0) THEN
			relia(icat) = contab(icat,2) / (contab(icat,1) + contab(icat,2))
			PRINT 203,float(icat*5),relia(icat)*100.,frequse(icat)
		ELSE
			PRINT 203,float(icat*5),relia(icat),frequse(icat)
		ENDIF
  	    203 format(f5.1,3x,2(f8.3,3x))
	END DO  !icat
	
	! ---- perform a resampling to generate the confidence intervals for reliability
	!      sampling the days with replacement
	
	idum = -12345
	relia_resa(:,:) = -99.9999
	DO iresa = 1, nresa
        PRINT *,'iresa = ',iresa
		contab_sum = 0.
		DO idate = 1, ndates
			2345 cran = ran3(idum)
			idate2 = MIN(1+NINT(cran*REAL(ndates)),ndates)
            PRINT *,'cran,idate2 = ',cran,idate2
            IF (bs_MME_daily(idate2) .gt. 0.0) THEN
			    contab_sum(:,:) = contab_sum(:,:) + contab_daily(idate2,:,:)
            ELSE
                GOTO 2345
            ENDIF
		END DO
		DO icat = 0,20
			IF ((contab_sum(icat,1) + contab_sum(icat,2)) .gt. 0) THEN
				relia_resa(iresa,icat) = contab_sum(icat,2) / (contab_sum(icat,1) + contab_sum(icat,2))
			ENDIF
		END DO ! icat
	END DO  !iresa
	
    ! ---- now find the 5th and 95th percentiles of the resampled distribution, accounting
	!      for the possibility that at high probabilities there may be no samples in some 
	!      cases.

	DO i = 0, nclasses-1
		rsamps(:) = relia_resa(:,i)
		CALL sort(nresa, rsamps)
		ibegin = 1
		DO iresa = 1, nresa
			IF (rsamps(iresa) .gt. -99.) THEN
				ibegin = iresa
				goto 3212
			ENDIF
		END DO ! iresa
3212	isamp = MIN(ibegin + NINT(0.05* (nresa-ibegin+1) ),nresa)
		relia_05(i) = rsamps(isamp)
		isamp = MIN(ibegin + NINT(0.95* (nresa-ibegin+1) ),nresa)
		relia_95(i) = rsamps(isamp)
	END DO ! i 
	
	! ---- copy to output arrays
		
	IF (itype .eq. 1) THEN
		relia_raw = relia
		relia_raw_05 = relia_05
		relia_raw_95 = relia_95
		frequse_raw = frequse
		bss_raw = bss
	ELSE IF (itype .eq. 2) THEN
		relia_qmapped = relia
		relia_qmapped_05 = relia_05
		relia_qmapped_95 = relia_95
		frequse_qmapped = frequse
		bss_qmapped = bss
	ELSE IF (itype .eq. 3) THEN
		relia_dressed = relia
		relia_dressed_05 = relia_05
		relia_dressed_95 = relia_95
		frequse_dressed = frequse
		bss_dressed = bss
	ENDIF
	
END DO  ! itype

RETURN
END SUBROUTINE verify_relia_bss_mme

