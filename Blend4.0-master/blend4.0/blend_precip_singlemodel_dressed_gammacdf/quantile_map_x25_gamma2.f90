
SUBROUTINE quantile_map_x25_gamma2(nxa, nya, nstride, nstencil, &
    nmultiply, nens_qmap, iadim, ipdim, conusmask, gamma_shape_qmap_forecast, &
    gamma_scale_qmap_forecast, fraction_zero_qmap_forecast, &
    gamma_shape_qmap_analysis, gamma_scale_qmap_analysis, &
    fraction_zero_qmap_analysis, forecast, &
    quantile_table, alpha_values, cumprob_values, &
    forecast_x25)
    
    ! blend version

! --- Perform the quantile mapping; here we use 25 (24 + original) 
!     surrounding grid points forecasts as well to increase sample  
!     size and account for position error

INTEGER, INTENT(IN) :: nxa, nya, nstride, nstencil, nmultiply, nens_qmap
INTEGER, INTENT(IN) :: iadim, ipdim
INTEGER*2, INTENT(IN), DIMENSION(nxa,nya) :: conusmask
REAL, INTENT(IN), DIMENSION(nxa, nya) :: gamma_shape_qmap_forecast
REAL, INTENT(IN), DIMENSION(nxa, nya) :: gamma_scale_qmap_forecast
REAL, INTENT(IN), DIMENSION(nxa, nya) :: fraction_zero_qmap_forecast
REAL, INTENT(IN), DIMENSION(nxa, nya) :: gamma_shape_qmap_analysis
REAL, INTENT(IN), DIMENSION(nxa, nya) :: gamma_scale_qmap_analysis
REAL, INTENT(IN), DIMENSION(nxa, nya) :: fraction_zero_qmap_analysis

REAL*8, INTENT(IN), DIMENSION(iadim, ipdim) :: quantile_table
REAL*8, INTENT(IN), DIMENSION(iadim) :: alpha_values
REAL*8, INTENT(IN), DIMENSION(ipdim) :: cumprob_values
REAL, INTENT(IN), DIMENSION(nxa,nya) :: forecast

REAL, INTENT(OUT), DIMENSION(nstencil,nxa,nya) :: forecast_x25

REAL :: weight, weight2, rmean
DOUBLE PRECISION alpha_fcst, beta_fcst, fz_fcst
DOUBLE PRECISION alpha_anal, beta_anal, fz_anal
DOUBLE PRECISION fcst
DOUBLE PRECISION precip_qmapped
DOUBLE PRECISION qgamma
INTEGER :: idum,itest,jtest, jyn, ixn

REAL, DIMENSION(5) :: a95_to_a99, f95_to_f99

rmean = 1000.*SUM(forecast)/REAL(nxa*nya)
idum = -1*INT(rmean)
forecast_x25(:,:,:) = -99.99 ! set to missing

! ---- Not for points inside conus mask, apply the procedure of doing quantile
!      mapping using the forecast at (i,j), but also surrounding locations.
!      also add random number to the input forecast quantile to also account for the
!      tendency of the ensemble to be over-certain of its amount.

weight = 0.0
DO jya = 1, nya
	DO ixa = 1, nxa
!DO jya = nya/2, nya/2
!    DO ixa = 388, 388

        !PRINT *,'processing ixa, jya, conusmask ',ixa,jya,conusmask(ixa,jya)
		IF (conusmask(ixa,jya) .le. 0) THEN
         	! If outside the conus, then make the output forecast array simply replicates
         	! of the input forecast array
         	forecast_x25(:,ixa,jya) = forecast(ixa,jya)
      	ELSE
        
         	ktr = 0

         	! ---- Loop thru the nstencil*nstencil grid points with (ixa,jya) in center

         	DO jyn = jya-nmultiply*nstride, jya+nmultiply*nstride, nstride
            	DO ixn = ixa-nmultiply*nstride, ixa+nmultiply*nstride, nstride
					
                    ! ---- deal with the possibility that the stencil point is outside 
                    !      domain.  In that case, use the value on the nearest border point
                    
                    IF (ixn .lt. 1) THEN
                        ixs = 1
                    ELSE IF (ixn .gt. nxa) THEN
                        ixs = nxa
                    ELSE
                        ixs = ixn
                    ENDIF
                    IF (jyn .lt. 1) THEN
                        jys = 1
                    ELSE IF (jyn .gt. nya) THEN
                        jys = nya
                    ELSE
                        jys = jyn
                    ENDIF
                    
                    !PRINT *,'ixs, jys = ',ixs,jys, conusmask(ixs,jys)
                    
                    ktr = ktr + 1
                    IF (conusmask(ixs,jys) .gt. 0) THEN

                        fcst = forecast(ixs,jys)
                        !PRINT *,'fcst = ', fcst
                        
						IF (forecast(ixs,jys) .eq. 0.0) THEN
							forecast_x25(ktr,ixa,jya) = 0.0
                        ELSE

                     		! ---- Perform quantile mapping.  First determine the quantile of 
                            !      forecast distribution associated with today's forecast value.
                            
                            alpha_fcst = gamma_shape_qmap_forecast(ixs,jys)
                            beta_fcst = gamma_scale_qmap_forecast(ixs,jys)
                            fz_fcst = fraction_zero_qmap_forecast(ixs,jys)
                            alpha_anal = gamma_shape_qmap_analysis(ixa,jya)
                            beta_anal = gamma_scale_qmap_analysis(ixa,jya)
                            fz_anal = fraction_zero_qmap_analysis(ixa,jya)
                            !PRINT *,'alpha_fcst, beta_fcst, fz_fcst = ', &
                            !    alpha_fcst, beta_fcst, fz_fcst
                            !PRINT *,'alpha_anal, beta_anal, fz_anal = ', &
                            !    alpha_anal, beta_anal, fz_anal  
                            IF (alpha_fcst .ne. alpha_fcst .or. alpha_anal .ne. alpha_anal) THEN
                                forecast_x25(ktr,ixa,jya) = fcst  
                            ELSE   
                                CALL gamma_quantile_map(iadim, ipdim, fcst, alpha_fcst, &
                                    beta_fcst, fz_fcst, alpha_anal, beta_anal, &
                                    fz_anal, quantile_table, alpha_values, cumprob_values, & 
                                    precip_qmapped)       
                                forecast_x25(ktr,ixa,jya) = precip_qmapped
                            ENDIF
                            forecast_x25(ktr,ixa,jya) = precip_qmapped 
						ENDIF
                        !PRINT *,'Forecast before, after = ', fcst, forecast_x25(ktr,ixa,jya)
					
                    ELSE 
                        
                        ! ---- stencil point is either outside the area where quantile mapping possible.
                        !      As a substitute, take the quantile-mapped value at center point and
                        !      add a small amount of random noise.
                        
                        fcst = forecast(ixa,jya) 
                        !print *,'fcst = ',fcst
						IF (forecast(ixa,jya) .eq. 0.0) THEN
							forecast_x25(ktr,ixa,jya) = 0.0
                        ELSE
                        
                            alpha_fcst = gamma_shape_qmap_forecast(ixa,jya)
                            beta_fcst = gamma_scale_qmap_forecast(ixa,jya)
                            fz_fcst = fraction_zero_qmap_forecast(ixa,jya)
                            alpha_anal = gamma_shape_qmap_analysis(ixa,jya)
                            beta_anal = gamma_scale_qmap_analysis(ixa,jya)
                            fz_anal = fraction_zero_qmap_analysis(ixa,jya)
                            
                            !PRINT *,'alpha_fcst, beta_fcst, fz_fcst = ', &
                            !    alpha_fcst, beta_fcst, fz_fcst
                            !PRINT *,'alpha_anal, beta_anal, fz_anal = ', &
                            !    alpha_anal, beta_anal, fz_anal   
                            IF (alpha_fcst .ne. alpha_fcst .or. alpha_anal .ne. alpha_anal) THEN
                                forecast_x25(ktr,ixa,jya) = fcst  
                            ELSE 
                                CALL gamma_quantile_map(iadim, ipdim, fcst, alpha_fcst, &
                                    beta_fcst, fz_fcst, alpha_anal, beta_anal, &
                                    fz_anal, quantile_table, alpha_values, cumprob_values, &
                                    precip_qmapped)
                                forecast_x25(ktr,ixa,jya) = precip_qmapped
                            ENDIF
                            r = 1 + 0.1*(ran1(idum) - 0.5)
                            forecast_x25(ktr,ixa,jya) = forecast_x25(ktr,ixa,jya)*r 
                        ENDIF  
                        
                 	END IF ! jyn>1, jyn < nya, etc.
                    
                    !PRINT *,'forecast_x25(ktr,ixa,jya) = ', forecast_x25(ktr,ixa,jya)
                    
            	END DO ! ixn
         	END DO ! jyn
      	END IF ! conusmask
   	END DO ! ixa
END DO ! jya

RETURN
END SUBROUTINE quantile_map_x25_gamma2

! ==================================================================

SUBROUTINE gamma_quantile_map (iadim, ipdim, fcst, &
    alpha_fcst, beta_fcst, fz_fcst, alpha_anal, beta_anal, &
    fz_anal, quantile_table, alpha_values, cumprob_values, & 
    qmapped_fcst)

INTEGER, INTENT(IN) :: iadim, ipdim
DOUBLE PRECISION, INTENT(IN) :: fcst, alpha_fcst, beta_fcst, fz_fcst
DOUBLE PRECISION, INTENT(IN) :: alpha_anal, beta_anal, fz_anal
REAL*8, INTENT(IN), DIMENSION(iadim, ipdim) :: quantile_table
REAL*8, INTENT(IN), DIMENSION(iadim) :: alpha_values
REAL*8, INTENT(IN), DIMENSION(ipdim) :: cumprob_values
DOUBLE PRECISION, INTENT(OUT) :: qmapped_fcst 
    
DOUBLE PRECISION ksi, cum, ccum, nonexceedance_prob, qgamma
LOGICAL tootrue
LOGICAL toofalse

REAL, DIMENSION(9) :: a90_to_a99, f90_to_f99
    
! ---- determine the cumulative probability of being less than 
!      or equal to the forecast amount using the CDF defined by
!      the forecast gamma distribution parameters and the forecast
!      fraction zero  

tootrue = .TRUE.
toofalse = .FALSE.

ksi = fcst / beta_fcst
CALL cumgam(ksi, alpha_fcst, cum, ccum)
nonexceedance_prob = fz_fcst + (1.-fz_fcst)*cum
!PRINT *,'ksi, fcst, beta_fcst, alpha_fcst, fz_fcst, nonexceedance_prob = ', &
!    ksi, fcst, beta_fcst, alpha_fcst, fz_fcst, nonexceedance_prob

! ---- Next, using Michael Scheuerer's code, use the quantile 
!      function to retrieve the analyzed value associates with
!      the same cumulative probability.

! ---- do not apply quantile mapping with zeros

IF (nonexceedance_prob .gt. fz_anal) THEN
    cum = (nonexceedance_prob - fz_anal) / (1.-fz_anal)
    !cum = nonexceedance_prob 
    !PRINT *,'fcst, cum, fz_fcst = ',fcst, cum, fz_fcst 
    !IF (cum .gt. 0.95) THEN 
    IF (nonexceedance_prob .gt. 0.90) THEN     

        !PRINT *,'cum = ',cum,' so we are in the top 10 percentile'
        ! ---- with minimal training data, estimation of the correction
        !      at the tails of the distribution may be error prone.
        !      Accordingly, if we are above the 90th percentile of the 
        !      distribution, we will apply an alternative procedure to
        !      direct quantile mapping.  We will estimate forecast and 
        !      analyzed values associated with the 91st to 
        !      99th percentiles of the distribution, and use this
        !      data to form a regression relationship.  The actual
        !      correction applied will be a regression correction based
        !      on this data.   First step here is to get the 90-99th
        !      percentiles, analyzed and forecast        

        CALL gamma_get_90_to_99(iadim, ipdim, alpha_fcst, beta_fcst, fz_fcst, &
            alpha_anal, beta_anal, fz_anal, quantile_table, &
            alpha_values, cumprob_values, a90_to_a99, f90_to_f99)

        !      Determine the regression slope associated with the correction
        !      following Schuerer's method discussed in the article
        !      http://journals.ametsoc.org/doi/pdf/10.1175/MWR-D-15-0061.1
        !      Apply that assuming regression intercept
        !      is set by the analyzed value at 90th percentile. 

        slope = SUM((a90_to_a99(2:9) - a90_to_a99(1)) * &
        	(f90_to_f99(2:9) - f90_to_f99(1))) / &
        	SUM((f90_to_f99(2:9) - f90_to_f99(1))**2)
        qmapped_fcst = a90_to_a99(1) + (fcst-f90_to_f99(1))*slope
        
        ! ---- If above the 99th percentile, the amount that the forecast
        !      exceeds the 99th forecast percentile is added to the regressed
        !      amount. 
        
        IF (fcst .gt. f90_to_f99(9)) qmapped_fcst = a90_to_a99(1) + &
            slope*(f90_to_f99(9) - f90_to_f99(1)) + (fcst - f90_to_f99(9))
            
            !PRINT *,'a90_to_a99(1), slope,f90_to_f99(9-1), fcst =  ',&
            !    a90_to_a99(1), slope, f90_to_f99(9), f90_to_f99(1), fcst 
    ELSE
 
        CALL get_alpha_index(iadim, alpha_values, alpha_anal, ia_anal)
        CALL get_cum_index(ipdim, cumprob_values, cum, ic_anal)
        
        !PRINT *,'alpha_anal = ', alpha_anal,' cum = ',cum,&
        !    'ia_anal, ic_anal = ', ia_anal, ic_anal
        IF (ia_anal .gt. 1 .and. ic_anal .gt. 1) THEN
            qmapped_fcst = beta_anal * quantile_table(ia_anal, ic_anal)
        ELSE
            IF (cum .lt. 0) THEN
                PRINT *,' cum < 0 ',cum
                STOP
            ENDIF
            qmapped_fcst = REAL(qgamma(cum,alpha_anal,beta_anal,tootrue,toofalse))  
        ENDIF
        !PRINT *,'fcst, qmapped_fcst, beta_anal, quantile_table(ia_anal, ic_anal) = ',&
        !    fcst, qmapped_fcst, beta_anal, quantile_table(ia_anal, ic_anal)
        
    ENDIF
    IF (qmapped_fcst .lt. 0.0) qmapped_fcst = 0.0
    IF (qmapped_fcst .ne. qmapped_fcst) qmapped_fcst = fcst

ELSE ! map to zero
    !PRINT *, 'quantile mapping to zero'
    qmapped_fcst = 0.0

ENDIF
!print *,'qmapped_fcst = ', qmapped_fcst


RETURN
END SUBROUTINE gamma_quantile_map

! ==================================================================================

SUBROUTINE gamma_get_90_to_99(iadim, ipdim, alpha_fcst, beta_fcst, fz_fcst, &
    alpha_anal, beta_anal, fz_anal, quantile_table, alpha_values, cumprob_values, &    
    a90_to_a99, f90_to_f99)
    
INTEGER, INTENT(IN) :: iadim, ipdim 
DOUBLE PRECISION, INTENT(IN) :: alpha_fcst, beta_fcst, fz_fcst
DOUBLE PRECISION, INTENT(IN) :: alpha_anal, beta_anal, fz_anal 
REAL*8, INTENT(IN), DIMENSION(iadim, ipdim) :: quantile_table
REAL*8, INTENT(IN), DIMENSION(iadim) :: alpha_values
REAL*8, INTENT(IN), DIMENSION(ipdim) :: cumprob_values

DOUBLE PRECISION ksi, cum, cuma, cumf, ccum, nonexceedance_prob, qgamma, qstart, qstride
LOGICAL tootrue
LOGICAL toofalse

REAL, DIMENSION(9) :: a90_to_a99, f90_to_f99
    
tootrue = .TRUE.
toofalse = .FALSE.
    
!PRINT *,'max, min quantile_table = ', maxval(quantile_table), minval(quantile_table)
DO i = 1,9
    qstart = MAX(0.9, fz_anal)
    qstride = (1.0 - qstart)/10.
    ccum = qstart + DBLE(i)*qstride   
    cuma = 1.0 - (1. - ccum) / (1.-fz_anal)
    
    qstart = MAX(0.9, fz_fcst)
    qstride = (1.0 - qstart)/10.
    ccum = qstart + DBLE(i)*qstride     
    cumf = 1.0 - (1. - ccum) / (1.-fz_fcst)
    IF (cuma .lt. 0) THEN
        PRINT *,'cuma, ccum, fz_anal = ',cuma, ccum, fz_anal
        stop
    ENDIF
    
    IF (cumf .lt. 0) THEN
        PRINT *,'cumf, ccum, fz_anal = ',cumf, ccum, fz_fcst
        stop
    ENDIF
    
    !a90_to_a99(i) = REAL(qgamma(cuma,alpha_anal,beta_anal,tootrue,toofalse))  
    CALL get_alpha_index(iadim, alpha_values, alpha_anal, ia_anal)
    CALL get_cum_index(ipdim, cumprob_values, cuma, ic_anal)
    !!PRINT *,'i = ',i,' beta_anal, cuma = ',beta_anal,cuma,'ia_anal, ic_anal = ', ia_anal, ic_anal
    IF (ia_anal .gt. 1 .and. ic_anal .gt. 1) THEN
        a90_to_a99(i) = beta_anal * quantile_table(ia_anal, ic_anal)
    ELSE
        a90_to_a99(i) = REAL(qgamma(cuma,alpha_anal,beta_anal,tootrue,toofalse))  
    ENDIF
    
    !f90_to_f99(i) = REAL(qgamma(cumf,alpha_fcst,beta_fcst,tootrue,toofalse))
    CALL get_alpha_index(iadim, alpha_values, alpha_fcst, ia_fcst)
    CALL get_cum_index(ipdim, cumprob_values, cumf, ic_fcst)
    !PRINT *,'i = ',i,' beta_fcst, cumf = ',beta_fcst, cumf,'ia_fcst, ic_fcst = ', ia_fcst, ic_fcst
    IF (ia_fcst .gt. 1 .and. ic_fcst .gt. 1) THEN
        f90_to_f99(i) = beta_fcst * quantile_table(ia_fcst, ic_fcst)  
    ELSE
        f90_to_f99(i) = REAL(qgamma(cumf,alpha_fcst,beta_fcst,tootrue,toofalse))
    ENDIF
    !PRINT *,'a90_to_a99(i), f90_to_f99(i) = ',a90_to_a99(i), f90_to_f99(i)
END DO

RETURN
END SUBROUTINE gamma_get_90_to_99

! ============================================================================

SUBROUTINE get_alpha_index(iadim,alpha_values, alpha, ia)

! return the index of alpha values for which the value of alpha is closest
INTEGER, INTENT(IN) :: iadim
REAL*8, INTENT(IN), DIMENSION(iadim) :: alpha_values
DOUBLE PRECISION, INTENT(IN) :: alpha
INTEGER, INTENT(OUT) :: ia

REAL*8 alpha_begin, delta_alpha

alpha_begin = alpha_values(1)
delta_alpha = alpha_values(2) - alpha_values(1)

ia = NINT( (1.0/delta_alpha) * (alpha - alpha_begin))
IF (ia .lt. 1 .or. ia .gt. iadim) ia = -99

RETURN
END SUBROUTINE get_alpha_index


! ============================================================================

SUBROUTINE get_cum_index(ipdim, cum_values, cum, ic)

! return the index of cum_values for which the value of cum is closest

INTEGER, INTENT(IN) :: ipdim
REAL*8, INTENT(IN), DIMENSION(ipdim) :: cum_values
DOUBLE PRECISION, INTENT(IN) :: cum
INTEGER, INTENT(OUT) :: ic

REAL*8 cum_begin, cum_end, delta_cum

cum_begin = cum_values(1)
delta_cum = cum_values(2) - cum_values(1)

ic = NINT( (1.0/delta_cum) * (cum - cum_begin))
IF (ic .lt. 1 .or. ic .gt. ipdim) ic = -99

RETURN
END SUBROUTINE get_cum_index




    

