SUBROUTINE raw_ensemble_probs_singlemodel(nxa, nya, nens, nthreshes, pthreshes, &
    ensemble_ccpa, conusmask, prob_forecast, ensemble_mean)

INTEGER, INTENT(IN) :: nxa, nya, nens, nthreshes

REAL, INTENT(IN), DIMENSION(nthreshes) :: pthreshes  ! event threshold amount
REAL, INTENT(IN), DIMENSION(nxa,nya,nens) :: ensemble_ccpa
INTEGER*2, INTENT(IN), DIMENSION(nxa,nya) :: conusmask
REAL, INTENT(OUT), DIMENSION(nxa,nya,nthreshes) :: prob_forecast
REAL, INTENT(OUT), DIMENSION(nxa,nya) :: ensemble_mean

REAL, DIMENSION(nens) :: ensemble

emax = MAXVAL(ensemble_ccpa)

IF (emax .le. 0.0) THEN
    PRINT *,'bad data for this day.  Setting all outpust to missing.'
    prob_forecast(:,:,:) = -99.99
ELSE
    DO ithresh = 1, nthreshes 
        rthresh = pthreshes(ithresh)
        DO jya = 1, nya
            DO ixa = 1, nxa
                ensemble(:) = ensemble_ccpa(ixa,jya,:)
                sume = 0.0
                ncount = 0
                ndenom = 0
                DO imem = 1, nens
                    IF (ensemble_ccpa(ixa,jya,imem) .ge. rthresh) ncount = ncount + 1
                    IF (ensemble_ccpa(ixa,jya,imem) .ge. -98.) ndenom = ndenom + 1
                    IF (ensemble_ccpa(ixa,jya,imem) .ge. -98. .and. ithresh .eq. 1) &
                        sume = sume + ensemble_ccpa(ixa,jya,imem) 
                END DO
                IF (ndenom .gt. 0.0) THEN
                    prob_forecast(ixa,jya, ithresh) = REAL(ncount) / REAL(ndenom)
                    IF (ithresh .eq. 1) ensemble_mean(ixa,jya) = sume / REAL(ndenom)
                ELSE
                    prob_forecast(ixa,jya, ithresh) = -99.99
                    IF (ithresh .eq. 1) ensemble_mean(ixa,jya) = -99.99 
                ENDIF
                    
            END DO  ! ixa
        END DO  ! jya
        write (6,fmt='(A40,f8.4)') 'Raw probabilities for thresh = ',rthresh
        write (6,fmt='(3(A10,1X))')'MIN','MAX','MEAN'
        write (6,fmt='(1X,3(F10.5,1X))') &
            MINVAL(conusmask(:,:)*prob_forecast(:,:,ithresh)), &
            MAXVAL(conusmask(:,:)*prob_forecast(:,:,ithresh)), &
            SUM(conusmask(:,:)*prob_forecast(:,:,ithresh))/&
            REAL(SUM(INT(conusmask)))   
    END DO   ! ithresh
    
ENDIF

RETURN
END SUBROUTINE raw_ensemble_probs_singlemodel
