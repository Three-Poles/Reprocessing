SUBROUTINE ensemble_probs_x25_singlemodel(nxa, nya, nens, nthreshes, &
    n25, pthreshes, ensemble_x25, conusmask, prob_forecast)

! --- purpose:  generate probabilities from the quantile-mapped ensemble, here using 
!     also data from not only the current grid point but also 25 surrounding grid
!     points  (Hamill July 2017 modification)

INTEGER, INTENT(IN) :: nxa, nya, nens, nthreshes, n25

REAL, INTENT(IN), DIMENSION(nthreshes) :: pthreshes  ! event threshold amounts
REAL, INTENT(IN), DIMENSION(n25,nxa,nya,nens) :: ensemble_x25
INTEGER*2, INTENT(IN), DIMENSION(nxa,nya) :: conusmask
REAL, INTENT(OUT), DIMENSION(nxa,nya, nthreshes) :: prob_forecast
REAL, DIMENSION(n25,nens) :: ensemble

REAL rthresh
INTEGER, DIMENSION(nxa,nya) :: iktr

DO ithresh = 1, nthreshes
    rthresh = pthreshes(ithresh)
    DO jya = 1, nya
        DO ixa = 1, nxa
            ensemble(:,:) = ensemble_x25(:,ixa,jya,:)
            ncount = 0
            ndenom = 0
            DO i25 = 1, n25
                DO imem = 1, nens
                    IF (ensemble(i25,imem) .ge. rthresh) ncount = ncount + 1
                    IF (ensemble(i25,imem) .ge. 0.0) ndenom = ndenom + 1
                END DO
                IF (ndenom .gt. 0.0) THEN
                    prob_forecast(ixa,jya,ithresh) = REAL(ncount) / REAL(ndenom)
                ELSE
                    prob_forecast(ixa,jya,ithresh) = -99.99
                ENDIF
            END DO  ! i9
        END DO  ! ixa
    END DO  ! jya
    
    write (6,fmt='(A40,f8.4)') 'Quantile-mapped probabilities for thresh = ',rthresh
    write (6,fmt='(3(A10,1X))')'MIN','MAX','MEAN'
    write (6,fmt='(1X,3(F10.5,1X))')&
        MINVAL(conusmask(:,:)*prob_forecast(:,:,ithresh)), &
        MAXVAL(conusmask(:,:)*prob_forecast(:,:,ithresh)), &
        SUM(conusmask(:,:)*prob_forecast(:,:,ithresh))/&
        REAL(SUM(INT(conusmask)))      
    
END DO  ! ithresh

RETURN
END SUBROUTINE ensemble_probs_x25_singlemodel
