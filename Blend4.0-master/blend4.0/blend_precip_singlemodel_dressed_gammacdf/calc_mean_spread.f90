SUBROUTINE calc_mean_spread(n25, nxa, nya, nmembers, &
    ensemble_ccpa_x25, conusmask, ensmean, stddev, POP)
    
USE netcdf

INTEGER, INTENT(IN) :: n25, nxa, nya, nmembers
REAL, INTENT(IN), DIMENSION(n25, nxa, nya, nmembers) :: ensemble_ccpa_x25
INTEGER*2, INTENT(IN), DIMENSION(nxa,nya) :: conusmask
REAL, INTENT(OUT), DIMENSION(nxa, nya) :: ensmean, stddev
REAL, INTENT(OUT), DIMENSION(nxa, nya) :: POP

REAL*8 :: sumxi, sumxi2, sn

ensmean(:,:) = 0.0
stddev(:,:) = 0.0
POP(:,:) = 0.0

DO jya = 1, nya
    DO ixa = 1, nxa
        ntot = 0
        npos = 0
        sumxi = 0.0
        sumxi2 = 0.0
        sn = 0.0
        DO i25 = 1, n25
            DO imem = 1, nmembers
                e = ensemble_ccpa_x25(i25, ixa, jya, imem)
                IF (e .ge. 0.0) THEN
                    sumxi = sumxi + e
                    sumxi2 = sumxi2 + e**2
                    IF (e .gt. 0.254) npos = npos+1
                    ntot = ntot + 1
                    sn = sn + 1.0
                ENDIF
            END DO
        END DO
        
        IF (sn .gt. 0.0) THEN
            ensmean(ixa,jya) = sumxi / sn
            IF (sn .gt. 1.0) THEN
                stddev(ixa,jya) = (sumxi2 - sn*ensmean(ixa,jya)**2)/(sn-1.0)
                stddev(ixa,jya) = SQRT(stddev(ixa,jya))
            ELSE
                stddev(ixa,jya) = 0.0
            ENDIF
            POP(ixa,jya) = REAL(npos) / REAL(ntot)
        ELSE
            ensmean(ixa,jya) = -99.99
            stddev(ixa,jya) = -99.99
            POP(ixa,jya) = -99.99
        ENDIF
            
    END DO
END DO

RETURN 
END SUBROUTINE calc_mean_spread
        