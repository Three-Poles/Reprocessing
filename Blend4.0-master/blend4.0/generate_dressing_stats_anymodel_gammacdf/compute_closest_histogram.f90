SUBROUTINE compute_closest_histogram (n25, nxa, nya, nmembers, &
    thresh_light, thresh_mod, thresh_high, &
    ensemble_x25, analysis, conusmask, outfilename_nc, &
    istat)

    ! Here is where we actually compare ensemble forecasts for a given prediction 
    ! system to the verifying analyses, tallying the information we need
    ! to generate (later, over many case days) the closest-histogram statistics
    ! There is also a dependence for both
    ! the closest histogram and the gamma-distribution statistics on the 
    ! ensemble-mean precipitation amount.   Data is saved to a netcdf file
    ! for this particular model and forecast initial time and forecast lead.
    
    ! coded by: Tom Hamill, ESRL/PSD, tom.hamill@noaa.gov (303) 497-3060

USE netcdf
    
INTEGER, INTENT(IN) :: n25, nxa, nya, nmembers
REAL, INTENT(IN) :: thresh_light, thresh_mod, thresh_high 
    ! for breaking up closest histogram etc
REAL, INTENT(IN), DIMENSION (n25, nxa, nya, nmembers) :: ensemble_x25 ! qmapped ensemble
REAL, INTENT(IN), DIMENSION (nxa, nya) :: analysis ! analyzed precip amt.
INTEGER*2, INTENT(IN), DIMENSION(nxa,nya) :: conusmask 
CHARACTER*(*), INTENT(IN) :: outfilename_nc

REAL, INTENT(OUT) :: istat

! ---- Here are the variables we intend to populate with the forecast, analysis information.
!      dimension 2 below is for light, mod, heavy precip, with thresholds between
!      them set by inputs thresh_light, thresh_mod, thresh_high.
!      dimension 3 below is for lowest sorted member(1), intermediate(2), highest(3) 

INTEGER, DIMENSION(nmembers*n25, 3) :: closest_histogram


REAL rdum 
INTEGER idum

integer :: dimid_2d(2)

print *, 'subroutine compute_closest_histogram'

istat = 1

! ----- initialize to zero

closest_histogram(:,:) = 0.0

! ---- set a random number seed for the closest member histogram
!      to use to allocate closest member in the case of ties.

rdum = 0.0
DO ixa = 1, nxa
	rdum = rdum + ensemble_x25(1,ixa,nya/2,1)
END DO
idum = NINT(rdum)

! ---- only process if input analysis, ensemble have valid data ....

rm = MINVAL(ensemble_x25)
rma = MAXVAL(analysis)
IF (rm .lt. -98. .and. rma .lt. -98.) THEN
	istat = -1
	PRINT *,'identified bad ensemble or analysis data.  rm, rma = ',rm, rma
	PRINT *,'setting istat to -1 in tally_gamma_stats_full_x25, exiting subroutine.'
ELSE
    
    ! --- loop thru and process all points in the CONUS
    
	DO jya = 1, nya
		DO ixa = 1, nxa
			
			IF (conusmask(ixa, jya) .eq. 1 .and. analysis(ixa,jya) .ge. 0.0) THEN
            
				! --- determine which member is the closest to the
				!     analyzed and how many members have values lower 
				!     than or equal to the analyzed value
			
				rsum = 0.0
				iktr = 0
				DO i25 = 1, n25
					DO imem = 1, nmembers
						IF (ensemble_x25(i25,ixa,jya,imem) .ge. 0.0) THEN
							rsum = rsum + ensemble_x25(i25,ixa,jya,imem)
							iktr = iktr+1
						ENDIF
					END DO
				END DO
				emean = rsum / REAL(iktr)
                emin = minval(ensemble_x25(:,ixa,jya,:))
                emax = maxval(ensemble_x25(:,ixa,jya,:))
                a = analysis(ixa, jya)

                ! ---- we will tally stats one of two ways below.  For situations
                !      where the precipitation is heavier than a threshold that
                !      is slightly greater than zero, we keep track of closest
                !      histogram, fraction zero, and Gamma alpha and beta parameters.
                !      In situations where the mean precip is zero or really 
                !      close to zero, we will tally the fraction zero, and the 
                !      Gamma alpha and beta parameters as a function of the 
                !      climatological probability of precipitation.
    
                IF (emean .ge. thresh_light) THEN ! nonzero precip
                    
                    ! ---- determine what the category is associated with this 
                    !      (quantile-mapped) ensemble-mean precipitation amount
                    
                    IF (emean .ge. thresh_light .and. emean .lt. thresh_mod) THEN
                        ipcat = 1
                    ELSE IF (emean .ge. thresh_mod .and. emean .lt. thresh_high) THEN
                        ipcat = 2
                    ELSE IF (emean .ge. thresh_high) THEN
                        ipcat = 3
                    ELSE
                        ipcat = -99
                        PRINT *,'invalid value for ipcat'
                        PRINT *,'tally_gamma_stats_full_n25.'
                        PRINT *,'Stopping.  ipcat, emean = ',ipcat, emean
                        STOP
                    ENDIF  
                    
                    ! ----- Increment closest histogram ----------------
                    
    				! ---- find the ensemble forecast member that is the closest
                    !      to the analyzed

    				i25closest = 1
    				imemclosest = 1
    				rclosest = 9999.
                    eclosest = 0.
    				DO i25 = 1, n25
    					DO imem = 1, nmembers
    						e = ensemble_x25(i25,ixa,jya,imem)
    						diff = ABS(a - e)
    						IF (diff .lt. rclosest .and. e .gt. -99) THEN
    							rclosest = diff
    							eclosest = e
    							i25closest = i25
    							imemclosest = imem
    						ENDIF
    					END DO
    				END DO                    

    				! ---- determine how many other members are lower than the
    				!      closest member, and how many are equal,  in order
                    !      to determine the rank of the closest in the sorted
                    !      ensemble.
				
    				ibelow = 0
    				iequal = 0
    				DO i25 = 1, n25
    					DO imem = 1, nmembers
    						e = ensemble_x25(i25,ixa,jya,imem)
    						IF (i25 .eq. i25closest .and. imem .eq. imemclosest) THEN
                                CONTINUE
                            ELSE
    							IF (e .lt. eclosest) ibelow = ibelow + 1
    							IF (e .eq. eclosest) iequal = iequal + 1
    						ENDIF
    					END DO
    				END DO

    				! --- determine the closest_histogram rank, + a randomization procedure
                    !     in the case of ties of ensemble and analysis (common when both = 0)
				
    				IF (iequal .eq. 0) THEN
    					iclosest = ibelow + 1			
    				ELSE
    					r = ran3(idum) * REAL(iequal)
    					ir = INT(r)
    					IF (ir .gt. iequal) ir = iequal
    					iclosest = ibelow + ir + 1
    				ENDIF                
                    
                    ! ---- increment closest histogram for nonzero ens mean precip 
                    
                    closest_histogram(iclosest,ipcat) = &
                        closest_histogram(iclosest,ipcat) + 1
                                    
                ENDIF  ! emean .ge. thresh_light	
            ENDIF ! conusmask
		END DO ! ixa
	END DO ! jya

    ! ---- print out the closest histogram sample

    PRINT *,'closest_histogram near 0.1-2.0 mm precip'
    DO i25 = 1,n25*nmembers, 20
        istart = i25
        iend = istart + 19
        IF (iend .le. n25*nmembers) PRINT 200, closest_histogram(istart:iend,1)
        200 FORMAT(20(i4,1x))
    END DO
    
    PRINT *,'closest_histogram 2-6 mm precip'
    DO i25 = 1,n25*nmembers, 20
        istart = i25
        iend = istart + 19
        IF (iend .le. n25*nmembers) PRINT 200, closest_histogram(istart:iend,2)
    END DO
    
    PRINT *,'closest_histogram> 6 mm precip'
    DO i25 = 1,n25*nmembers, 20
        istart = i25
        iend = istart + 19
        IF (iend .le. n25*nmembers) PRINT 200, closest_histogram(istart:iend,3)
    END DO    
    
    ! ====================== WRITE DAILY OUTPUT TO NETCDF FILE ==========================
    
    ! ---- Create the netCDF file.
    
    print *,'writing to ',TRIM(outfilename_nc)
    CALL check( nf90_create(TRIM(outfilename_nc), NF90_CLOBBER, ncid) )
    
    ! ---- Define the array dimensions. NetCDF will hand back an ID for each.
  
    PRINT *,'array dimensions'
    CALL check( nf90_def_dim(ncid, "nhist", nmembers*n25, nhist_dimid) )
    CALL check( nf90_def_dim(ncid, "nscalar", 1, nscalar_dimid) )
    CALL check( nf90_def_dim(ncid, "nprecipcats", 3, nprecipcats_dimid) )
    dimid_2d =  (/ nhist_dimid, nprecipcats_dimid /)

    ! ---- Define the variables and associated IDs

    PRINT *,'defining variables'
    CALL check( nf90_def_var(ncid, "closest_histogram", &
        NF90_INT, dimid_2d, nhist_varid) ) 

    CALL check( nf90_def_var(ncid, "thresh_light", &
        NF90_FLOAT, nscalar_dimid, nthresh_light_varid) ) 

    CALL check( nf90_def_var(ncid, "thresh_mod", &
        NF90_FLOAT, nscalar_dimid, nthresh_mod_varid) ) 

    CALL check( nf90_def_var(ncid, "thresh_high", &
        NF90_FLOAT, nscalar_dimid, nthresh_high_varid) ) 

      
    ! --- End define mode. This tells netCDF we are done defining metadata.

    CALL check( nf90_enddef(ncid) )

    ! ---- write the data.  

    PRINT *,'writing data'
    CALL check( nf90_put_var(ncid, nhist_varid, closest_histogram))
    CALL check( nf90_put_var(ncid, nthresh_light_varid, thresh_light))
    CALL check( nf90_put_var(ncid, nthresh_mod_varid, thresh_mod))    
    CALL check( nf90_put_var(ncid, nthresh_high_varid, thresh_high))
    
    ! ---- Close the file. This frees up any internal netCDF resources
    !      associated with the file, and flushes any buffers.
    
    CALL check( nf90_close(ncid) )
    PRINT *, "*** SUCCESS writing netcdf file "
    
ENDIF
				
RETURN
END SUBROUTINE compute_closest_histogram

