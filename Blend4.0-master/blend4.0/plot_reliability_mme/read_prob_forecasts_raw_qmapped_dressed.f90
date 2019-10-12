! ============================================================================

SUBROUTINE read_prob_forecasts_raw_qmapped_dressed(nxa, nya, cyyyymmddhh, &
    cmodel, cftype, clead_use, rthresh, nthreshes, pthreshes, conusmask, &
    rlonsa, rlatsa, prob_forecast_raw, prob_forecast_qmapped, &
    prob_forecast_dressed, climo_prob)

USE netcdf

INTEGER, INTENT(IN) :: nxa, nya
CHARACTER*(*), INTENT(IN) :: cyyyymmddhh
CHARACTER*(*), INTENT(IN) :: cmodel 
CHARACTER*(*), INTENT(IN) :: cftype
CHARACTER*(*), INTENT(IN) :: clead_use
REAL, INTENT(IN) :: rthresh
REAL, INTENT(OUT), DIMENSION(nthreshes) :: pthreshes

CHARACTER*120 infile
CHARACTER*25 cfield
INTEGER*2, INTENT(OUT), DIMENSION(nxa,nya) :: conusmask
REAL, INTENT(OUT), DIMENSION(nxa,nya) :: rlonsa, rlatsa, &
    prob_forecast_dressed, climo_prob, prob_forecast_raw, prob_forecast_qmapped
LOGICAL iex


infile = '/Users/thamill/precip/ecmwf_data/' // TRIM(cmodel) // '/'// &
    TRIM (cftype) // '/' // TRIM(cmodel)//'_'//TRIM(clead_use) // &
    'h_IC'//cyyyymmddhh//'.nc'
    

PRINT *,TRIM(infile)
INQUIRE (file=infile,exist=iex)
IF (iex) THEN
    
    ! ---- open the file

    netid = 0
    CALL check (nf90_open(infile,NF90_NOWRITE,netid))

    ! ---- read in the list of dates/times in yyyymmddhh format associated with
    !      each time index stored in the netcdf file

    cfield ='pthreshes'
    CALL check(nf90_inq_varid(netid,trim(adjustl(cfield)),ivar))
    CALL check(nf90_get_var(netid,ivar,pthreshes,&
    	start=(/1/),count=(/nthreshes/)))
        
    DO ithresh = 1, nthreshes
        IF (ABS(pthreshes(ithresh) - rthresh) .lt. 0.01) GOTO 1000
    END DO
    PRINT *,'Did not find threshold = ', rthresh,' in subroutine verify_relia_bss_anymodel.f90'
    PRINT *,'Stopping.'
    STOP
    1000 CONTINUE    
    !PRINT *,'using thresh = ',ithresh,rthresh
        
    cfield ='conusmask'
    CALL check(nf90_inq_varid(netid,trim(adjustl(cfield)),ivar))
    CALL check(nf90_get_var(netid,ivar,conusmask,&
        start=(/1,1/),count=(/nxa,nya/)))
    
    cfield ='rlonsa'
    CALL check(nf90_inq_varid(netid,trim(adjustl(cfield)),ivar))
    CALL check(nf90_get_var(netid,ivar,rlonsa,&
        start=(/1,1/),count=(/nxa,nya/)))

    cfield ='rlatsa'
    CALL check(nf90_inq_varid(netid,trim(adjustl(cfield)),ivar))
    CALL check(nf90_get_var(netid,ivar,rlatsa,&
        start=(/1,1/),count=(/nxa,nya/)))
    	
    cfield ='prob_forecast'
    CALL check(nf90_inq_varid(netid,trim(adjustl(cfield)),ivar))
    CALL check(nf90_get_var(netid,ivar,prob_forecast_dressed,&
        start=(/1,1,ithresh/),count=(/nxa,nya,1/)))
        
    cfield ='prob_forecast_raw'
    CALL check(nf90_inq_varid(netid,trim(adjustl(cfield)),ivar))
    CALL check(nf90_get_var(netid,ivar,prob_forecast_raw,&
        start=(/1,1,ithresh/),count=(/nxa,nya,1/)))
            
    cfield ='prob_forecast_qmapped'
    CALL check(nf90_inq_varid(netid,trim(adjustl(cfield)),ivar))
    CALL check(nf90_get_var(netid,ivar,prob_forecast_qmapped,&
        start=(/1,1,ithresh/),count=(/nxa,nya,1/)))
        
    pmax1 = maxval(prob_forecast_dressed)
    pmax2 = maxval(prob_forecast_raw)
    pmax3 = maxval(prob_forecast_qmapped)  
    IF (pmax1 .lt. 0.0 .or. pmax2 .lt. 0.0 .and. pmax3 .lt. 0.0) THEN  
        prob_forecast_dressed(:,:) = -99.99
        prob_forecast_raw(:,:)= -99.99
        prob_forecast_qmapped(:,:) = -99.99   
    ENDIF        
    	
    cfield ='climo_prob'
    CALL check(nf90_inq_varid(netid,trim(adjustl(cfield)),ivar))
    CALL check(nf90_get_var(netid,ivar,climo_prob,&
        start=(/1,1,ithresh/),count=(/nxa,nya,1/)))

    CALL check(nf90_close(netid))
    
ELSE
    PRINT *,'******* Missing data for ',TRIM(infile)
    prob_forecast_raw(:,:) = -99.99
    prob_forecast_qmapped(:,:) = -99.99
    prob_forecast_dressed(:,:) = -99.99
    climo_prob(:,:) = -99.99
ENDIF
RETURN
END SUBROUTINE read_prob_forecasts_raw_qmapped_dressed
