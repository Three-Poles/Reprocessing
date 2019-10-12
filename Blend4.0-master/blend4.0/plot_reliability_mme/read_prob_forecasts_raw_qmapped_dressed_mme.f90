! ============================================================================

SUBROUTINE read_prob_forecasts_raw_qmapped_dressed_mme(nxa, nya, cyyyymmddhh, &
    cftype, clead_use, rthresh, nthreshes, pthreshes, conusmask, rlonsa, rlatsa, &
    prob_forecast_raw, prob_forecast_qmapped, prob_forecast_dressed, climo_prob)

USE netcdf

INTEGER, INTENT(IN) :: nxa, nya
CHARACTER*(*), INTENT(IN) :: cyyyymmddhh
CHARACTER*(*), INTENT(IN) :: clead_use, cftype
REAL, INTENT(IN) :: rthresh
REAL, INTENT(OUT), DIMENSION(nthreshes) :: pthreshes

CHARACTER*120 infile1, infile2, infile3
CHARACTER*25 cfield
INTEGER*2, INTENT(OUT), DIMENSION(nxa,nya) :: conusmask
REAL, INTENT(OUT), DIMENSION(nxa,nya) :: rlonsa, rlatsa, &
    prob_forecast_dressed, climo_prob, prob_forecast_raw, &
    prob_forecast_qmapped
LOGICAL iex1, iex2, iex3

REAL, DIMENSION(nxa,nya) :: prob_forecast1, prob_forecast2, &
    prob_forecast3

infile1 = '/Users/thamill/precip/ecmwf_data/NCEP/'// TRIM(cftype) // &
    '/NCEP_'//TRIM(clead_use)//'h_IC'//cyyyymmddhh//'.nc'
infile2 = '/Users/thamill/precip/ecmwf_data/CMC/'// TRIM(cftype) // &
    '/CMC_'//TRIM(clead_use)//'h_IC'//cyyyymmddhh//'.nc'
infile3 = '/Users/thamill/precip/ecmwf_data/ECMWF/'// TRIM(cftype) // &
    '/ECMWF_'//TRIM(clead_use)//'h_IC'//cyyyymmddhh//'.nc'

    
PRINT *,TRIM(infile1)
PRINT *,TRIM(infile2)
PRINT *,TRIM(infile3)

INQUIRE (file=infile1,exist=iex1)
INQUIRE (file=infile2,exist=iex2)
INQUIRE (file=infile3,exist=iex3)
PRINT *,'iex1, iex2, iex3 = ', iex1, iex2, iex3 

IF (iex1 .and. iex2 .and. iex3) THEN
    
    ! ---- open the file

    netid = 0
    CALL check (nf90_open(infile1,NF90_NOWRITE,netid1))
    CALL check (nf90_open(infile2,NF90_NOWRITE,netid2))
    CALL check (nf90_open(infile3,NF90_NOWRITE,netid3))

    ! ---- read in the list of dates/times in yyyymmddhh format associated with
    !      each time index stored in the netcdf file

    cfield ='pthreshes'
    CALL check(nf90_inq_varid(netid1,trim(adjustl(cfield)),ivar))
    CALL check(nf90_get_var(netid1,ivar,pthreshes,&
    	start=(/1/),count=(/nthreshes/)))
        
    DO ithresh = 1, nthreshes
        IF (ABS(pthreshes(ithresh) - rthresh) .lt. 0.01) GOTO 1000
    END DO
    PRINT *,'Did not find threshold = ', rthresh,' in subroutine verify_relia_bss_anymodel.f90'
    PRINT *,'Stopping.'
    STOP
    1000 CONTINUE    
    PRINT *,'using thresh = ',ithresh,rthresh
        
    cfield ='conusmask'
    CALL check(nf90_inq_varid(netid1,trim(adjustl(cfield)),ivar))
    CALL check(nf90_get_var(netid1,ivar,conusmask,&
        start=(/1,1/),count=(/nxa,nya/)))
    
    cfield ='rlonsa'
    CALL check(nf90_inq_varid(netid1,trim(adjustl(cfield)),ivar))
    CALL check(nf90_get_var(netid1,ivar,rlonsa,&
        start=(/1,1/),count=(/nxa,nya/)))

    cfield ='rlatsa'
    CALL check(nf90_inq_varid(netid1,trim(adjustl(cfield)),ivar))
    CALL check(nf90_get_var(netid1,ivar,rlatsa,&
        start=(/1,1/),count=(/nxa,nya/)))
    	
    cfield ='prob_forecast'
    CALL check(nf90_inq_varid(netid1,trim(adjustl(cfield)),ivar1))
    CALL check(nf90_get_var(netid1,ivar1,prob_forecast1,&
        start=(/1,1,ithresh/),count=(/nxa,nya,1/)))
    CALL check(nf90_inq_varid(netid2,trim(adjustl(cfield)),ivar2))
    CALL check(nf90_get_var(netid2,ivar2,prob_forecast2,&
        start=(/1,1,ithresh/),count=(/nxa,nya,1/)))        
    CALL check(nf90_inq_varid(netid3,trim(adjustl(cfield)),ivar3))
    CALL check(nf90_get_var(netid3,ivar3,prob_forecast3,&
            start=(/1,1,ithresh/),count=(/nxa,nya,1/)))     
    pmax1 = maxval(prob_forecast1)
    pmax2 = maxval(prob_forecast2)
    pmax3 = maxval(prob_forecast3)  
    IF (pmax1 .gt. 0.0 .and. pmax2 .gt. 0.0 .and. pmax3 .gt. 0.0) THEN  
        prob_forecast_dressed = 0.25 * prob_forecast1 + 0.25 * prob_forecast2 + &
            0.5 * prob_forecast3
    ELSE
        prob_forecast_dressed(:,:) = -99.99
    ENDIF
        
    cfield ='prob_forecast_raw'
    CALL check(nf90_inq_varid(netid1,trim(adjustl(cfield)),ivar1))
    CALL check(nf90_get_var(netid1,ivar1,prob_forecast1,&
        start=(/1,1,ithresh/),count=(/nxa,nya,1/)))
    CALL check(nf90_inq_varid(netid2,trim(adjustl(cfield)),ivar2))
    CALL check(nf90_get_var(netid2,ivar2,prob_forecast2,&
        start=(/1,1,ithresh/),count=(/nxa,nya,1/)))        
    CALL check(nf90_inq_varid(netid3,trim(adjustl(cfield)),ivar3))
    CALL check(nf90_get_var(netid3,ivar3,prob_forecast3,&
            start=(/1,1,ithresh/),count=(/nxa,nya,1/)))        
    prob_forecast_raw = 0.25 * prob_forecast1 + 0.25 * prob_forecast2 + &
        0.5 * prob_forecast3
        
    IF (pmax1 .gt. 0.0 .and. pmax2 .gt. 0.0 .and. pmax3 .gt. 0.0) THEN  
        prob_forecast_raw = 0.25 * prob_forecast1 + 0.25 * prob_forecast2 + &
            0.5 * prob_forecast3
    ELSE
        prob_forecast_raw(:,:) = -99.99
    ENDIF
            
    cfield ='prob_forecast_qmapped'
    CALL check(nf90_inq_varid(netid1,trim(adjustl(cfield)),ivar1))
    CALL check(nf90_get_var(netid1,ivar1,prob_forecast1,&
        start=(/1,1,ithresh/),count=(/nxa,nya,1/)))
    CALL check(nf90_inq_varid(netid2,trim(adjustl(cfield)),ivar2))
    CALL check(nf90_get_var(netid2,ivar2,prob_forecast2,&
        start=(/1,1,ithresh/),count=(/nxa,nya,1/)))        
    CALL check(nf90_inq_varid(netid3,trim(adjustl(cfield)),ivar3))
    CALL check(nf90_get_var(netid3,ivar3,prob_forecast3,&
            start=(/1,1,ithresh/),count=(/nxa,nya,1/)))
            
    IF (pmax1 .gt. 0.0 .and. pmax2 .gt. 0.0 .and. pmax3 .gt. 0.0) THEN  
        prob_forecast_qmapped = 0.25 * prob_forecast1 + 0.25 * prob_forecast2 + &
            0.5 * prob_forecast3
    ELSE
        prob_forecast_qmapped(:,:) = -99.99
    ENDIF          

    !print *,'prob_forecast_qmapped(1:nxa:5,nya/2) = ',prob_forecast_qmapped(1:nxa:5,nya/2)
    	
    cfield ='climo_prob'
    CALL check(nf90_inq_varid(netid1,trim(adjustl(cfield)),ivar))
    CALL check(nf90_get_var(netid1,ivar,climo_prob,&
        start=(/1,1,ithresh/),count=(/nxa,nya,1/)))

    CALL check(nf90_close(netid1))
    CALL check(nf90_close(netid2))
    CALL check(nf90_close(netid3))  
ELSE
    prob_forecast_raw(:,:) = -99.99
    prob_forecast_qmapped(:,:) = -99.99
    prob_forecast_dressed(:,:) = -99.99
    climo_prob(:,:) = -99.99
ENDIF

RETURN
END SUBROUTINE read_prob_forecasts_raw_qmapped_dressed_mme
