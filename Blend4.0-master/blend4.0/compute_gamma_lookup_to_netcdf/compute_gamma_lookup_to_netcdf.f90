PROGRAM compute_gamma_lookup_to_netcdf
    
! intended application is to produce a lookup table to compute the quantile function
! of the gamma distribution for a range of cumulative probabilities and 
! for a range of alpha parameter values.   The resultant values of the lookup 
! table will be written to a netCDF file.   Beta value is fixed at 1.0, but
! resultant quantile mapping code can adjust to different beta values (a simple
! multiple of the value for beta=1)
!
USE netcdf

INTEGER, PARAMETER :: ipdim = 999 ! cumulative probability, 0.001 to 0.999 by 0.001
INTEGER, PARAMETER :: iadim = 50000 ! alpha, 0.0001 to 5.0 by 0.0001
DOUBLE PRECISION, PARAMETER :: beta = 1.0
DOUBLE PRECISION, DIMENSION(iadim, ipdim) :: quantile
DOUBLE PRECISION, DIMENSION(iadim) :: alpha
DOUBLE PRECISION, DIMENSION(ipdim) :: cumprob
DOUBLE PRECISION :: qgamma
INTEGER cumprob_varid, alpha_varid, quantile_varid

CHARACTER*256, PARAMETER :: data_directory = '/Users/thamill/precip/ecmwf_data/'
CHARACTER*256 :: outfile

LOGICAL tootrue
LOGICAL toofalse
    
INTEGER :: dimid_2d(2)

tootrue = .TRUE.
toofalse = .FALSE.

! ---- initialize cumprob vector, 0.001 to 0.999 by 0.001

DO ip = 1, ipdim
    cumprob(ip) = DBLE(ip) / 1000.
END DO

! ---- initialize alpha vector, 0.0001 to 5 by 0.0001

DO ia = 1, iadim
    alpha(ia) = DBLE(ia) / 10000.
END DO

! ---- generate lookup table with calls to qgamma

DO ip = 1, ipdim
    IF (MOD(ip,50) .eq. 0) PRINT *,'processing ip = ',ip,' of ',ipdim
    DO ia = 1, iadim
        quantile(ia,ip) = qgamma(cumprob(ip),alpha(ia),beta,tootrue,toofalse)
    END DO
END DO

! =========================================================================
! write the netCDF file
! =========================================================================

outfile = TRIM(data_directory) // 'gamma_quantile_function_lookuptable.nc'
PRINT *,'writing to ', TRIM(outfile)

! ---- Create the netCDF file.

CALL check( nf90_create(TRIM(outfile), NF90_CLOBBER, ncid) )

! ---- Define the array dimensions. NetCDF will hand back an ID for each.

CALL check( nf90_def_dim(ncid, "ipdim", ipdim, ip_dimid) )
CALL check( nf90_def_dim(ncid, "iadim", iadim, ia_dimid) )
PRINT *,'defined array dimensions'
dimid_2d =  (/ ia_dimid, ip_dimid /)

! ---- Define the variables and associated IDs

CALL check( nf90_def_var(ncid, "cumprob", &
    NF90_DOUBLE, ip_dimid, cumprob_varid) )
    PRINT *,'cumprob'
CALL check( nf90_def_var(ncid, "alpha", &
    NF90_DOUBLE, ia_dimid, alpha_varid) ) 
    PRINT *,'alpha'
CALL check( nf90_def_var(ncid, "quantile", &
    NF90_DOUBLE, dimid_2d, quantile_varid) )    
PRINT *,'end define mode'
 
! --- End define mode. This tells netCDF we are done defining metadata.

CALL check( nf90_enddef(ncid) )

! ---- write the data.

CALL check( nf90_put_var(ncid, cumprob_varid, cumprob))
PRINT *,'wrote cumprob'
CALL check( nf90_put_var(ncid, alpha_varid, alpha))
PRINT *,'wrote alpha'
CALL check( nf90_put_var(ncid, quantile_varid, quantile))
PRINT *,'wrote quantile'
PRINT *,'max, min quantile = ', maxval(quantile), minval(quantile)

! ---- Close the file. This frees up any internal netCDF resources
!      associated with the file, and flushes any buffers.

CALL check( nf90_close(ncid) )

RETURN 
END PROGRAM compute_gamma_lookup_to_netcdf
        