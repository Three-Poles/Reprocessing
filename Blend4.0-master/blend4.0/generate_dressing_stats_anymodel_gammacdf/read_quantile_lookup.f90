SUBROUTINE read_quantile_lookup(iadim, ipdim, data_directory, &
    quantile_table, alpha_values, cumprob_values)
    
    ! read in a lookup table for defining the precipitation 
    ! amount (the quantile) associated with a particular alpha
    ! value of the Gamma distribution and a particular
    ! cumulative probability.  The lookup table was 
    ! previously generated (see directory 
    ! compute_gamma_lookup_to_netcdf)

USE netcdf
INTEGER, INTENT(IN) :: iadim, ipdim
CHARACTER*(*), INTENT(IN) :: data_directory

REAL*8, INTENT(OUT), DIMENSION(iadim, ipdim) :: quantile_table
REAL*8, INTENT(OUT), DIMENSION(iadim) :: alpha_values
REAL*8, INTENT(OUT), DIMENSION(ipdim) :: cumprob_values
      
CHARACTER*256 infilename
CHARACTER*20 cfield


infilename = TRIM(data_directory) // &
    'gamma_quantile_function_lookuptable.nc'
PRINT *, 'quantile lookup table read from ', TRIM(infilename)

! --- read cdf and precip analysis cdf

CALL check (nf90_open(infilename,NF90_NOWRITE,netid))

cfield='cumprob'
CALL check(nf90_inq_varid(netid,trim(adjustl(cfield)),ivar))
CALL check(nf90_get_var(netid,ivar,cumprob_values,&
    start=(/1/),count=(/ipdim/)))
    
cfield='alpha'
CALL check(nf90_inq_varid(netid,trim(adjustl(cfield)),ivar))
CALL check(nf90_get_var(netid,ivar,alpha_values,&
    start=(/1/),count=(/iadim/)))
    
cfield='quantile'
CALL check(nf90_inq_varid(netid,trim(adjustl(cfield)),ivar))
CALL check(nf90_get_var(netid,ivar,quantile_table,&
    start=(/1,1/),count=(/iadim,ipdim/)))
    
CALL check(nf90_close(netid))
PRINT *, 'done reading lookup table'

RETURN
END SUBROUTINE read_quantile_lookup
