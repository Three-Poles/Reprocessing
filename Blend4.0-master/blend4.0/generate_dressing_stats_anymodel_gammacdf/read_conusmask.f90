
SUBROUTINE read_conusmask(nxa, nya, pclimo_infile, &
    rlonsa, rlatsa, conusmask)

! ---- read in the CONUS mask and the associated lat/lons

USE netcdf

INTEGER, INTENT(IN) :: nxa, nya
CHARACTER*(*), INTENT(IN) :: pclimo_infile
REAL, INTENT(OUT), DIMENSION(nxa,nya) :: rlonsa, rlatsa
INTEGER*2, INTENT(OUT), DIMENSION(nxa,nya) :: conusmask

CHARACTER*20 cfield

! ---- Initialize

conusmask(:,:) = 0
rlonsa(:,:) = 0.
rlatsa(:,:) = 0.

! ---- Open the file, use the number of times in the file later to allocate a date array

netid=0
PRINT *,'netid, reading from ',netid, TRIM(pclimo_infile)
CALL check (nf90_open(pclimo_infile, NF90_NOWRITE, netid))

cfield='conusmask'
CALL check(nf90_inq_varid(netid, trim(adjustl(cfield)), ivar))
CALL check(nf90_get_var(netid, ivar, conusmask, &
    start=(/1,1/), count=(/nxa,nya/)))

cfield='lonsa'
CALL check(nf90_inq_varid(netid, trim(adjustl(cfield)), ivar))
CALL check(nf90_get_var(netid, ivar, rlonsa,&
    start=(/1,1/), count=(/nxa,nya/)))

cfield='latsa'
CALL check(nf90_inq_varid(netid,trim(adjustl(cfield)),ivar))
CALL check(nf90_get_var(netid, ivar, rlatsa,&
    start=(/1,1/), count=(/nxa,nya/)))

! ---- Close netcdf file.

CALL check(nf90_close(netid))
PRINT *,'done reading conus mask'

! ---- check data values

rminlon = MINVAL(rlonsa)
IF (rminlon .gt. 0.0) THEN
    PRINT *,'minimum longitude input does not have longitudes below zero.   Fixing by subtracting 360.'
    rlonsa = rlonsa-360.
ENDIF

rminlat = MINVAL(rlatsa)
rmaxlat = MAXVAL(rlatsa)
IF (rminlat .lt. -90.0 .or. rmaxlat .gt. 90.0) THEN
    PRINT *, 'Error in read_precip_climatology_local.  Stopping. '
    PRINT *, 'Latitudes out of bounds; min, max lat = ', rminlat, rmaxlat
    STOP
ENDIF

RETURN
END subroutine read_conusmask
