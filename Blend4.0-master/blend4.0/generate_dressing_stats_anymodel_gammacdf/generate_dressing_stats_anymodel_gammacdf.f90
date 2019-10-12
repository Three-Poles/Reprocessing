PROGRAM generate_dressing_stats_anymodel_gammacdf
 
! This fortran90 program generates closest-member 
! histogram dressing statistics for precipitation
! forecasts used in the National Blend.  User inputs the the date/time in 
! yyyymmddhh format (UTC), name of the modeling system (currently CMC for 
! Canadian Center ensemble, NCEP for GEFS ensemble, ECMWF for their system's 
! ensemble), and the ending lead time, e.g., 24 for 12 to 24-h forecast.  

! The output is a netCDF file with the synthesized closest-member histogram 
! information, such that at a later step when combined with information 
! from other training days, it is possible to estimate the closest-member 
! histograms with large-enough samples.
! 
! coded by Tom Hamill, 2017-2018. tom.hamill@noaa.gov, (303) 497-3060
  
USE netcdf

INTEGER, PARAMETER :: nxa = 464  ! number of grid pts in x-dir for 1/8-deg analysis grid
INTEGER, PARAMETER :: nya = 224  ! number of grid pts in y-dir for 1/8-deg analysis grid 
INTEGER, PARAMETER :: n25 = 25
INTEGER, PARAMETER :: nens_ecmwf = 50 ! number of ECMWF perturbed ensemble members
INTEGER, PARAMETER :: nens_cmc   = 20 ! number of CMC perturbed ensemble members
INTEGER, PARAMETER :: nens_ncep  = 20 ! number of NCEP perturbed ensemble members
INTEGER, PARAMETER :: ifcstint = 12 ! number of hours for accumulated precip
INTEGER, PARAMETER :: ipdim = 999 ! cumulative probability, 0.001 to 0.999 by 0.001
INTEGER, PARAMETER :: iadim = 50000 ! alpha, 0.0001 to 5.0 by 0.0001
REAL, PARAMETER :: thresh_light = 0.01 ! divider amount between closest histogram none and low
REAL, PARAMETER :: thresh_mod = 2.0 ! divider amount between closest histogram low and mod
REAL, PARAMETER :: thresh_high = 6.0 ! divider amount between closest histogram mod and high
  
CHARACTER*256, PARAMETER :: data_directory = '/Users/thamill/precip/ecmwf_data/'
    ! data_directory = where I kept all the data 

INTEGER :: nstride ! number of grid pts between samples in 3x3 stencil


CHARACTER*2 chh, cmm  ! character version of UTC hour and month
CHARACTER*3, DIMENSION(12) :: cmonths ! 3-letter version of months of year
CHARACTER*5 :: center ! the prediction center
CHARACTER*256 pclimo_infile ! name of file with climo probs
CHARACTER*256 infile_early ! ecmwf forecast data for today
CHARACTER*256 infile_late ! ecmwf forecast data for today
CHARACTER*256 infile

CHARACTER*256 outfile      ! name of flat fortran file with output prob forecasts
CHARACTER*256 outfile_nc   ! output file name for netcdf data
CHARACTER*10 cyyyymmddhh   ! year,month,day,hour of initial time of forecast
CHARACTER*3 cleade         ! ending hour of precip forecast accumulation, 3 digits, e.g., '024'
CHARACTER*3 cleadb         ! beginning hour of precip forecast accumulation, 3 digits, e.g., '012'
CHARACTER*5 cprefix        ! prefix for data of interest e.g., 'ECMWF'
CHARACTER*5 cmodel
CHARACTER*2 cmembers ! number of members for this system (<100)
INTEGER*2, DIMENSION(nxa,nya) :: conusmask  ! inherited from CCPA data set

! ---- 1/8 deg. Lat/Lon arrays (i.e. CCPA grid)

REAL, ALLOCATABLE, DIMENSION(:,:,:) :: ensemble_ccpa ! input ens precip forecast on 1/8-deg ccpa grid
REAL, ALLOCATABLE, DIMENSION(:,:) :: analysis ! precipitation analysis
REAL, ALLOCATABLE, DIMENSION(:,:) :: rlonsa ! precip analysis grid lat/lons
REAL, ALLOCATABLE, DIMENSION(:,:) :: rlatsa ! precip analysis grid lat/lons
REAL, ALLOCATABLE, DIMENSION(:,:) :: ensemble_mean ! precip analysis grid lat/lons
REAL, ALLOCATABLE, DIMENSION(:,:,:) :: gamma_shape_qmap_forecast ! contains the shape parameter of climatological
!   distribution of Gamma forecasts, used in quantile mapping.
REAL, ALLOCATABLE, DIMENSION(:,:,:) :: gamma_scale_qmap_forecast ! same but Gamma scale parameter
REAL, ALLOCATABLE, DIMENSION(:,:,:) :: fraction_zero_qmap_forecast ! same but fraction of samples with zero precip
REAL, ALLOCATABLE, DIMENSION(:,:) :: gamma_shape_qmap_analysis ! same but for analyzed data
REAL, ALLOCATABLE, DIMENSION(:,:) :: gamma_scale_qmap_analysis
REAL, ALLOCATABLE, DIMENSION(:,:) :: fraction_zero_qmap_analysis
REAL, ALLOCATABLE, DIMENSION(:,:,:) :: precip_anal_cdf
REAL, ALLOCATABLE, DIMENSION(:,:,:,:) :: ensemble_cdf

REAL*8, ALLOCATABLE, DIMENSION(:,:) :: quantile_table
REAL*8, ALLOCATABLE, DIMENSION(:) :: alpha_values
REAL*8, ALLOCATABLE, DIMENSION(:) :: cumprob_values

! ---- x25 array

REAL, ALLOCATABLE, DIMENSION(:,:,:,:) :: ensemble_ccpa_x25 ! ecmwf ens quantile mapped 
    ! precip forecast on 1/8-deg ccpa grid, with 5x5 stencil of surrounding 
    ! grid points in first dimension

INTEGER :: ierr      ! return variable for BAOPEN
INTEGER :: ios       ! return variable for Fortran I/O, Allocation statements

INTEGER :: iyyyymmddhh,jyyyymmddhh
INTEGER :: iyear,imo,iday,ihour,idoy ! Parsed date variables from iyyyymmddhh
INTEGER :: jyear,jmo,jday,jhour,jdoy ! Parsed date variables for valid date

LOGICAL exchangeable

! ---- Initialize

DATA cmonths /'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'/

! =====================================================================================
! START PROGRAM EXECUTION BY READING INFO FROM COMMAND LINE AND SET VARIOUS VARIABLES
! =====================================================================================

! --- Via command line, read in the input year/mo/day/hr and the forecast resolution 
!     we're working with.  Process date to determine the day of the month as an integer

CALL getarg(1, cyyyymmddhh) ! input year month day hour of initial condition, 'yyyymmddhh' format
CALL getarg(2, cmodel)  ! name of the model "ECMWF" "NCEP" "CMC" 
CALL getarg(3, cleade) ! forecast lead time for beginning of precip accum period, hours, e.g.'060'

cmm = cyyyymmddhh(5:6)

! ---- Convert character based variables from command line to integers

READ (cmm,'(i2)') imonth
READ (cyyyymmddhh, '(i10)') iyyyymmddhh
READ (cleade, '(i3)') ileade
PRINT *,'ileade = ',ileade

IF (TRIM(cmodel) .eq. 'ECMWF' .or. TRIM(cmodel) .eq. 'ecmwf') THEN
    cmodel = 'ECMWF' ! make sure not lowercase
    nens = nens_ecmwf
    nens_qmap = 1  ! exchangeable, CDFs same for all members
    exchangeable = .TRUE.
ELSE IF (TRIM(cmodel) .eq. 'NCEP' .or. TRIM(cmodel) .eq. 'ncep') THEN
    cmodel = 'NCEP' 
    nens = nens_ncep
    nens_qmap = 1 
    exchangeable = .TRUE.
ELSE IF (TRIM(cmodel) .eq. 'CMC' .or. TRIM(cmodel) .eq. 'cmc') THEN
    cmodel = 'CMC'
    nens = nens_cmc
    nens_qmap = nens_cmc 
    exchangeable = .FALSE.
ELSE
    PRINT *,'Invalid input model input to generate_dressing_stats_anymodel_gammacdf.x ',TRIM(cmodel)
    PRINT *,'Stopping.'
    STOP
ENDIF    

nmembersx25 = nens*n25
PRINT *,'nens, nmembersx25 = ', nens, nmembersx25

! ---- Output the input information.

PRINT *, '-------------------------------------------------------------'
PRINT *, 'generate_dressing_stats_anymodel_gammacdf.x '
PRINT *, '-------------------------------------------------------------'
PRINT *,' '
write(6,*)' ******* Command line arguments: ******'
write(6,110)  cmodel, cyyyymmddhh, cleade
110 format(1x,'model: ',A/1x,'yyyymmddhh: ',A/1x,'lead in hours: ',A/)

ileadb = ileade - ifcstint
WRITE (cleadb,'(i3)') ileadb

! ---- Set nstride as a function of ileade.  This controls how far away from 
!      each other the grid point are when making a 5x5 stencil of points
!      around the grid point of interest.

nstride = nint(1.+4.*ileade/168.)

! ---- Parse the initializtion date; determine the valid hour.  
!      This is dependent on the precip variable (cpcpvar), 
!      the model initialization (iyyyymmddhh), and forecast 
!      ending hour (ileade).

iendhour=0
CALL doy(iyyyymmddhh,iyear,imo,iday,ihour,idoy)
CALL updat(iyyyymmddhh,ileade,jyyyymmddhh)
CALL doy(jyyyymmddhh,jyear,jmo,jday,jhour,jdoy)
iendhour = jhour

WRITE (6,*)'Model Initialization: ',iyyyymmddhh
WRITE (6,*)'Forecast Projection Ending: ',ileade
WRITE (6,*)'Forecast Valid Date/Hour: ',jyyyymmddhh, iendhour

! ================================================================
! Allocate dynamic arrays
! ================================================================

ALLOCATE(analysis(nxa,nya))
ALLOCATE(ensemble_ccpa_x25(25,nxa,nya,nens))
ALLOCATE(ensemble_ccpa(nxa,nya,nens))
ALLOCATE(ensemble_mean(nxa,nya))
ALLOCATE(fraction_zero_qmap_analysis(nxa,nya))
ALLOCATE(fraction_zero_qmap_forecast(nxa,nya,nens_qmap))
ALLOCATE(gamma_shape_qmap_forecast(nxa,nya,nens_qmap))
ALLOCATE(gamma_scale_qmap_forecast(nxa,nya,nens_qmap))
ALLOCATE(gamma_shape_qmap_analysis(nxa,nya))
ALLOCATE(gamma_scale_qmap_analysis(nxa,nya))
ALLOCATE(rlonsa(nxa,nya))
ALLOCATE(rlatsa(nxa,nya))

! ================================================================
! READ IN DATA
! ================================================================

! ---- use this file to read in the conus mask, lon/lat array. These
!      were previously calculated by the python script 
!      compute_climatology_ppn_multithresh.py in the
!      compute_climatology_ppn_multithresh directory.

write(6,*) 'Calling read_precip_climatology_local'
IF (ihour .eq. 0) THEN 
    pclimo_infile = TRIM(data_directory) // &
        'apcp_climatologies_12_to_00UTC_'//cmonths(imonth)//'_2002_to_2016.nc'
ELSE
    pclimo_infile = TRIM(data_directory) // &
        'apcp_climatologies_00_to_12UTC_'//cmonths(imonth)//'_2002_to_2016.nc'
ENDIF 
   
CALL read_conusmask(nxa, nya, pclimo_infile, rlonsa, rlatsa, conusmask)   

! ---- read the precipitation analysis valid for this lead time.
!      the input file was created by the python script ccpa_to_netcdf.py

infile = TRIM(data_directory)//&    
    'precip_analyses_ccpa_v1_2002010100_to_2016123100.nc'
PRINT *,'reading precipitation analysis from ', TRIM(infile)
CALL read_precipitation_analysis(nxa, nya, jyyyymmddhh,&
    infile, analysis, istat)
PRINT *,'min, max of analyzed precip = ', &
    minval(analysis), maxval(analysis)

! ---- Read precipitation forecasts for all models/ensemble members 

IF (TRIM(cmodel) .eq. 'ECMWF') THEN
    cprefix = 'ECMWF'
ELSE IF (TRIM(cmodel) .eq. 'NCEP') THEN
    cprefix = 'NCEP'
ELSE IF (TRIM(cmodel) .eq. 'CMC') THEN
    cprefix = 'CMC'
ELSE
    PRINT *,'invalid model name ',cprefix
    STOP
ENDIF

infile_late = TRIM(data_directory) // TRIM(cmodel) // &
    '_' // cyyyymmddhh // '_leadtime' // TRIM(ADJUSTL(cleade)) // 'h.nc'
infile_early = TRIM(data_directory) // TRIM(cmodel) // &
    '_' // cyyyymmddhh // '_leadtime' // TRIM(ADJUSTL(cleadb)) // 'h.nc'
PRINT *,'reading forecasts from ', TRIM(infile_early)
PRINT *,'reading forecasts from ', TRIM(infile_late)
CALL read_forecasts_on_CCPA (nxa, nya, nens, infile_early, &
    infile_late, ensemble_ccpa, ensemble_mean)
PRINT *,'min, maxval(ensemble_ccpa) in driver = ', &
    MINVAL(ensemble_ccpa), MAXVAL(ensemble_ccpa)

! ===================================================================
! DATA PROCESSING COMMENCES HERE
! ===================================================================
    
! --- only process if ensemble mean is greater than zero, which will 
!     only happen if valid data has been read in

IF (MAXVAL(ensemble_mean) .gt. 0.0) THEN  ! only if so is data for this date ok 
  
    ALLOCATE(quantile_table(iadim, ipdim))
    ALLOCATE(alpha_values(iadim))
    ALLOCATE(cumprob_values(ipdim))
    PRINT *,'calling read_quantile_lookup'
    
    ! --- read in a lookup table that will be used to efficiently estimate
    !     CDFs for the Gamma distribution
    
    CALL read_quantile_lookup(iadim, ipdim, data_directory, &
        quantile_table, alpha_values, cumprob_values)

    ! --- this subroutine will read in the previous 60 days of synthesized
    !     information on Gamma distributions as well as supplemental location 
    !     information, calculating for each grid point an estimated fraction
    !     zero (fraction of the samples with zero precip amount) and Gamma 
    !     distribution parameters shape (alpha) and scale (beta).  This is done
    !     for both forecast and analyzed data.
  
    PRINT *,'calling determine_gamma_parameters_for_quantile_mapping'
    CALL determine_gamma_parameters_for_quantile_mapping(nxa, nya, &
        iyyyymmddhh, imonth, nens_qmap, data_directory, cmodel, cleade, &
        cmonths, conusmask, gamma_shape_qmap_forecast, gamma_scale_qmap_forecast, &
        fraction_zero_qmap_forecast, gamma_shape_qmap_analysis, &
        gamma_scale_qmap_analysis, fraction_zero_qmap_analysis)
    
    ! ---- perform the quantile mapping using the data at a stencil of 
    !      5 x 5 grid points surrounding the point in question (which is
    !      at the center of the array).  Populate the _x25 arrays with 
    !      the 5 x 5 stencil of quantile-mapped values

    PRINT *,'calling control_quantile_mapping_singlemodel_gamma'
    CALL control_quantile_mapping_singlemodel_gamma(nxa, nya, &
        nstride, nens, nens_qmap, n25, iadim, ipdim, &
        exchangeable, conusmask, ensemble_ccpa, &
        gamma_shape_qmap_forecast, gamma_scale_qmap_forecast, &
        fraction_zero_qmap_forecast, gamma_shape_qmap_analysis, &
        gamma_scale_qmap_analysis, fraction_zero_qmap_analysis, &
        quantile_table, alpha_values, cumprob_values, &
        data_directory, cyyyymmddhh, cmodel, cleade, &
        ensemble_ccpa_x25)

    ! ---- compute closest-member histogram on the quantile-mapped 
    !      forecast values, and save information to file.

    outfile_nc = TRIM(data_directory) // 'closest_hist_stats_'//&
        TRIM(cprefix)//'_'//cyyyymmddhh//'_fhour'//&
        TRIM(cleade)//'.nc'
    
    CALL compute_closest_histogram (n25, nxa, nya, nens, &
        thresh_light, thresh_mod, thresh_high, &
        ensemble_ccpa_x25, analysis, conusmask, outfile_nc, &
        istat)

ELSE 
    PRINT *,'unable to find valid ensemble data for this day.  Stopping.'
    STOP
ENDIF  ! maxval(ensemble_mean) .gt. 0.0 (good data)

! ===================================================================
! CLEAN UP AND QUIT
! ===================================================================
 		
DEALLOCATE(analysis, ensemble_ccpa_x25, ensemble_ccpa, &
    ensemble_mean, fraction_zero_qmap_analysis, &
    fraction_zero_qmap_forecast, gamma_shape_qmap_forecast, &
    gamma_scale_qmap_forecast, gamma_shape_qmap_analysis, &
    gamma_scale_qmap_analysis, rlonsa, rlatsa, stat=ios)

write(6,*) 'Deallocation Status = ',ios
write(6,*) 'Done!'

END PROGRAM generate_dressing_stats_anymodel_gammacdf
