module PrecStru !precision and struct
!-------------------------------------------------------------------------------
! Purpose:
!  Define the structure 
!  throughout the model.
!-------------------------------------------------------------------------------
implicit none 
  integer, parameter :: r4 = selected_real_kind(5)
  integer, parameter :: r8 = selected_real_kind(8)
  integer, parameter :: i8 = selected_int_kind(13)
  integer, parameter :: c50 = 50
  integer, parameter :: c100 =220
  integer, parameter :: c1000 =4000
  !integer, parameter :: Tmin =240

  type :: SEPP
    !input================================================================= 
    character(c100)  EpFile       !  observed data file for ensemble members(shaake shuffle)  
    character(c100)  ObFile       !  observed data file   
    character(c100)  SiFile       !  Simulated data file   
    character(c100)  SiFile1      !  Simulated data file   
    character(c100)  PaFile       !  parameter Result data file   
    character(c100)  DoFile       !  domain file, merge    
    character(c100)  ReFile       !  output file, ensemble forecast    
    character(c100)  CoFile       !  Control file   
    character(c100)  EName        !  Name of Forcing 'tmp2m' 'prate'

    real(r8),pointer ::  ep(:,:) ! Observed data for ensemble members [date, lead time]
    real(r8),pointer ::  ob(:,:) ! Observed data [date, lead time]
    real(r8),pointer ::  si(:,:) ! Simulation [date, lead time]
    integer,pointer ::  Estart(:)  ! event start time in lead time
    integer,pointer ::  Estop(:)   ! event stop time in lead time

    real(r8),pointer ::   pthresh_obs(:,:)  ! precipitation threshold [365,Events]
    real(r8),pointer ::   pthresh_fcst(:,:) ! forecast threshold
    real(r8),pointer ::   popobs(:,:)   ! climatological probability of obs > pthreshold 
    real(r8),pointer ::   cavgobs(:,:)  ! climatological mean of obs when obs > pthreshold
    real(r8),pointer ::   ccvobs(:,:)   ! climatological cv of obs when obs > pthreshold
    real(r8),pointer ::   popfcst(:,:)  ! climatological probability of fcst > pthreshold
    real(r8),pointer ::   cavgfcst(:,:) ! climatological mean of fcst when fcst > pthreshold
    real(r8),pointer ::   ccvfcst(:,:)  ! climatological cv of fcst when fcst > pthreshold
    real(r8),pointer ::   rho_est(:,:)  ! recommended coefficient of correlation between z(fcst) and z(obs)
    !Parameter=============================================================
    integer x         ! column  number
    integer y         ! row number 
    integer x1,x2       !column  number
    integer y1,y2       !row number 
    integer flagOs    !select Operating system, 1-windows  2-linux 
    integer BeginYear
    integer EndYear 
    integer BeginYearEns
    integer EndYearEns

    integer iYear 
    integer imonth 
    integer iday 
    integer ncols  !column  total
    integer nrows  !row total 

    integer flagCG  !select read GFS data, 1-yes  0-not   
    integer flagT   !select read temperature data, 1-yes  0-not   
    integer leadT   ! Lead time 
    integer leadT1  ! Lead time 
	integer nparint !interval of calculation in 365 days
	integer nparint1  !interval of calculation in 365 days
    integer ndays     ! total time
    integer ndays1    ! total time for ensemble member of history observed data
    integer Events    ! Canonical Events numbers
    integer iobstran  ! observations probability distribution
    integer ifcsttran ! forecast probability distribution
    integer nems      ! number of ensemble members 

    !output================================================================
    real(r8),pointer ::   Eep(:,:)    !Events Observed data for ensemble members [date, Events]
    real(r8),pointer ::   Eob(:,:)    !Events Observed data [date, Events]
    real(r8),pointer ::   Esi(:,:)    !Events Simulation [date, Events]
    real(r8),pointer ::   Eens(:,:,:) !Events Ensemble [date,Events,members]
    real(r8),pointer ::   Pens(:,:,:) !Ensemble [date,lead time,members]
    real    ,pointer ::   PensArea(:,:,:,:) !Ensemble [Lon,Lat,lead time,members]
    integer ,pointer ::   indEep(:,:,:) !order of ensemble members [365,Events,members]
    integer ,pointer ::   indrho(:,:) !order of coefficient of correlation between z(fcst) and z(obs)
	                                    ![365,Events]

  end type SEPP


  type ::  CanonicalE  !Calculate the Canonical Events Paramater
    !input================================================================= 
    character(c100)  ObFile       !  observed data file   
    character(c100)  SiFile       !  Simulated data file  GFS  
    character(c100)  SiFile1      !  Simulated data file  CFS 
    character(c100)  ReFile       !  Result data file   
    character(c100)  DoFile       !  domain file, merge    
    character(c100)  CoFile       !  Control File   
    character(c50),pointer :: model(:)  !models name
    character(c50),pointer :: basin(:)  !basins number

    integer Events                !Canonical Events numbers
    integer,pointer ::  Estart(:)  ! event start time in lead time
    integer,pointer ::  Estop(:)   ! event stop time in lead time
    real(r8),pointer ::   ob(:,:)     !Observed data [date, lead time]
    real(r8),pointer ::   si(:,:)     !Simulation [date, lead time]
    !  real(r8),pointer ::   ob(:,:,:)       !Observed data [year,month,day]
    !  real(r8),pointer ::   si(:,:,:,:)     !Simulation [year,month,day, lead time]
    !Parameter=============================================================
    integer BeginYear
    integer EndYear 
    integer x       !column  number
    integer y       !row number 
    integer x1,x2       !column  number
    integer y1,y2       !row number 
    integer flagOs  !select Operating system, 1-windows  2-linux 
    integer flagPo  !select compute daily statistic, 1-yes  0-not   
    !integer PTmin   !Minum temperature in order to increase CV    
    integer flagSE  !select simulated/ensemble data, 1-simulated  0-ensemble   
    integer flagCG  !select read GFS data, 1-yes  0-not   
    integer flagT   !select read temperature data, 1-yes  0-not   
    integer leadT   !Lead time 
    integer leadT1   !Lead time 
    integer ndays   !total time
    integer ndays1   !total time  CFS 5 days
	integer nparint !interval of calculation in 365 days
    integer win     !windows wide of data
    integer maxwin  !max windows wide of data
    integer minwin  !min windows wide of data
    integer minall  !minimum number of all data 150
    integer minpos  !minimum number of data > pthresh_fcst 50
    integer nems    !  number of ensemble members 
    real(r8) pop_fraction !fraction of total precip > 0
    integer iobstran  ! observations probability distribution
    integer ifcsttran ! forecast probability distribution
    integer iopt_rho !correlation coefficient selection     
	integer verstats !1-check calibration 0-no check                     
    real(r8) cor_weight
    real(r8) cavgobsx  ! maximum conditional mean precipitation value (mm/6hr)50.8
    real(r8) cavgfcstx ! maximum conditional mean precipitation value (mm/6hr)50.8

    !output================================================================
    real(r8),pointer ::   Eob(:,:)    !Observed data [date, Events]
    real(r8),pointer ::   Esi(:,:)    !Simulation [date, Events]

    !integer    Nout          !number of out in a day
    real(r8),pointer ::   r(:,:)            !correlation coefficient [365,Events]
    real(r8),pointer ::   DC(:,:)           ! Nash-Sutcliffe efficiency
    real(r8),pointer ::   Bias(:,:)      !Balance
    real(r8),pointer ::   RMSE(:,:)         !Root-mean-square deviation

    integer,pointer ::    nall(:,:)         !number in window [365,Events] 
    real(r8),pointer ::   pthresh_obs(:,:)  ! precipitation threshold
    real(r8),pointer ::   pthresh_fcst(:,:) ! forecast threshold
    real(r8),pointer ::   avgobs(:,:)   ! average obervation
    real(r8),pointer ::   popobs(:,:)   ! climatological probability of obs > pthreshold 
    real(r8),pointer ::   cavgobs(:,:)  ! climatological mean of obs when obs > pthreshold
    real(r8),pointer ::   ccvobs(:,:)   ! climatological cv of obs when obs > pthreshold
    real(r8),pointer ::   avgfcst(:,:)  ! average forecast 
    real(r8),pointer ::   popfcst(:,:)  ! climatological probability of fcst > pthreshold
    real(r8),pointer ::   cavgfcst(:,:) ! climatological mean of fcst when fcst > pthreshold
    real(r8),pointer ::   ccvfcst(:,:)  ! climatological cv of fcst when fcst > pthreshold
    real(r8),pointer ::   rho_est(:,:)  ! recommended coefficient of correlation between z(fcst) and z(obs)
    real(r8),pointer ::   rmsfcst(:,:)  ! rms difference between fcst and obs
    real(r8),pointer ::   effnsfcst(:,:) ! Nash-Sutcliffe efficiency of forecast
    real(r8),pointer ::   ratio(:,:)     ! avgobs/avgfcst
    real(r8),pointer ::   eqts_fcst(:,:) ! equitable threat score for wet day forecast
    real(r8),pointer ::   cor(:,:)       ! correlation coefficient between obs and fcst (includes zeros)
    real(r8),pointer ::   ccor(:,:)      ! correlation coefficient between wet obs and fcsts
    real(r8),pointer ::   trccor(:,:)    ! correlation coefficient between transformed wet obs and fcsts
    real(r8),pointer ::   rhoopt(:,:)    ! correlation coefficient from transformed forecasts and obs
    integer,pointer ::    nvalcal(:,:,:) ! number [365,Events,4] 1-(<,<) 2-(<,>) 3-(>,<) 4-(>,>) 

    real(r8),pointer ::   avgobscal(:,:,:)   !average obervation
    real(r8),pointer ::   stdobscal(:,:,:)   !standard deviation of obervation
    real(r8),pointer ::   avgfcstcal(:,:,:)  !average forecast
    real(r8),pointer ::   stdfcstcal(:,:,:)  !standard deviation   forecast
        
    real(r8),pointer ::   coro(:,:,:) ![1-365, Events,Events]
	                                  !coeffient of correlation among Events of observed
    real(r8),pointer ::   corf(:,:,:) ![1-365, Events,Events]
	                                  !coeffient of correlation among Events of forecast
    real(r8),pointer ::   corX(:,:,:) ![1-365, Events,Events]
	                                  !coeffient of correlation among Events of forecast-observed
    real(r8),pointer ::   xm(:,:,:)   ![1-365, Events,Events] cross mean, mean(x*y)

    real(r8),pointer :: crps_avgobs(:,:) 
    real(r8),pointer :: crps_clim(:,:) 
    real(r8),pointer :: crps_fcst(:,:) 
    real(r8),pointer :: crps_ensavg(:,:) 
    real(r8),pointer :: crps_ens(:,:) 
    real(r8),pointer :: crps_ensmed(:,:) 
    real(r8),pointer :: crps_cavgobs(:,:) 
    real(r8),pointer :: crps_cclim(:,:) 
    real(r8),pointer :: crps_cfcst(:,:) 
    real(r8),pointer :: crps_censavg(:,:) 
    real(r8),pointer :: crps_cens(:,:) 
    real(r8),pointer :: crps_censmed(:,:) 
    real(r8),pointer :: roc_area(:,:) 

    real(r8),pointer ::biasfcst(:,:) 
    real(r8),pointer ::biasens(:,:) 
    real(r8),pointer ::corfcst(:,:) 
    real(r8),pointer ::corens(:,:) 
    real(r8),pointer ::corensfcst(:,:) 
    !real(r8),pointer ::rmsfcst(:,:) 
    real(r8),pointer ::rmsens(:,:) 
    !real(r8),pointer ::effnsfcst(:,:) 
    real(r8),pointer ::effnsens(:,:) 
    real(r8),pointer ::bss(:,:) 
    real(r8),pointer ::rmspop(:,:) 
    real(r8),pointer ::bsscpc1(:,:) 
    real(r8),pointer ::bsscpc2(:,:) 
    real(r8),pointer ::bsscpc3(:,:) 
    real(r8),pointer ::cdfdiffmax(:,:) 
    real(r8),pointer ::rmscdf(:,:) 
        

  end type CanonicalE

  type :: SgenerateCaE

    !input================================================================= 
    integer leadT   !Lead time 
    integer ndays   !total time
    integer Events  !Canonical Events numbers
	integer nparint !interval of calculation in 365 days
    real(r8),pointer ::   ob(:,:)     !climate data [date, lead time]
    real(r8),pointer ::   hob(:,:)    !historical climate data [date, lead time]
    integer,pointer ::  Estart(:)  ! event start time in lead time
    integer,pointer ::  Estop(:)   ! event stop time in lead time
    !output================================================================
    real(r8),pointer ::   Eob(:,:)    !Canonical Events data [date, Events]

  end type SgenerateCaE

  type :: SreadEppPara
    !input================================================================= 
    character(c100)  PaFile  !  Parameter file   
    integer Events              !Canonical Events numbers
    integer x       ! column  number
    integer y       ! row number 
    !output================================================================
    real(r8),pointer ::   pthresh_obs(:,:)  ! precipitation threshold
    real(r8),pointer ::   pthresh_fcst(:,:) ! forecast threshold
    real(r8),pointer ::   popobs(:,:)   ! climatological probability of obs > pthreshold 
    real(r8),pointer ::   cavgobs(:,:)  ! climatological mean of obs when obs > pthreshold
    real(r8),pointer ::   ccvobs(:,:)   ! climatological cv of obs when obs > pthreshold
    real(r8),pointer ::   popfcst(:,:)  ! climatological probability of fcst > pthreshold
    real(r8),pointer ::   cavgfcst(:,:) ! climatological mean of fcst when fcst > pthreshold
    real(r8),pointer ::   ccvfcst(:,:)  ! climatological cv of fcst when fcst > pthreshold
    real(r8),pointer ::   rho_est(:,:)  ! recommended coefficient of correlation between z(fcst) and z(obs)

  end type SreadEppPara

  type :: SreadData
    !input================================================================= 
    character(c100)  EpFile       !  observed data file for ensemble members(shaake shuffle)  
    character(c100)  ObFile       !  observed data file   
    character(c100)  SiFile       !  Simulated data file   
    character(c50),pointer :: model(:)  !models name
    character(c50),pointer :: basin(:)  !basins number

    !Parameter=============================================================
    integer flagOs  !select Operating system, 1-windows  2-linux 
    integer flagT   !select read temperature data, 1-Tmax  2-Tmin  0-not   
    integer x       !column  number
    integer y       !row number 
    integer BeginYear
    integer EndYear 
    integer BeginYearEns
    integer EndYearEns
	integer Bday
	integer Bmonth
    integer leadT   !Lead time 
    integer ndays   !total time
    integer ndays1   !total time CFS 5 days
	integer nparint !interval of calculation in 365 days
    !output================================================================
    real(r8),pointer ::   ob(:,:)    ! Observed data [date, lead time]
    real(r8),pointer ::   si(:,:)    ! Simulated [date, lead time]
    real(r8),pointer ::   ep(:,:)    ! mean data for ensemble members [date, lead time]
    !real(r8),pointer ::   obflow(:)! Observed stream flow [date]
    !real(r8),pointer ::   siflow(:,:)! Simulated stream flow [date, lead time]
    !real(r8),pointer ::   poflow(:,:)! post processing stream flow [date, lead time]

  end type SreadData

  type :: Swindow
    !input================================================================= 
    real(r8),pointer ::   Eob(:,:)    !Observed data [date, Events]
    real(r8),pointer ::   Esi(:,:)    !Simulation [date, Events]
    !Parameter=============================================================
    integer BeginYear
    integer EndYear 
    integer ndays    !total time
    integer ne       !Canonical Events number
    integer day      !day 1-365
    integer maxwin   !max windows wide of data
    integer minwin   !min windows wide of data
    integer minall   !minimum number of all data 150
    integer minpos   !minimum number of data > pthresh_fcst 50
    real(r8) pop_fraction !fraction of total precip > 0
    !output================================================================
	real(r8) pthresh_obs  ! precipitation threshold
	real(r8) pthresh_fcst ! forecast threshold
    integer  win          ! windows wide of data
    integer  nall         ! number in window
    real(r8),pointer ::   obs(:) !Observed data [date] of an Events
    real(r8),pointer ::   fcst(:)!Simulation [date] of an Events

    integer ,pointer ::   Enall(:,:) !number in window [1-365, Events] of an Events
    real(r8),pointer ::   Eobs(:,:,:) !Observed data [1-365, Events,nall] of an Events
    real(r8),pointer ::   Efcst(:,:,:)!Simulation [1-365, Events,nall] of an Events
    real(r8),pointer ::   corf(:,:,:) ![1-365, Events,Events]
	                                  !coeffient of correlation among Events of forecast
    real(r8),pointer ::   corX(:,:,:) ![1-365, Events,Events]
	                                  !coeffient of correlation among Events of forecast-observed
    
  end type Swindow

  type :: Sensverify  !ensemble verify
    !input================================================================= 
     
  
    integer   nobs ! number of observations/forecasts  
	integer iobstran  ! observations probability distribution
	integer ifcsttran ! forecast probability distribution
	                  ! 0  no  transformation   1 GAMM = Gamma
                      ! 2 LOGN = Log Normal     3 EXPO = Exponetial
	                  ! 4  normal  distribution 5 WEIB = Weibull
    real(r8)  pthresh ! threshold value to define conditional event (obs>pthresh or xppx>pthresh)
                
    real(r8),pointer :: obs(:) ! vector of observed values, obs(i), i!1,nobs
    real(r8),pointer :: fcst(:)!Simulation [date] of an Events
    integer  nmem  ! vector specifying number of ensemble members in each forecast,
                            ! nmem(i), i=1,nobs
    real(r8),pointer ::  xppx(:,:) ! array of forecast member values, xppx(j,i), j=1,nmem(i), i=1,nobs
    !output================================================================
    real(r8),pointer :: avgpx(:) ! ensemble average for forecast i (unconditional)
    real(r8),pointer :: stdpx(:) ! ensemble standard deviation of forecast i (unconditional)
    real(r8),pointer :: cavgpx(:) ! conditional average of forecast i (> pthresh)
    real(r8),pointer :: cstdpx(:) ! conditional stand deviation of forecast i (> pthresh)
    real(r8),pointer :: ccvpx(:) ! conditional coefficient of variation of forecast i  (> pthresh)
                 
    real(r8),pointer :: poppx(:) ! propability of precipitation of forecast i
 
    real(r8) avgobs ! average value of obs(i)
    real(r8) stdobs ! standard deviation of obs(i)
    real(r8) cvobs ! coefficient of variation of obs(i)
    real(r8) cavgobs ! conditional average value of obs(i) (obs>pthresh)
    real(r8) cstdobs ! conditional stanedard deviation of obs(i) (obs>pthresh)
    real(r8) ccvobs ! conditional coefficient of variation of obs(i) (obs>pthresh)
    real(r8) popobs ! observed probability of precipitation

    real(r8) avgens ! average of the ensemble forecast averages
    real(r8) stdens ! standard deviation of the ensemble forecast averages
    real(r8) cvens ! coefficient of variation of the ensemble forecast averages
    real(r8) popens ! fraction of avgpx(i)>pthresh
    real(r8) cavgens ! average of the conditional ensemble forecast averages
    real(r8) cstdens ! standard deviation of the conditional ensemble forecast averages              
    real(r8) ccvens ! coefficient of variation of the conditional ensemble forecast averages             
    real(r8) avgpopens ! average of the ensemble forecast pop
    real(r8) stdpopens ! standard deviation of ensemble forecast pop
    real(r8) avgcavgens ! average of the ensemble conditional average
    real(r8) stdcavgens ! standard deviation of the ensemble conditional average
    real(r8) avgccvens ! average of the ensemble forecast ccv
    real(r8) stdccvens ! standard deviation of the ensemble forecast ccv
    real(r8) biasens ! (avgens - avgobs)/avgobs + 1
    real(r8) cbiasens ! (cavgens - cavgobs)/cavgobs + 1
    real(r8) corens ! coeffient of correlation between ensemble mean and observed value
    real(r8) corensfcst ! coeffient of correlation between ensemble mean and forecast value
               
    real(r8) ccorens ! coefficient of correlation between conditional ensemble 
                 ! mean and observed value, given obs>pthresh
    real(r8) tccorens ! coefficient of correlation between transformed 
                  !  conditional ensemble mean and transformed observed value,
                  !  given obs>pthresh
    real(r8) rhoopt ! contingency estimate of bivariate normal correlation 
                ! parameter for avgpx(i)
    real(r8) effens ! Nash-Sutcliffe efficiency of ensemble mean forecast
    real(r8) ceffens ! Nash-Sutcliffe efficiency of ensemble conditional mean 
                 !forecast, given obs>pthresh
    real(r8) rmsens ! rms error of ensemble mean
    real(r8) crmsens ! real(r8)rms error of conditional ensemble mean, given obs>pthresh

    real(r8) avgensmem ! average value of all ensemble members
    real(r8) stdensmem ! standard deviation of all ensemble members
    real(r8) cvensmem ! coefficient of variation of all ensemble members

    real(r8) cavgensmem ! conditional average value of all ensemble members
    real(r8) cstdensmem ! conditional standard deviation of all ensemble members
    real(r8) ccvensmem ! conditional coefficient of variation of all ensemble members
                    
    real(r8) popensmem ! pop of all ensemble memmbers
    real(r8) biasmem ! (avgensmem - avgobs)/avgobs + 1
    real(r8) cbiasmem ! (cavgensmem - cavgobs)/cavgobs + 1

    real(r8) bs ! brier score of ensemble pop
    real(r8) bss ! brier skill score of ensemble pop
    real(r8) rmspop ! rms difference between avgpclass and cobspop
    real(r8) popclass(4) ! limits of pop categories (increasing from 0 to 1) k=1,4
    integer  npclass(4) ! number of forecasts in category k
    real(r8) avgpclass(4) ! average forecast pop in class k
    integer  nposobs(4) ! number of obs>pthresh given pop fcst in category k
    real(r8) cobspop(4) ! obs avg pop in category k, (nposobs(:)/npclass(:))
    real(r8) expclass(4) ! class limit of conditional cdf category k, k=1,4
    real(r8) ctrpclass(4) ! centroid value of category k
    integer nexpclass(4) ! number of forecasts when obs fell in conditional 
                             ! forecast proability category k
    real(r8) avgexpclass(4) ! average value of forecast probability in category k
    real(r8) fracpclass(4) ! fraction of forecasts falling in category k
    real(r8) cdfpclass(4) ! cumulative fraction of forecasts falling in category k
    real(r8) cdfdiffmax ! maximum difference between cdfpclass and uniform 
                  ! distribution (! Kolmogov statistic)
    real(r8) rmscdf ! rms difference between cdfdiffmax and uniform distribution
    real(r8) eq_ts ! equitable threat score for avgpx(:)>pthresh
    real(r8),pointer ::  cumestprob(:) ! est probability of event being less than or equal to obs(:)
                   
    integer npopobs ! number of obs > pthresh
    real(r8),pointer ::  xobs(:) ! list of obs > pthresh (i=1,npopobs)
    real(r8),pointer ::  cdfobs(:) ! cum dist fn of xobs(:) (i=1,npopobs)
    real(r8) bsscpc(3) ! brier skill score for cpc tercile conditional 
                  ! prob fcsts (k!1,3)
    real(r8),pointer ::  icatobs(:) ! tercile category (1,2,3) of the i-th observation
    integer nhits33 ! number of hits of tercile categorical forecasts
    integer nfcsts33 ! number of hits of tercile categorical forecasts

    real(r8),pointer ::  prhit33(:,:) ! forecast prob of i-th positive obs being in k-th tercile
                      
    real(r8) heidke ! heidke skill score
    real(r8) sskuipers ! Kuiper's skill score for categorical preditiction of 
                   !  positive obs
    real(r8) prob_err(3) ! transformed rank histogram quartile value
                      !  k=1 => bias error
                      !  k=2 => spread error
                      !  k=3 => random error

    real(r8) pval(2,2) ! fraction of events, class = i,j
    real(r8) avgobsval(2,2) ! average observation, class = i,j
    real(r8) stdobsval(2,2) ! std dev observation, class = i,j
    real(r8) avgensval(2,2) ! average ens mean, class = i,j
    real(r8) stdensval(2,2) ! std dev ens mean, class = i,j
  end type Sensverify

  type :: Sroc
    !input=================================================================
    integer nfcst          !  number of ensemble forecasts w/ obs
    integer   nmem         !  actual number of members in array xx
    integer  nthresh       !  number of observations thresholds  to be used
    real(r8),pointer :: xx(:,:)!  xx(j,i) = j-th member value for ensemble forecast i 
    real(r8),pointer :: obs(:) !  obs(i) = value of observation for i-th forecast
    !output================================================================
    real(r8),pointer ::obsthresh(:) ! obsthresh(i) = observation threshold (i=1,nthresh)
    real(r8)   roc_area          ! area of skill for relative operating characteristics curve 
    real(r8)   roc_score         ! relative roc skill = roc_area/0.5
    integer  istatus         ! error status (=0, no errors)
    real(r8),pointer :: hr(:)    ! hr(i) = hit rate for i-th observation threshold 
    real(r8),pointer :: far(:)   ! far(i) = false alarm rate for i-th observation threshold 
  end type Sroc

  type :: Scrpseval
    !input=================================================================
    integer nobs
    real(r8),pointer :: fcst(:)
    real(r8),pointer :: obs(:)
    real(r8),pointer :: ensem(:,:)
	integer  nmem
	real(r8) pthresh 
    !output================================================================
    real(r8) crps_avgobs,crps_clim,crps_fcst,crps_ensavg,crps_ens,crps_ensmed,&
         crps_cavgobs,crps_cclim,crps_cfcst,crps_censavg,crps_cens,crps_censmed                                      
                           
  end type Scrpseval

  type :: Scrps
    real(r8) x0
    integer n
    real(r8),pointer :: xx(:)
  end type Scrps

  type :: Sextractp
    !input=================================================================
    integer         nv       ! number of values of input distribution of p(.)
    integer         npp      ! number of values of pp(.)
    real(r8),pointer :: p(:)     ! input conditional mean variate in interval <j-1,j>
    real(r8),pointer :: ccdfv(:) ! cdf of p (cdf at end of interval j)
	real(r8)        obspz    ! probability that p ! 0.
    !output================================================================
    real(r8),pointer :: pp(:)    ! output conditional mean variate in interval i:i!1,npp

  
  end type Sextractp

  type :: Sfcst2ensem
    !input=================================================================

    real(r8) fcst    ! deterministic forecast
    integer  ifcst   ! climatological forecast distribution
    real(r8) fpop    ! climatological probability of fcst > pthreshold
    real(r8) fcavg   ! climatological mean of fcst when fcst > pthreshold
    real(r8) fccv    ! climatological cv of fcst when fcst > pthreshold
    integer  iobs ! climatological observations distribution
    real(r8) obspop  ! climatological probability of obs > pthreshold
    real(r8) obscavg ! climatological mean of obs when obs > pthreshold
    real(r8) obsccv  ! climatological cv of obs when obs > pthreshold
    real(r8) rho     ! coefficient of correlation between z(fcst) and z(obs)
    integer  npp     ! number of values of pp (number of ensemble members)
	real(r8) pthresh_obs  ! precipitation threshold
	real(r8) pthresh_fcst  ! forecast threshold

    !output================================================================
    real(r8) genavg  ! average value of forecast ensemble distribution (internal)
    real(r8) genpop  ! pop of forecast ensemble distribution (internal)
    real(r8) gencavg ! conditional mean of fcst ensemble distribution (internal)
    real(r8) genccv  ! conditional cv of forecast ensemble distribution (internal)
    real(r8) vmax    ! maximum value of forecast ensemble distribution (internal)
    real(r8) exprobmin ! excedence probability associated with vmax
                   !note:  The output ensemble distribution is pp and has npp members.  
                   !This is computed from a more detailed internal distribution 
                   !with nvmax members that has a maximum value of vmax
    real(r8),pointer :: pp(:) ! vector of precipitation values sorted with pp(1) = smallest value
    real(r8) pmin  ! value climatological cdf probability of observing a value 
                   !corresponding to the first (lowest) point on the internal 
                   !ensemble forecast distribution
    real(r8) pmax  ! value climatological cdf probability of observing a value 
                   !corresponding to the last (highest) point on the internal 
                   !ensemble forecast distribution                                 

  end type Sfcst2ensem


  type :: Sfcstparm
    !input=================================================================
    !verstats = verification statistics switch (compute if 'yes') 
    integer nobs   ! number of observations = number of forecasts
	real(r8) pthresh1  ! precipitation threshold
	real(r8) pthresh2  ! forecast threshold
	integer iobstran  ! observations probability distribution
	integer ifcsttran ! forecast probability distribution
	                  ! 0  no  transformation   1 GAMM = Gamma
                      ! 2 LOGN = Log Normal     3 EXPO = Exponetial
	                  ! 4  normal  distribution 5 WEIB = Weibull
    real(r8),pointer :: obsx(:) ! observed values  
    real(r8),pointer :: fcstx(:)! deterministic forecasts
    real(r8)  cor_weight
	integer iopt_rho !correlation coefficient selection                          
    !output================================================================
    real(r8) cavgobsx
    real(r8) cavgfcstx
    real(r8) avgobs   ! average obervation
    real(r8) popobs   ! climatological probability of obs > pthreshold 
    real(r8) cavgobs  ! climatological mean of obs when obs > pthreshold
    real(r8) ccvobs   ! climatological cv of obs when obs > pthreshold
    real(r8) avgfcst  ! average forecast 
    real(r8) popfcst  ! climatological probability of fcst > pthreshold
    real(r8) cavgfcst ! climatological mean of fcst when fcst > pthreshold
    real(r8) ccvfcst  ! climatological cv of fcst when fcst > pthreshold
    real(r8) rho_est  ! recommended coefficient of correlation between z(fcst) and z(obs)
    real(r8) rmsfcst  ! rms difference between fcst and obs
    real(r8) effnsfcst ! Nash-Sutcliffe efficiency of forecast
    real(r8) ratio     ! avgobs/avgfcst
    real(r8) eqts_fcst ! equitable threat score for wet day forecast
    real(r8) cor       ! correlation coefficient between obs and fcst (includes zeros)
    real(r8) ccor      ! correlation coefficient between wet obs and fcsts
    real(r8) trccor    ! correlation coefficient between transformed wet obs and fcsts
    real(r8) rhoopt    ! correlation coefficient from transformed forecasts and obs
    integer nvalcal(2,2)
    real(r8) avgobscal(2,2)  
    real(r8) stdobscal(2,2)  
    real(r8) avgfcstcal(2,2)  
    real(r8) stdfcstcal(2,2)     
  end type Sfcstparm


  type :: Scgauss
    !input=================================================================
    real(r8) u     != given value of correlated standard normal deviate
    real(r8) u0    != value of u corresponding to g(u) = 0, u<=u0
    real(r8) rho   != coefficient of correlation between u and v
    integer nv != number of values of vval
    real(r8) vmin  != minimum value of vval range
    real(r8) vmax  != maximum value of vval range
    !output================================================================
    real(r8),pointer :: vval(:)  != values of standard normal deviate v
    real(r8),pointer :: cpdfv(:) != probability density function of v
    real(r8),pointer :: ccdfv(:) != cumulative distribution function of v
  end type Scgauss

  type :: Sbivar
    !input=================================================================
	real(r8) ax,bx,cx,tol,rho22
    real(r8) p11,p12,p21,p22,umin,u0,umax,vmin,v0,vmax
	!     p11 = prob v<=0 if u<=0
    !     p12 = prob v>0 if u<=0
    !     p21 = prob v<=0 if u>0
    !     p22 = prob v>0 if u>0
    !output================================================================
    real(r8) rho
  end type Sbivar

  type :: Svtrans
    !input=================================================================
    real(r8) avg,cv,pz,obs
	integer ivopt
    !output================================================================
    !real(r8) vtrans
  end type Svtrans

  type :: Sbeta
    !input=================================================================
    real(r8) a,b,f
    !output================================================================
    real(r8) beta
  end type Sbeta

  type :: Sthreshold
    !input=================================================================
    real(r8),pointer :: obs(:)  !observations array      
    !para==================================================================
    integer nobs            !number of observations
    real(r8) pop_fraction       !fraction of total precip > 0
    !output================================================================
    real(r8) pthresh            !precipitation threshold
  end type Sthreshold


  type :: Sgammcf  ! gammcf,a,x,gln
    !input=================================================================
    real(r8)   a
    !para==================================================================
    real(r8) x
    !output================================================================
    real(r8) gammcf,gln
  end type Sgammcf


  type :: Ssort  !sort a
    !input=================================================================
    real(r8),pointer :: a(:)
    !para==================================================================
    integer n
    !output================================================================
    integer,pointer :: indx(:) 
  end type Ssort

  type :: Scaldat   
    !input=================================================================
    integer julian
    !para==================================================================
    !output================================================================
   integer mm,id,iyyy
  end type Scaldat

  type ::  CRelEffP  !Calculate the Nash-Sutcliffe efficiency  and correlation coefficient 
    real(r8),pointer ::   yy(:)       !Observed
    real(r8),pointer ::   yc(:)       !Simulation
    integer     TimeSum           !total time
    integer     aSum              !Accumulation time
    real(r8)   r                             !correlation coefficient
    real(r8)   DC                            ! Nash-Sutcliffe efficiency
    real(r8)   Bias                          !Balance
    real(r8)   xm                            !cross mean, mean(yy*yc)

    real(r8)   RMSE                          !Root-mean-square deviation
    real(r8)   sr                            !Sliding correlation coefficient
    real(r8)   sDC                           !Sliding Nash-Sutcliffe efficiency
    real(r8)   sBias                     !Sliding Balance
    real(r8)   sRMSE                         !Sliding Root-mean-square deviation
    real(r8)   ar                            !Cumulative correlation coefficient
    real(r8)   aDC                           !Cumulative Nash-Sutcliffe efficiency
    real(r8)   aBias                     !Cumulative Balance
    real(r8)   aRMSE                         !Cumulative Root-mean-square deviation
  end type CRelEffP

  type :: stat24
    !input=================================================================
    real,pointer :: Obdata(:,:,:)       !matrix of observed data [year,month,day]  
    !para==================================================================
    integer    BeginYear
    integer    EndYear 
	real       pthresh
    !output================================================================
    real       avg(12),pop(12),cavg(12),ccv(12)
    integer    npos(12),nobs(12)
  end type stat24

  type :: stat06
    !input=================================================================
   real,pointer :: Obdata06(:,:,:,:)       !matrix of observed data [year,month,day,4]    
    !para==================================================================
    integer    BeginYear
    integer    EndYear 
	real       pthresh
    !output================================================================
    real  avg(4,12),pop(4,12),cavg(4,12),ccv(4,12)
    integer*4 npos(4,12),nobs(4,12)
  end type stat06

 
  type ::  Obsdata      
    character(c100)  BName       !Name of file   
    character(c100)  obw         !Observed revise weight of file  
    character(c100)  BName06       !Name of 6 hours file 
    character(c100)  CoFile       !  Control File 
    integer    flagOs      ! select Operating system, 1-windows  2-linux 
    !integer    flagF       !select forcing 1-precipitation 2-Tmax 3-Tmin
    integer    SumN             !days  total 
    integer    BeginYear
    integer    EndYear 
    integer    ncols            !column  total
    integer    nrows            !row total 
    integer    NVara            !Numbers of variable
    real       undef            !no data vale   
    integer    x                !column  number
    integer    y                !row number 
    integer    Tim              !read total times
    integer    CuVa             !Current  variable 
    integer    k                !read times

    real,pointer :: Obdata(:,:,:,:,:)       !matrix of observed data [x,y,year,month,day]  
    real,pointer :: Obdata06(:,:,:,:)       !matrix of observed data [year,month,day,4]     
       
    real  avg(5,12),pop(5,12),cavg(5,12),ccv(5,12) 
      !average, >0.25   Probability, >0.25 average,Coefficient of variation
    integer*4 npos(5,12),nobs(5,12)
    ! day number of precipitaion>0.25, day number of precipitaion>0
  end type Obsdata 

 
  type ::  GFSdata      
    character(c100)  BName       !Name of input data file
    character(c100)  Weightfile  !Weight file for between two data
    character(c100)  MName       !Name of merge file
    character(c100)  CoFile       !  Control File 
    integer    flagOs      ! select Operating system, 1-windows  2-linux 
    integer    SumN             !grid  total 
    integer    ncols            !column  total
    integer    nrows            !row total 
    integer    leadT            !Lead time 
    real(r8) cellsize          
    real(r8) NODATA_value     
    real(r8) xllcorner       
    real(r8) yllcorner   
    integer    BeginYear
    integer    EndYear 
    integer    ndays
    real(r8),pointer :: Gdata(:,:,:,:)       !matrix of GFS [year,month,day,leadT]        
    real(r8),pointer :: Gdatat(:,:,:,:,:,:)       !matrix of GFS [Lon,lat,year,month,day,leadT]        
  end type GFSdata 



  type ::  Stack      
    integer    X                !X coordinate
    integer    Y                !Y coordinate
  end type Stack 

  type ::  MWnet     
    type(Stack),pointer :: net(:)
  end type MWnet 

 

  type ::  ASCFile    !ASCII file
    integer    flagAF           !judge whether record exist, 1-Yes ,0-No
    integer    SumN             !Grid total
    integer    ncols            !column  total
    integer    nrows            !row total 
    real(r8) cellsize        
    real(r8) NODATA_value  
    real(r8) xllcorner    
    real(r8) yllcorner        
    real(r8) ,pointer :: MValue(:,:)         !matrix of value
  end type ASCFile 


  type ::  CFSData    !CFS data
    character(c100)  ObFile       !  observed data file  
    character(c100)  CFile       !Name of coordinate file   
    character(c100)  DFile       !Name of data file   
    character(c100)  WFile       !Name of weight file   
    character(c100)  OutFile     !Name of output data file   
    character(c100)  MName       !Name of merge file
    character(c100)  EName       !Name of Forcing 'tmp2m' 'prate'
    character(c100)  CoFile       !  Control File 
    integer    ifile            !file no.
    logical    alive            !file exist
    integer    flagOs           !select Operating system, 1-windows  2-linux 
    integer    flagAF           !judge whether record exist, 1-Yes ,0-No
    integer    ntime            !time total
    integer    nstep            !time step /hours
    integer    ncols            !column  total
    integer    nrows            !row total     
	integer    ix              !point longitude   
	integer    iy              !point latitude
	integer    Crow            !current row in binary file
	integer    Crow1            !current row in binary file
    integer    x1,x2       !column  number
    integer    y1,y2       !row number 

    integer    leadT           !Lead time 
    integer    BeginYear
    integer    EndYear 
    integer    iYear 
    integer    imonth 
    integer    iday 
    integer    cday           !current day 1-365
    integer    ensm           !current member 1-4
    integer    ihour 
    integer    monthday(12,32)

    integer    ndays
	integer    nparint !interval of calculation in 365 days
	  
    real(r8) NODATA_value  
    real(r8) ,pointer :: x(:)    
    real(r8) ,pointer :: y(:)       
	type(Stack), pointer :: position(:,:) !(ix,iy) point position in raw CFS data 
    type(ASCFile)  ::  obw
    integer  ,pointer :: mask(:,:) !(384,190)!mask of each CFS grid
    real  ,pointer :: w(:,:,:) !(464,224,4)!weight of least 4 points (lon, lat,4) real   w(464,224,4)
    real  ,pointer :: MValue(:,:,:) !CFS  raw value (lead time, lon, lat)
    real  ,pointer :: MValue1(:,:,:) !CFS  raw value (lead time, lon, lat)
    real  ,pointer :: MValue00(:,:,:) !CFS  raw value (lead time, lon, lat)
    real  ,pointer :: MValue06(:,:,:) !CFS  raw value (lead time, lon, lat)
    real  ,pointer :: MValue12(:,:,:) !CFS  raw value (lead time, lon, lat)
    real  ,pointer :: MValue18(:,:,:) !CFS  raw value (lead time, lon, lat)
    real  ,pointer :: PValue(:)   !CFS new point value ( lead time)
    real  ,pointer :: MValueN(:,:,:)   !CFS new point value ( lead time, lon, lat)
    real  ,pointer :: MValueN1(:,:,:)   !CFS new point value ( lead time, lon, lat)

    integer  ,pointer :: PLValue(:,:,:)   !CFS file point value (time,lat,lead time) (time,4,lead time)
    character(10),pointer :: PIValue(:)   !CFS file point value ( lead time)
    real(r8),pointer ::   cor(:,:,:)       ! correlation coefficient between obs and fcst members (includes zeros)

  end type CFSData 



end module PrecStru
