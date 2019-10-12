SUBROUTINE control_quantile_mapping_singlemodel_gamma(nxa, nya, &
    nstride, nens, nens_qmap, n25, iadim, ipdim, &
    exchangeable, conusmask, ensemble_ccpa, &
    gamma_shape_qmap_forecast, gamma_scale_qmap_forecast, &
    fraction_zero_qmap_forecast, gamma_shape_qmap_analysis, &
    gamma_scale_qmap_analysis, fraction_zero_qmap_analysis, &
    quantile_table, alpha_values, cumprob_values, &
    data_directory, cyyyymmddhh, cmodel, cleade, &
    ensemble_ccpa_x25)
!   
INTEGER, INTENT(IN) :: nxa, nya ! grid dimensions 
INTEGER, INTENT(IN) :: nstride ! stride length when skipping grid pts
INTEGER, INTENT(IN) :: nens
INTEGER, INTENT(IN) :: nens_qmap
INTEGER, INTENT(IN) :: n25
INTEGER, INTENT(IN) :: iadim
INTEGER, INTENT(IN) :: ipdim
LOGICAL, INTENT(IN) :: exchangeable

INTEGER*2, INTENT(IN), DIMENSION(nxa,nya) :: conusmask ! conus mask

REAL, INTENT(IN), DIMENSION(nxa,nya,nens) :: ensemble_ccpa

REAL, INTENT(IN), DIMENSION(nxa, nya, nens_qmap) :: gamma_shape_qmap_forecast
REAL, INTENT(IN), DIMENSION(nxa, nya, nens_qmap) :: gamma_scale_qmap_forecast
REAL, INTENT(IN), DIMENSION(nxa, nya, nens_qmap) :: fraction_zero_qmap_forecast
REAL, INTENT(IN), DIMENSION(nxa, nya) :: gamma_shape_qmap_analysis
REAL, INTENT(IN), DIMENSION(nxa, nya) :: gamma_scale_qmap_analysis
REAL, INTENT(IN), DIMENSION(nxa, nya) :: fraction_zero_qmap_analysis
REAL*8, INTENT(IN), DIMENSION(iadim, ipdim) :: quantile_table
REAL*8, INTENT(IN), DIMENSION(iadim) :: alpha_values
REAL*8, INTENT(IN), DIMENSION(ipdim) :: cumprob_values
CHARACTER*(*), INTENT(IN) :: data_directory
CHARACTER*(*), INTENT(IN) :: cmodel
CHARACTER*(*), INTENT(IN) :: cleade
CHARACTER*(*), INTENT(IN) :: cyyyymmddhh

REAL, INTENT(OUT), DIMENSION(n25,nxa,nya,nens)  :: ensemble_ccpa_x25

REAL, DIMENSION(n25,nxa,nya) :: forecast_x25 ! work array

! --- and now the quantile mappings for the ensemble forecasts

PRINT *,'Subroutine control_quantile_mapping_singlemodel_gamma'
PRINT *,'Quantile mapping and dressing of ensemble, # members = ',nens
DO imem = 1, nens
    IF (exchangeable .eqv. .TRUE.) THEN
        imemout = 1
    ELSE 
        imemout = imem
    ENDIF
    PRINT *,'****  processing imem, imemout = ',imem

    CALL quantile_map_x25_gamma2(nxa, nya, nstride, n25, &
        nens_qmap, iadim, ipdim, conusmask, &
        gamma_shape_qmap_forecast(1,1,imemout), &
        gamma_scale_qmap_forecast(1,1,imemout), &
        fraction_zero_qmap_forecast(1,1,imemout), &
        gamma_shape_qmap_analysis(1,1), &
        gamma_scale_qmap_analysis(1,1), &
        fraction_zero_qmap_analysis(1,1), ensemble_ccpa(1,1,imem), &
        quantile_table, alpha_values, cumprob_values, &
        forecast_x25) 

    ensemble_ccpa_x25(:,:,:,imem) = forecast_x25(:,:,:)
END DO

RETURN
END SUBROUTINE control_quantile_mapping_singlemodel_gamma
