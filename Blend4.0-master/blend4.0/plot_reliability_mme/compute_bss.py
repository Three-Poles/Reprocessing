""" 
    compute_bss.py 
    this routine will compute Brier skill scores
    for raw, multi-model, quantile-mapped forecasts.  
    Data saved to cPickle files
"""


import numpy as np
from numpy import ma
import os, sys
from netCDF4 import Dataset
print 'importing netCDF4'
from verify_bss_raw_qmapped_dressed import verify_bss_raw_qmapped_dressed
print verify_bss_raw_qmapped_dressed.__doc__
import pygrib
from dateutils import daterange, dateshift
import cPickle

# ---- read in precipitation analyses

cftype = ''
    
cyyyymmddhh_start = '2016040100'
cyyyymmddhh_end = '20160707000'
date_list = daterange(cyyyymmddhh_start, cyyyymmddhh_end, 24)


nxa = 464 # CCPA 1/8 degree grid over CONUS
nya = 224 #

date_list_anal = daterange(cyyyymmddhh_start, cyyyymmddhh_end, 24)
print 'date_list_anal = ',date_list_anal
ndates = len(date_list_anal)
apcp_anal_t = np.zeros((nxa,nya,ndates), dtype=np.float32)
mninetynine = -99.99*np.ones((nya,nxa), dtype=np.float32)

for cleade in ['24','48','72','96','120','144','168']:
    ileade = int(cleade)
    ileadb = ileade - 12
    cleadb = str(ileadb)
    
    for idate, date in zip(range(ndates), date_list_anal):

        print '------------------------- getting precipitation for idate = ',idate,date
        date_fearly = dateshift(date, ileade-6)
        date_flate = dateshift(date, ileade)

        # --- read in the first of the two precip analysis files, 6 hourly
    
        infile = '/data/thamill/Rf2_tests/ccpa_v1/0.125d/ccpa.'+date_fearly[0:8]+\
            '/18/ccpa.t18z.06h.0p125.conus.gb2'
        print infile
        fexist1 = os.path.exists(infile)
        if fexist1:
            afile1 = pygrib.open(infile)
            grb1 = afile1.select()[0]    # --- read in the first of the two precip analysis files, 6 hourly
    
        infile2 = '/data/thamill/Rf2_tests/ccpa_v1/0.125d/ccpa.'+date_flate[0:8]+\
            '/00/ccpa.t00z.06h.0p125.conus.gb2'
        print infile2
        fexist2 = os.path.exists(infile2)
        if fexist2:
            afile2 = pygrib.open(infile2)
            grb2 = afile2.select()[0]
        if fexist1 and fexist2:
            apcp_anal = grb1.values + grb2.values
            apcp_anal = np.where(apcp_anal > 500., mninetynine, apcp_anal)
            afile1.close()
            afile2.close()
        else:
            print 'no analysis data for this date'
            apcp_anal = -99.99*np.ones((nya,nxa),dtype=np.float32)
        apcp_anal_t[:,:,idate] = np.transpose(apcp_anal[:,:])
    
    print '+++++++++++ Processing lead = ',cleade
    
    for cthresh in ['POP','10.0']:
                    
        print '   ++++++  Processing cthresh = ', cthresh
        if cthresh == 'POP':
            rthresh = 0.254
            cthresh_plot = 'POP'
        else:
            rthresh = float(cthresh)
            cthresh_plot = r'$\geq$ '+cthresh+' mm'

        # ---- call fortran routine to do the reading in of forecasts and the generation
        #      of reliability, frequency of use, and bss information.

        validdays = np.ones((ndates),dtype=np.int32)

        # NCEP declarations

        bs_NCEP_raw_daily = np.zeros(ndates,dtype=np.float32)
        bs_NCEP_qmapped_daily = np.zeros(ndates,dtype=np.float64)
        bs_NCEP_dressed_daily = np.zeros(ndates,dtype=np.float64)
        bs_climo_daily = np.zeros(ndates,dtype=np.float64)

        bss = 0.
        bss_NCEP_raw = 0.
        bss_NCEP_qmapped = 0.
        bss_NCEP_dressed = 0.

        # CMC declarations

        bs_CMC_raw_daily = np.zeros(ndates,dtype=np.float64)
        bs_CMC_qmapped_daily = np.zeros(ndates,dtype=np.float64)
        bs_CMC_dressed_daily = np.zeros(ndates,dtype=np.float64)

        bss_CMC_raw = 0.
        bss_CMC_qmapped = 0.
        bss_CMC_dressed = 0.

        # ECMWF declarations

        bs_ECMWF_raw_daily = np.zeros(ndates,dtype=np.float64)
        bs_ECMWF_qmapped_daily = np.zeros(ndates,dtype=np.float64)
        bs_ECMWF_dressed_daily = np.zeros(ndates,dtype=np.float64)

        bss_ECMWF_raw = 0.
        bss_ECMWF_qmapped = 0.
        bss_ECMWF_dressed = 0.

        # MME declarations

        bs_MME_raw_daily = np.zeros(ndates,dtype=np.float64)
        bs_MME_qmapped_daily = np.zeros(ndates,dtype=np.float64)
        bs_MME_dressed_daily = np.zeros(ndates,dtype=np.float64)

        bss_MME_raw = 0.
        bss_MME_qmapped = 0.
        bss_MME_dressed = 0.
        
        cmodel = 'MME'
        print 'Processing MME'
        print verify_bss_raw_qmapped_dressed.__doc__
        bss_MME_raw, bss_MME_qmapped, bss_MME_dressed, \
            bs_MME_raw_daily, bs_MME_qmapped_daily, bs_MME_dressed_daily, bs_climo_daily = \
            verify_bss_raw_qmapped_dressed(cleade, cmodel, cftype, rthresh, \
            date_list_anal, apcp_anal_t, validdays, nxa, nya, ndates)   
    
        cmodel = 'NCEP'
        print 'Processing NCEP'
        bss_NCEP_raw, bss_NCEP_qmapped, bss_NCEP_dressed, \
            bs_NCEP_raw_daily, bs_NCEP_qmapped_daily, bs_NCEP_dressed_daily, bs_climo_daily  = \
            verify_bss_raw_qmapped_dressed(cleade, cmodel, cftype, rthresh, \
            date_list_anal, apcp_anal_t, validdays, nxa, nya, ndates)       

        cmodel = 'CMC'
        print 'Processing CMC'
        bss_CMC_raw, bss_CMC_qmapped, bss_CMC_dressed, \
            bs_CMC_raw_daily, bs_CMC_qmapped_daily, bs_CMC_dressed_daily, bs_climo_daily  = \
            verify_bss_raw_qmapped_dressed(cleade, cmodel, cftype, rthresh, \
            date_list_anal, apcp_anal_t, validdays, nxa, nya, ndates)
    
        cmodel = 'ECMWF'
        print 'Processing ECMWF'
        bss_ECMWF_raw, bss_ECMWF_qmapped, bss_ECMWF_dressed, \
            bs_ECMWF_raw_daily, bs_ECMWF_qmapped_daily, bs_ECMWF_dressed_daily, bs_climo_daily = \
            verify_bss_raw_qmapped_dressed(cleade, cmodel, cftype, rthresh, \
            date_list_anal, apcp_anal_t, validdays, nxa, nya, ndates)
            
        filename = 'bss_MME_'+cleade+'h'+cthresh+'.pydat' 
        print 'writing data to ',filename 
        ouf = open(filename,'w') 
        cPickle.dump(bss_MME_raw,ouf) 
        cPickle.dump(bss_MME_qmapped,ouf) 
        cPickle.dump(bss_MME_dressed,ouf)   
        cPickle.dump(bs_MME_raw_daily,ouf) 
        cPickle.dump(bs_MME_qmapped_daily,ouf) 
        cPickle.dump(bs_MME_dressed_daily,ouf)
        cPickle.dump(bs_climo_daily,ouf)
        ouf.close()
        
        print 'writing data to ',filename 
        filename = 'bss_NCEP_'+cleade+'h'+cthresh+'.pydat'  
        ouf = open(filename,'w') 
        cPickle.dump(bss_NCEP_raw,ouf) 
        cPickle.dump(bss_NCEP_qmapped,ouf) 
        cPickle.dump(bss_NCEP_dressed,ouf)   
        cPickle.dump(bs_NCEP_raw_daily,ouf) 
        cPickle.dump(bs_NCEP_qmapped_daily,ouf) 
        cPickle.dump(bs_NCEP_dressed_daily,ouf)
        cPickle.dump(bs_climo_daily,ouf)
        ouf.close()
        
        print 'writing data to ',filename 
        filename = 'bss_CMC_'+cleade+'h'+cthresh+'.pydat'  
        ouf = open(filename,'w') 
        cPickle.dump(bss_CMC_raw,ouf) 
        cPickle.dump(bss_CMC_qmapped,ouf) 
        cPickle.dump(bss_CMC_dressed,ouf)   
        cPickle.dump(bs_CMC_raw_daily,ouf) 
        cPickle.dump(bs_CMC_qmapped_daily,ouf) 
        cPickle.dump(bs_CMC_dressed_daily,ouf)
        cPickle.dump(bs_climo_daily,ouf)
        ouf.close()
        
        print 'writing data to ',filename 
        filename = 'bss_ECMWF_'+cleade+'h'+cthresh+'.pydat'  
        ouf = open(filename,'w') 
        cPickle.dump(bss_ECMWF_raw,ouf) 
        cPickle.dump(bss_ECMWF_qmapped,ouf) 
        cPickle.dump(bss_ECMWF_dressed,ouf)   
        cPickle.dump(bs_ECMWF_raw_daily,ouf) 
        cPickle.dump(bs_ECMWF_qmapped_daily,ouf) 
        cPickle.dump(bs_ECMWF_dressed_daily,ouf)
        cPickle.dump(bs_climo_daily,ouf)
        ouf.close()
            




