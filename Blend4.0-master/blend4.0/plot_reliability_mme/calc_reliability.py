""" 
    calc_reliability.py 
    this routine will plot reliability diagrams as well as maps of Brier skill scores
    for raw multi-model, quantile-mapped, and weighted, dressed forecasts.  
    Plots are saved to pdf files.  Now with bootstrapped confidence intervals.
"""

from mpl_toolkits.basemap import Basemap
import matplotlib
print 'importing matplotlib'
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import rcParams
import numpy as np
from numpy import ma
import os, sys
from netCDF4 import Dataset
print 'importing netCDF4'
from verify_relia_bss_raw_qmapped_dressed import verify_relia_bss_raw_qmapped_dressed
print verify_relia_bss_raw_qmapped_dressed.__doc__
import pygrib
from dateutils import daterange, dateshift

rcParams['legend.fontsize']='small'
rcParams['legend.fancybox']=True
rcParams['xtick.labelsize']='small'

nclasses = 21 # 0 to 100% probability by 5%
nxa = 464 # CCPA 1/8 degree grid over CONUS
nya = 224 #


cyyyymmddhh_start = '2016040100'
cyyyymmddhh_end = '20160706000'
date_list = daterange(cyyyymmddhh_start, cyyyymmddhh_end, 24)
date_list_anal = daterange(cyyyymmddhh_start, cyyyymmddhh_end, 24)
ndates = len(date_list_anal)

for ilead, cleade in zip([24,48,72,96,120,144,168],['24','48','72','96','120','144','168']):

    ileade = int(cleade)
    ileadb = ileade - 12
    cleadb = str(ileadb)
    
    for ithresh, cthresh, cthresh_title in zip([0,1], ['POP','10.0'],['POP',' > 10mm']):
        
        if cthresh == 'POP':
            rthresh = 0.254
        else:
            rthresh = float(cthresh)
            
        # ---- read in precipitation analyses

        apcp_anal_t = np.zeros((nxa,nya,ndates), dtype=np.float32)
        mninetynine = -99.99*np.ones((nya,nxa), dtype=np.float32)

        for idate, date in zip(range(ndates), date_list_anal):

            print '----- getting precipitation for idate = ',idate,date
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
            
        # ---- call fortran routine to do the reading in of forecasts and the generation
        #      of reliability, frequency of use, and bss information.

        relia = np.zeros(21, dtype=np.float32)
        relia_05 = np.zeros(21, dtype=np.float32)
        relia_95 = np.zeros(21, dtype=np.float32)
        frequse_95 = np.zeros(21, dtype=np.float32)

        validdays = np.ones((ndates),dtype=np.int32)

        # ---- NCEP declarations

        relia_NCEP_raw = np.zeros(21, dtype=np.float32)
        relia_NCEP_qmapped = np.zeros(21, dtype=np.float32)
        relia_NCEP_dressed = np.zeros(21, dtype=np.float32)

        relia_NCEP_raw_05 = np.zeros(21, dtype=np.float32)
        relia_NCEP_qmapped_05 = np.zeros(21, dtype=np.float32)
        relia_NCEP_dressed_05 = np.zeros(21, dtype=np.float32)

        relia_NCEP_raw_95 = np.zeros(21, dtype=np.float32)
        relia_NCEP_qmapped_95 = np.zeros(21, dtype=np.float32)
        relia_NCEP_dressed_95 = np.zeros(21, dtype=np.float32)

        frequse_NCEP_raw_95 = np.zeros(21, dtype=np.float32)
        frequse_NCEP_qmapped_95 = np.zeros(21, dtype=np.float32)
        frequse_NCEP_dressed_95 = np.zeros(21, dtype=np.float32)

        bs_NCEP_raw_daily = np.zeros(ndates,dtype=np.float64)
        bs_NCEP_qmapped_daily = np.zeros(ndates,dtype=np.float64)
        bs_NCEP_dressed_daily = np.zeros(ndates,dtype=np.float64)
        bs_climo_daily = np.zeros(ndates,dtype=np.float64)

        bss = 0.
        bss_NCEP_raw = 0.
        bss_NCEP_qmapped = 0.
        bss_NCEP_dressed = 0.

        # ---- CMC declarations

        relia_CMC_raw = np.zeros(21, dtype=np.float32)
        relia_CMC_qmapped = np.zeros(21, dtype=np.float32)
        relia_CMC_dressed = np.zeros(21, dtype=np.float32)

        relia_CMC_raw_05 = np.zeros(21, dtype=np.float32)
        relia_CMC_qmapped_05 = np.zeros(21, dtype=np.float32)
        relia_CMC_dressed_05 = np.zeros(21, dtype=np.float32)

        relia_CMC_raw_95 = np.zeros(21, dtype=np.float32)
        relia_CMC_qmapped_95 = np.zeros(21, dtype=np.float32)
        relia_CMC_dressed_95 = np.zeros(21, dtype=np.float32)

        frequse_CMC_raw_95 = np.zeros(21, dtype=np.float32)
        frequse_CMC_qmapped_95 = np.zeros(21, dtype=np.float32)
        frequse_CMC_dressed_95 = np.zeros(21, dtype=np.float32)

        bs_CMC_raw_daily = np.zeros(ndates,dtype=np.float64)
        bs_CMC_qmapped_daily = np.zeros(ndates,dtype=np.float64)
        bs_CMC_dressed_daily = np.zeros(ndates,dtype=np.float64)

        bss_CMC_raw = 0.
        bss_CMC_qmapped = 0.
        bss_CMC_dressed = 0.

        # ----- ECMWF declarations

        relia_ECMWF_raw = np.zeros(21, dtype=np.float32)
        relia_ECMWF_qmapped = np.zeros(21, dtype=np.float32)
        relia_ECMWF_dressed = np.zeros(21, dtype=np.float32)

        relia_ECMWF_raw_05 = np.zeros(21, dtype=np.float32)
        relia_ECMWF_qmapped_05 = np.zeros(21, dtype=np.float32)
        relia_ECMWF_dressed_05 = np.zeros(21, dtype=np.float32)

        relia_ECMWF_raw_95 = np.zeros(21, dtype=np.float32)
        relia_ECMWF_qmapped_95 = np.zeros(21, dtype=np.float32)
        relia_ECMWF_dressed_95 = np.zeros(21, dtype=np.float32)

        frequse_ECMWF_raw_95 = np.zeros(21, dtype=np.float32)
        frequse_ECMWF_qmapped_95 = np.zeros(21, dtype=np.float32)
        frequse_ECMWF_dressed_95 = np.zeros(21, dtype=np.float32)

        bs_ECMWF_raw_daily = np.zeros(ndates,dtype=np.float64)
        bs_ECMWF_qmapped_daily = np.zeros(ndates,dtype=np.float64)
        bs_ECMWF_dressed_daily = np.zeros(ndates,dtype=np.float64)

        bss_ECMWF_raw = 0.
        bss_ECMWF_qmapped = 0.
        bss_ECMWF_dressed = 0.

        # ----- MME declarations

        relia_MME_raw = np.zeros(21, dtype=np.float32)
        relia_MME_qmapped = np.zeros(21, dtype=np.float32)
        relia_MME_dressed = np.zeros(21, dtype=np.float32)

        relia_MME_raw_05 = np.zeros(21, dtype=np.float32)
        relia_MME_qmapped_05 = np.zeros(21, dtype=np.float32)
        relia_MME_dressed_05 = np.zeros(21, dtype=np.float32)

        relia_MME_raw_95 = np.zeros(21, dtype=np.float32)
        relia_MME_qmapped_95 = np.zeros(21, dtype=np.float32)
        relia_MME_dressed_95 = np.zeros(21, dtype=np.float32)

        frequse_MME_raw_95 = np.zeros(21, dtype=np.float32)
        frequse_MME_qmapped_95 = np.zeros(21, dtype=np.float32)
        frequse_MME_dressed_95 = np.zeros(21, dtype=np.float32)

        bs_MME_raw_daily = np.zeros(ndates,dtype=np.float64)
        bs_MME_qmapped_daily = np.zeros(ndates,dtype=np.float64)
        bs_MME_dressed_daily = np.zeros(ndates,dtype=np.float64)

        bss_MME_raw = 0.
        bss_MME_qmapped = 0.
        bss_MME_dressed = 0.

        nxa = 464 # CCPA 1/8 degree grid over CONUS
        nya = 224 #

        cftype = ''
        cmodel = 'MME'
        relia_MME_raw, relia_MME_raw_05, relia_MME_raw_95, frequse_MME_raw, \
            bss_MME_raw, bs_MME_raw_daily,  relia_MME_qmapped, relia_MME_qmapped_05, \
            relia_MME_qmapped_95, frequse_MME_qmapped, bss_MME_qmapped, \
            bs_MME_qmapped_daily, relia_MME_dressed, relia_MME_dressed_05, \
            relia_MME_dressed_95, frequse_MME_dressed, bss_MME_dressed, \
            bs_MME_dressed_daily, bs_climo_daily = \
            verify_relia_bss_raw_qmapped_dressed(cleade, cmodel, cftype, nclasses, \
            rthresh, date_list_anal, apcp_anal_t, validdays, nxa, nya, ndates)    

        cmodel = 'NCEP '
        relia_NCEP_raw, relia_NCEP_raw_05, relia_NCEP_raw_95, frequse_NCEP_raw, \
            bss_NCEP_raw, bs_NCEP_raw_daily, relia_NCEP_qmapped, relia_NCEP_qmapped_05, \
            relia_NCEP_qmapped_95, frequse_NCEP_qmapped, bss_NCEP_qmapped, bs_NCEP_qmapped_daily, \
            relia_NCEP_dressed, relia_NCEP_dressed_05, relia_NCEP_dressed_95, \
            frequse_NCEP_dressed, bss_NCEP_dressed, bs_NCEP_dressed_daily, \
            bs_climo_daily = verify_relia_bss_raw_qmapped_dressed(cleade, \
            cmodel, cftype, nclasses, rthresh, date_list_anal, apcp_anal_t, \
            validdays, nxa, nya, ndates)

        cmodel = 'CMC  '
        relia_CMC_raw, relia_CMC_raw_05, relia_CMC_raw_95, frequse_CMC_raw, \
            bss_CMC_raw, bs_CMC_raw_daily, relia_CMC_qmapped, relia_CMC_qmapped_05, \
            relia_CMC_qmapped_95, frequse_CMC_qmapped, bss_CMC_qmapped, bs_CMC_qmapped_daily, \
            relia_CMC_dressed, relia_CMC_dressed_05, relia_CMC_dressed_95, \
            frequse_CMC_dressed, bss_CMC_dressed, bs_CMC_dressed_daily, bs_climo_daily = \
            verify_relia_bss_raw_qmapped_dressed(cleade, cmodel, cftype, nclasses, \
            rthresh, date_list_anal, apcp_anal_t, validdays, nxa, nya, ndates)
    
        cmodel = 'ECMWF'
        relia_ECMWF_raw, relia_ECMWF_raw_05, relia_ECMWF_raw_95, frequse_ECMWF_raw, \
            bss_ECMWF_raw, bs_ECMWF_raw_daily, relia_ECMWF_qmapped, \
            relia_ECMWF_qmapped_05, relia_ECMWF_qmapped_95, frequse_ECMWF_qmapped, \
            bss_ECMWF_qmapped, bs_ECMWF_qmapped_daily, relia_ECMWF_dressed, \
            relia_ECMWF_dressed_05, relia_ECMWF_dressed_95, frequse_ECMWF_dressed, \
            bss_ECMWF_dressed, bs_ECMWF_dressed_daily, bs_climo_daily = \
            verify_relia_bss_raw_qmapped_dressed(cleade, cmodel, cftype, nclasses, \
            rthresh, date_list_anal, apcp_anal_t, validdays, nxa, nya, ndates) 


        # ---- now make multi-panel reliability diagrams for each of the forecasts

        fig = plt.figure(figsize=(9.,3.8))
        plt.suptitle(cthresh +' reliability diagrams for +'+cleadb+' to +'+cleade+' hour forecasts',fontsize=15)

        rcParams['xtick.labelsize']='x-small'
        rcParams['ytick.labelsize']='x-small'
        probs = np.arange(nclasses) * 100./np.real(nclasses-1)

        for itype in range(3):
            if itype == 0:
                ov_ax = [0.06,0.14,0.25,0.85]
                ctitle = '(a) Raw ensemble'
                frequse_NCEP = frequse_NCEP_raw
                frequse_CMC = frequse_CMC_raw
                frequse_ECMWF = frequse_ECMWF_raw
                frequse_MME = frequse_MME_raw
                relia_NCEP = relia_NCEP_raw
                relia_CMC = relia_CMC_raw
                relia_ECMWF = relia_ECMWF_raw
                relia_MME = relia_MME_raw
                bss_NCEP = bss_NCEP_raw
                bss_CMC = bss_CMC_raw
                bss_ECMWF = bss_ECMWF_raw
                bss_MME = bss_MME_raw
            elif itype == 1:
                ov_ax = [0.39,0.14,0.25,0.85]
                ctitle = '(b) Quantile-mapped, dressed'
                frequse_NCEP = frequse_NCEP_qmapped
                frequse_CMC = frequse_CMC_qmapped
                frequse_ECMWF = frequse_ECMWF_qmapped
                frequse_MME = frequse_MME_qmapped
                relia_NCEP = relia_NCEP_qmapped
                relia_CMC = relia_CMC_qmapped
                relia_ECMWF = relia_ECMWF_qmapped
                relia_MME = relia_MME_qmapped
                bss_NCEP = bss_NCEP_qmapped
                bss_CMC = bss_CMC_qmapped
                bss_ECMWF = bss_ECMWF_qmapped
                bss_MME = bss_MME_qmapped
            elif itype == 2:
                ov_ax = [0.72,0.14,0.25,0.85]  
                ctitle = '(c) ... and weighted'
                frequse_NCEP = frequse_NCEP_dressed
                frequse_CMC = frequse_CMC_dressed
                frequse_ECMWF = frequse_ECMWF_dressed
                frequse_MME = frequse_MME_dressed  
                relia_NCEP = relia_NCEP_dressed
                relia_CMC = relia_CMC_dressed
                relia_ECMWF = relia_ECMWF_dressed
                relia_MME = relia_MME_dressed 
                bss_NCEP = bss_NCEP_dressed 
                bss_CMC = bss_CMC_dressed 
                bss_ECMWF = bss_ECMWF_dressed 
                bss_MME = bss_MME_dressed         
       
            # --- Frequency of usage inset diagram

            a2 = fig.add_axes(ov_ax,zorder=3)
            a2.yaxis.tick_right()
            a2.yaxis.set_label_position("right")
            a2.set_frame_on(False)
        
            a2.semilogy(probs, frequse_NCEP,color='Red',linewidth=1, alpha=1,linestyle='--')
            a2.semilogy(probs, frequse_CMC,color='LimeGreen',linewidth=1, alpha=1,linestyle='--')
            a2.semilogy(probs, frequse_ECMWF,color='RoyalBlue',linewidth=1, alpha=1,linestyle='--')
            a2.semilogy(probs, frequse_MME,color='DarkGray',linewidth=1, alpha=1,linestyle='--')    
                    
            a2.set_ylabel('Frequency of Usage',fontsize=8,verticalalignment='center')
            a2.set_yticks([0.0001,0.001,0.01,0.1,1.0])
            a2.set_yticklabels(['0.0001','0.001','0.01','0.1','1.0'])
            a2.set_ylim(0.001,1.)
            a2.set_xlim(-1,101)
            a2.set_xticks([0,100])
            a2.set_xticklabels([' ',' '],zorder=0)
            a2.set_xlabel('',zorder=0)

            # -- make reliability diagram 

            a1 = fig.add_axes(ov_ax,zorder=2)
            a1.set_title(ctitle,fontsize=9)

            strbss_NCEP = 'NCEP BSS = %0.2f' %(bss_NCEP)
            strbss_CMC = 'CMC BSS = %0.2f' %(bss_CMC)
            strbss_ECMWF = 'ECMWF BSS = %0.2f' %(bss_ECMWF)
            strbss_MME = 'MME BSS = %0.2f' %(bss_MME)
            relia_NCEP_m = ma.array(relia_NCEP)
            relia_NCEP_m = ma.masked_where(relia_NCEP_m < 0.0, relia_NCEP_m)
            relia_CMC_m = ma.array(relia_CMC)
            relia_CMC_m = ma.masked_where(relia_CMC_m < 0.0, relia_CMC_m)
            relia_ECMWF_m = ma.array(relia_ECMWF)
            relia_ECMWF_m = ma.masked_where(relia_ECMWF_m < 0.0, relia_ECMWF_m)
            relia_MME_m = ma.array(relia_MME)
            relia_MME_m = ma.masked_where(relia_MME_m < 0.0, relia_MME_m)
    
            a1.plot([0,100],[0,100],'-',linewidth=0.5,color='LightGray')
            a1.plot(probs,100.*relia_NCEP_m,'o',color='r',markersize=1.2,zorder=2)
            a1.plot(probs,100.*relia_NCEP_m,'-',color='r',linewidth=1.2,zorder=2)
            a1.plot(probs,100.*relia_CMC_m,'o',color='LimeGreen',markersize=1.2,zorder=2)
            a1.plot(probs,100.*relia_CMC_m,'-',color='LimeGreen',linewidth=1.2,zorder=2)
            a1.plot(probs,100.*relia_ECMWF_m,'o',color='RoyalBlue',markersize=1.2,zorder=2)
            a1.plot(probs,100.*relia_ECMWF_m,'-',color='RoyalBlue',linewidth=1.2,zorder=2)
            a1.plot(probs,100.*relia_MME_m,'o',color='DarkGray',markersize=1.2,zorder=2)
            a1.plot(probs,100.*relia_MME_m,'-',color='DarkGray',linewidth=1.2,zorder=2)
    
            if itype == 0:  
                a1.set_ylabel('Observed relative\nfrequency (%)',fontsize=8)
            else:
                a1.set_ylabel('',fontsize=8)
            a1.set_xlabel('Forecast probability (%)',fontsize=8)
            a1.set_ylim(-1,101)
            a1.set_xlim(-1,101)
            a1.set_xticks([0,20,40,60,80,100])
            a1.set_yticks([0,20,40,60,80,100])

            # -- BSS inserted here

            a1.text(8,93,strbss_NCEP,fontsize=6,color='Red')
            a1.text(8,87,strbss_CMC,fontsize=6,color='LimeGreen')        
            a1.text(8,81,strbss_ECMWF,fontsize=6,color='RoyalBlue')
            a1.text(8,75,strbss_MME,fontsize=6,color='DarkGray') 
                   
        plot_title = 'relia_'+cthresh+'_hour'+cleade+'.pdf'
        print 'saving plot to file = ',plot_title
        plt.savefig(plot_title)
        print 'Plot done'




