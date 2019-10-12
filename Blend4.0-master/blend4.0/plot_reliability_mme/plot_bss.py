""" 
    plot_reliability_raw_qmapped.py 
    this routine will plot reliability diagrams as well as maps of Brier skill scores
    for raw multi-model, quantile-mapped forecasts.  
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
import cPickle
rcParams['legend.fontsize']='small'
rcParams['legend.fancybox']=True
rcParams['xtick.labelsize']='small'

def paired_bootstrap(x,y,z):
    """ given forecast brier score vectors x and y and climatology z, \]
    return a paired bootstrap estimate of the 5th and 95th percentiles
    of the resampled bss distribution"""
    import numpy as np

    print 'x = ',x
    print 'y = ',y
    print 'z = ',z
    xcopy = np.copy(x)
    ycopy = np.copy(y)
    zcopy = np.copy(z)
    nshape = np.shape(xcopy)
    nelts = nshape[0]
    for ielt in range(nelts):
        if xcopy[ielt] < 0.0 or ycopy[ielt] < 0.0 or zcopy[ielt] < 0.0:
            xcopy[ielt] = 0.0
            ycopy[ielt] = 0.0
            zcopy[ielt] = 0.0

    nresa = 1000   # 1000 resamplings
    sum1 = np.zeros((nresa),dtype=np.float64)
    sum2 = np.zeros((nresa),dtype=np.float64)
    bss1 = np.zeros((nresa),dtype=np.float32)
    bss2 = np.zeros((nresa),dtype=np.float32)
    ones = np.ones((nelts),dtype=np.float)
    zeros = np.zeros((nelts),dtype=np.float)
    for i in range(nresa):
        x0 = np.random.rand(nelts)
        iusex = np.where(x0 < 0.5, ones, zeros)
        iusey = 1 - iusex
        sum1[i] = np.mean(xcopy*iusey + ycopy*iusex)
        sum2[i] = np.mean(xcopy*iusex + ycopy*iusey)
        bss1[i] = 1.0 - sum1[i]/np.sum(zcopy)
        bss2[i] = 1.0 - sum2[i]/np.sum(zcopy)

    bssdiff = bss1 - bss2
    dsort = np.sort(bssdiff)
    d05 = dsort[np.int(.05*np.float(nresa))]
    d95 = dsort[np.int(.95*np.float(nresa))]
    print 'd05, d95 = ',d05,d95
    return d05, d95

# ---- get forecast type

cftype = ''

# ---- call fortran routine to do the reading in of forecasts and the generation
#      of reliability, frequency of use, and bss information.

filename = 'bss_MME_24hPOP.pydat'
inf = open(filename,'r')
bss_MME_raw_in = cPickle.load(inf)
bss_MME_qmapped_in = cPickle.load(inf)
bss_MME_dressed_in = cPickle.load(inf)
bs_MME_raw_daily = cPickle.load(inf)
bs_MME_qmapped_daily = cPickle.load(inf)
inf.close()
ndates = len(bs_MME_qmapped_daily)


# NCEP declarations

bs_NCEP_raw_daily = np.zeros(ndates,dtype=np.float64)
bs_NCEP_qmapped_daily = np.zeros(ndates,dtype=np.float64)
bs_NCEP_dressed_daily = np.zeros(ndates,dtype=np.float64)
bs_climo_daily = np.zeros(ndates,dtype=np.float64)

bss = 0.
bss_NCEP_raw = 0.
bss_NCEP_qmapped = 0.
bss_NCEP_dressed = 0.

bss_NCEP_raw_flead = np.zeros((7,2), dtype=np.float32)
bss_NCEP_qmapped_flead = np.zeros((7,2), dtype=np.float32)
bss_NCEP_dressed_flead = np.zeros((7,2), dtype=np.float32)

# CMC declarations

bs_CMC_raw_daily = np.zeros(ndates,dtype=np.float64)
bs_CMC_qmapped_daily = np.zeros(ndates,dtype=np.float64)
bs_CMC_dressed_daily = np.zeros(ndates,dtype=np.float64)

bss_CMC_raw = 0.
bss_CMC_qmapped = 0.
bss_CMC_dressed = 0.

bss_CMC_raw_flead = np.zeros((7,2), dtype=np.float32)
bss_CMC_qmapped_flead = np.zeros((7,2), dtype=np.float32)
bss_CMC_dressed_flead = np.zeros((7,2), dtype=np.float32)

# ECMWF declarations

bs_ECMWF_raw_daily = np.zeros(ndates,dtype=np.float64)
bs_ECMWF_qmapped_daily = np.zeros(ndates,dtype=np.float64)
bs_ECMWF_dressed_daily = np.zeros(ndates,dtype=np.float64)

bss_ECMWF_raw = 0.
bss_ECMWF_qmapped = 0.
bss_ECMWF_dressed = 0.

bss_ECMWF_raw_flead = np.zeros((7,2), dtype=np.float32)
bss_ECMWF_qmapped_flead = np.zeros((7,2), dtype=np.float32)
bss_ECMWF_dressed_flead = np.zeros((7,2), dtype=np.float32)

# MME declarations

bs_MME_raw_daily = np.zeros(ndates,dtype=np.float64)
bs_MME_qmapped_daily = np.zeros(ndates,dtype=np.float64)
bs_MME_dressed_daily = np.zeros(ndates,dtype=np.float64)

bss_MME_raw = 0.
bss_MME_qmapped = 0.
bss_MME_dressed = 0.

bss_MME_raw_flead = np.zeros((7,2), dtype=np.float32)
bss_MME_qmapped_flead = np.zeros((7,2), dtype=np.float32)
bss_MME_dressed_flead = np.zeros((7,2), dtype=np.float32)

bss_05 = np.zeros((7,2), dtype=np.float32)
bss_95 = np.zeros((7,2), dtype=np.float32)

bss_MME_05 = np.zeros((7,2), dtype=np.float32)
bss_MME_95 = np.zeros((7,2), dtype=np.float32)

bss_ECMWF_05 = np.zeros((7,2), dtype=np.float32)
bss_ECMWF_95 = np.zeros((7,2), dtype=np.float32)

errorbars = np.zeros((7,2), dtype=np.float32)


for ilead, cleade in zip(range(7),['24','48','72','96','120','144','168']):

    ileade = int(cleade)
    ileadb = ileade - 12
    cleadb = str(ileadb)
    
    for ithresh, cthresh in zip(range(2), ['POP', '10.0']):

        filename = 'bss_MME_'+cleade+'h'+cthresh+'.pydat'
        inf = open(filename,'r')
        bss_MME_raw_in = cPickle.load(inf)
        bss_MME_qmapped_in = cPickle.load(inf)
        bss_MME_dressed_in = cPickle.load(inf)
        bs_MME_raw_daily = cPickle.load(inf)
        bs_MME_qmapped_daily = cPickle.load(inf)
        bs_MME_dressed_daily = cPickle.load(inf)
        bs_climo_daily = cPickle.load(inf)
        inf.close()
        
        filename = 'bss_NCEP_'+cleade+'h'+cthresh+'.pydat'
        inf = open(filename,'r')
        bss_NCEP_raw_in = cPickle.load(inf)
        bss_NCEP_qmapped_in = cPickle.load(inf)
        bss_NCEP_dressed_in = cPickle.load(inf)
        bs_NCEP_raw_daily = cPickle.load(inf)
        bs_NCEP_qmapped_daily = cPickle.load(inf)
        bs_NCEP_dressed_daily = cPickle.load(inf)
        inf.close()
                      
        filename = 'bss_CMC_'+cleade+'h'+cthresh+'.pydat'
        inf = open(filename,'r')
        bss_CMC_raw_in = cPickle.load(inf)
        bss_CMC_qmapped_in = cPickle.load(inf)
        bss_CMC_dressed_in = cPickle.load(inf)
        bs_CMC_raw_daily = cPickle.load(inf)
        bs_CMC_qmapped_daily = cPickle.load(inf)
        bs_CMC_dressed_daily = cPickle.load(inf)
        inf.close()
               
        filename = 'bss_ECMWF_'+cleade+'h'+cthresh+'.pydat'
        inf = open(filename,'r')
        bss_ECMWF_raw_in = cPickle.load(inf)
        bss_ECMWF_qmapped_in = cPickle.load(inf)
        bss_ECMWF_dressed_in = cPickle.load(inf)
        bs_ECMWF_raw_daily = cPickle.load(inf)
        bs_ECMWF_qmapped_daily = cPickle.load(inf)
        bs_ECMWF_dressed_daily = cPickle.load(inf)
        inf.close()
        
        # --- get block bootstrap confidence intervals around MME forecast 
       
        bss_05[ilead,ithresh], bss_95[ilead,ithresh] = paired_bootstrap \
            (bs_ECMWF_dressed_daily, bs_MME_dressed_daily, bs_climo_daily)
        errorbars[ilead,ithresh] = bss_95[ilead,ithresh] - bss_05[ilead,ithresh]

        # ---- now set the final output array values
        
        bss_NCEP_raw_flead[ilead,ithresh] = 1. - np.sum(bs_NCEP_raw_daily) / np.sum(bs_climo_daily)
        bss_NCEP_qmapped_flead[ilead,ithresh] = 1. - np.sum(bs_NCEP_qmapped_daily) / np.sum(bs_climo_daily)
        bss_NCEP_dressed_flead[ilead,ithresh] = 1. - np.sum(bs_NCEP_dressed_daily) / np.sum(bs_climo_daily)

        bss_CMC_raw_flead[ilead,ithresh] = 1. - np.sum(bs_CMC_raw_daily) / np.sum(bs_climo_daily)
        bss_CMC_qmapped_flead[ilead,ithresh] = 1. - np.sum(bs_CMC_qmapped_daily) / np.sum(bs_climo_daily)
        bss_CMC_dressed_flead[ilead,ithresh] = 1. - np.sum(bs_CMC_dressed_daily) / np.sum(bs_climo_daily)

        bss_ECMWF_raw_flead[ilead,ithresh] = 1. - np.sum(bs_ECMWF_raw_daily) / np.sum(bs_climo_daily)
        bss_ECMWF_qmapped_flead[ilead,ithresh] = 1. - np.sum(bs_ECMWF_qmapped_daily) / np.sum(bs_climo_daily)
        bss_ECMWF_dressed_flead[ilead,ithresh] = 1. - np.sum(bs_ECMWF_dressed_daily) / np.sum(bs_climo_daily)

        bss_MME_raw_flead[ilead,ithresh] = 1. - np.sum(bs_MME_raw_daily) / np.sum(bs_climo_daily)
        bss_MME_qmapped_flead[ilead,ithresh] = 1. - np.sum(bs_MME_qmapped_daily) / np.sum(bs_climo_daily)
        bss_MME_dressed_flead[ilead,ithresh] = 1. - np.sum(bs_MME_dressed_daily) / np.sum(bs_climo_daily)
        
# ---- now make multi-panel reliability diagrams for each of the forecasts

fig = plt.figure(figsize=(9,3.5))
plt.suptitle('Brier Skill Scores',fontsize=16)

rcParams['xtick.labelsize']='medium'
rcParams['ytick.labelsize']='medium'

# ---- panel a

axes = [0.08,0.13,0.4,0.7]
a1 = fig.add_axes(axes,)
a1.set_title('(a) 12-h POP',fontsize=14)    
a1.plot(range(1,8),bss_NCEP_dressed_flead[:,0],linestyle='-',color='Red',label='NCEP weighted, dressed')
a1.plot(range(1,8),bss_CMC_dressed_flead[:,0],linestyle='-',color='LimeGreen',label='CMC weighted, dressed')
a1.plot(range(1,8),bss_ECMWF_dressed_flead[:,0],\
    linestyle='-',color='RoyalBlue',label='ECMWF weighted, dressed')
a1.plot(range(1,8),bss_MME_dressed_flead[:,0],\
    linestyle='-',color='Gray',label='MME weighted, dressed')
    
a1.plot(range(1,8),bss_NCEP_dressed_flead[:,0],linestyle=':',color='Red',label='NCEP quantile-mapped')
a1.plot(range(1,8),bss_CMC_dressed_flead[:,0],linestyle=':',color='LimeGreen',label='CMC quantile-mapped')
a1.plot(range(1,8),bss_ECMWF_dressed_flead[:,0],\
    linestyle=':',color='RoyalBlue',label='ECMWF quantile-mapped')
a1.plot(range(1,8),bss_MME_dressed_flead[:,0],\
    linestyle=':',color='Gray',label='MME quantile-mapped')    
    
a1.plot(range(1,8),bss_NCEP_raw_flead[:,0],linestyle='--',color='Red',label='NCEP raw')
a1.plot(range(1,8),bss_CMC_raw_flead[:,0],linestyle='--',color='LimeGreen',label='CMC raw')
a1.plot(range(1,8),bss_ECMWF_raw_flead[:,0],\
    linestyle='--',color='RoyalBlue',label='ECMWF raw')
a1.plot(range(1,8),bss_MME_raw_flead[:,0],\
    linestyle='--',color='Gray',label='MME raw')

a1.set_xlabel('Lead time (days)',fontsize=11)
a1.set_ylabel('Brier Skill Score',fontsize=11)
a1.set_ylim(-0.05,0.65)
a1.set_xlim(1.5,7.5)
a1.set_xticks([1,2,3,4,5,6,7])
a1.set_yticks([-0.1,0,0.1,0.2,0.3,0.4,0.5,0.6])
a1.grid(color='LightGray',linestyle='-.',linewidth=0.25)
a1.plot([0.5,7.5],[0,0], linewidth=1,color='Black')

# ---- panel b

axes = [0.58,0.13,0.4,0.7]
a1 = fig.add_axes(axes,)
a1.set_title('(b) > 10 mm / 12 h',fontsize=14)    
a1.plot(range(1,8),bss_NCEP_dressed_flead[:,1],linestyle='-',color='Red',label='NCEP dressed')
a1.plot(range(1,8),bss_CMC_dressed_flead[:,1],linestyle='-',color='LimeGreen',label='CMC dressed')
a1.plot(range(1,8),bss_ECMWF_dressed_flead[:,1],\
    linestyle='-',color='RoyalBlue',label='ECMWF dressed')
a1.plot(range(1,8),bss_MME_dressed_flead[:,1],\
    linestyle='-',color='Gray',label='MME dressed')
    
a1.plot(range(1,8),bss_NCEP_dressed_flead[:,1],linestyle=':',color='Red',label='NCEP quantile-mapped')
a1.plot(range(1,8),bss_CMC_dressed_flead[:,1],linestyle=':',color='LimeGreen',label='CMC quantile-mapped')
a1.plot(range(1,8),bss_ECMWF_dressed_flead[:,1],\
    linestyle=':',color='RoyalBlue',label='ECMWF quantile-mapped')
a1.plot(range(1,8),bss_MME_dressed_flead[:,1],\
    linestyle=':',color='Gray',label='MME quantile-mapped')    

a1.plot(range(1,8),bss_NCEP_raw_flead[:,1],linestyle='--',color='Red',label='NCEP raw')
a1.plot(range(1,8),bss_CMC_raw_flead[:,1],linestyle='--',color='LimeGreen',label='CMC raw')
a1.plot(range(1,8),bss_ECMWF_raw_flead[:,1],\
    linestyle='--',color='RoyalBlue',label='ECMWF raw')
a1.plot(range(1,8),bss_MME_raw_flead[:,1],\
    linestyle='--',color='Gray',label='MME raw')

a1.errorbar(range(1,8),bss_ECMWF_dressed_flead[:,1],\
    yerr=errorbars[:,1],linestyle='-',color='RoyalBlue')
a1.errorbar(range(1,8),bss_MME_dressed_flead[:,1],\
    yerr=errorbars[:,1],linestyle='-',color='Gray')
a1.set_xlabel('Lead time (days)',fontsize=11)
a1.set_ylabel('Brier Skill Score',fontsize=11)
a1.set_ylim(-0.1,0.4)
a1.set_xlim(0.5,7.5)
a1.set_xticks([1,2,3,4,5,6,7])
a1.set_yticks([-0.1,0,0.1,0.2,0.3,0.4])
a1.grid(color='LightGray',linestyle='-.',linewidth=0.25)
a1.plot([0.5,7.5],[0,0], linewidth=1,color='Black')
a1.legend(loc=0)
    
plot_title = 'BSS.pdf'
print 'saving plot to file = ',plot_title
plt.savefig(plot_title)
print 'Plot done'

