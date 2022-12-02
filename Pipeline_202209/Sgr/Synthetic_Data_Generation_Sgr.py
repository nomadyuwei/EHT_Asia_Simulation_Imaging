# Note: this is an example sequence of commands to run in ipython
# The matplotlib windows may not open/close properly if you run this directly as a script

from __future__ import division
from __future__ import print_function

import matplotlib.pyplot as plt
import numpy as np
import ehtim as eh
from   ehtim.calibrating import self_cal as sc
#from  ehtim.plotting import self_cal as sc




from ehtim.const_def import *
from ehtim.statistics.dataframes import *
from ehtim.statistics.stats import *
from ehtim.observing.obs_helpers import *
import itertools as it


from astropy.utils import iers
iers.conf.auto_download = False 

#import ehtim.observing.obs_helpers as obsh
#import ehtim.const_def as ehc

#added by yuwei
import os
os.environ['KMP_DUPLICATE_LIB_OK']='True'



# Load the image and the array
im = eh.image.load_txt('../../models/avery_sgra_eofn.txt')
eht = eh.array.load_txt('../../arrays_new_version_more_sites/EHT_2022_sg.txt')
eht2 = eh.array.load_txt('../../arrays_new_version_more_sites/EHT_2022_sg_reduced.txt')



# Look at the image
im.display()

# Observe the image
# tint_sec is the integration time in seconds, and tadv_sec is the advance time between scans
# tstart_hr is the GMST time of the start of the observation and tstop_hr is the GMST time of the end
# bw_hz is the  bandwidth in Hz
# sgrscat=True blurs the visibilities with the Sgr A* scattering kernel for the appropriate image frequency
# ampcal and phasecal determine if gain variations and phase errors are included
tint_sec = 10
tadv_sec = 600


tstart_hr = 8
tstop_hr = 16
bw_hz = 8e9


obs = im.observe(eht, tint_sec, tadv_sec, tstart_hr, tstop_hr, bw_hz,timetype='GMST',elevmin=10,
                 sgrscat=False, ampcal=True, phasecal=False)

obs.save_txt('simulated_Sgr_with_SG.txt') # exports a text file with the visibilities
obs.save_uvfits('simulated_Sgr_with_SG.uvfits') # exports a UVFITS file modeled on template.UVP
beamparams_with_SG = obs.fit_beam() # fitted beam parameters (fwhm_maj, fwhm_min, theta) in radians


obs_no_SG=obs.flag_sites('SG')
obs_no_SG.save_txt('simulated_Sgr_no_SG.txt') # exports a text file with the visibilities
obs_no_SG.save_uvfits('simulated_Sgr_no_SG.uvfits') # exports a UVFITS file modeled on template.UVP
beamparams_no_SG = obs_no_SG.fit_beam() # fitted beam parameters (fwhm_maj, fwhm_min, theta) in radians


#added by yuwei    
##################################Elevations VS GMST time
times_sid=np.array([])
#for ii in range(0, 24, 0.1):
for ii in range(0, 240, 1):
    iii=ii/10
    times_sid=np.append(times_sid,[iii])

sourcefile = open('../../models/avery_sgra_eofn.txt')
sourcesrc = ' '.join(sourcefile.readline().split()[2:])
sourcera = sourcefile.readline().split()
sourcera2 = float(sourcera[2]) + float(sourcera[4]) / 60.0 + float(sourcera[6]) / 3600.0
thetas = np.mod((times_sid - sourcera2) * HOUR, 2 * np.pi)



times_sid_new=times_sid


color_station=['crimson','green','blue','cyan','fuchsia','gold','orange','limegreen','dodgerblue','purple','red','grey','lightgrey','brown']
color_i=0

fig1, ax = plt.subplots()


for ii in range(len(eht2.tarr)):
    coords = np.array([eht2.tarr[ii]['x'], eht2.tarr[ii]['y'], eht2.tarr[ii]['z']])
    elevation=elev(earthrot(coords, thetas), obs.sourcevec())/DEGREE
    elevation_new=elevation
    
    
    if np.max(elevation_new)>10:
        plt.ylim(10, 90)
        ax.plot(elevation_new,label=eht2.tarr[ii]['site'],color=color_station[color_i])
        color_i=color_i+1
        #ax.legend(loc='center left', bbox_to_anchor=(0.04, -0.37),ncol=3,fontsize=16)
        ax.legend(loc='center left', bbox_to_anchor=(-0.07, -0.55), ncol=3, fontsize=16)
        xticks=list(range(0,len(times_sid_new),40)) # 这里设置的是x轴点的位置（40设置的就是间隔了）
        xlabels=[times_sid_new[x] for x in xticks] #这里设置X轴上的点对应在数据集中的值（这里用的数据为totalSeed）
        ax.set_xticks(xticks)
        ax.set_xticklabels(xlabels)
plt.xlabel('GMST (Hour)', size=16)
plt.ylabel('EL (Deg)', size=16)
plt.title("Sgr A*", size=16)
plt.savefig('./ELE_GMST_Sgr.pdf',format='pdf', dpi=400, bbox_inches = 'tight')
plt.show()
##################################




##added by yuwei    
##################################UV Coverage with SG
timetype=False

field1='u'
field2='v'
conj=True
debias=True
tag_bl=False
ang_unit='deg'
timetype=False
axis=False
rangex=False
rangey=False
snrcut=0.
color=SCOLORS[0]
marker='o'
markersize=MARKERSIZE
label=None
grid=True
ebar=True
axislabels=True
legend=False
show=True
export_pdf=""

# Label individual baselines
# ANDREW TODO this is way too slow, make  it faster??


clist0 = SCOLORS
clist = SCOLORS[1:len(SCOLORS)]

print(clist[0])

# make a color coding dictionary
cdict = {}
ii = 0
baselines = list(it.combinations(obs.tarr['site'],2))
for baseline in baselines:
    cdict[(baseline[0],baseline[1])] = clist[ii%len(clist)]
    cdict[(baseline[1],baseline[0])] = clist[ii%len(clist)]
    ii+=1

# get unique baselines -- TODO easier way? separate function?
alldata = []
allsigx = []
allsigy = []
bllist = []
colors = []
bldata = obs.bllist(conj=True)
for bl in bldata:
    t1 = bl['t1'][0]
    t2 = bl['t2'][0]
    bllist.append((t1,t2))
        

    if t1=='SG' or t2=='SG':
        colors.append(clist0[1])
    else:
        colors.append(clist0[0])
        
    # Unpack data
    dat = obs.unpack_dat(bl, [field1,field2], ang_unit=ang_unit, debias=debias, timetype=timetype)
    alldata.append(dat)
    
    #print(obs.timetype)

    # X error bars
    if sigtype(field1):
        allsigx.append(self.unpack_dat(bl,[sigtype(field1)], ang_unit=ang_unit)[sigtype(field1)])
    else:
        allsigx.append(None)

    # Y error bars
    if sigtype(field2):
        allsigy.append(self.unpack_dat(bl,[sigtype(field2)], ang_unit=ang_unit)[sigtype(field2)])
    else:
        allsigy.append(None)
        
        
        
        
factor=1e9    
        
# make plot(s)

fig=plt.figure()
x = fig.add_subplot(111)

xmins = []
xmaxes = []
ymins = []
ymaxes = []
for i in range(len(alldata)):
    data = alldata[i]
    sigy = allsigy[i]
    sigx = allsigx[i]
    color = colors[i]
    bl = bllist[i]

    # Flag out nans (to avoid problems determining plotting limits)
    mask = ~(np.isnan(data[field1]) + np.isnan(data[field2]))

    data = data[mask]
    if not sigy is None: sigy = sigy[mask]
    if not sigx is None: sigx = sigx[mask]
    if len(data) == 0:
        continue

    xmins.append(np.min(data[field1])/factor)
    xmaxes.append(np.max(data[field1])/factor)
    
    
    ymins.append(np.min(data[field2])/factor)
    ymaxes.append(np.max(data[field2])/factor)
    

    # Plot the data
    tolerance = len(data[field2])

    if label is None:
        labelstr="%s-%s"%((str(bl[0]),str(bl[1])))

    else:
        labelstr=str(label)

    if ebar and (np.any(sigy) or np.any(sigx)):
        x.errorbar(data[field1], data[field2], xerr=sigx, yerr=sigy, label=labelstr,
                   fmt=marker, markersize=markersize, color=color,picker=tolerance,linewidth=1)
    else:
        #x.plot(data[field1], data[field2],
        #       color=color,label=labelstr, linewidth=1.5)
        x.plot(data[field1]/factor, data[field2]/factor, marker, markersize=1.5, color=color)
        #x.plot(data[field1], data[field2], marker, markersize=1.5, color="b")

        
        
min_new=[np.min(xmins),np.min(ymins)]
max_new=[np.max(xmins),np.max(ymins)]

if not rangex:
    rangex = [np.min(min_new) - 0.2 * np.abs(np.min(min_new)),
              np.max(max_new) + 0.2 * np.abs(np.max(max_new))]
    if np.any(np.isnan(np.array(rangex))):
        print("Warning: NaN in data x range: specifying rangex to default")
        rangex = [-100,100]

if not rangey:
    rangey = [np.min(min_new) - 0.2 * np.abs(np.min(min_new)),
              np.max(max_new) + 0.2 * np.abs(np.max(max_new))]
    if np.any(np.isnan(np.array(rangey))):
        print("Warning: NaN in data y range: specifying rangey to default")
        rangey = [-100,100]


x.set_xlim(rangex)
x.set_ylim(rangey)

x.invert_xaxis()


x.set_xlabel(r'u (G$\lambda$)', size=16)
x.set_ylabel(r'v (G$\lambda$)', size=16)


plt.axes().set_aspect('equal')



if legend and tag_bl:
    plt.legend(ncol=2)
elif legend:
    plt.legend()
    
    
    
plt.savefig('./UV_coverage_with_SG.pdf',format='pdf', dpi=400, bbox_inches = 'tight')
    
if export_pdf != "" and not axis:
    fig.savefig(export_pdf, bbox_inches='tight')
if show:
    plt.show(block=False)
         
####################added by yuwei
    
    
##added by yuwei    
##################################UV Coverage without SG
#timetype=False
#
#field1='u'
#field2='v'
#conj=True
#debias=True
#tag_bl=False
#ang_unit='deg'
#timetype=False
#axis=False
#rangex=False
#rangey=False
#snrcut=0.
#color=SCOLORS[0]
#marker='o'
#markersize=MARKERSIZE
#label=None
#grid=True
#ebar=True
#axislabels=True
#legend=False
#show=True
#export_pdf=""

# Label individual baselines
# ANDREW TODO this is way too slow, make  it faster??


clist0 = SCOLORS
clist = SCOLORS[1:len(SCOLORS)]

print(clist[0])

# make a color coding dictionary
cdict = {}
ii = 0
baselines = list(it.combinations(obs_no_SG.tarr['site'],2))
for baseline in baselines:
    cdict[(baseline[0],baseline[1])] = clist[ii%len(clist)]
    cdict[(baseline[1],baseline[0])] = clist[ii%len(clist)]
    ii+=1

# get unique baselines -- TODO easier way? separate function?
alldata = []
allsigx = []
allsigy = []
bllist = []
colors = []
bldata = obs_no_SG.bllist(conj=True)
for bl in bldata:
    t1 = bl['t1'][0]
    t2 = bl['t2'][0]
    bllist.append((t1,t2))
        

    if t1=='SG' or t2=='SG':
        colors.append(clist0[1])
    else:
        colors.append(clist0[0])
        
    # Unpack data
    dat = obs_no_SG.unpack_dat(bl, [field1,field2], ang_unit=ang_unit, debias=debias, timetype=timetype)
    alldata.append(dat)
    

    # X error bars
    if sigtype(field1):
        allsigx.append(self.unpack_dat(bl,[sigtype(field1)], ang_unit=ang_unit)[sigtype(field1)])
    else:
        allsigx.append(None)

    # Y error bars
    if sigtype(field2):
        allsigy.append(self.unpack_dat(bl,[sigtype(field2)], ang_unit=ang_unit)[sigtype(field2)])
    else:
        allsigy.append(None)
        
        
        
        
factor=1e9    
        
# make plot(s)

fig=plt.figure()
x = fig.add_subplot(111)

xmins = []
xmaxes = []
ymins = []
ymaxes = []
for i in range(len(alldata)):
    data = alldata[i]
    sigy = allsigy[i]
    sigx = allsigx[i]
    color = colors[i]
    bl = bllist[i]

    # Flag out nans (to avoid problems determining plotting limits)
    mask = ~(np.isnan(data[field1]) + np.isnan(data[field2]))

    data = data[mask]
    if not sigy is None: sigy = sigy[mask]
    if not sigx is None: sigx = sigx[mask]
    if len(data) == 0:
        continue

    xmins.append(np.min(data[field1])/factor)
    xmaxes.append(np.max(data[field1])/factor)
    
    
    ymins.append(np.min(data[field2])/factor)
    ymaxes.append(np.max(data[field2])/factor)
    

    # Plot the data
    tolerance = len(data[field2])

    if label is None:
        labelstr="%s-%s"%((str(bl[0]),str(bl[1])))

    else:
        labelstr=str(label)

    if ebar and (np.any(sigy) or np.any(sigx)):
        x.errorbar(data[field1], data[field2], xerr=sigx, yerr=sigy, label=labelstr,
                   fmt=marker, markersize=markersize, color=color,picker=tolerance,linewidth=1)
    else:
        #x.plot(data[field1], data[field2],
        #       color=color,label=labelstr, linewidth=1.5)
        x.plot(data[field1]/factor, data[field2]/factor, marker, markersize=1.5, color=color)
        #x.plot(data[field1], data[field2], marker, markersize=1.5, color="b")

        
        
min_new=[np.min(xmins),np.min(ymins)]
max_new=[np.max(xmins),np.max(ymins)]

if not rangex:
    rangex = [np.min(min_new) - 0.2 * np.abs(np.min(min_new)),
              np.max(max_new) + 0.2 * np.abs(np.max(max_new))]
    if np.any(np.isnan(np.array(rangex))):
        print("Warning: NaN in data x range: specifying rangex to default")
        rangex = [-100,100]

if not rangey:
    rangey = [np.min(min_new) - 0.2 * np.abs(np.min(min_new)),
              np.max(max_new) + 0.2 * np.abs(np.max(max_new))]
    if np.any(np.isnan(np.array(rangey))):
        print("Warning: NaN in data y range: specifying rangey to default")
        rangey = [-100,100]


x.set_xlim(rangex)
x.set_ylim(rangey)

x.invert_xaxis()


x.set_xlabel(r'u (G$\lambda$)', size=16)
x.set_ylabel(r'v (G$\lambda$)', size=16)


plt.axes().set_aspect('equal')



if legend and tag_bl:
    plt.legend(ncol=2)
elif legend:
    plt.legend()
    
    
    
plt.savefig('./UV_coverage_no_SG.pdf',format='pdf', dpi=400, bbox_inches = 'tight')
    
if export_pdf != "" and not axis:
    fig.savefig(export_pdf, bbox_inches='tight')
if show:
    plt.show(block=False)
         
####################added by yuwei    

##added by yuwei    
##################################Amp_vs_Baseline 
timetype=False

field1='uvdist'
field2='amp'
conj=False
debias=True
tag_bl=False
ang_unit='deg'
timetype=False
axis=False
rangex=False
rangey=False
snrcut=0.
color=SCOLORS[0]
marker='o'
markersize=MARKERSIZE
label=None
grid=True
ebar=True
axislabels=True
legend=False
show=True
export_pdf=""

# Label individual baselines
# ANDREW TODO this is way too slow, make  it faster??


clist0 = SCOLORS
clist = SCOLORS[1:len(SCOLORS)]

print(clist[0])

# make a color coding dictionary
cdict = {}
ii = 0
baselines = list(it.combinations(obs.tarr['site'],2))
for baseline in baselines:
    cdict[(baseline[0],baseline[1])] = clist[ii%len(clist)]
    cdict[(baseline[1],baseline[0])] = clist[ii%len(clist)]
    ii+=1


##################################with the SG site
alldata = []
allsigx = []
allsigy = []
bllist = []
colors = []
bldata = obs.bllist(conj=True)
for bl in bldata:
    t1 = bl['t1'][0]
    t2 = bl['t2'][0]
    bllist.append((t1,t2))

    if t1=='SG' or t2=='SG':
        colors.append(clist0[1])
    else:
        colors.append(clist0[0])

    # Unpack data
    dat = obs.unpack_dat(bl, [field1,field2], ang_unit=ang_unit, debias=debias, timetype=timetype)
    alldata.append(dat)
    

    # X error bars
    if sigtype(field1):
        allsigx.append(obs.unpack_dat(bl,[sigtype(field1)], ang_unit=ang_unit)[sigtype(field1)])
    else:
        allsigx.append(None)

    # Y error bars
    if sigtype(field2):
        allsigy.append(obs.unpack_dat(bl,[sigtype(field2)], ang_unit=ang_unit)[sigtype(field2)])
    else:
        allsigy.append(None)
        
        
        
factor=1e9  

# make plot(s)

fig=plt.figure()
x = fig.add_subplot(111)

xmins = []
xmaxes = []
ymins = []
ymaxes = []
for i in range(len(alldata)):
    data = alldata[i]
    sigy = allsigy[i]
    sigx = allsigx[i]
    color = colors[i]
    bl = bllist[i]

    # Flag out nans (to avoid problems determining plotting limits)
    mask = ~(np.isnan(data[field1]) + np.isnan(data[field2]))

    data = data[mask]
    if not sigy is None: sigy = sigy[mask]
    if not sigx is None: sigx = sigx[mask]
    if len(data) == 0:
        continue

    xmins.append(np.min(data[field1]))
    xmaxes.append(np.max(data[field1]))
    ymins.append(np.min(data[field2]))
    ymaxes.append(np.max(data[field2]))

    # Plot the data
    tolerance = len(data[field2])

    if label is None:
        labelstr="%s-%s"%((str(bl[0]),str(bl[1])))

    else:
        labelstr=str(label)

    if ebar and (np.any(sigy) or np.any(sigx)):
        x.errorbar(data[field1], data[field2], xerr=sigx, yerr=sigy, label=labelstr,
                   fmt=marker, markersize=markersize, color=color,picker=tolerance,linewidth=1)
    else:
        x.plot(data[field1], data[field2], marker, markersize=1.5, color=color)
        

# Data ranges
if not rangex:
    rangex = [np.min(xmins) - 0.2 * np.abs(np.min(xmins)),
              np.max(xmaxes) + 0.2 * np.abs(np.max(xmaxes))]
    if np.any(np.isnan(np.array(rangex))):
        print("Warning: NaN in data x range: specifying rangex to default")
        rangex = [-100,100]

if not rangey:
    rangey = [np.min(ymins) - 0.2 * np.abs(np.min(ymins)),
              np.max(ymaxes) + 0.2 * np.abs(np.max(ymaxes))]
    if np.any(np.isnan(np.array(rangey))):
        print("Warning: NaN in data y range: specifying rangey to default")
        rangey = [-100,100]
        
# label and save
        
if axislabels:
    try:
        x.set_xlabel(FIELD_LABELS[field1])
        x.set_ylabel(FIELD_LABELS[field2])
    except:
        x.set_xlabel(field1.capitalize())
        x.set_ylabel(field2.capitalize())
     
x.set_xlabel(r'u-v Distance ($\lambda$)', fontsize=16)        
x.set_ylabel(r'Correlated Flux Density (Jy)', fontsize=16)
#
#plt.axes().set_aspect('equal')

x.set_yscale("log")



if legend and tag_bl:
    plt.legend(ncol=2)
elif legend:
    plt.legend()
if grid:
    x.grid()
    
plt.savefig('./Amp_Baseline_Sgr_with_SG.pdf',format='pdf', dpi=400, bbox_inches = 'tight')

if export_pdf != "" and not axis:
    fig.savefig(export_pdf, bbox_inches='tight')
if show:
    plt.show(block=False)
    
    
    
################################## AMP_with_Baseline no SG site
alldata = []
allsigx = []
allsigy = []
bllist = []
colors = []
bldata = obs_no_SG.bllist(conj=True)
for bl in bldata:
    t1 = bl['t1'][0]
    t2 = bl['t2'][0]
    bllist.append((t1,t2))

    if t1=='SG' or t2=='SG':
        colors.append(clist0[1])
    else:
        colors.append(clist0[0])

    # Unpack data
    dat = obs_no_SG.unpack_dat(bl, [field1,field2], ang_unit=ang_unit, debias=debias, timetype=timetype)
    alldata.append(dat)
    

    # X error bars
    if sigtype(field1):
        allsigx.append(obs_no_SG.unpack_dat(bl,[sigtype(field1)], ang_unit=ang_unit)[sigtype(field1)])
    else:
        allsigx.append(None)

    # Y error bars
    if sigtype(field2):
        allsigy.append(obs_no_SG.unpack_dat(bl,[sigtype(field2)], ang_unit=ang_unit)[sigtype(field2)])
    else:
        allsigy.append(None)
        
        
        
factor=1e9  

# make plot(s)

fig=plt.figure()
x = fig.add_subplot(111)

xmins = []
xmaxes = []
ymins = []
ymaxes = []
for i in range(len(alldata)):
    data = alldata[i]
    sigy = allsigy[i]
    sigx = allsigx[i]
    color = colors[i]
    bl = bllist[i]

    # Flag out nans (to avoid problems determining plotting limits)
    mask = ~(np.isnan(data[field1]) + np.isnan(data[field2]))

    data = data[mask]
    if not sigy is None: sigy = sigy[mask]
    if not sigx is None: sigx = sigx[mask]
    if len(data) == 0:
        continue

    xmins.append(np.min(data[field1]))
    xmaxes.append(np.max(data[field1]))
    ymins.append(np.min(data[field2]))
    ymaxes.append(np.max(data[field2]))

    # Plot the data
    tolerance = len(data[field2])

    if label is None:
        labelstr="%s-%s"%((str(bl[0]),str(bl[1])))

    else:
        labelstr=str(label)

    if ebar and (np.any(sigy) or np.any(sigx)):
        x.errorbar(data[field1], data[field2], xerr=sigx, yerr=sigy, label=labelstr,
                   fmt=marker, markersize=markersize, color=color,picker=tolerance,linewidth=1)
    else:
        x.plot(data[field1], data[field2], marker, markersize=1.5, color=color)
        

# Data ranges
if not rangex:
    rangex = [np.min(xmins) - 0.2 * np.abs(np.min(xmins)),
              np.max(xmaxes) + 0.2 * np.abs(np.max(xmaxes))]
    if np.any(np.isnan(np.array(rangex))):
        print("Warning: NaN in data x range: specifying rangex to default")
        rangex = [-100,100]

if not rangey:
    rangey = [np.min(ymins) - 0.2 * np.abs(np.min(ymins)),
              np.max(ymaxes) + 0.2 * np.abs(np.max(ymaxes))]
    if np.any(np.isnan(np.array(rangey))):
        print("Warning: NaN in data y range: specifying rangey to default")
        rangey = [-100,100]
        
# label and save
        
if axislabels:
    try:
        x.set_xlabel(FIELD_LABELS[field1])
        x.set_ylabel(FIELD_LABELS[field2])
    except:
        x.set_xlabel(field1.capitalize())
        x.set_ylabel(field2.capitalize())
     
x.set_xlabel(r'u-v Distance ($\lambda$)', fontsize=16)        
x.set_ylabel(r'Correlated Flux Density (Jy)', fontsize=16)
#
#plt.axes().set_aspect('equal')

x.set_yscale("log")



if legend and tag_bl:
    plt.legend(ncol=2)
elif legend:
    plt.legend()
if grid:
    x.grid()
    
plt.savefig('./Amp_Baseline_Sgr_no_SG.pdf',format='pdf', dpi=400, bbox_inches = 'tight')

if export_pdf != "" and not axis:
    fig.savefig(export_pdf, bbox_inches='tight')
if show:
    plt.show(block=False)
    
###########################cphase


site1='SG'
site2='SPART'
#site3='DELINGHA'
site3=['Yonsei','AMT','SPT']


site1_ii='SG'
site2_ii='SPART'
#site3='DELINGHA'
site3_ii=['Yonsei','AMT','SPT']



vtype='vis'
cphases=[]
force_recompute=False
ang_unit = 'deg'
timetype = 'GMST'
snrcut=0.
axis=False
rangex=False
rangey=False

marker='o'
markersize=ehc.MARKERSIZE
#label=None



if ang_unit == 'deg':
    angle = 1.0
else:
    angle = ehc.DEGREE


# Plot the data
if axis:
    x = axis
else:
    fig = plt.figure()
    x = fig.add_subplot(1, 1, 1)

# Data ranges
if not rangex:
    rangex = [obs.tstart, obs.tstop]
    if np.any(np.isnan(np.array(rangex))):
        print("Warning: NaN in data x range: specifying rangex to default")
        rangex = [8, 12]

if not rangey:
    if ang_unit == 'deg':
        rangey = [-190, 190]
    else:
        rangey = [-1.1 * np.pi, 1.1 * np.pi]

x.set_xlim(rangex)
x.set_ylim(rangey)



for ii in range(len(site3)):
    
    color=ehc.SCOLORS[ii]
    label_ii='%s-%s-%s' % (site1_ii, site2_ii, site3_ii[ii])
    
    if (len(cphases) == 0) and (obs.cphase is not None) and not force_recompute:
        cphases = obs.cphase
    
    cpdata = obs.cphase_tri(site1, site2, site3[ii], vtype=vtype, timetype=timetype,
                             cphases=cphases, force_recompute=force_recompute, snrcut=snrcut)
    plotdata = np.array([[obs['time'], obs['cphase'] * angle, obs['sigmacp']]
                         for obs in cpdata])
    
    nan_mask = np.isnan(plotdata[:, 1])
    plotdata = plotdata[~nan_mask]
    

    
    if ebar and np.any(plotdata[:, 2]):

        x.errorbar(plotdata[:, 0], plotdata[:, 1], yerr=plotdata[:, 2],
                   fmt=marker, markersize=markersize,
                   color=color, label=label_ii)
    else:
        x.plot(plotdata[:, 0], plotdata[:, 1], marker, markersize=markersize, color=color, label=label_ii)  
    

    x.legend(loc='upper right', fontsize=13)
 
    
    

if axislabels:
    x.set_xlabel(obs.timetype + ' (Hour)', fontsize=16)
    if ang_unit == 'deg':
        x.set_ylabel(r'Closure Phase (deg)', fontsize=16)
#        x.set_ylabel(r'Closure Phase $(^\circ)$')
    else:
        x.set_ylabel(r'Closure Phase (radian)', fontsize=16)


if grid:
    x.grid()
if legend:
    plt.legend()
    
plt.savefig('./cphase_Sgr.pdf',format='pdf', dpi=400, bbox_inches = 'tight')

    
if export_pdf != "" and not axis:
    fig.savefig(export_pdf, bbox_inches='tight')
if export_pdf != "" and axis:
    fig = plt.gcf()
    fig.savefig(export_pdf, bbox_inches='tight')
if show:
    plt.show(block=False)