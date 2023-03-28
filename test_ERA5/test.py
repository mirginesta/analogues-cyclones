#!/usr/bin/env python
# coding: utf-8

# In[2]:


from netCDF4 import Dataset as ncdf
from datetime import date, timedelta, datetime
import numpy as np
import xarray as xr
import pandas as pd

# In[3]:


n_analogues = 30

"""

We define the following variables needed to find the n analogues (n_analogues):

time_cyclone : time step of the minimum sea level pressure of the cyclone 

lim_lon0,lim_lon1,lim_lat0,lim_lat1 : limits of which the analogues of the slp field will be found

----------

Important for ERA5 data:
    the limits need to be specified in the following order:
    minimum longitude, maximum longitude, maximum latitude, minimum latitude

    and the ranges are: 

    longitude --> from -180 to 180
    latitude --> from -90 to 90
    
"""
    
time_cyclone = '2020-10-02-06'

lim_lon0,lim_lon1,lim_lat0,lim_lat1 = [-20., 20.,65.,35.]


# In[5]:


# READ DATA:

years=np.arange(1950,2022,1)

geopt_open_list=[]
msl_open_list=[]
for y in years:
    geopt_open_list.append('/home/estimr3/bdd/ERA5/North_Atlantic_22p5N_70N_80W_50E/6hourly/z500/geopotential_500hPa_6hourly_'+str(y)+'.nc')
    msl_open_list.append('/home/estimr3/bdd/ERA5/North_Atlantic_22p5N_70N_80W_50E/6hourly/msl/mean_sea_level_pressure_6hourly_'+str(y)+'.nc')
    
    
mean_sea_level_pressure=xr.open_mfdataset(msl_open_list, combine='by_coords', parallel=True)    
geopotential_500=xr.open_mfdataset(geopt_open_list, parallel=True)

msl=mean_sea_level_pressure.msl
geopt=(geopotential_500.z)/9.8 

time=msl.time.dt.strftime("%Y-%m-%d-%H") 

lon=mean_sea_level_pressure.longitude
lat=mean_sea_level_pressure.latitude


# In[9]:


# Selection of the time of the cyclone 

slp_cycl=msl.sel(time=str(time_cyclone))
zg_cycl=geopt.sel(time=str(time_cyclone))
cycl_date=msl.time.sel(time=str(time_cyclone))

# and spatial (longitud--latitude) region to find the analogues

lonA=lon.sel(longitude=slice(-lim_lon0, lim_lon1))
latA=lat.sel(latitude=slice(lim_lat0,lim_lat1))
slp_cycl_A=slp_cycl.sel(longitude=slice(lim_lon0, lim_lon1)).sel(latitude=slice(lim_lat0,lim_lat1))
zg_cycl_A=zg_cycl.sel(longitude=slice(lim_lon0, lim_lon1)).sel(latitude=slice(lim_lat0,lim_lat1))


# In[11]:


# SELECTING TWO PERIODS PAST AND PRESENT

nyr=35 # number of years for both periods

msl_past = msl[0:365*nyr*4,:,:]
msl_present =  msl[-365*nyr*4:len(time),:,:]

geopt_past = geopt[0:365*nyr*4,:,:]
geopt_present = geopt[-365*nyr*4:len(time),:,:]

time_past=time[0:365*26*4]
time_present=time[-365*nyr*4:len(time)]

msl_past_A = msl[0:365*nyr*4,:,:].sel(longitude=slice(lim_lon0, lim_lon1)).sel(latitude=slice(lim_lat0,lim_lat1))#.groupby('time.year').mean('time')
msl_present_A =  msl[-365*nyr*4:len(time),:,:].sel(longitude=slice(lim_lon0, lim_lon1)).sel(latitude=slice(lim_lat0,lim_lat1))


# EUCLIDEAN DISTANCES PAST

dist_matrix_past=np.sqrt((slp_cycl_A-msl_past_A)**2)
dist_past=((dist_matrix_past.mean(dim='longitude')).mean(dim='latitude'))

# EUCLIDEAN DISTANCES PRESENT

dist_matrix_present=np.sqrt((slp_cycl_A-msl_present_A)**2)
dist_present=((dist_matrix_present.mean(dim='longitude')).mean(dim='latitude'))


time_past=(msl_past.time).dt.strftime("%Y-%m-%d-%H")
time_present=(msl_present.time).dt.strftime("%Y-%m-%d-%H")


# In[34]:


def analogues(Y,N): 
    
    """
    Find the N analogues to a reference event based on their Euclidean distances.

    Parameters
    ----------
    Y : array_like
        The timeseries of the spatially-averaged distances.
    N : int
        The number of closest events to find.

    Returns
    -------
    dist_analogues : ndarray
        The Euclidean distances of the analogues.
    index_analogue : ndarray
        The indices of the position of the analogues.

    """
    
    # Find the values below a quantile p; that is, the (p*100)% lowest Euclidean distances.
    
    p=0.001
    u=np.quantile(Y, p, axis=0)
    Li=np.where(Y<u)[0]
    
    # Aggrupate the values below the quantile that are consecutive and split them if they are
    # separated a maximum of 7 days (that is, 7*4 timesteps, if we use 6-hourly data).
    # These aggrupations are termed as "clusters"
    
    clusters_index=np.split(Li, np.where(np.diff(Li) >= 7*4)[0]+1) # clusters is an array of indices 
    
    # Adjust the quantile p in order to have at least N clusters.
    
    n_clus=len(clusters_index)
    while n_clus < N:
        p=p+0.0005
        u=np.quantile(Y, p, axis=0)
        Li=np.where(Y<u)[0]
        
        clusters_index=np.split(Li, np.where(np.diff(Li) >= 7*4)[0]+1) 
        list_index_stepi=np.zeros((len(clusters_index)))
        
        n_clus=len(clusters_index)
        if n_clus >= N:
            break
        
    # for each cluster, find the time step that has the minimum euclidean distance (i.e. the analogue)
    
    index_analogue=np.zeros((len(clusters_index)))
    for clus in range(n_clus):
        current_cluster=clusters_index[clus] # indices 
        cycl=Y[current_cluster]
        index_analogue_cluster=np.where(cycl==min(cycl)) # position inside the cluster of the minimum distance
        
        index_analogue[clus]=(current_cluster[index_analogue_cluster]).astype(int) # index in Y of the minimum distance timestep in each cluster
       
    dist_analogues=Y[(index_analogue.astype(int))] 
    
    # sort out the distances and select the N lowest (those are the N analogues)
    
    N_analogues=np.argsort(np.array(dist_analogues))[0:N] 
    
    dist_N_analogues=dist_analogues[N_analogues] 
    
    index_N_analogues=index_analogue[N_analogues]
    
    return dist_N_analogues, index_N_analogues.astype(int) 


# In[37]:


# PAST ANALOGUES:

analogues_past=analogues(dist_past, n_analogues)

dist_analogues_past=analogues_past[0] 

index_i_past=analogues_past[1] 

print('Num of ANALOGUES PAST')
print(len(index_i_past))


# In[38]:


# PRESENT ANALOGUES:

analogues_present=analogues(dist_present, n_analogues+1) 

dist_analogues_present=analogues_present[0][1:] # [1:] to exclude the cyclone itself, which is in position 0

index_i_present=analogues_present[1][1:] 

print('Num of ANALOGUES PRESENT')
print(len(index_i_present))


# In[39]:


# ANALOGUES TWO PERIODS

slp_analogues_past=msl_past[index_i_past]
zg_analogues_past=geopt_past[index_i_past]
dates_analogues_past=(time_past[index_i_past]).values

slp_analogues_present=msl_present[index_i_present]
zg_analogues_present=geopt_present[index_i_present]
dates_analogues_present=(time_present[index_i_present]).values

