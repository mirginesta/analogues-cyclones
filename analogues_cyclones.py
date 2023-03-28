#!/usr/bin/env python
# coding: utf-8


from datetime import date, timedelta, datetime
import numpy as np


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



