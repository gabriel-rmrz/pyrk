#Crea i plot dai file .h5 

import pandas as pd
import numpy as np
from glob import glob
import matplotlib
import matplotlib.pyplot as plt
import uproot_methods
from uproot_methods.classes.TH1 import from_numpy
import uproot
#funzione che fa istogrammi con numpy

def histo(arr, bins=10, rangex=None, normed=None, weights=None, density=None,
            **kwargs):
    
    content, edges = np.histogram(
        arr, bins = bins, range = rangex, normed = normed, 
        weights = weights, density = density
    )
    if 'fmt' not in kwargs:
        kwargs['fmt'] = 'o'
    centers = (edges[:-1] + edges[1:]) / 2
    errs = np.sqrt(content) # TODO: add frequentist asymmetric errors    
    plt.errorbar(centers, content, yerr = errs, **kwargs)
    return content, edges

#funzione che crea un unico plot con i tre istogrammi (dati, mc1 e mc 2) (va usato stack in realtà)
def plot_and_save(name, pf_t_var, pf_m_var, pf_data_var, bins, rangex):
    
    #create numpy histo
    plt.figure(figsize=(10, 10))

    
    plt.hist([pf_t_var,pf_m_var],bins = bins, 
             label = 'MC mu', 
             range = rangex,  color = [colors[1],colors[2]],histtype='barstacked',stacked =True)

    
    histo(pf_data_var, bins = bins, 
             label = 'Data', 
             rangex = rangex, fmt = 'o', c = colors[0])
    
    plt.xlabel(name)
    plt.ylabel('Counts')
    plt.legend(loc = 'best')
    plt.title(name)
    plt.grid()
    
    plt.savefig(name+'_plt.png')
    print('Saved plot '+name+'_plt.png')
    
    plt.figure(figsize=(10, 10))
    #save as TH1 histos (uproot)
    out_data = from_numpy(histo(pf_data_var, bins = bins, label = 'data_obs', rangex = rangex))
    out_m = from_numpy(histo(pf_m_var, bins = bins, label = 'mc_mu', rangex = rangex))
    out_t = from_numpy(histo(pf_t_var, bins = bins, label = 'mc_tau', rangex = rangex))
    
    #open a .root file for the datacard
    file = uproot.recreate(name + '.root', compression=None)
    file["mc_tau"]=out_t
    file["mc_mu"]=out_m
    file["data_obs"]=out_data
    file.close()
    print('Saved file '+name+'.root')


#Apre i file .h5 salvati
pf_t = pd.read_hdf('BcToJpsiTauNu.h5', 'pf')
pf_m = pd.read_hdf('BcToJpsiMuNu_JpsiToMuMu.h5', 'pf')
pf_data = pd.read_hdf('Run2018D.h5', 'pf')

#opt matplot
plt.rcParams.update({'font.size': 18})
colors = matplotlib.rcParams['axes.prop_cycle'].by_key()['color']

#crea e salva plot e file.root dei vari branch
plot_and_save("B_mass", pf_t.Bmass, pf_m.Bmass, pf_data.Bmass,20,(0,8))