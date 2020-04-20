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

    print("Creating original histos...")
    mc_mu_info=plt.hist(pf_m_var,bins = bins, 
             label = 'MC mu', 
             range = rangex,  color = colors[1])
    
    mc_tau_info=plt.hist(pf_t_var,bins = bins, 
             label = 'MC tau', 
             range = rangex,  color = colors[2])

    plt.figure()
    
    data_info=histo(pf_data_var, bins = bins, 
             label = 'Data', 
             rangex = rangex, fmt = 'o', c = colors[0])

    print("Scaling mc histos values...")
    #scale of mc_mu
    first_bin=next(i for i,v in enumerate(mc_mu_info[0]) if v >= 3)
    scale_factor_mu=((data_info[0][first_bin]+data_info[0][first_bin+1])/(mc_mu_info[0][first_bin]+mc_mu_info[0][first_bin+1]))
    mc_mu_scaled=[(scale_factor_mu * item) for item in mc_mu_info[0]]
   
    #scale of mc_tau
    last_bin=next(i for i,v in reversed(list(enumerate(mc_tau_info[0]))) if v>=3)
    scale_factor_tau=(data_info[0][last_bin]+data_info[0][last_bin-1])/(mc_tau_info[0][last_bin]+mc_tau_info[0][last_bin-1])
    mc_tau_scaled=[(scale_factor_tau*item) for item in mc_tau_info[0]]
    
    
    
    print("Creating final plot...")
    fig2=plt.figure(figsize=(10, 10))
    #counts, bins = np.histogram(data)
    plt.hist([mc_mu_info[1][:-1],mc_tau_info[1][:-1]], mc_mu_info[1], label=['MC mu','MC tau'],weights=[mc_mu_scaled,mc_tau_scaled],color= [colors[1],colors[2]],histtype='barstacked',stacked =True)

    histo(pf_data_var, bins = bins, 
             label = 'Data', 
             rangex = rangex, fmt = 'o', c = colors[0])

    plt.xlabel(name)
    plt.ylabel('Counts')
    plt.legend(loc = 'best')
    plt.title(name)
    plt.grid()
    
    plt.savefig('plt/'+name+'_plt.png')
    print('Saved plot '+'plt/'+name+'_plt.png')
    
    plt.figure(figsize=(10, 10))

    print("Creating root histos...")
    #save as TH1 histos (uproot)
    out_data = from_numpy(histo(pf_data_var, bins = bins, label = 'data_obs', rangex = rangex))
    out_m = from_numpy(histo(mc_mu_info[1][:-1], mc_mu_info[1], weights=mc_mu_scaled,c=colors[1],label = 'MC mu',rangex = rangex))
    out_t = from_numpy(histo(mc_tau_info[1][:-1], mc_tau_info[1], weights=mc_tau_scaled, c= colors[2],label = 'MC tau',rangex = rangex))
    print("Creating root file...")
    #open a .root file for the datacard
    file = uproot.recreate('root_files/'+name + '.root', compression=None)
    file["mc_tau"]=out_t
    file["mc_mu"]=out_m
    file["data_obs"]=out_data
    file.close()
    print('Saved file '+'root_files/'+name+'.root')


print("Opening files...")
#Apre i file .h5 salvati
pf_t = pd.read_hdf('hd5_files/BcToJpsiTauNu.h5', 'pf')
pf_m = pd.read_hdf('hd5_files/BcToJpsiMuNu_JpsiToMuMu.h5', 'pf')
pf_data = pd.read_hdf('hd5_files/data_files_path.h5', 'pf')

#opt matplot
plt.rcParams.update({'font.size': 18})
colors = matplotlib.rcParams['axes.prop_cycle'].by_key()['color']

print("Creating plot for branch B_m_miss_sq")
#crea e salva plot e file.root dei vari branch
plot_and_save("B_m_miss_sq", pf_t.Bm_miss_sq, pf_m.Bm_miss_sq, pf_data.Bm_miss_sq,20,(0,10))