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

#compute the scale factor of mc histo comparing last or first bin of data histo
def compute_scale_factor(pf_data,pf_mc,bins,rangex,bin="first"):
    print("Computing the scale factor...")
    mc_info=plt.hist(pf_mc,bins = bins,
                          label = 'MC',
                          range = rangex,  color = colors[3])

    data_info=histo(pf_data, bins = bins,
             label = 'Data',
             rangex = rangex, fmt = 'o', c = colors[0])


    if(bin=="first"):
        first_bin=next(i for i,v in enumerate(mc_info[0]) if v >= 3)
        scale_factor=((data_info[0][first_bin]+data_info[0][first_bin+1])/(mc_info[0][first_bin]+mc_info[0][first_bin+1]))
        
    if(bin=="last"):
        last_bin=next(i for i,v in reversed(list(enumerate(mc_info[0]))) if v>=3)
        scale_factor=(data_info[0][last_bin]+data_info[0][last_bin-1])/(mc_info[0][last_bin]+mc_info[0][last_bin-1])
    else:
        print("Mistake")
    mc_scaled=[(scale_factor * item) for item in mc_info[0]]
    return mc_scaled,mc_info[1]    


def plot_and_save(name,file_name_add, pf_t, pf_m, pf_o, pf_data, bins, rangex, last_bin=False, comp_scale_factor=False):
    print("Creating plot for "+name)
  
    if(last_bin==True):
        print("Overflow option chosen.")
        last=rangex[-1]-1
        pf_m_var=[last if item>=last else item for item in pf_m]
        pf_t_var=[last if item>=last else item for item in pf_t]
        pf_o_var=[last if item>=last else item for item in pf_o]
        pf_data_var=[last if item>=last else item for item in pf_data]
        
    else:
        pf_m_var=pf_m[:]
        pf_t_var=pf_t[:]
        pf_o_var=pf_o[:]
        pf_data_var=pf_data[:]

    if(comp_scale_factor==True):
        print("Computing the scale factors...")
        mc_mu,bin_mu=compute_scale_factor(pf_data_var,pf_m_var,bins,rangex,"first")
        mc_tau,bin_tau=compute_scale_factor(pf_data_var,pf_t_var,bins,rangex,"last")
        mc_onia,bin_onia=compute_scale_factor(pf_data_var,pf_o_var,bins,rangex,"last")
        print("Plotting histo...")
        fig2, ax=plt.subplots(figsize=(10, 10))
        stack_histo=plt.hist([bin_mu[:-1],bin_tau[:-1],bin_onia[:-1]], bin_onia, label=['MC mu','MC tau','MC X'],weights=[mc_mu,mc_tau,mc_onia],color= [colors[1],colors[2],colors[3]],histtype='barstacked',stacked =True)

    elif(comp_scale_factor==False):
        scale_factor_mu= 49.04117647058823 * 0.89402  #values from combine
        scale_mu=[scale_factor_mu for item in pf_m_var]
        scale_factor_tau= 29.035714285714285 * 0.27791   #values from combine
        scale_tau=[scale_factor_tau for item in pf_t_var]
        scale_factor_onia=22.933333333333334 * 0.13395    #values from combine
        scale_onia=[scale_factor_onia for item in pf_o_var]
        print("Plotting histo...")
        fig2, ax=plt.subplots(figsize=(10, 10))
        stack_histo=plt.hist([pf_m_var,pf_t_var,pf_o_var], bins=bins, range=rangex, weights=[scale_mu,scale_tau,scale_onia],label=['MC mu','MC tau','MC X'],color= [colors[1],colors[2],colors[3]],histtype='barstacked',stacked =True)
        

    data_info=histo(pf_data_var, bins = bins, 
             label = 'Data', 
             rangex = rangex, fmt = 'o', c = colors[0])
    
    if(last_bin==True):
        xlabels = [n for n in range(-1,last+1)]
        xlabels=np.array(xlabels).astype(str)
        xlabels[-1] += '+'
        ax.set_xticklabels(xlabels)
    
    plt.xlabel(name)
    plt.ylabel('Counts')
    plt.title(name)
    plt.legend(loc = 'best')
    plt.grid()
    
    plt.savefig('plt/'+name+file_name_add+'_plt.png')
    print('Saved plot '+'plt/'+name+file_name_add+'_plt.png')
    
    tau_int=sum(stack_histo[0][1]-stack_histo[0][0])
    mu_int=sum(stack_histo[0][0])
    onia_int= sum(stack_histo[0][2]-stack_histo[0][1])
    data_int=sum(data_info[0])
    print("The MC tau integral is ",tau_int)
    print("The MC mu integral is ",mu_int)
    print("The MC X integral is ",onia_int)
    print("The data integrale is ",data_int)

    print("Creating root histos...")
    out_data = from_numpy(histo(pf_data_var, bins = bins, label = 'data_obs', rangex = rangex))
    out_m = from_numpy(histo(stack_histo[1][:-1], stack_histo[1], weights=stack_histo[0][0],c=colors[1],label = 'MC mu',rangex = rangex))
    out_t = from_numpy(histo(stack_histo[1][:-1], stack_histo[1], weights=stack_histo[0][1]-stack_histo[0][0],c=colors[2],label = 'MC tau',rangex = rangex))
    out_x = from_numpy(histo(stack_histo[1][:-1], stack_histo[1], weights=stack_histo[0][2]-stack_histo[0][1],c=colors[3],label = 'MC x',rangex = rangex)) 
    
    print("Creating root file...")
    file = uproot.recreate('root_files/'+name +file_name_add+'.root', compression=None)
    file["mc_tau"]=out_t
    file["mc_mu"]=out_m
    file["mc_x"]=out_x
    file["data_obs"]=out_data
    file.close()
    print('Saved file '+'root_files/'+name+file_name_add+'.root')

#function that returns the data with new selections
def selection(df):
#    return df[ (df.Blxy_sig > 10)]
    return df[(df.Blxy_sig > 10) & (df.Bsvprob > 0.05) & (df.Bcos2D > 0.95) & \
              (df.Bmass <6.3) & (df.kpt >4) & (df.Bpt >15) & \
              (df.Bmll_raw < 3.1969) & (df.Bmll_raw>2.9969) & \
              ((df.Bb_iso03/df.Bpt)<0.2)\
#              (bcands.k.tightId == 1) & (bcands.mu1.mediumId == 1) & (bcands.mu2.mediumId == 1) & \
              #(df.k.tightId == 1) & (df.mu1.mediumId == 1) & (df.mu2.mediumId == 1) & \ already applied in the .h5 files      
          ]


print("Opening files...")
#Apre i file .h5 salvati
pf_t = pd.read_hdf('hd5_files/OLD3/BcToJpsiTauNu.h5', 'pf')
pf_m = pd.read_hdf('hd5_files/OLD3/BcToJpsiMuNu.h5', 'pf')
pf_o = pd.read_hdf('hd5_files/OLD3/OniaX.h5', 'pf')
pf_data = pd.read_hdf('hd5_files/OLD3/data_files_path.h5', 'pf')

#apply new selections
pf_t  = selection(pf_t)
pf_m = selection(pf_m)
pf_o = selection(pf_o)
pf_data = selection(pf_data)

#opt matplot
plt.rcParams.update({'font.size': 18})
colors = matplotlib.rcParams['axes.prop_cycle'].by_key()['color']

plot_and_save("Bm_miss_sq", "2", pf_t.Bm_miss_sq, pf_m.Bm_miss_sq, pf_o.Bm_miss_sq, pf_data.Bm_miss_sq, 6, (0,7),comp_scale_factor=False,last_bin=True)

#plot_and_save("B_lxy_sig", "", pf_t.Blxy_sig, pf_m.Blxy_sig, pf_data.Blxy_sig, 20, (0,20))
#plot_and_save("B_pt_var", "", pf_t.Bpt_var, pf_m.Bpt_var, pf_data.Bpt_var, 20, (-20,45))
#plot_and_save("B_pt_miss", "", pf_t.Bpt_miss_vec, pf_m.Bpt_miss_vec, pf_data.Bpt_miss_vec, 20, (0,30))
#plot_and_save("B_DR", "", pf_t.BDR, pf_m.BDR, pf_data.BDR, 20, (0,0.9))
#plot_and_save("B_Q_sq", "", pf_t.BQ_sq, pf_m.BQ_sq, pf_data.BQ_sq, 15, (0,11))


#plot_and_save("B_E_mu_star", "", pf_t.BE_mu_star, pf_m.BE_mu_star, pf_data.BE_mu_star, 10, (0,3))
#plot_and_save("B_E_mu_canc", "", pf_t.BE_mu_canc, pf_m.BE_mu_canc, pf_data.BE_mu_canc, 10, (0,6))
#plot_and_save("B_mass", "", pf_t.Bmass, pf_m.Bmass, pf_data.Bmass, 10, (1,8))
#plot_and_save("B_pt", "", pf_t.Bpt, pf_m.Bpt, pf_data.Bpt, 20, (14,50))
#plot_and_save("B_m_jpsi", "", pf_t.Bmll_raw, pf_m.Bmll_raw, pf_data.Bmll_raw, 20, (2.8,3.4))




#RICORDA DI TOGLIERE IL TAGLIO IN ISO
'''
plot_and_save("B_iso03", "", pf_t.Bb_iso03/pf_t.Bpt, pf_m.Bb_iso03/pf_t.Bpt, pf_data.Bb_iso03/pf_t.Bpt, 20, (0,1))
plot_and_save("B_iso04", "", pf_t.Bb_iso04/pf_t.Bpt, pf_m.Bb_iso04/pf_t.Bpt, pf_data.Bb_iso04/pf_t.Bpt, 20, (0,1))

plot_and_save("l1_iso03", "", pf_t.Bl1_iso03/pf_t.mu1pt, pf_m.Bl1_iso03/pf_t.mu1pt, pf_data.Bl1_iso03/pf_t.mu1pt, 20, (0,1))
plot_and_save("l1_iso04", "", pf_t.Bl1_iso04/pf_t.mu1pt, pf_m.Bl1_iso04/pf_t.mu1pt, pf_data.Bl1_iso04/pf_t.mu1pt, 20, (0,1))

plot_and_save("l2_iso03", "", pf_t.Bl2_iso03/pf_t.mu2pt, pf_m.Bl2_iso03/pf_t.mu2pt, pf_data.Bl2_iso03/pf_t.mu2pt, 20, (0,1))
plot_and_save("l2_iso04", "", pf_t.Bl2_iso04/pf_t.mu2pt, pf_m.Bl2_iso04/pf_t.mu2pt, pf_data.Bl2_iso04/pf_t.mu2pt, 20, (0,1))

plot_and_save("k_iso03", "", pf_t.Bk_iso03/pf_t.kpt, pf_m.Bk_iso03/pf_t.kpt, pf_data.Bk_iso03/pf_t.kpt, 20, (0,1))
plot_and_save("k_iso04", "", pf_t.Bk_iso04/pf_t.kpt, pf_m.Bk_iso04/pf_t.kpt, pf_data.Bk_iso04/pf_t.kpt, 20, (0,1))
'''
