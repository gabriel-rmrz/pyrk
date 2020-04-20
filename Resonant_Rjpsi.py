#!/usr/bin/env python
# coding: utf-8

from mybatch import *

import awkward as awk
import numpy as np
import uproot
from nanoframe import NanoFrame
import os
import particle
import pandas as pd
import uproot_methods


final_dfs = {
    'pf' : pd.DataFrame(),
    #'lpt' : pd.DataFrame(),
    #'mix' : pd.DataFrame(),
}

nprocessed = 0

from pdb import set_trace

#loop sui dataset

for dataset in [args.f_1,args.f_2,args.f_3]: 
    if(dataset==''):
        continue
    print("Opening file", dataset)
    f=open(dataset,"r")
    paths = f.readlines()
    final_dfs['pf'] = None
    #print(final_dfs['pf'])
    for i,fname in enumerate(paths):
        fname= fname.strip('\n')
        if(i%1==0):
            print("Processing file ", fname)
        #nf=NanoFrame('root://t3dcachedb03.psi.ch:1094///pnfs/psi.ch/cms/trivcat/store/user/friti/crab_job_2020Apr06/Charmonium/crab_data_Run2018D_pr_v2/200406_134114/0000/RJpsi_data_2020Apr06_535.root',)
        nf = NanoFrame(fname, )#branches = branches)
        # Load the needed collections, NanoFrame is just an empty shell until we call the collections
        evt = nf['event']
        muons = nf['Muon']
        bcands = nf['BTommm']
        hlt = nf['HLT']
        nprocessed += hlt.shape[0]

        # Attach the objects to the candidates
        bcands['mu1'] = muons[bcands['l1Idx']]
        bcands['mu2'] = muons[bcands['l2Idx']]
        bcands['k'] = muons[bcands['kIdx']]

        bcands['event'] = nf['event']
        bcands['run'] = nf['run']
        bcands['luminosityBlock'] = nf['luminosityBlock']    
        bcands['l_xy_sig'] = bcands.l_xy / np.sqrt(bcands.l_xy_unc)

        # Candidate selection, cut-based for the moment

        '''
        b_selection = (bcands.k.tightId == 1) & (bcands.mu1.mediumId == 1) & (bcands.mu2.mediumId == 1) & \
                  ((bcands.l_xy / bcands.l_xy_unc) > 10) & (bcands.svprob > 0.05) & (bcands.cos2D > 0.95) & \
                  (bcands.mass <6.3) & (bcands.k.p4.pt >4) & (bcands.p4.pt >15) & \
                  ((bcands.k_iso03/bcands.p4.pt)<0.2) & (bcands.m_jpsi < 3.1969) & (bcands.m_jpsi>2.9969) 

        b_pf = bcands.mu1.isPFcand & bcands.mu2.isPFcand & bcands.k.isPFcand
        '''
        b_selection = (bcands.k.tightId == 1) & (bcands.mu1.mediumId == 1) & (bcands.mu2.mediumId == 1)
        best_pf_cand = bcands[b_selection].svprob.argmax()
        bcands_pf = (bcands[b_selection][best_pf_cand]).flatten()

        dfs = {}

        for name, tab, sel in [
                ('pf', bcands_pf, b_selection), 
        ]:
            dfs[name] = pd.DataFrame()
            df = dfs[name]
            df['event'] = tab['event']
            df['run'] = tab['run']
            df['luminosityBlock'] = tab['luminosityBlock']
            df['mu1pt'] = tab.mu1.p4.pt
            df['mu2pt'] = tab.mu2.p4.pt
            df['kpt'] = tab.k.p4.pt
            df['mu1mass'] = tab.mu1.p4.mass
            df['mu2mass'] = tab.mu2.p4.mass
            df['kmass'] = tab.k.p4.mass
            df['mu1_isPF'] = tab.mu1.isPFcand
            df['mu2_isPF'] = tab.mu2.isPFcand
            df['k_isPF'] = tab.k.isPFcand
            df['Bcharge'] = tab.charge
            df['Bpt'] = tab.p4.pt
            df['Beta'] = tab.p4.eta
            df['Bsvprob'] = tab.svprob
            df['Bcos2D'] = tab.cos2D
            df['Blxy_sig'] = (tab.l_xy / tab.l_xy_unc)
            df['Bmll_raw'] = tab.m_jpsi
            df['Bm_miss_sq'] = tab.m_miss_sq
            df['Bmass'] = tab.p4.mass
            df['Bb_iso03'] = tab.b_iso03
            df['Bb_iso04'] = tab.b_iso04
            df['Bk_iso03'] = tab.k_iso03
            df['Bk_iso04'] = tab.k_iso04
            df['Bl1_iso03'] = tab.l1_iso03
            df['Bl1_iso04'] = tab.l1_iso04
            df['Bl2_iso03'] = tab.l2_iso03
            df['Bl2_iso04'] = tab.l2_iso04
            df['BE_mu_star']=tab.E_mu_star
            df['BE_mu_canc']=tab.E_mu_canc
            df['BQ_sq']=tab.Q_sq
            df['Bpt_var']=tab.pt_var
            df['Bpt_miss_vec']=tab.pt_miss_vec
            df['BDR']=tab.DR
            df['nB'] = sel.sum()[sel.sum() != 0]

        #final_dfs['lpt'] = pd.concat((final_dfs['lpt'], dfs['lpt']))
        #print(final_dfs['pf'])
        final_dfs['pf'] = pd.concat((final_dfs['pf'], dfs['pf']))
        #print(final_dfs['pf'])
        #final_dfs['mix'] = pd.concat((final_dfs['mix'], dfs['mix']))
        #print(len(final_dfs['pf']))
        #final_dfs['lpt'].to_hdf(args.f_out, 'lpt', mode = 'w')
    dataset=dataset.strip('.txt')
    name=dataset.split('/')
    final_dfs['pf'].to_hdf('hd5_files/'+name[len(name)-1]+'.h5', 'pf', mode = 'w')
    print("Saved file "+ 'hd5_files/'+ name[len(name)-1]+'.h5')
    print("")
#final_dfs['mix'].to_hdf(args.f_out, 'mix', mode = 'a')
print('DONE! Processed events: ', nprocessed)
print('PF size:', dfs['pf'].shape[0])
''', 'LPT size:', 
      final_dfs['lpt'].shape[0], 'MIX size:', 
      final_dfs['mix'].shape[0])
'''
