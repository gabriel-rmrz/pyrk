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

for dataset in [args.f_data,args.f_mc_mu,args.f_mc_tau]: 
    print("Opening file", dataset)
    f=open(dataset,"r")
    paths = f.readlines()
    final_dfs['pf'] = None
    #print(final_dfs['pf'])
    for fname in paths:
        fname= fname.strip('\n')
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

        # Attach the trigger muon, identified as the closest in dz to the lead electron
        '''
        muon_trg_mask = (muons.isTriggering == 1)
        for path, pt_thr in [('Mu8_IP5', 8), ('Mu10p5_IP3p5', 10), ('Mu8_IP3', 8), ('Mu8p5_IP3p5', 8.5), 
                             ('Mu9_IP5', 9), ('Mu7_IP4', 7), ('Mu9_IP4', 9), ('Mu9_IP6', 9), 
                             ('Mu8_IP6', 8), ('Mu12_IP6', 12)]:
            if not any(path in i for i in hlt.columns): # the trigger is not here
                continue
            else:
                #merge all the parts and compute an or
                hlt_fired = np.hstack(
                    [hlt[i].reshape((hlt[i].shape[0], 1)) for i in hlt.columns if path in i]
                ).any(axis = 1)
                muon_trg_mask = muon_trg_mask | (hlt_fired & (muons.p4.pt > pt_thr))

        one_trg_muon = (muon_trg_mask.sum() != 0)
        trig_mu = muons[muon_trg_mask][one_trg_muon]
        # bcands = bcands[one_trg_muon]

        # e1z, muz = bcands.e1.vz.cross(trig_mu.vz, nested = True).unzip()
        # closest_mu = np.abs(e1z - muz).argmin().flatten(axis = 1)
        # bcands['trg_mu'] = trig_mu[closest_mu]
        '''

        # Candidate selection, cut-based for the moment

        b_selection = (bcands.k.tightId == 1) & (bcands.mu1.mediumId == 1) & (bcands.mu2.mediumId == 1) & \
                  ((bcands.l_xy / bcands.l_xy_unc) > 10) & (bcands.svprob > 0.05) & (bcands.cos2D > 0.95) & \
                  (bcands.mass <6.3) & (bcands.k.p4.pt >4) & (bcands.p4.pt >15) & \
                  ((bcands.k_iso03/bcands.p4.pt)<0.2) & (bcands.m_jpsi < 3.1969) & (bcands.m_jpsi>2.9969) 

        b_pf = bcands.mu1.isPFcand & bcands.mu2.isPFcand & bcands.k.isPFcand
    #    b_lpt = bcands.e1.isLowPt & bcands.e2.isLowPt & (bcands.e1.mvaId > 3.96) & (bcands.e2.mvaId > 3.96)
    #    b_mix = (bcands.e1.isLowPt & (bcands.e1.mvaId > 3.96) & bcands.e2.isPF) | \
     #           (bcands.e2.isLowPt & (bcands.e2.mvaId > 3.96) & bcands.e1.isPF)


        best_pf_cand = bcands[b_selection & b_pf].svprob.argmax()
        bcands_pf = (bcands[b_selection & b_pf][best_pf_cand]).flatten()
        ''' 
        best_lpt_cand = bcands[b_selection & b_lpt].svprob.argmax()
        bcands_lpt = (bcands[b_selection & b_lpt][best_lpt_cand]).flatten()

        best_mix_cand = bcands[b_selection & b_mix].svprob.argmax()
        bcands_mix = (bcands[b_selection & b_mix][best_mix_cand]).flatten()
        '''
        dfs = {}

        for name, tab, sel in [
                ('pf', bcands_pf, b_selection & b_pf), 
           #     ('lpt', bcands_lpt, b_selection & b_lpt),
           #     ('mix', bcands_mix, b_selection & b_mix),
        ]:
            dfs[name] = pd.DataFrame()
            df = dfs[name]
            df['event'] = tab['event']
            df['run'] = tab['run']
            df['luminosityBlock'] = tab['luminosityBlock']
            df['mu1pt'] = tab.mu1.p4.pt
            df['mu2pt'] = tab.mu2.p4.pt
            df['kpt'] = tab.k.p4.pt
            df['Bcharge'] = tab.charge
            df['Bpt'] = tab.p4.pt
            df['Beta'] = tab.p4.eta
            df['Bsvprob'] = tab.svprob
            df['Bcos2D'] = tab.cos2D
            df['Blxy_sig'] = (tab.l_xy / tab.l_xy_unc)
            df['Bmll_raw'] = tab.m_jpsi
            df['Bm_miss_sq'] = tab.m_miss_sq
            df['Bmass'] = tab.p4.mass
            df['nB'] = sel.sum()[sel.sum() != 0]

        #final_dfs['lpt'] = pd.concat((final_dfs['lpt'], dfs['lpt']))
        #print(final_dfs['pf'])
        final_dfs['pf'] = pd.concat((final_dfs['pf'], dfs['pf']))
        #print(final_dfs['pf'])
        #final_dfs['mix'] = pd.concat((final_dfs['mix'], dfs['mix']))

#final_dfs['lpt'].to_hdf(args.f_out, 'lpt', mode = 'w')
    final_dfs['pf'].to_hdf(dataset.strip('.txt')+'.h5', 'pf', mode = 'w')
    print("Saved file "+dataset.strip('.txt')+'.h5')
    print("")
#final_dfs['mix'].to_hdf(args.f_out, 'mix', mode = 'a')
print('DONE! Processed events: ', nprocessed)
print('PF size:', dfs['pf'].shape[0])
''', 'LPT size:', 
      final_dfs['lpt'].shape[0], 'MIX size:', 
      final_dfs['mix'].shape[0])
'''
