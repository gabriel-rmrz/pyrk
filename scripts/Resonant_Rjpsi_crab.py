#from nanoAOD root files, to 1D hd5 files

from mybatch import *
from coffea.analysis_objects import JaggedCandidateArray
import awkward as awk
import numpy as np
import uproot
from nanoframe import NanoFrame
import os
import particle
import pandas as pd
import uproot_methods
import ROOT
from pdb import set_trace
from root_pandas import to_root
from uproot_methods import TLorentzVectorArray
from uproot_methods import TVector3Array
from uproot_methods import TLorentzVector
from uproot_methods import TVector3
from scipy.constants import c as speed_of_light

maxEvents = -1
checkDoubles = True

nMaxFiles = HOOK_MAX_FILES
skipFiles = HOOK_SKIP_FILES


## lifetime weights ##
def weight_to_new_ctau(old_ctau, new_ctau, ct):
    '''
    Returns an event weight based on the ratio of the normalised lifetime distributions.
    old_ctau: ctau used for the sample production
    new_ctau: target ctau
    ct      : per-event lifetime
    '''
    weight = old_ctau/new_ctau * np.exp( (1./old_ctau - 1./new_ctau) * ct )
    return weight
    
def lifetime_weight(pf, fake = True):
    print("Adding lifetime weight branch...")
    if fake:
        ctau_weight_central = np.ones(len(pf))
        ctau_weight_up = np.ones(len(pf))
        ctau_weight_down = np.ones(len(pf))
        pf['ctau_weight_central'] = ctau_weight_central
        pf['ctau_weight_up'] = ctau_weight_up
        pf['ctau_weight_down'] = ctau_weight_down
        return pf
    else:
        Bc_mass = 6.274
        ctau_pdg    = 0.510e-12 * speed_of_light * 1000. # in mm
        ctau_actual = 0.1358
        ctau_up     = (0.510+0.009)*1e-12 * speed_of_light * 1000. # in mm
        ctau_down   = (0.510-0.009)*1e-12 * speed_of_light * 1000. # in mm
        
        ctau_weight_central = []
        ctau_weight_up = []
        ctau_weight_down = []

        for i in range(len(pf)):
            flag = 0
            #jpsi vertex
            if( abs(pf.mu1_mother_pdgId[i]) == 443 ):
                jpsi_vertex = TVector3(pf.mu1_mother_vx[i],pf.mu1_mother_vy[i],pf.mu1_mother_vz[i])
            elif( abs(pf.mu2_mother_pdgId[i]) == 443 ):
                jpsi_vertex = TVector3(pf.mu2_mother_vx[i],pf.mu2_mother_vy[i],pf.mu2_mother_vz[i])
             
            else: 
                flag = 1
        
            #Bc vertex
            if(abs(pf.mu1_grandmother_pdgId[i]) == 541):
                Bc_vertex = TVector3(pf.mu1_grandmother_vx[i],pf.mu1_grandmother_vy[i],pf.mu1_grandmother_vz[i])
                Bc_p4 = TLorentzVector.from_ptetaphim(pf.mu1_grandmother_pt[i],pf.mu1_grandmother_eta[i],pf.mu1_grandmother_phi[i],Bc_mass)
            elif(abs(pf.mu2_grandmother_pdgId[i]) == 541):
                Bc_vertex = TVector3(pf.mu2_grandmother_vx[i],pf.mu2_grandmother_vy[i],pf.mu2_grandmother_vz[i])
                Bc_p4 = TLorentzVector.from_ptetaphim(pf.mu2_grandmother_pt[i],pf.mu2_grandmother_eta[i],pf.mu2_grandmother_phi[i],Bc_mass)

            else:
                flag = 1
        
            if(flag == 1):
                ctau_weight_central.append(1)
                ctau_weight_up.append (1)
                ctau_weight_down.append(1)
       
            else:
                # distance
                lxyz = (jpsi_vertex - Bc_vertex).mag
                beta = Bc_p4.beta
                gamma = Bc_p4.gamma
                ct = lxyz/(beta * gamma)
                #print(lxyz,beta,gamma,ct)
                ctau_weight_central.append( weight_to_new_ctau(ctau_actual, ctau_pdg , ct*10.))
                ctau_weight_up.append (weight_to_new_ctau(ctau_actual, ctau_up  , ct*10.))
                ctau_weight_down.append(weight_to_new_ctau(ctau_actual, ctau_down, ct*10.))

        pf['ctau_weight_central'] = ctau_weight_central
        pf['ctau_weight_up'] = ctau_weight_up
        pf['ctau_weight_down'] = ctau_weight_down
        return pf
## end lifetime weights ##

def DR_jpsimu(pf):
    print("Adding DR between jpsi and mu branch...")
    mu1_p4 = TLorentzVectorArray.from_ptetaphim(pf.mu1pt,pf.mu1eta,pf.mu1phi,pf.mu1mass)
    mu2_p4 = TLorentzVectorArray.from_ptetaphim(pf.mu2pt,pf.mu2eta,pf.mu2phi,pf.mu2mass)

    jpsi_p4= mu1_p4 + mu2_p4 
    #    jpsi_p4 = TLorentzVectorArray.from_ptetaphim((pf.mu1pt+pf.mu2pt),(pf.mu1eta,+pf.mu2eta),(pf.mu1phi+pf.mu2phi),(pf.mu1mass+pf.mu2mass))
    mu_p4 = TLorentzVectorArray.from_ptetaphim(pf.kpt,pf.keta,pf.kphi,pf.kmass)
    #    print(jpsi_p4.delta_r(mu_p4))
    #    pf.copy()
    pf['DR_jpsimu'] = jpsi_p4.delta_r(mu_p4)
    return pf

def jpsi_branches(pf):
    print("Adding jpsi four momentum branches...")
    mu1_p4 = TLorentzVectorArray.from_ptetaphim(pf.mu1pt,pf.mu1eta,pf.mu1phi,pf.mu1mass)
    mu2_p4 = TLorentzVectorArray.from_ptetaphim(pf.mu2pt,pf.mu2eta,pf.mu2phi,pf.mu2mass)
    jpsi_p4= mu1_p4 + mu2_p4
    
    #pf = pf.copy()

    pf['jpsi_pt'] = jpsi_p4.pt
    pf['jpsi_eta'] = jpsi_p4.eta
    pf['jpsi_phi'] = jpsi_p4.phi
    pf['jpsi_mass'] = jpsi_p4.mass
    return pf
    
# two decay time branches... which one to choose?
def decaytime(pf):
    print("Adding decay time branch...")
    PV_pos = TVector3Array(pf.pv_x,pf.pv_y,pf.pv_z)
    jpsiVertex_pos1 = TVector3Array(pf.mu1_vx,pf.mu1_vy,pf.mu1_vz)
    jpsiVertex_pos2 = TVector3Array(pf.mu2_vx,pf.mu2_vy,pf.mu2_vz)

    dist1 = (PV_pos - jpsiVertex_pos1).mag
    dist2 = (PV_pos - jpsiVertex_pos2).mag

    decay_time1 = dist1 * 6.276 / (pf.Bpt_reco * 2.998e+10)
    decay_time2 = dist2 * 6.276 / (pf.Bpt_reco * 2.998e+10)
    #pf.copy()
    pf['decay_time1'] = decay_time1
    pf['decay_time2'] = decay_time2
    return pf


nprocessedAll = 0
#loop on datasets
for dataset in [args.data,args.mc_mu,args.mc_tau,args.mc_x,args.mc_onia,args.mc_gen]: 
    if(dataset==''):
        continue
    print(" ")
    
    #    if (dataset == args.data or dataset == args.mc_x ):
    channels = ['BTommm','BTopmm','BTokmm']
    
    '''else:
        channels =['BTommm']
    '''
    print("Opening file", dataset)
    f=open(dataset,"r")
    paths = f.readlines()

    #################################################################################
    # MC BcToXToJpsi #
    #################################################################################
    if(dataset == args.mc_x):

        final_dfs_mmm = {
            'is_jpsi_mu' : pd.DataFrame(),
            'is_jpsi_tau' : pd.DataFrame(),
            'is_jpsi_pi' : pd.DataFrame(),
            'is_psi2s_mu' :pd.DataFrame(),
            'is_chic0_mu' : pd.DataFrame(),
            'is_chic1_mu' : pd.DataFrame(),
            'is_chic2_mu' : pd.DataFrame(),
            'is_hc_mu' : pd.DataFrame(),
            'is_psi2s_tau' : pd.DataFrame(),
            'is_jpsi_3pi' : pd.DataFrame(),
            'is_jpsi_hc' : pd.DataFrame(),
        }
        
        final_dfs_pmm = {
            'is_jpsi_mu' : pd.DataFrame(),
            'is_jpsi_tau' : pd.DataFrame(),
            'is_jpsi_pi' : pd.DataFrame(),
            'is_psi2s_mu' :pd.DataFrame(),
            'is_chic0_mu' : pd.DataFrame(),
            'is_chic1_mu' : pd.DataFrame(),
            'is_chic2_mu' : pd.DataFrame(),
            'is_hc_mu' : pd.DataFrame(),
            'is_psi2s_tau' : pd.DataFrame(),
            'is_jpsi_3pi' : pd.DataFrame(),
            'is_jpsi_hc' : pd.DataFrame(),
        }
        
        final_dfs_kmm = {
            'is_jpsi_mu' : pd.DataFrame(),
            'is_jpsi_tau' : pd.DataFrame(),
            'is_jpsi_pi' : pd.DataFrame(),
            'is_psi2s_mu' :pd.DataFrame(),
            'is_chic0_mu' : pd.DataFrame(),
            'is_chic1_mu' : pd.DataFrame(),
            'is_chic2_mu' : pd.DataFrame(),
            'is_hc_mu' : pd.DataFrame(),
            'is_psi2s_tau' : pd.DataFrame(),
            'is_jpsi_3pi' : pd.DataFrame(),
            'is_jpsi_hc' : pd.DataFrame(),
        }
        
        '''
        final_dfs_mmm['is_jpsi_mu']=None
        final_dfs_mmm['is_jpsi_tau']=None
        final_dfs_mmm['is_jpsi_pi']=None
        final_dfs_mmm['is_bcbkg']=None
        
        final_dfs_pmm['is_jpsi_mu']=None
        final_dfs_pmm['is_jpsi_tau']=None
        final_dfs_pmm['is_jpsi_pi']=None
        final_dfs_pmm['is_bcbkg']=None
        
        final_dfs_kmm['is_jpsi_mu']=None
        final_dfs_kmm['is_jpsi_tau']=None
        final_dfs_kmm['is_jpsi_pi']=None
        final_dfs_kmm['is_bcbkg']=None
        '''
    else:
        final_dfs_mmm = {
            'ptmax' : pd.DataFrame(),
        }

        final_dfs_pmm = {
            'ptmax' : pd.DataFrame(),
        }

        final_dfs_kmm = {
            'ptmax' : pd.DataFrame(),
        }

    
    nprocessedDataset = 0
    nFiles = 0
    for i,fname in enumerate(paths):
        fname= fname.strip('\n')
        if(i%1==0):
            print("Processing file ", fname)
        
        
        if(i < skipFiles): # if I want to skip 1 file, I want to skip i=0 -> i+1
            print("Skipping the file...")
            continue
           
        if(nFiles >= nMaxFiles and nMaxFiles != -1):
            print("No! Max number of files reached!")
            break
        nFiles +=1
       

        for channel in channels:
            print("In channel "+channel)
            # Load the needed collections, NanoFrame is just an empty shell until we call the collections
            nf = NanoFrame(fname, )#branches = branches)
            evt = nf['event']
            muons = nf['Muon']
            bcands = nf[channel]
            hlt = nf['HLT']
            gen= nf['GenPart']
            bcands['event'] = nf['event']
            bcands['run'] = nf['run']
            bcands['luminosityBlock'] = nf['luminosityBlock']    
            if(len(nf[channel]) == 0):
                print("Empty file! Check why?")
                continue
            bcands['l_xy_sig'] = bcands.l_xy / np.sqrt(bcands.l_xy_unc)
            bcands['pv'] = nf['PV']
            bcands['PV_npvsGood'] = nf['PV_npvsGood']
            
            #number of events processed
            nprocessedDataset += hlt.shape[0]
            nprocessedAll+=hlt.shape[0]

            #add muon infos        
            mu1 = JaggedCandidateArray.zip({n: muons[bcands['l1Idx']][n] for n in muons[bcands['l1Idx']].columns})
            mu2 = JaggedCandidateArray.zip({n: muons[bcands['l2Idx']][n] for n in muons[bcands['l2Idx']].columns})
            if channel == 'BTommm':
                k = JaggedCandidateArray.zip({n: muons[bcands['kIdx']][n] for n in muons[bcands['kIdx']].columns})
            else:
                tracks = nf['ProbeTracks']
                k = JaggedCandidateArray.zip({n: tracks[bcands['kIdx']][n] for n in tracks[bcands['kIdx']].columns})
                

            # the MC ONia needs a specific treatment for gen, because there are cases in which thegen collection is zero! And in these cases we can not directly add the gen branch to the muons branch. So, consdiering that to use the MC Onia we need an additional requirement: that the third muon iha pdgId==+-13, we can mask the events in which k.genPrtIdx==-1
            if(dataset==args.mc_onia):
                mask = (k.genPartIdx != -1)
                mu1_new = mu1[mask]
                mu2_new = mu2[mask]
                k_new = k[mask]
                
                mu1 = JaggedCandidateArray.zip({n: mu1_new[n] for n in mu1_new.columns})
                mu2 = JaggedCandidateArray.zip({n: mu2_new[n] for n in mu2_new.columns})
                k = JaggedCandidateArray.zip({n: k_new[n] for n in k_new.columns})
            
                bcands = JaggedCandidateArray.zip({n: bcands[mask][n] for n in bcands[mask].columns})

            
            # add gen info as a column of the muon
            if (dataset!=args.data):
                #pile up weights only for mc
                
                bcands['puWeight'] = nf['puWeight']
                bcands['puWeightUp'] = nf['puWeightUp']
                bcands['puWeightDown'] = nf['puWeightDown']
                
                mu1['gen'] = gen[mu1.genPartIdx]
                mu2['gen'] = gen[mu2.genPartIdx]
                k['gen'] = gen[k.genPartIdx]
                
                mu1['mother'] = gen[gen[mu1.genPartIdx].genPartIdxMother]
                mu2['mother'] = gen[gen[mu2.genPartIdx].genPartIdxMother]
                k['mother'] = gen[gen[k.genPartIdx].genPartIdxMother]
                
                mu1['grandmother'] = gen[gen[gen[mu1.genPartIdx].genPartIdxMother].genPartIdxMother]
                mu2['grandmother'] = gen[gen[gen[mu2.genPartIdx].genPartIdxMother].genPartIdxMother]
                k['grandmother'] = gen[gen[gen[k.genPartIdx].genPartIdxMother].genPartIdxMother]
                
                if (dataset!=args.mc_mu):
                    mu1['grandgrandmother'] =gen[gen[gen[gen[mu1.genPartIdx].genPartIdxMother].genPartIdxMother].genPartIdxMother]
                    mu2['grandgrandmother'] =gen[gen[gen[gen[mu2.genPartIdx].genPartIdxMother].genPartIdxMother].genPartIdxMother]
                    k['grandgrandmother'] =gen[gen[gen[gen[k.genPartIdx].genPartIdxMother].genPartIdxMother].genPartIdxMother]                

            bcands['mu1']= mu1
            bcands['mu2'] = mu2
            bcands['k'] = k
            
            b_selection = ((bcands.p4.mass < 10 ) & (bcands.svprob > 1e-7))
            x_selection= (bcands.k.p4.pt > -99)
        
        
            #Delete the signal from the JpsiX MC
            if (dataset==args.mc_onia):
                x_selection= ~ ((bcands.k.genPartIdx>=0) & ( bcands.mu1.genPartIdx>=0) & (bcands.mu2.genPartIdx>=0) & (abs(bcands.mu1.mother.pdgId) == 443) & (abs(bcands.mu2.mother.pdgId) == 443) & (abs(bcands.mu1.grandmother.pdgId) == 541) & (abs(bcands.mu2.grandmother.pdgId) == 541) & ( (abs(bcands.k.mother.pdgId)==541) | ( (abs(bcands.k.mother.pdgId)==15) & (abs(bcands.k.grandmother.pdgId)== 541))))

            if(dataset == args.mc_x):
            
                jpsi_tau_sel = (bcands.is_jpsi_tau == 1)
                jpsi_mu_sel = (bcands.is_jpsi_mu == 1)
                jpsi_pi_sel = (bcands.is_jpsi_pi == 1)
                psi2s_mu_sel = (bcands.is_psi2s_mu == 1)
                chic0_mu_sel = (bcands.is_chic0_mu == 1)
                chic1_mu_sel = (bcands.is_chic1_mu == 1)
                chic2_mu_sel = (bcands.is_chic2_mu == 1)
                hc_mu_sel = (bcands.is_hc_mu == 1)
                psi2s_tau_sel = (bcands.is_psi2s_tau == 1)
                jpsi_3pi_sel = (bcands.is_jpsi_3pi == 1)
                jpsi_hc_sel = (bcands.is_jpsi_hc == 1)

                flag_selection= [jpsi_tau_sel,jpsi_mu_sel,jpsi_pi_sel,psi2s_mu_sel,chic0_mu_sel,chic1_mu_sel,chic2_mu_sel,hc_mu_sel,psi2s_tau_sel,jpsi_3pi_sel,jpsi_hc_sel]
                flag_names = ['is_jpsi_tau','is_jpsi_mu','is_jpsi_pi','is_psi2s_mu','is_chic0_mu','is_chic1_mu','is_chic2_mu','is_hc_mu','is_psi2s_tau','is_jpsi_3pi','is_jpsi_hc']
            else:
                flag_selection = [(bcands.p4.pt>-99)]
                flag_names = ['ptmax']
            for selection,name in zip(flag_selection, flag_names):
                best_pf_cand_pt = bcands[b_selection & x_selection & selection ].p4.pt.argmax() #B con pt massimo
                bcands_flag = (bcands[b_selection & x_selection & selection][best_pf_cand_pt]).flatten()
                dfs = {}
            
                for chan, tab, sel in [
                        (channel, bcands_flag, b_selection & x_selection & selection), 
                ]:
                    dfs[name] = pd.DataFrame()
                    df = dfs[name]
                    df['event'] = tab['event']
                    df['run'] = tab['run']
                    df['luminosityBlock'] = tab['luminosityBlock']
                    # Number of vertices
                
                    #Kinematic
                    df['mu1pt'] = tab.mu1.p4.pt
                    df['mu2pt'] = tab.mu2.p4.pt
                    df['kpt'] = tab.k.p4.pt
                    df['mu1mass'] = tab.mu1.p4.mass
                    df['mu2mass'] = tab.mu2.p4.mass
                    df['kmass'] = tab.k.p4.mass
                    df['mu1phi'] = tab.mu1.p4.phi
                    df['mu2phi'] = tab.mu2.p4.phi
                    df['kphi'] = tab.k.p4.phi
                    df['mu1eta'] = tab.mu1.p4.eta
                    df['mu2eta'] = tab.mu2.p4.eta
                    df['keta'] = tab.k.p4.eta
                    
                    df['Bpt'] = tab.p4.pt
                    df['Bmass'] = tab.p4.mass
                    df['Beta'] = tab.p4.eta
                    df['Bphi'] = tab.p4.phi
                    df['Bpt_reco'] = (tab.p4.pt * 6.275 / tab.p4.mass)
                    
                    df['mu1_dxy'] = tab.mu1.dxy
                    df['mu2_dxy'] = tab.mu2.dxy
                    df['k_dxy'] = tab.k.dxy
                    df['mu1_dxy_sig'] = tab.mu1.dxyErr
                    df['mu2_dxy_sig'] = tab.mu2.dxyErr
                
                    df['mu1_dz'] = tab.mu1.dz
                    df['mu2_dz'] = tab.mu2.dz
                    df['k_dz'] = tab.k.dz
                    df['mu1_dz_sig'] = tab.mu1.dzErr
                    df['mu2_dz_sig'] = tab.mu2.dzErr             
                    if(chan == 'BTommm'):
                        df['k_dxy_sig'] = tab.k.dxyErr
                        df['k_dz_sig'] = tab.k.dzErr
                        
                    else:
                        df['k_dxy_sig'] = tab.k.dxyS
                        df['k_dz_sig'] = tab.k.dzS
                        
                    #B Vertex properties
                    df['Blxy_sig'] = (tab.l_xy / tab.l_xy_unc)
                    df['Blxy'] = tab.l_xy
                    df['Blxy_unc'] = tab.l_xy_unc
                    df['Bsvprob'] = tab.svprob
                    df['Bcos2D'] = tab.cos2D
                    if (dataset== args.data or dataset== args.mc_x): #only because others dont have it for now
                        df['Bchi2'] = tab.chi2

                        #Jpsi vertex properties
                        '''
                        if(chan == 'BTommm'):
                            df['jpsi_chi2'] = tab.jpsi_chi2
                            df['jpsi_svprob'] = tab.jpsi_svprob
                        '''
                    
                    df['mu1_vx'] = tab.mu1.vx
                    df['mu2_vx'] = tab.mu2.vx
                    df['k_vx'] = tab.k.vx
                    
                    df['mu1_vy'] = tab.mu1.vy
                    df['mu2_vy'] = tab.mu2.vy
                    df['k_vy'] = tab.k.vy
                    
                    df['mu1_vz'] = tab.mu1.vz
                    df['mu2_vz'] = tab.mu2.vz
                    df['k_vz'] = tab.k.vz
                    
                    #PV position
                    df['pv_x'] = tab.pv.x
                    df['pv_y'] = tab.pv.y
                    df['pv_z'] = tab.pv.z
                    
                    df['npv_good'] = tab.PV_npvsGood
                    #our variables
                    df['m_miss_sq'] = tab.m_miss_sq
                    df['Q_sq']=tab.Q_sq
                    df['pt_var']=tab.pt_var
                    df['pt_miss_vec']=tab.pt_miss_vec
                    df['pt_miss_scal'] = tab.pt_miss
                    df['DR_mu1mu2']=tab.DR 
                    #            df['DR_jpsi_mu'] =   # next step, here it is difficult
                    df['E_mu_star']=tab.E_mu_star
                    df['E_mu_canc']=tab.E_mu_canc
                    
                    
                    #IDs
                    df['mu1_mediumID']= tab.mu1.mediumId
                    df['mu2_mediumID']= tab.mu2.mediumId
                    df['mu1_tightID']= tab.mu1.tightId
                    df['mu2_tightID']= tab.mu2.tightId
                    df['mu1_softID']= tab.mu1.softId
                    df['mu2_softID']= tab.mu2.softId
                    
                    if(chan == 'BTommm'):
                        df['k_tightID']= tab.k.tightId
                        df['k_mediumID']=tab.k.mediumId
                        df['k_softID']=tab.k.softId
                    
                    #is PF ?
                    df['mu1_isPF'] = tab.mu1.isPFcand
                    df['mu2_isPF'] = tab.mu2.isPFcand
                    if(chan == 'BTommm'):
                        df['k_isPF'] = tab.k.isPFcand
                        
                    #iso
                    df['b_iso03'] = tab.b_iso03
                    df['b_iso04'] = tab.b_iso04
                    df['k_iso03'] = tab.k_iso03
                    df['k_iso04'] = tab.k_iso04
                    df['l1_iso03'] = tab.l1_iso03
                    df['l1_iso04'] = tab.l1_iso04
                    df['l2_iso03'] = tab.l2_iso03
                    df['l2_iso04'] = tab.l2_iso04
                    
                    #others
                    df['Bcharge'] = tab.charge
                    df['mll_raw'] = tab.m_jpsi
                    
                    #B vertex (we don't use this)
                    df['Bvtx_x'] = tab.vtx_x
                    df['Bvtx_y'] = tab.vtx_y
                    df['Bvtx_z'] = tab.vtx_z
                    df['Bvtx_ex'] = tab.vtx_ex
                    df['Bvtx_ey'] = tab.vtx_ey
                    df['Bvtx_ez'] = tab.vtx_ez

                    #fit values
                    df['fit_Beta'] = tab.fit_eta
                    df['fit_Bphi'] = tab.fit_phi
                    df['fit_Bpt'] = tab.fit_pt
                    df['fit_Bmass'] = tab.fit_mass

                    df['fit_mu1_eta'] = tab.fit_l1_eta
                    df['fit_mu1_phi'] = tab.fit_l1_phi
                    df['fit_mu1_pt'] = tab.fit_l1_pt

                    df['fit_mu2_eta'] = tab.fit_l2_eta
                    df['fit_mu2_phi'] = tab.fit_l2_phi
                    df['fit_mu2_pt'] = tab.fit_l2_pt

                    df['fit_k_eta'] = tab.fit_k_eta
                    df['fit_k_phi'] = tab.fit_k_phi
                    df['fit_k_pt'] = tab.fit_k_pt

                    df['fit_Bmass_err'] = tab.fit_massErr
                    df['fit_Bcos2D'] = tab.fit_cos2D
                    
                    df['nB'] = sel.sum()[sel.sum() != 0]
                    
                    #Gen variables
                    if(dataset!=args.data):
                        
                        #PU weight
                        
                        df['puWeight'] = tab.puWeight
                        df['puWeightUp'] = tab.puWeightUp
                        df['puWeightDown'] = tab.puWeightDown
                        
                        #gen Part Flavour e gen Part Idx  -> if I need to access the gen info, this values tell me is it is a valid info or not
                        df['mu1_genPartFlav'] = tab.mu1.genPartFlav
                        df['mu2_genPartFlav'] = tab.mu2.genPartFlav
                        df['k_genPartFlav'] = tab.k.genPartFlav
                        
                        df['mu1_genPartIdx'] = tab.mu1.genPartIdx
                        df['mu2_genPartIdx'] = tab.mu2.genPartIdx
                        df['k_genPartIdx'] = tab.k.genPartIdx
                        
                        #lifetime (gen info)
                        df['mu1_gen_vx'] = tab.mu1.gen.vx
                        df['mu2_gen_vx'] = tab.mu2.gen.vx
                        df['k_gen_vx'] = tab.k.gen.vx
                        
                        df['mu1_gen_vy'] = tab.mu1.gen.vy
                        df['mu2_gen_vy'] = tab.mu2.gen.vy
                        df['k_gen_vy'] = tab.k.gen.vy
                    
                        df['mu1_gen_vz'] = tab.mu1.gen.vz
                        df['mu2_gen_vz'] = tab.mu2.gen.vz
                        df['k_gen_vz'] = tab.k.gen.vz
                        
                        #pdgId
                        df['mu1_pdgId'] = tab.mu1.pdgId
                        df['mu2_pdgId'] = tab.mu2.pdgId
                        df['k_pdgId'] = tab.k.pdgId
                        
                        df['mu1_genpdgId'] = tab.mu1.gen.pdgId
                        df['mu2_genpdgId'] = tab.mu2.gen.pdgId
                        df['k_genpdgId'] = tab.k.gen.pdgId
                        
                        #mother info
                        df['mu1_mother_pdgId'] = tab.mu1.mother.pdgId
                        df['mu2_mother_pdgId'] = tab.mu2.mother.pdgId
                        df['k_mother_pdgId'] = tab.k.mother.pdgId
                        
                        df['mu1_mother_pt'] = tab.mu1.mother.p4.pt
                        df['mu2_mother_pt'] = tab.mu2.mother.p4.pt
                        df['k_mother_pt'] = tab.k.mother.p4.pt
                        
                        df['mu1_mother_eta'] = tab.mu1.mother.p4.eta
                        df['mu2_mother_eta'] = tab.mu2.mother.p4.eta
                        df['k_mother_eta'] = tab.k.mother.p4.eta
                        
                        df['mu1_mother_phi'] = tab.mu1.mother.p4.phi
                        df['mu2_mother_phi'] = tab.mu2.mother.p4.phi
                        df['k_mother_phi'] = tab.k.mother.p4.phi
                        
                        df['mu1_mother_vx'] = tab.mu1.mother.vx
                        df['mu2_mother_vx'] = tab.mu2.mother.vx
                        df['k_mother_vx'] = tab.k.mother.vx
                        df['mu1_mother_vy'] = tab.mu1.mother.vy
                        df['mu2_mother_vy'] = tab.mu2.mother.vy
                        df['k_mother_vy'] = tab.k.mother.vy
                        df['mu1_mother_vz'] = tab.mu1.mother.vz
                        df['mu2_mother_vz'] = tab.mu2.mother.vz
                        df['k_mother_vz'] = tab.k.mother.vz
                        
                        #grandmother info
                        df['mu1_grandmother_pdgId'] = tab.mu1.grandmother.pdgId
                        df['mu2_grandmother_pdgId'] = tab.mu2.grandmother.pdgId
                        df['k_grandmother_pdgId'] = tab.k.grandmother.pdgId
                        
                        df['mu1_grandmother_pt'] = tab.mu1.grandmother.p4.pt
                        df['mu2_grandmother_pt'] = tab.mu2.grandmother.p4.pt
                        df['k_grandmother_pt'] = tab.k.grandmother.p4.pt
                        
                        df['mu1_grandmother_eta'] = tab.mu1.grandmother.p4.eta
                        df['mu2_grandmother_eta'] = tab.mu2.grandmother.p4.eta
                        df['k_grandmother_eta'] = tab.k.grandmother.p4.eta
                        
                        df['mu1_grandmother_phi'] = tab.mu1.grandmother.p4.phi
                        df['mu2_grandmother_phi'] = tab.mu2.grandmother.p4.phi
                        df['k_grandmother_phi'] = tab.k.grandmother.p4.phi
                        
                        df['mu1_grandmother_vx'] = tab.mu1.grandmother.vx
                        df['mu2_grandmother_vx'] = tab.mu2.grandmother.vx
                        df['k_grandmother_vx'] = tab.k.grandmother.vx
                        df['mu1_grandmother_vy'] = tab.mu1.grandmother.vy
                        df['mu2_grandmother_vy'] = tab.mu2.grandmother.vy
                        df['k_grandmother_vy'] = tab.k.grandmother.vy
                        df['mu1_grandmother_vz'] = tab.mu1.grandmother.vz
                        df['mu2_grandmother_vz'] = tab.mu2.grandmother.vz
                        df['k_grandmother_vz'] = tab.k.grandmother.vz
                        
                        if (dataset!=args.mc_mu):
                            #grand grand mother info
                            df['mu1_grandgrandmother_pdgId'] = tab.mu1.grandgrandmother.pdgId
                            df['mu2_grandgrandmother_pdgId'] = tab.mu2.grandgrandmother.pdgId
                            df['k_grandgrandmother_pdgId'] = tab.k.grandgrandmother.pdgId
                            
                            df['mu1_grandgrandmother_pt'] = tab.mu1.grandgrandmother.p4.pt
                            df['mu2_grandgrandmother_pt'] = tab.mu2.grandgrandmother.p4.pt
                            df['k_grandgrandmother_pt'] = tab.k.grandgrandmother.p4.pt
                            
                            df['mu1_grandgrandmother_eta'] = tab.mu1.grandgrandmother.p4.eta
                            df['mu2_grandgrandmother_eta'] = tab.mu2.grandgrandmother.p4.eta
                            df['k_grandgrandmother_eta'] = tab.k.grandgrandmother.p4.eta
                            
                            df['mu1_grandgrandmother_phi'] = tab.mu1.grandgrandmother.p4.phi
                            df['mu2_grandgrandmother_phi'] = tab.mu2.grandgrandmother.p4.phi
                            df['k_grandgrandmother_phi'] = tab.k.grandgrandmother.p4.phi
                            
                            df['mu1_grandgrandmother_vx'] = tab.mu1.grandgrandmother.vx
                            df['mu2_grandgrandmother_vx'] = tab.mu2.grandgrandmother.vx
                            df['k_grandgrandmother_vx'] = tab.k.grandgrandmother.vx
                            df['mu1_grandgrandmother_vy'] = tab.mu1.grandgrandmother.vy
                            df['mu2_grandgrandmother_vy'] = tab.mu2.grandgrandmother.vy
                            df['k_grandgrandmother_vy'] = tab.k.grandgrandmother.vy
                            df['mu1_grandgrandmother_vz'] = tab.mu1.grandgrandmother.vz
                            df['mu2_grandgrandmother_vz'] = tab.mu2.grandgrandmother.vz
                            df['k_grandgrandmother_vz'] = tab.k.grandgrandmother.vz
                    
                        if(dataset == args.mc_x):
                            df['weightGen'] = tab.weightGen
                            '''
                            if (name == 'is_bcbkg'):
                                df['is_psi2s_mu'] = tab.is_psi2s_mu
                                df['is_chic0_mu'] = tab.is_chic0_mu
                                df['is_chic1_mu'] = tab.is_chic1_mu
                                df['is_chic2_mu'] = tab.is_chic2_mu
                                df['is_hc_mu'] = tab.is_hc_mu
                                df['is_psi2s_tau'] = tab.is_psi2s_tau
                                df['is_jpsi_3pi'] = tab.is_jpsi_3pi
                                df['is_jpsi_hc'] = tab.is_jpsi_hc
                            '''
                    if(dataset == args.mc_mu or dataset == args.mc_tau or dataset == args.mc_x):
                        df = lifetime_weight(df, fake = False)
                    else:
                        df = lifetime_weight(df)
                    df = jpsi_branches(df)
                    df = DR_jpsimu(df)
                    df = decaytime(df)
                        
                    #print("Finito di processare la flag ",name, " la concateno a quella totale.")
                    
                    if(channel=='BTommm'):
                        final_dfs_mmm[name] = pd.concat((final_dfs_mmm[name], dfs[name])) # da scrievre esplicito?
                    elif(channel == 'BTopmm'):
                        final_dfs_pmm[name] = pd.concat((final_dfs_pmm[name], dfs[name])) # da scrievre esplicito
                    elif(channel == 'BTokmm'):
                        final_dfs_kmm[name] = pd.concat((final_dfs_kmm[name], dfs[name])) # da scrievre esplicito
                    if(nprocessedDataset > maxEvents and maxEvents != -1):
                        break

    dataset=dataset.strip('.txt')
    name=dataset.split('/')
    d=name[len(name)-1].split('_')
    adj=''
    if(dataset != args.mc_onia):
        adj = '_UL_flags'
    
    for flag in flag_names:
        for channel in channels:
            if channel == 'BTommm':
                final_dfs_mmm[flag].to_root('HOOK_FILE_OUT'+'_'+flag+'.root', key=channel)
            elif (channel == 'BTopmm'):
                final_dfs_pmm[flag].to_root('HOOK_FILE_OUT'+'_'+flag+'.root', key=channel, mode = 'a')
            elif (channel == 'BTokmm'):
                final_dfs_kmm[flag].to_root('HOOK_FILE_OUT'+'_'+flag+'.root', key=channel, mode = 'a')

        print("Saved file "+ 'HOOK_FILE_OUT'+'_'+flag+'.root')


print('DONE! Processed events: ', nprocessedAll)
