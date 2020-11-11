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

maxEvents = -1
checkDoubles = True

nMaxFiles = 2
skipFiles = 0

channels = ['BTommm','BTopmm','BTokmm']


nprocessedAll = 0
#loop on datasets
for dataset in [args.data,args.mc_mu,args.mc_tau,args.mc_x,args.mc_onia,args.mc_gen]: 
    if(dataset==''):
        continue
    print(" ")
    print("Opening file", dataset)
    f=open(dataset,"r")
    paths = f.readlines()

    #################################################################################
    # MC BcToXToJpsi #
    #################################################################################
    if(dataset == args.mc_x):

        final_dfs_mmm = {
            'is_jpsi_mu' : pd.DataFrame(),
            'is_bcbkg' : pd.DataFrame(),
            'is_jpsi_tau' : pd.DataFrame(),
            'is_jpsi_pi' : pd.DataFrame(),
        }
        
        final_dfs_pmm = {
            'is_jpsi_mu' : pd.DataFrame(),
            'is_bcbkg' : pd.DataFrame(),
            'is_jpsi_tau' : pd.DataFrame(),
            'is_jpsi_pi' : pd.DataFrame(),
        }
        
        final_dfs_kmm = {
            'is_jpsi_mu' : pd.DataFrame(),
            'is_bcbkg' : pd.DataFrame(),
            'is_jpsi_tau' : pd.DataFrame(),
            'is_jpsi_pi' : pd.DataFrame(),
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
       

        nf = NanoFrame(fname, )#branches = branches)
        for channel in channels:
            print("In channel "+channel)
            # Load the needed collections, NanoFrame is just an empty shell until we call the collections
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
                '''
                bcands['puWeight'] = nf['puWeight']
                bcands['puWeightUp'] = nf['puWeightUp']
                bcands['puWeightDown'] = nf['puWeightDown']
                '''
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

            jpsi_tau_sel = (bcands.is_jpsi_tau == 1)
            jpsi_mu_sel = (bcands.is_jpsi_mu == 1)
            jpsi_pi_sel = (bcands.is_jpsi_pi == 1)
            bcbkg_sel = ((bcands.is_psi2s_mu == 1) | (bcands.is_chic0_mu == 1) | (bcands.is_chic1_mu == 1) | (bcands.is_chic2_mu == 1) | (bcands.is_hc_mu ==1 ) | (bcands.is_psi2s_tau ==1) | (bcands.is_jpsi_3pi==1) | (bcands.is_jpsi_hc ==1))

            
            if(dataset == args.mc_x):
                flag_selection= [jpsi_tau_sel,jpsi_mu_sel,jpsi_pi_sel,bcbkg_sel]
                flag_names = ['is_jpsi_tau','is_jpsi_mu','is_jpsi_pi','is_bcbkg']
            else:
                flag_selection = [(bcands.pt>-99)]
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
                        
                    #Vertex properties
                    df['Blxy_sig'] = (tab.l_xy / tab.l_xy_unc)
                    df['Blxy'] = tab.l_xy
                    df['Blxy_unc'] = tab.l_xy_unc
                    df['Bsvprob'] = tab.svprob
                    df['Bcos2D'] = tab.cos2D
                    df['Bchi2'] = tab.chi2
                    # decay time! Muon vtrx info and PV info!!
                    
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
                    
                    
                    df['nB'] = sel.sum()[sel.sum() != 0]
                    
                    #Gen variables
                    if(dataset!=args.data):
                        
                        #PU weight
                        '''
                        df['puWeight'] = tab.puWeight
                        df['puWeightUp'] = tab.puWeightUp
                        df['puWeightDown'] = tab.puWeightDown
                        '''
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
                final_dfs_mmm[flag].to_root('dataframes/'+d[0]+'_'+flag+adj+'.root', key=channel)
            elif (channel == 'BTopmm'):
                final_dfs_pmm[flag].to_root('dataframes/'+d[0]+'_'+flag+adj+'.root', key=channel, mode = 'a')
            elif (channel == 'BTokmm'):
                final_dfs_kmm[flag].to_root('dataframes/'+d[0]+'_'+flag+adj+'.root', key=channel, mode = 'a')

        print("Saved file dataframes/"+ d[0]+'_'+flag+adj+'.root')


print('DONE! Processed events: ', nprocessedAll)
