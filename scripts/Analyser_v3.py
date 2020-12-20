import sys
import math
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.model_selection import train_test_split
from xgboost import XGBClassifier, plot_importance
import xgboost as xgb
from collections import OrderedDict
from itertools import product
from sklearn.metrics         import roc_curve, roc_auc_score 
import seaborn as sns
from sklearn.model_selection import GridSearchCV
import sklearn
from uproot_methods.classes.TH1 import from_numpy
import uproot
from create_datacard_v2 import create_datacard_signal
from create_datacard_v2 import create_datacard_control
import pickle
from uproot_methods import TLorentzVectorArray
from uproot_methods import TVector3Array
from uproot_methods import TLorentzVector
from uproot_methods import TVector3
import root_numpy as rp
import ROOT
from variable import variable
from sample import sample_dic
#from sample import sample
from cmsstyle import CMS_lumi
from makePlot import makeSinglePlot
from makePlot import makeComparison
from makePlot import makeStack
import os
from scipy.constants import c as speed_of_light
from root_pandas import to_root
from root_pandas import read_root
import random

overflow_underflow = False   #overflow underflow option in plots
fac = False   #factor to the signal yield for the fit
flag = "v3" #flag for plot directory and dataset names

# fit options (only one at the time True for now)
fit = False
feature_importance = False
conf_matr = False

computeFR = False
prob_shapes_plot  = False

fitPassFail  = True

saveonwww = True
saveonwww_all = False
#BDT cut
cut = 1
cutStr = "1"

fr = 0.774791533164465
misId_norm=fr/(1-fr)

#no pop up canvases
ROOT.gROOT.SetBatch()
ROOT.gStyle.SetOptStat(0)


#branches for training MVA
'''
branches=[

    'Q_sq',

    'pt_miss_vec',

    'E_mu_star',

    'm_miss_sq',

    'Bmass',

    'kpt',

    'Bpt',
 
    'bvtx_cos2D',
 
    'DR_mu1mu2',
    
    'Beta',

    'pt_var',
    
    'mu1pt',

    'mu2pt',

    'mu1eta',

    'mu2eta',

    'keta',    

]
'''
branches = [
    'Bpt'         ,
    'Bmass'       ,
    'kpt'         ,
    'Bpt_reco'    ,
    'bvtx_lxy_sig'    ,
    'bvtx_svprob'     ,
#     'log10_svprob',
    'm_miss_sq'   ,
    'Q_sq'        ,
    'pt_var'      ,
    'pt_miss_vec' ,
    'E_mu_star'   ,
    'E_mu_canc'   ,
    'b_iso03' ,
    #    'dr13'        ,
    #    'dr23'        ,
    'DR_jpsimu'  ,
]
#BDT fit features function
def BDT_fit(X_train, y_train,X_valid,y_valid):
    
    clf = XGBClassifier(
        #max_depth        = 3,
        learning_rate    = 1e-2, # 0.01
        n_estimators     = 10000, # 400
        #silent           = False,
        #subsample        = 0.6,
        #colsample_bytree = 0.8,
        #min_child_weight = 1,
        #gamma            = 0, # 20
        #seed             = 123,
        #scale_pos_weight = 1,
        #reg_alpha        = 0.01,
        #reg_lambda       = 1.2,
        objective        = 'multi:softprob',  # error evaluation for multiclass training
        num_class       = 3,
    )

    clf.fit(
        X_train, 
        y_train,
        eval_set              = [(X_train, y_train), (X_valid, y_valid)],
        early_stopping_rounds = 1000,
        eval_metric           = 'mlogloss',
        verbose               = True,
        #sample_weight         = train['weight'],
    )
    
    pickle.dump(clf, open("bdtModel/BDT_Model_" +flag+ ".dat", "wb"))
    print("Model saved.")
    return clf


def createHisto(sample, var, norm=1, over = True):
    histo = ROOT.TH1F(sample.histoName,"" ,var.nbins,var.xmin,var.xmax)
    for item,pileup,ctau in zip(sample.df[var.name],sample.df['puWeight'],sample.df['ctau_weight_central']):
            histo.Fill(item,pileup*ctau*norm)

    if(over == True):
        histo.SetBinContent(1, histo.GetBinContent(0) + histo.GetBinContent(1))
        histo.SetBinError(1, math.sqrt(pow(histo.GetBinError(0),2) + pow(histo.GetBinError(1),2)))
        histo.SetBinContent(var.nbins, histo.GetBinContent(var.nbins) + histo.GetBinContent(var.nbins+1))
        histo.SetBinError(var.nbins, math.sqrt(pow(histo.GetBinError(var.nbins),2) + pow(histo.GetBinError(var.nbins+1),2)))

    
    #    his.Scale(norm)
    his_int=histo.Integral()
    return histo, his_int

def addShapeUncHisto(sample, var, upweight,downweight, addFileName):
    fopen = ROOT.TFile.Open("rootFiles/"+ var.name + addFileName +  ".root", "update")
    
    his = ROOT.TH1F(sample.histoName,"" ,var.nbins,var.xmin,var.xmax)
    for item,pileup,ctau in zip(sample.df[var.name],sample.df['puWeight'], upweight):
            his.Fill(item,pileup*ctau)

    his.Scale(sample.norm)
    his.SetName(his.GetName() + '_ctau_UP')
    fopen.cd()
    his.Write()

    hisDown = ROOT.TH1F(sample.histoName,"" ,var.nbins,var.xmin,var.xmax)
    for item,pileup,ctau in zip(sample.df[var.name],sample.df['puWeight'], downweight):
            hisDown.Fill(item,pileup*ctau)

    hisDown.Scale(sample.norm)  #for mu and tau works
    hisDown.SetName(hisDown.GetName() + '_ctau_DOWN')
    fopen.cd()
    hisDown.Write()

    fopen.Close()

def inv_mass(pf, k= True):
    print("Adding the new invariant mass column...")
    if (k==True):
        k_p4= TLorentzVectorArray.from_ptetaphim(pf.kpt,pf.keta,pf.kphi,0.493677)                 
    else:   
        k_p4= TLorentzVectorArray.from_ptetaphim(pf.kpt,pf.keta,pf.kphi,pf.kmass)
    
    mu1_p4= TLorentzVectorArray.from_ptetaphim(pf.mu1pt,pf.mu1eta,pf.mu1phi,pf.mu1mass)
    mu2_p4= TLorentzVectorArray.from_ptetaphim(pf.mu2pt,pf.mu2eta,pf.mu2phi,pf.mu2mass)
    sum2 = mu1_p4 + mu2_p4
    pf['inv_mass']= (k_p4+sum2).mass 


#Gaus fit function; ps are the fit parameters
def gaus_fit(histo,var,plt_name="",path="pngPlots",ps=[None,None,None,None,None,None],pol="pol1"):
        func= ROOT.TF1("func","gaus(0) +"+ pol+"(3)")
        for i,p in enumerate(ps):
            if(p!=None):
                func.SetParameter(i,p)

        c = ROOT.TCanvas("","",700, 700) 
        histo.SetTitle( '; Events / bin;' + var.xlabel+var.unit+';Counts')
        histo.SetMarkerStyle(20)
        histo.SetMarkerSize(0.9)


        for i in range(1,histo.GetNbinsX()+1):
            histo.SetBinError(i, math.sqrt(histo.GetBinContent(i)))
        fit_result=histo.Fit(func,"S")
        fit_gaus=ROOT.TF1("gaus","gaus",var.xmin,var.xmax)
        fit_gaus.SetParameter(0,fit_result.Parameter(0))
        fit_gaus.SetParameter(1,fit_result.Parameter(1))
        fit_gaus.SetParameter(2,fit_result.Parameter(2))

        #legend
        legend = ROOT.TLegend(0.65,0.70,0.98,0.86)
        legend.AddEntry(histo, "data", "lp")
        legend.AddEntry(fit_gaus, "Gauss fit", "l")
        legend.SetTextFont(43)
        legend.SetTextSize(15)
        legend.SetBorderSize(0)
        legend.SetFillColor(0)

        c.SetTicks(True)
        c.SetBottomMargin(2)
        c.SetLeftMargin(0.15)
        c.SetRightMargin(0.15)
        c.Draw()
        c.cd()

        histo.Draw("e")
        #fit_gaus.Draw("lsame")
        legend.Draw("SAME")

        histo.GetYaxis().SetTitle("Events / bin")
        histo.GetYaxis().SetLabelOffset(0.01)
        histo.GetYaxis().SetTitleOffset(2)
        histo.GetYaxis().SetLabelFont(43)
        histo.GetYaxis().SetLabelSize(15)
        histo.GetYaxis().SetTitleFont(43)
        histo.GetYaxis().SetTitleSize(18)
        histo.SetMaximum(histo.GetMaximum()*0.7)
        histo.GetXaxis().SetLabelOffset(0.01)
        histo.GetXaxis().SetTitleOffset(1.6)
        histo.GetXaxis().SetLabelFont(43)
        histo.GetXaxis().SetLabelSize(15)
        histo.GetXaxis().SetTitleFont(43)
        histo.GetXaxis().SetTitleSize(18)
        histo.GetXaxis().SetTitle(var.xlabel + " " + var.unit)
        c.RedrawAxis()
        CMS_lumi(c,4,1)

        c.SaveAs(path+"/"+var.name+"_gausFit_"+plt_name+".png")
        print("Saved '"+path+"/"+var.name+"_gausFit_"+plt_name+".png'")

        c.SaveAs(path+"/"+var.name+"_gausFit_"+plt_name+".pdf")
        print("Saved '"+path+"/"+var.name+"_gausFit_"+plt_name+".pdf'")

        c2 = ROOT.TCanvas()
        fit_gaus.Draw("l")
        c2.SaveAs("pngPlots/"+var.name+"_Onlygaus_"+plt_name+".png")
        print("Saved 'pngPlots/"+var.name+"_Onlygaus_"+plt_name+".png'")

        f=ROOT.TFile.Open("rootPlots/"+var.name+"_gausFit_"+plt_name+".root","RECREATE")
        histo.Write()
        fit_gaus.Write()
        f.Write()
        f.Close()
        print("Saved 'rootPlots/"+var.name+"_gausFit_"+plt_name+".root'")
        
        integral = fit_gaus.Integral(var.xmin,var.xmax)
        return integral/((var.xmax-var.xmin)/var.nbins)

#fake pileup branch for data (to speed the compuation)
def pileupFake(pf):
    print("Adding fake pile up branches to data ...")
    puWeight = np.ones(len(pf))
    puWeightUp = np.ones(len(pf))
    puWeightDown = np.ones(len(pf))
    
    pf['puWeight'] = puWeight
    pf['puWeightUp'] = puWeightUp
    pf['puWeightDown'] = puWeightDown
    return pf

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
def dr13(pf):
    mu1_p4 = TLorentzVectorArray.from_ptetaphim(pf.mu1pt,pf.mu1eta,pf.mu1phi,pf.mu1mass)
    mu_p4 = TLorentzVectorArray.from_ptetaphim(pf.kpt,pf.keta,pf.kphi,pf.kmass)
    pf.copy()
    pf['dr13'] = mu_p4.delta_r(mu1_p4)
    return pf
def dr23(pf):
    mu2_p4 = TLorentzVectorArray.from_ptetaphim(pf.mu2pt,pf.mu2eta,pf.mu2phi,pf.mu2mass)
    mu_p4 = TLorentzVectorArray.from_ptetaphim(pf.kpt,pf.keta,pf.kphi,pf.kmass)
    pf.copy()
    pf['dr13'] = mu_p4.delta_r(mu2_p4)
    return pf


def DR_jpsimu(pf):
    print("Adding DR between jpsi and mu branch...")
    mu1_p4 = TLorentzVectorArray.from_ptetaphim(pf.mu1pt,pf.mu1eta,pf.mu1phi,pf.mu1mass)
    mu2_p4 = TLorentzVectorArray.from_ptetaphim(pf.mu2pt,pf.mu2eta,pf.mu2phi,pf.mu2mass)

    jpsi_p4= mu1_p4 + mu2_p4 
    #    jpsi_p4 = TLorentzVectorArray.from_ptetaphim((pf.mu1pt+pf.mu2pt),(pf.mu1eta,+pf.mu2eta),(pf.mu1phi+pf.mu2phi),(pf.mu1mass+pf.mu2mass))
    mu_p4 = TLorentzVectorArray.from_ptetaphim(pf.kpt,pf.keta,pf.kphi,pf.kmass)
    #    print(jpsi_p4.delta_r(mu_p4))
    pf.copy()
    pf['DR_jpsimu'] = jpsi_p4.delta_r(mu_p4)
    return pf

def jpsi_branches(pf):
    print("Adding jpsi four momentum branches...")
    mu1_p4 = TLorentzVectorArray.from_ptetaphim(pf.mu1pt,pf.mu1eta,pf.mu1phi,pf.mu1mass)
    mu2_p4 = TLorentzVectorArray.from_ptetaphim(pf.mu2pt,pf.mu2eta,pf.mu2phi,pf.mu2mass)
    jpsi_p4= mu1_p4 + mu2_p4
    
    pf = pf.copy()

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
    pf.copy()
    pf['decay_time1'] = decay_time1
    pf['decay_time2'] = decay_time2
    return pf

#Do it!
def mcor(pf):
    print("Adding mcor branch...")
    

#Various third muon ID selections
def selMediumID(df,region):
    if region == 1:
        return df[(df.k_mediumID == 1)]
    if region == 0:
        return df[(df.k_mediumID == 0)]
def selSoftID(df,region):
    if region == 1:
        return df[(df.k_softID == 1)]
    if region == 0:
        return df[(df.k_softID == 0)]
def selTightID(df,region):
    if region == 1:
        return df[(df.k_tightID == 1)]
    if region == 0:
        return df[(df.k_tightID == 0)]


#preselection cuts
def preselection(df):
    return df[(df.mu1pt >3) & (df.mu2pt>3) & (df.kpt > 2.5) & (abs(df.mu1eta) < 2.5) & (abs(df.mu2eta) < 2.5) & (abs(df.keta) < 2.5) & (df.bvtx_svprob > 1e-7) & (abs(df.k_dxy) <0.2) & (abs(df.mu1_dxy)<0.2) & (abs(df.mu2_dxy) < 0.2)   & (df.bvtx_cos2D > 0.95)  & (df.Bmass < 6.3) & (df.k_isPF == 1) & (df.mu1_mediumID == 1) & (df.mu2_mediumID == 1) & (df.Bpt_reco > 15) & (abs(df.mu1_dz-df.mu2_dz)<0.4) & (abs(df.mu1_dz-df.k_dz)<0.4) & (abs(df.mu2_dz-df.k_dz)<0.4) ]


print("Flag = ", flag)

#Training MVA and save samples with BDT branches and other added branches
if(fit):
    print("Computing the BDT test arrays from scratch...")
    print("Opening files...")
    pf_file=[]
    pf=[]
    root_files = []

    path1 = '/pnfs/psi.ch/cms/trivcat/store/user/friti/dataframes_2020Dec09/'
    path2 = '/pnfs/psi.ch/cms/trivcat/store/user/friti/dataframes_2020Dec10/'

    
    root_files.append(path1 + 'BcToJpsiTauNu_ptmax_merged.root') 
    root_files.append(path1 + 'BcToJpsiMuNu_ptmax_merged.root')
    root_files.append(path2 + 'OniaX_ptmax_merged.root') #UL samples
    
    for fil in root_files:
        pf_file.append(read_root(fil, 'BTommm')) #only this channel needed
    
    #apply selection on ID 
    print("Selecting the Signal region...")
    pf.append(selMediumID(pf_file[0],1)) 
    pf.append(selMediumID(pf_file[1],1))
    pf.append(selMediumID(pf_file[2],1))



    #preselection
    for i in range(len(pf)):
        #        pf[i] = jpsi_branches(pf[i])
        pf[i] = preselection(pf[i])
        pf[i] = dr13(pf[i])
        pf[i] = dr23(pf[i])

    sig_all = pf[0].copy()
    mu_all = pf[1].copy()
    bkg_all = pf[2].copy()
    
    sig_all['target'] = np.ones(sig_all.shape[0]).astype(np.int)
    bkg_all['target'] = np.zeros(bkg_all.shape[0]).astype(np.int)
    mu_all['target'] = np.full((mu_all.shape[0]),2)


    #adapt all at the minimum number of events 
    len_min=min(len(sig_all),len(bkg_all),len(mu_all))
    sig=sig_all[:len_min]
    bkg=bkg_all[:len_min]
    mu=mu_all[:len_min]

    print("Splitting into train, validation and test...")
    data = pd.concat([bkg, sig, mu])
    data['id'] = np.arange(len(data))
    train, test = train_test_split(data, test_size=0.4, random_state=123)
    split_value= int(len_min/3 *2)
    X_train, X_valid = train[:split_value][branches], train[split_value:][branches]
    y_train, y_valid = train[:split_value]['target'], train[split_value:]['target']
    
    ##################################################################################################
    #fit
    
    print("Fitting BDT...")
    clf=BDT_fit(X_train, y_train,X_valid,y_valid)

    # plot feature importance    
    if feature_importance:
        fig=plt.figure()
        feat_imp = pd.Series(clf.get_booster().get_fscore()).sort_values(ascending=False)
        feat_imp.plot(kind='bar', title='Feature Importances')
        plt.ylabel('Feature Importance Score')
        plt.tight_layout()
        fig.savefig("features.png")
        print("Saved figure of features importance score as features.png")

    print("Testing on BDT...")
    p_test = clf.predict_proba(test[branches])  
    test.index = [i for i in range(len(test))]
    xgb_test_pred=pd.DataFrame(data=p_test,columns=['bkg','tau','mu'])
    xgb_test_pred.loc[:,'max']=xgb_test_pred.idxmax(axis=1)
    xgb_test_pred.loc[:,'label']=test['target'].to_numpy()
    xgb_test_pred= pd.concat([test,xgb_test_pred],axis=1,sort=False)
      
    #confusion matrix
    if conf_matr:
        confusion_matrix = pd.crosstab(xgb_test_pred['label'], xgb_test_pred['max'], rownames=['Actual'], colnames=['Predicted'])
        sns_plot=sns.heatmap(confusion_matrix, annot=True)
        fig = sns_plot.get_figure()
        fig.savefig('BDT_plts/conf_matrx.png') 
        print("Confusion matrix saved in BDT_plts/conf_matrx.png")

    print("Dividing test arrays into the original samples...")
    bkg_df = xgb_test_pred[(xgb_test_pred.label==0)]
    tau_df=xgb_test_pred[(xgb_test_pred.label==1)]
    mu_df=xgb_test_pred[(xgb_test_pred.label==2)]
    
    #Addition of parts f the samples eliminated before to do training and test with equal number of events
    if(len_min == len(bkg_all)):
        tmp_df_mu=pd.DataFrame(data=clf.predict_proba(mu_all[len_min:][branches]), columns=['bkg','tau','mu'])
        tmp_df_sig=pd.DataFrame(data=clf.predict_proba(sig_all[len_min:][branches]), columns=['bkg','tau','mu'])
        
        mu_all.index = [i for i in range(-len_min,len(mu_all[len_min:]))]
        tmp_df_mu = pd.concat([mu_all[len_min:],tmp_df_mu],axis=1,sort=False)
        mu_df = pd.concat([mu_df,tmp_df_mu])

        sig_all.index = [i for i in range(-len_min,len(sig_all[len_min:]))]
        tmp_df_sig = pd.concat([sig_all[len_min:],tmp_df_sig],axis=1,sort=False)
        tau_df = pd.concat([tau_df,tmp_df_sig])
    elif(len_min == len(sig_all)):
        tmp_df_mu=pd.DataFrame(data=clf.predict_proba(mu_all[len_min:][branches]), columns=['bkg','tau','mu'])
        tmp_df_bkg=pd.DataFrame(data=clf.predict_proba(bkg_all[len_min:][branches]), columns=['bkg','tau','mu'])

        mu_all.index = [i for i in range(-len_min,len(mu_all[len_min:]))]
        tmp_df_mu = pd.concat([mu_all[len_min:],tmp_df_mu],axis=1,sort=False)
        mu_df = pd.concat([mu_df,tmp_df_mu])

        bkg_all.index = [i for i in range(-len_min,len(bkg_all[len_min:]))]
        tmp_df_bkg = pd.concat([bkg_all[len_min:],tmp_df_bkg],axis=1,sort=False)
        bkg_df = pd.concat([bkg_df,tmp_df_bkg])
   
    elif(len_min == len(mu_all)):
        tmp_df_sig=pd.DataFrame(data=clf.predict_proba(sig_all[len_min:][branches]), columns=['bkg','tau','mu'])
        tmp_df_bkg=pd.DataFrame(data=clf.predict_proba(bkg_all[len_min:][branches]), columns=['bkg','tau','mu'])

        sig_all.index = [i for i in range(-len_min,len(sig_all[len_min:]))]
        tmp_df_sig = pd.concat([sig_all[len_min:],tmp_df_sig],axis=1,sort=False)
        tau_df = pd.concat([tau_df,tmp_df_sig])
   
        bkg_all.index = [i for i in range(-len_min,len(bkg_all[len_min:]))]
        tmp_df_bkg = pd.concat([bkg_all[len_min:],tmp_df_bkg],axis=1,sort=False)
        bkg_df = pd.concat([bkg_df,tmp_df_bkg])

    print("Appending rest of the samples...")

    mu_df.index = [i for i in range(0,len(mu_df))]
    tau_df.index = [i for i in range(0,len(tau_df))]
    bkg_df.index = [i for i in range(0,len(bkg_df))]

    #################################################################################################
    # Saving dataframes with all the necessary branches #
    #################################################################################################

    print("Saving complete hd5 files, with BDT results and new branches...")
    df_list = []
    df_list.append(read_root(path1 + 'BcToJpsiMuNu_ptmax_merged.root','BTommm' ))
    df_list.append( read_root(path1 + 'BcToJpsiTauNu_ptmax_merged.root', 'BTommm'))


    df_list.append( read_root(path2 + 'BcToXToJpsi_is_chic0_mu_merged.root', 'BTommm'))
    df_list.append( read_root(path2 + 'BcToXToJpsi_is_chic1_mu_merged.root', 'BTommm'))
    df_list.append( read_root(path2 + 'BcToXToJpsi_is_chic2_mu_merged.root', 'BTommm'))
    df_list.append( read_root(path2 + 'BcToXToJpsi_is_hc_mu_merged.root', 'BTommm'))
    #    df_list.append( read_root(path + 'BcToXToJpsi_is_jpsi_3pi_merged.root', 'BTommm'))
    df_list.append( read_root(path2 + 'BcToXToJpsi_is_jpsi_hc_merged.root', 'BTommm'))
    df_list.append( read_root(path2 + 'BcToXToJpsi_is_psi2s_mu_merged.root', 'BTommm'))
    #df_list.append( read_root(path + 'BcToXToJpsi_is_psi2s_tau_merged.root', 'BTommm'))
    df_list.append( read_root(path1 + 'data_ptmax_merged.root', 'BTommm'))

    df_list.append( read_root(path2 +'OniaX_ptmax_merged.root', 'BTommm')) 
    df_list.append( read_root(path2 +'OniaX_ptmax_merged.root', 'BTommm')) 

    #training with All Onia, but comb must be only comb
    df_list[-1] = df_list[-1][abs(df_list[-1].k_genpdgId)==13]

    #pre selection and new branches
    for j,sname in zip(range(len(df_list)),['mu','tau','chic0','chic1','chic2','hc_mu','jpsi_hc','psi2s_mu','data','oniaAll','comb']):
        print("Processing ",sname)
        df_list[j] = preselection(df_list[j])
        df_list[j].index = [ind for ind in range(len(df_list[j]))] # reorder the indexes

    # data does not have these branches -> to CORRECT (redo dataframes)
    #    df_list[-3] = jpsi_branches(df_list[-3])
    #    df_list[-3] = DR_jpsimu(df_list[-3])
    #    df_list[-3] = decaytime(df_list[-3])
    #    df_list[-3] = lifetime_weight(df_list[-3]) #fake lifetime weight
    df_list[-3] = pileupFake(df_list[-3]) #data has fake pu branch for data (only '1' values)

    for df,name in zip(df_list,['mu','tau','chic0','chic1','chic2','hc_mu','jpsi_hc','psi2s_mu','data','oniaAll','comb']):
        print("Adding new BDT branches...")

        df_pred = clf.predict_proba(df[branches])
        df['bkg'] = [row[0] for row in df_pred]
        df['tau'] = [row[1] for row in df_pred]
        df['mu'] = [row[2] for row in df_pred]

        df.to_root('dataframes_local/bdt/' + name + '_' + flag+ '.root', 'df' ,mode = 'w')
        print("Saved file 'dataframes_local/bdt/" + name + "_" + flag+ ".root")

    #To continue with the flow of the code
    '''
    mu = df_list[0]
    tau = df_list[1]
    data = df_list[2]
    oniaAll = df_list[3]
    onia = df_list[4]
    '''
else:
    #Loading the dataframes that already have the BDT branches
    print("Loading BDT dataframes...")
    

    sample_all = {}
    sample_names_wOA = ['mu','tau','data','chic0','chic1','chic2','hc_mu','jpsi_hc','psi2s_mu','comb','oniaAll']
    sample_names = ['mu','tau','chic0','chic1','chic2','hc_mu','jpsi_hc','psi2s_mu','comb','mis_id','data']

    for sname in sample_names_wOA:
        sample_all[sname] = read_root('dataframes_local/bdt/'+sname+'_' + flag+ '.root','df' )

    print("Loading BDT model...")
    clf = pickle.load(open("bdtModel/BDT_Model_"+ flag +".dat", "rb"))

sample_names = ['mu','tau','chic0','chic1','chic2','hc_mu','jpsi_hc','psi2s_mu','comb','mis_id','data']
sample_df={}
for sname in sample_names[:-2]:  #except mis-id bkg
    sample_df[sname] = sample_all[sname].copy()
    sample_df[sname] = selMediumID(sample_df[sname],1)
sample_df["data"] = sample_all["data"].copy()
sample_df["data"] = selMediumID(sample_df["data"],1)


if not os.path.exists("www"):
    os.makedirs("www")
pathDir = "www/" + flag + "/"
if not os.path.exists(pathDir+"FR_plots"):
            os.makedirs(pathDir+"FR_plots")

if not os.path.exists(pathDir+"BDT_output"):
            os.makedirs(pathDir+"BDT_output")

pathPassFail = "www/" + flag + "/"+flag+"_cut" + cutStr
if not os.path.exists(pathPassFail):
            # if there si already one it does not delete it
            os.makedirs(pathPassFail)

if not os.path.exists(pathDir+"postFit"):
    os.makedirs(pathDir+"postFit")

if not os.path.exists("rootFiles/"+flag+"_cut" + cutStr):
    # if there si already one it does not delete it
    os.makedirs("rootFiles/"+flag+"_cut" + cutStr)
        

# Compute the fake rate with ID method
if computeFR:
    print("")
    print("Computing fake rate with mediumId...")

    data_fr = sample_all["data"].copy()  #already with preselection
    inv_mass(data_fr)
    

    
    
    var=[]
    var.append(variable("inv_mass","m(3#mu) with mediumId=1","m(3#mu)","[GeV]", 25, 5., 5.45 ))
    var.append(variable("inv_mass","m(3#mu)","m(3#mu)", "[GeV]",30, 5, 5.45 ))
    sample_dic["data"].df = selMediumID(data_fr,1).copy()

    histo_medium1, t1_int = createHisto(sample_dic["data"],var[0],over = False)
    makeSinglePlot(histo_medium1,var[0],sample_dic["data"],"medium1_" + flag)
    medium1_gaus_int=gaus_fit(histo_medium1,var[0],"medium1_" + flag,pathDir+"FR_plots",[None,5.28,0.07,None,None,None],pol="pol2")

    sample_dic["data"].df = data_fr
    histo_mediumAll, tAll_int = createHisto(sample_dic["data"],var[1] ,over = False)
    makeSinglePlot(histo_mediumAll,var[1],sample_dic["data"],"mediumAll_" + flag)

    mediumall_gaus_int=gaus_fit(histo_mediumAll,var[1],"mediumAll_" + flag,pathDir+"FR_plots",[None,5.27,0.05,None,None,None])

    fr_mediumId = medium1_gaus_int/mediumall_gaus_int

    print("Fake rate computed with mediumId: %s"%fr_mediumId)
    if (saveonwww):
        print("Sinc with the website...")
        os.system("rsync -aP www/"+flag+"/FR_plots/ friti@lxplus.cern.ch:/afs/cern.ch/user/f/friti/eos/www/"+flag+"/FR_plots/.")


else:
    print("")
    print("Using fake rate = ",fr)
    print("MisId normalisation = ",misId_norm)

if prob_shapes_plot:
    print("Drawing the shapes of the BDt output...")
    
    var=[]
    var.append(variable("tau","prob_{#tau}","prob_{#tau}","", 10, 0., 1 ))
    var.append(variable("mu","prob_{#mu}","prob_{#mu}", "", 10, 0., 1 ))
    var.append(variable("bkg","prob_{comb}","prob_{comb}", "", 10, 0., 1 ))
    
    sample_dic['mu'].df = sample_all['mu']
    sample_dic['mu'].color = ROOT.kRed
    sample_dic['tau'].df = sample_all['tau']
    sample_dic['tau'].color = ROOT.kGreen
    sample_dic['comb'].df = sample_all['oniaAll']
    sample_dic['comb'].color = ROOT.kBlack
    #    sample_dic[''].color = ROOT.kBlack
    

    #fare gli histo che mi servono
    for v in var:
        histo = []
        for sname in ['mu','tau','comb']:
            histo.append(createHisto(sample_dic[sname],v)[0])
        makeComparison(histo, v, pathDir, addFile="" + flag)
        
    if (saveonwww):
        print("Sinc with the website...")
        os.system("rsync -aP www/"+flag+"/BDT_output/ friti@lxplus.cern.ch:/afs/cern.ch/user/f/friti/eos/www/"+flag+"/BDT_output/.")

    
norm_fact_onia = 184.65
print("Using normalization factor for combinatorial bkg = ",norm_fact_onia)

#produces plots, datacards and root files for the simultaneous Pass and Fail fit
if(fitPassFail):
    
        sample_0 = {}
        for sname in sample_all:
            sample_0[sname] = selMediumID(sample_all[sname],0).copy()

        #blind analysis
        random.seed(2)
        rand = random.randint(0, 10000)
        blind = rand/10000 *1.5 +0.5

        fromfit = 1

        #sistema
        sample_dic["mu"].norm = 1 * fromfit * 24433./38658.
        sample_dic["tau"].norm = 1 * fromfit * 290620./500805. * blind
        sample_dic["comb"].norm = norm_fact_onia * 0.34338 * 5.428866 *0.05 * 0.5
        sample_dic["data"].norm = 1
        sample_dic['psi2s_mu' ].norm = (0.336000000 + 0.177300000 + 0.032800000 + 0.001300000) * fromfit
        sample_dic['chic0' ].norm = 0.011600000 * fromfit
        sample_dic['chic1' ].norm = 0.344000000* fromfit
        sample_dic['chic2' ].norm = 0.195000000 * fromfit
        sample_dic['hc_mu'    ].norm = 0.01 * fromfit
        #sample_dic['psi2s_tau'].norm = (0.336000000 + 0.177300000 + 0.032800000 + 0.001300000) * fromfit
        #sample_dic['jpsi_3pi' ].norm = 1.
        sample_dic['jpsi_hc'  ].norm = 1. * fromfit
        
        variables=[]
        '''
        variables.append(variable("mu1pt","p_{T}^{#mu_1}","p_{T}^{#mu_1}","[GeV]",30,2,30))
        variables.append(variable("mu1eta","#eta_{#mu_1}","#eta_{#mu_1}","",30,-3,3))
        variables.append(variable("mu2pt","p_{T}^{#mu_2}","p_{T}^{#mu_2}","[GeV]",30,2,30))
        variables.append(variable("mu2eta","#eta_{#mu_2}","#eta_{#mu_2}","",30,-3,3))
        variables.append(variable("kpt","p_{T}^{#mu}","p_{T}^{#mu}","[GeV]",30,0,30))
        variables.append(variable("keta","#eta_{#mu}","#eta_{#mu}","",30,-3,3))
        
        #jpsi 
        variables.append(variable("jpsi_pt","p_T^{J/\psi}","p_T^{J/\psi}","",20,9,20))
        variables.append(variable("jpsi_eta","#eta^{J/\psi}","#eta^{J/\psi}","",10,-2.5,2.5))

        #Bc
        variables.append(variable("Bpt","p_{T}^{B_{vis}}","p_{T}^{B_{vis}}","[GeV]",30,9,50))
        variables.append(variable("Bpt_reco","p_{T}^{B_c}","p_{T}^{B_c}","[GeV]",30,15,50))
        variables.append(variable("Bmass","m(3#mu)","m(3#mu)","[GeV]",10,3,6.5))


        #dxy and dz
        variables.append(variable("mu1_dxy","d_{xy}^{\mu_1}","d_{xy}^{\mu_1}","[cms]",15,-0.2,0.2))
        variables.append(variable("mu1_dz","d_{z}^{\mu_1}","d_{z}^{\mu_1}","[cms]",40,-30,30))
        variables.append(variable("mu2_dxy","d_{xy}^{\mu_2}","d_{xy}^{\mu_2}","[cms]",15,-0.2,0.2))
        variables.append(variable("mu2_dz","d_{z}^{\mu_2}","d_{z}^{\mu_2}","[cms]",40,-30,30))
        variables.append(variable("k_dxy","d_{xy}^{\mu}","d_{xy}^{\mu}","[cms]",15,-0.2,0.2))
        variables.append(variable("k_dz","d_{z}^{\mu}","d_{z}^{\mu}","[cms]",40,-30,30))

        #VERTEX properties
        variables.append(variable("bvtx_cos2D","cos2D","cos2D","",30,0.99,1))
        variables.append(variable("bvtx_lxy_sig","#sigma_{L_{xy}}","#sigma_{L_{xy}}","[cm]",15,0,20))
        
        variables.append(variable("bvtx_svprob","vx_{prob}","vx_{prob}","",40,0,0.02))
        
        variables.append(variable("decay_time","#tau_{B_c}","#tau_{B_c}","[s]",100,0,10e-12))
        #        variables.append(variable("decay_time2","#tau2_{B_c}","#tau2_{B_c}","[s]",100,0,10e-12))



        #isolamenti per ogni mu e isolamenti per B
        variables.append(variable("b_iso03","B_{iso03}","B_{iso03}","[GeV]",15,0,4))
        variables.append(variable("b_iso04","B_{iso04}","B_{iso04}","[GeV]",15,0,4))
        variables.append(variable("l1_iso03","#mu_1^{iso03}","#mu_1^{iso03}","[GeV]",15,0,4))
        variables.append(variable("l1_iso04","#mu_1^{iso04}","#mu_1^{iso04}","[GeV]",15,0,4))
        variables.append(variable("l2_iso03","#mu_2^{iso03}","#mu_2^{iso03}","[GeV]",15,0,4))
        variables.append(variable("l2_iso04","#mu_2^{iso04}","#mu_2^{iso04}","[GeV]",15,0,4))
        variables.append(variable("k_iso03","#mu^{iso03}","#mu^{iso03}","[GeV]",15,0,4))
        variables.append(variable("k_iso04","#mu^{iso04}","#mu^{iso04}","[GeV]",15,0,4))

        '''
        variables.append(variable("pt_miss_vec","p_{T}^{miss-vec}","p_{T}^{miss-vec}","[GeV]",15,0,15))
        variables.append(variable("pt_miss_scal","p_{T}^{miss}","p_{T}^{miss}","[GeV]",15,0,15))

        variables.append(variable("Q_sq","Q^{2}","Q^{2}","[GeV^{2}]",15,3,11))

        variables.append(variable("m_miss_sq","m_{miss}^{2}","m_{miss}^{2}","[GeV^{2}]",10,0,8))
        variables.append(variable("E_mu_star","E_{#mu}^{*}","E_{#mu}^{*}","[GeV]",12,1,8))
        variables.append(variable("E_mu_canc","E_{#mu}^{#}","E_{#mu}^{#}","[GeV]",10,1,3))
        variables.append(variable("DR_jpsimu","DR_{J/#psi/#mu}","DR_{J/#psi/#mu}","",10,0,1))
        

        #muon ID per ogni muone
        
        variables.append(variable("tau","#tau_{prob}","#tau_{prob}","",10,0,1))
        variables.append(variable("mu","#mu_{prob}","#mu_{prob}","",10,0,1))

        variables.append(variable("bkg","bkg_{prob}","bkg_{prob}","",10,0,1))


        #plots path directory
        print("The plots are going to be saved in " + pathPassFail)
        
        
        for var in variables:
            #################################################################
            #FAIL region
            #################################################################
        
            
            histos0 = []
            for sname in sample_names[:-2]:  #except data and mis_id bkg
                sample_dic[sname].df = sample_0[sname][(sample_0[sname].bkg<cut)]
                histos0.append(createHisto(sample_dic[sname], var,  norm = sample_dic[sname].norm, over = overflow_underflow)[0])

            sample_dic["data"].df = sample_0["data"][(sample_0["data"].bkg<cut)]
            datah0 = createHisto(sample_dic["data"], var, norm = sample_dic["data"].norm, over = overflow_underflow)[0]
            
            #bkg mis-id in the fail region is data - everything else
            bkg_r = ROOT.TH1F("","", var.nbins, var.xmin, var.xmax)
            bkg_r = datah0.Clone()
            bkg_r.SetName("mis_id")
            sum1 = 0
            for i in range(len(histos0)):
                bkg_r.Add(histos0[i],-1)      
                sum1 += histos0[i].GetBinContent(12)  
            for j in range(bkg_r.GetNbinsX()):
                if(bkg_r.GetBinContent(j) < 0):
                    bkg_r.SetBinContent(j,0.0001)
            histos0.append(bkg_r.Clone())
            bkg_df = selMediumID(sample_all["data"],0).copy()  
            
            #################################################################
            #PASS region
            #################################################################
            histos = []
            for sname in sample_names[:-2]:  #except data and mis_id bkg
                sample_dic[sname].df = sample_df[sname][(sample_df[sname].bkg<cut)]
                histos.append(createHisto(sample_dic[sname], var,  norm = sample_dic[sname].norm, over = overflow_underflow)[0])

            bkg_r.Scale(misId_norm )
            histos.append(bkg_r) #using the subtracted bkg as mis id    

            sample_dic["data"].df = sample_df["data"][(sample_df["data"].bkg<cut)]
            datah = createHisto(sample_dic["data"], var,norm = sample_dic["data"].norm , over = overflow_underflow)[0]
      
            
            '''
            #factor obsolete
            if (over == True and fac == True):
                factor = (histos[0].Integral()+histos[0].GetBinContent(0)+histos[0].GetBinContent(var.nbins+1))/(histos[1].Integral()+histos[1].GetBinContent(0)+histos[1].GetBinContent(var.nbins+1))
            elif (over == False and fac == True):
                factor = histos[0].Integral() / histos[1].Integral()

            else:
                factor = 1
            '''
            makeStack(datah, histos, var, addFileName = "Pass_cut"+cutStr, rootFile = True, ratioPad = False,  path = pathPassFail, rootPath = "rootFiles/"+flag+"_cut" + cutStr)
            

            addShapeUncHisto(sample_dic["mu"], var, sample_dic["mu"].df["ctau_weight_up"], sample_dic["mu"].df["ctau_weight_down"], addFileName = "Pass_cut"+cutStr)
            addShapeUncHisto(sample_dic["tau"], var, sample_dic["tau"].df["ctau_weight_up"], sample_dic["tau"].df["ctau_weight_down"], addFileName = "Pass_cut"+cutStr)

            create_datacard_signal(datah.Integral(),histos, dat_name = var.name+"_Pass")
            
            makeStack(datah0, histos0, var, addFileName = "Fail_cut"+cutStr, rootFile = True, ratioPad = False,  path = pathPassFail, rootPath = "rootFiles/"+flag+"_cut" + cutStr)

            addShapeUncHisto(sample_dic["mu"], var, sample_dic["mu"].df["ctau_weight_up"], sample_dic["mu"].df["ctau_weight_down"], addFileName = "Fail_cut"+cutStr)
            addShapeUncHisto(sample_dic["tau"], var, sample_dic["tau"].df["ctau_weight_up"], sample_dic["tau"].df["ctau_weight_down"], addFileName = "Fail_cut"+cutStr)
            
            create_datacard_control(datah0.Integral(),histos0,dat_name = var.name+"_Fail")

            with open(pathPassFail+"/info.txt","w") as info:
                info.write(
'''
BDT branches:
{bdt_branches}

BDT cut:
{bdt_cut}

preselection:

FakeRate:
{FR}
'''.format(
    bdt_branches   = branches,
    bdt_cut        = cut,
    FR             = fr,
)
                )
                
            #Option to copy everything in lxplus and save on the website
        if (saveonwww):
            print("Sinc with the website...")
            os.system("rsync -aP "+pathPassFail+"/ friti@lxplus.cern.ch:/afs/cern.ch/user/f/friti/eos/"+pathPassFail+"/.")

if(saveonwww_all):
    print("Sinc with the website...")
    #    os.system("rsync -aP www/"+flag+"/ friti@lxplus.cern.ch:/afs/cern.ch/user/f/friti/eos/www/"+flag+"/.")
    os.system("rsync -aP www/ friti@lxplus.cern.ch:/afs/cern.ch/user/f/friti/eos/www/.")


