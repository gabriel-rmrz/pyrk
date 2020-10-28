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
from create_datacard import create_datacard_signal
from create_datacard import create_datacard_control
import pickle
from uproot_methods import TLorentzVectorArray
from uproot_methods import TVector3Array
from uproot_methods import TLorentzVector
from uproot_methods import TVector3
import root_numpy as rp
import ROOT
from variable import variable
from sample import sampleColl
from sample import sample
from cmsstyle import CMS_lumi
from makePlot import makeSinglePlot
from makePlot import makeComparison
from makePlot import makeStack
import os
from scipy.constants import c as speed_of_light

over = True   #overflow underflow option in plots
fac = False   #factor to the signal yield for the fit
flag = "UL_newctau"

# fit options
fit = False
feature_importance = False
conf_matr = False


computeFR = False 
prob_shapes_plot  = False
compCombNorm = False
compSignNorm = False

fitPassFail  = False
saveonwww = False

fr = 0.85373
misId_norm=fr/(1-fr)

#no pop up canvases
ROOT.gROOT.SetBatch()
ROOT.gStyle.SetOptStat(0)


#branches for training MVA
branches=[

    'Q_sq',

    'pt_miss_vec',

    'E_mu_star',

    'm_miss_sq',

    'Bmass',

    'kpt',

    'Bpt',
 
    'Bcos2D',
 
    'DR_mu1mu2',
    
    'Beta',

    'pt_var',
    
    'mu1pt',

    'mu2pt',

    'mu1eta',

    'mu2eta',

    'keta',    

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


def createHisto(sample, var, norm=1):
    his = ROOT.TH1F(sample.histoName,"" ,var.nbins,var.xmin,var.xmax)
    for item,pileup,ctau in zip(sample.df[var.name],sample.df['puWeight'],sample.df['ctau_weight_central']):
            his.Fill(item,pileup*ctau)

    his.Scale(norm)
    his_int=his.Integral()
    return his, his_int

#Gaus fit function; ps are the fit parameters
def gaus_fit(histo,var,plt_name="",ps=[None,None,None,None,None,None],pol="pol1"):
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

        c.SaveAs("pngPlots/"+var.name+"_gausFit_"+plt_name+".png")
        print("Saved 'pngPlots/"+var.name+"_gausFit_"+plt_name+".png'")

        c.SaveAs("pngPlots/"+var.name+"_gausFit_"+plt_name+".pdf")
        print("Saved 'pngPlots/"+var.name+"_gausFit_"+plt_name+".pdf'")

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

#fake pileup branch for data (to speed the computation)
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
    return df[(df.mu1pt >3) & (df.mu2pt>3) & (df.kpt > 2.5) & (df.jpsi_pt>6.5) & (df.Bmass < 6.275) & (df.Bpt > 15) & (df.k_isPF == 1) & (df.mu1_softID == 1) & (df.mu2_softID == 1) & (df.keta < 2.5) & (df.mu1eta < 2.5) & (df.mu2eta < 2.5) & (df.Bcos2D > 0.95)]


print("Flag = ", flag)

#Training MVA and save samples with BDT branches and other added branches
if(fit):
    print("Computing the BDT test arrays from scratch...")
    print("Opening hd5 files...")
    pf_file=[]
    pf=[]
    files=['hd5_files/OniaX.h5','hd5_files/BcToJpsiTauNu.h5','hd5_files/BcToJpsiMuNu.h5']

    for fil in files:
        pf_file.append(pd.read_hdf(fil,'pf'))
    

    #apply selection on ID 
    print("Selecting the Pass region...")
    pf.append(selSoftID(pf_file[1],1))
    pf.append(selSoftID(pf_file[2],1))
    pf.append(selSoftID(pf_file[0],1)) #OniaX

    #preselection
    for i in range(len(pf)):
        pf[i] = jpsi_branches(pf[i])
        pf[i] = preselection(pf[i])
        
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
    tmp_df_mu=pd.DataFrame(data=clf.predict_proba(mu_all[len_min:][branches]), columns=['bkg','tau','mu'])
    tmp_df_sig=pd.DataFrame(data=clf.predict_proba(sig_all[len_min:][branches]), columns=['bkg','tau','mu'])

    print("Appending rest of the samples...")
    mu_all.index = [i for i in range(-len_min,len(mu_all[len_min:]))]
    tmp_df_mu = pd.concat([mu_all[len_min:],tmp_df_mu],axis=1,sort=False)
    mu_df = pd.concat([mu_df,tmp_df_mu])

    sig_all.index = [i for i in range(-len_min,len(sig_all[len_min:]))]
    tmp_df_sig = pd.concat([sig_all[len_min:],tmp_df_sig],axis=1,sort=False)
    tau_df = pd.concat([tau_df,tmp_df_sig])

    mu_df.index = [i for i in range(0,len(mu_df))]
    tau_df.index = [i for i in range(0,len(tau_df))]
    bkg_df.index = [i for i in range(0,len(bkg_df))]

    #################################################################################################
    # Saving dataframes with all the necessary branches #
    #################################################################################################

    print("Saving complete hd5 files, with BDT results and new branches...")
    df_list = []
    df_list.append(pd.read_hdf('hd5_files/BcToJpsiMuNu.h5','pf' ))
    df_list.append( pd.read_hdf('hd5_files/BcToJpsiTauNu.h5', 'pf'))
    df_list.append( pd.read_hdf('hd5_files/data.h5', 'pf'))
    df_list.append( pd.read_hdf('hd5_files/OniaX.h5', 'pf'))
    df_list.append( pd.read_hdf('hd5_files/OniaX.h5', 'pf'))

    #training with All Onia, but comb must be only comb
    df_list[-1] = df_list[-1][abs(df_list[-1].k_genpdgId)==13]

    #pre selection and new branches
    for j in range(len(df_list)):
        df_list[j] = jpsi_branches(df_list[j])
        df_list[j] = preselection(df_list[j])
        df_list[j].index = [ind for ind in range(len(df_list[j]))] # reorder the indexes
        df_list[j] = DR_jpsimu(df_list[j])
        df_list[j] = decaytime(df_list[j])
        if (j<2): # only tau and mu need this weight
            df_list[j] = lifetime_weight(df_list[j], fake = False)  
        else: #others have a branch filled with 1 
            df_list[j] = lifetime_weight(df_list[j])
        
        df_list[2] = pileupFake(df_list[2]) #data has fake pu branch (only 1 values)

    for df,name in zip(df_list,['mu','tau','data','oniaAll','comb']):
        print("Adding new BDT branches...")

        df_pred = clf.predict_proba(df[branches])
        df['bkg'] = [row[0] for row in df_pred]
        df['tau'] = [row[1] for row in df_pred]
        df['mu'] = [row[2] for row in df_pred]

        df.to_hdf('hd5_files/bdt/' + name + '_' + flag+ '.h5','df' ,mode = 'w')
        print("Saved file 'hd5_files/bdt/" + name + "_" + flag+ ".h5")

    #To continue with the flow of the code
    mu = df_list[0]
    tau = df_list[1]
    data = df_list[2]
    oniaAll = df_list[3]
    onia = df_list[4]
    
else:
    #Loading the dataframes that already have the BDT branches
    print("Loading BDT dataframes...")
    mu = pd.read_hdf('hd5_files/bdt/mu_' + flag+ '.h5','df' )
    tau = pd.read_hdf('hd5_files/bdt/tau_' + flag+ '.h5','df')
    data = pd.read_hdf('hd5_files/bdt/data_' + flag+ '.h5','df' )
    onia = pd.read_hdf('hd5_files/bdt/comb_' + flag+ '.h5','df' )
    oniaAll = pd.read_hdf('hd5_files/bdt/oniaAll_' + flag+ '.h5','df' )
    
    print("Loading BDT model...")
    clf = pickle.load(open("bdtModel/BDT_Model_"+ flag +".dat", "rb"))

#_df are in the Pass region already
mu_df = mu.copy()
mu_df = selSoftID(mu_df,1)
tau_df = tau.copy()
tau_df = selSoftID(tau_df,1)
onia_df = onia.copy()
onia_df = selSoftID(onia_df,1)
data_df = data.copy()
data_df = selSoftID(data_df,1)

# Compute the fake rate with soft ID method
if computeFR:
    print("")
    print("Computing fake rate with softId...")
    data_fr = data.copy()  #already with ID selection on jpsi muons
    inv_mass(data_fr)
    
    var=[]
    var.append(variable("inv_mass","m(3#mu) with softId=1","m(3#mu)","[GeV]", 30, 5., 5.45 ))
    var.append(variable("inv_mass","m(3#mu)","m(3#mu)", "[GeV]",30, 5, 5.45 ))
    var.append(variable("inv_mass","m(3#mu) with softId=0","m(3#mu)", "[GeV]",60, 3., 8. ))
    var.append(variable("inv_mass","m(3#mu) with softId=0","m(3#mu)", "[GeV]",60, 5, 5.45 ))
    sampleColl[-1].df = selSoftID(data_fr,1).copy()

    histo_soft1, t1_int = createHisto(sampleColl[-1],var[0])
    makeSinglePlot(histo_soft1,var[0],sampleColl[-1],"soft1_" + flag)
    soft1_gaus_int=gaus_fit(histo_soft1,var[0],"soft1_" + flag,[None,5.28,0.07,None,None,None],pol="pol2")

    sampleColl[-1].df = data_fr
    histo_softAll, tAll_int = createHisto(sampleColl[-1],var[1] )
    makeSinglePlot(histo_softAll,var[1],sampleColl[-1],"softAll_" + flag)

    softall_gaus_int=gaus_fit(histo_softAll,var[1],"softAll_" + flag,[None,5.29,0.07,None,None,None])

    fr_softId = soft1_gaus_int/softall_gaus_int

    print("Fake rate computed with softId: %s"%fr_softId)


else:
    print("")
    print("Using fake rate = ",fr)
    print("MisId norm = ",misId_norm)

if prob_shapes_plot:

    var=[]
    var.append(variable("tau","prob_{#tau}","prob_{#tau}","", 10, 0., 1 ))
    var.append(variable("mu","prob_{#mu}","prob_{#mu}", "", 10, 0., 1 ))
    var.append(variable("bkg","prob_{comb}","prob_{comb}", "", 10, 0., 1 ))
    
    sampleColl[0].df = mu
    sampleColl[1].df = tau
    sampleColl[2].df = oniaAll
    sampleColl[2].color = ROOT.kBlack
    
    #fare gli histo che mi servono
    for v in var:
        histo = []
        for i in range(len(sampleColl)-2):
            histo.append(createHisto(sampleColl[i],v)[0])
        makeComparison(histo, v, addFile="" + flag)

#obsolete        
if(compCombNorm):

        #needed to compute the compCombNorm
        bkg_df = selMuonID(data,0).copy()  #jpsi selection

        print("Computing normalization for combine background in the region mediumId==0...")

        inv_mass(bkg_df)
        inv_mass(onia)

        cut=0.1
        
        var = variable("inv_mass","m(3#mu) with mediumId==0 & prob_{bkg}<0.1","m(3#mu)", "[GeV]", 50, 5., 5.5)
        #df = sample(bkg_df[(bkg_df.bkg<cut)],"bkg","bkgSample","misId bkg",ROOT.kBlue,misId_norm)
        bkg_smpl = sample("bkg","mis_id","misId",ROOT.kBlue)
        bkg_smpl.df = bkg_df[(bkg_df.bkg<cut)]

        histo_misid, misid_int = createHisto(bkg_smpl,var)
        makeSinglePlot(histo_misid, var, bkg_smpl, "medium_" + flag)
        misid_gaus_int=gaus_fit(histo_misid, var, "medium_" + flag,[None,5.28,0.05,None,None,None])

        
        #no BDT selection to compute this normalization! We are in the tightId=0 region!
#        histo_misid, misid_int = histo_plot(df,var,"bkgl01","e")
#        misid_gaus_int = gaus_fit(histo_misid,var,"bkgl01",[None,5.27,0.05,None,None,None])

        norm_fact_onia=misid_gaus_int/5.  #5 is the number of processes with K coming from a B+
 
        print("Onia normalization factor: ",norm_fact_onia)

else:
    #norm_fact_onia = 170 
    #    norm_fact_onia = 380.41  
    norm_fact_onia = 184.65
    print("Using normalization factor for combinatorial bkg = ",norm_fact_onia)

# obsolete
if(compSignNorm):
        print("Finding normalizations for signal and normalization processes...")
        #normalization factors mu e tau
        mu_totevents = 225666
        tau_totevents = 5952702
        eps_mu = len(mu)/mu_totevents
        eps_pi = 0.01672
        L_2018 = 59.7
        L_runII = 143
        pi_obs = 7629
        BR_mu_pi = 16
        norm_fact_mu = (L_2018 /(0.01672 * 59.7 + 0.01032 * 41.86 + 0.01049 * 39) * eps_mu * pi_obs * BR_mu_pi ) / len(mu_df)
        norm_fact_tau = 0.28 * mu_totevents/tau_totevents * norm_fact_mu

        print("norm mu ",norm_fact_mu)
        print("norm_tau",norm_fact_tau)
   
else:
        norm_fact_mu = 18.068365882579833
        norm_fact_tau = 0.1917912973759333
        print("Using normalization factor for signal = ", norm_fact_tau)
        print("Using normalization factor for normalization = ", norm_fact_mu)

#produces plots, datacards and root files for the simultaneous Pass and Fail fit
if(fitPassFail):
        data_0 = selSoftID(data,0).copy()
        tau_0 = selSoftID(tau,0).copy()
        mu_0 = selSoftID(mu,0).copy()
        onia_0 = selSoftID(onia,0).copy()


        cut = 0.4
        cutStr = "04"
        
        sampleColl[0].norm = norm_fact_mu *0.044319
        sampleColl[1].norm = norm_fact_tau *18.442375
        sampleColl[2].norm = norm_fact_onia * 0.34338 
        sampleColl[-1].norm = 1


        variables=[]
        variables.append(variable("mu1pt","p_{T}^{#mu_1}","p_{T}^{#mu_1}","[GeV]",30,2,30))
        variables.append(variable("mu1eta","#eta_{#mu_1}","#eta_{#mu_1}","",30,-3,3))
        variables.append(variable("mu2pt","p_{T}^{#mu_2}","p_{T}^{#mu_2}","[GeV]",30,2,30))
        variables.append(variable("mu2eta","#eta_{#mu_2}","#eta_{#mu_2}","",30,-3,3))
        variables.append(variable("kpt","p_{T}^{#mu}","p_{T}^{#mu}","[GeV]",30,0,30))
        variables.append(variable("keta","#eta_{#mu}","#eta_{#mu}","",30,-3,3))

        #jpsi 
        variables.append(variable("jpsi_pt","p_T^{J/\psi}","p_T^{J/\psi}","",20,9,20))
        variables.append(variable("jpsi_eta","#eta^{J/\psi}","#eta^{J/\psi}","",10,-3,3))

        #Bc
        variables.append(variable("Bpt","p_{T}^{B_{vis}}","p_{T}^{B_{vis}}","[GeV]",30,8,50))
        variables.append(variable("Bpt_reco","p_{T}^{B_c}","p_{T}^{B_c}","[GeV]",30,8,50))
        variables.append(variable("Bmass","m(3#mu)","m(3#mu)","[GeV]",10,3.6,6))


        #dxy and dz
        variables.append(variable("mu1_dxy","d_{xy}^{\mu_1}","d_{xy}^{\mu_1}","[cms]",15,-0.2,0.2))
        variables.append(variable("mu1_dz","d_{z}^{\mu_1}","d_{z}^{\mu_1}","[cms]",40,-30,30))
        variables.append(variable("mu2_dxy","d_{xy}^{\mu_2}","d_{xy}^{\mu_2}","[cms]",15,-0.2,0.2))
        variables.append(variable("mu2_dz","d_{z}^{\mu_2}","d_{z}^{\mu_2}","[cms]",40,-30,30))
        variables.append(variable("k_dxy","d_{xy}^{\mu}","d_{xy}^{\mu}","[cms]",15,-0.2,0.2))
        variables.append(variable("k_dz","d_{z}^{\mu}","d_{z}^{\mu}","[cms]",40,-30,30))

        #VERTEX properties
        variables.append(variable("Bcos2D","cos2D","cos2D","",40,0.95,1))
        variables.append(variable("Blxy_sig","#sigma_{L_{xy}}","#sigma_{L_{xy}}","[cm]",15,0,20))
        
        variables.append(variable("Bsvprob","vx_{prob}","vx_{prob}","",50,0,0.1))
        
        variables.append(variable("decay_time1","#tau1_{B_c}","#tau1_{B_c}","[s]",100,0,10e-12))
        variables.append(variable("decay_time2","#tau2_{B_c}","#tau2_{B_c}","[s]",100,0,10e-12))

        
        variables.append(variable("pt_miss_vec","p_{T}^{miss-vec}","p_{T}^{miss-vec}","[GeV]",15,0,20))
        variables.append(variable("pt_miss_scal","p_{T}^{miss}","p_{T}^{miss}","[GeV]",15,0,20))
        variables.append(variable("Q_sq","Q^{2}","Q^{2}","[GeV^{2}]",12,1,10))
        variables.append(variable("m_miss_sq","m_{miss}^{2}","m_{miss}^{2}","[GeV^{2}]",10,0,9))
        variables.append(variable("E_mu_star","E_{#mu}^{*}","E_{#mu}^{*}","[GeV]",12,0.2,2.2))
        variables.append(variable("E_mu_canc","E_{#mu}^{#}","E_{#mu}^{#}","[GeV]",20,0.2,6))
        variables.append(variable("DR_jpsimu","DR_{J/#psi/#mu}","DR_{J/#psi/#mu}","",10,0,1))


        #muon ID per ogni muone

        #isolamenti per ogni mu e isolamenti per B
        variables.append(variable("b_iso03","B_{iso03}","B_{iso03}","[GeV]",15,0,20))
        variables.append(variable("b_iso04","B_{iso04}","B_{iso04}","[GeV]",15,0,30))
        variables.append(variable("l1_iso03","#mu_1^{iso03}","#mu_1^{iso03}","[GeV]",15,0,20))
        variables.append(variable("l1_iso04","#mu_1^{iso04}","#mu_1^{iso04}","[GeV]",15,0,30))
        variables.append(variable("l2_iso03","#mu_2^{iso03}","#mu_2^{iso03}","[GeV]",15,0,20))
        variables.append(variable("l2_iso04","#mu_2^{iso04}","#mu_2^{iso04}","[GeV]",15,0,30))
        variables.append(variable("k_iso03","#mu^{iso03}","#mu^{iso03}","[GeV]",15,0,20))
        variables.append(variable("k_iso04","#mu^{iso04}","#mu^{iso04}","[GeV]",15,0,30))

        variables.append(variable("tau","#tau_{prob}","#tau_{prob}","",10,0,0.7))
        variables.append(variable("mu","#mu_{prob}","#mu_{prob}","",10,0,1))
        variables.append(variable("bkg","bkg_{prob}","bkg_{prob}","",10,0,1))

        #plots path directory
        if not os.path.exists("stackPlt/www"):
            os.makedirs("stackPlt/www")
        pathDir = "stackPlt/www/" + flag + "_cut" + cutStr
        print("The plots are going to be saved in " + pathDir)
        if not os.path.exists(pathDir):
            # if there si already one it does not delete it
            os.makedirs(pathDir)
        
        for var in variables:
            #################################################################
            #FAIL region
            #################################################################
            
            sampleColl[0].df = mu_0[(mu_0.bkg<cut)]
            sampleColl[1].df = tau_0[(tau_0.bkg<cut)]
            sampleColl[2].df = onia_0[(onia_0.bkg<cut)]
            sampleColl[-1].df = data_0[(data_0.bkg<cut)]

            histos0 = []
            for i in range(len(sampleColl)-1): #data not in this histo list
                if i == 3: #bkg mis_Id
                    continue
                histos0.append(createHisto(sampleColl[i], var,  norm = sampleColl[i].norm)[0])
            datah0 = createHisto(sampleColl[-1], var, norm = sampleColl[-1].norm)[0]

            #bkg mis-id in the fail region is data - everything else
            bkg_r = ROOT.TH1F("","", var.nbins, var.xmin, var.xmax)
            bkg_r = datah0.Clone()
            bkg_r.SetName("mis_id")

            for i in range(len(histos0)):
                bkg_r.Add(histos0[i],-1)      
            for j in range(bkg_r.GetNbinsX()):
                if(bkg_r.GetBinContent(j) < 0):
                    bkg_r.SetBinContent(j,0.1)
            histos0.append(bkg_r.Clone())
            
            bkg_df = selSoftID(data,0).copy()  

            #################################################################
            #PASS region
            #################################################################
            sampleColl[0].df = mu_df[(mu_df.bkg<cut)]
            sampleColl[1].df = tau_df[(tau_df.bkg<cut)]
            sampleColl[2].df = onia_df[(onia_df.bkg<cut)]
            sampleColl[3].df = bkg_df[(bkg_df.bkg<cut)]
            sampleColl[-1].df = data_df[(data_df.bkg<cut)]
        
            histos = []
            histos_int = []
            for i in range(len(sampleColl)-1): #data not in this histo list
                if i == 3:
                    continue
                histos.append(createHisto(sampleColl[i], var, norm = sampleColl[i].norm)[0])
            bkg_r.Scale(misId_norm )
            histos.append(bkg_r) #using the subtracted bkg as mis id    
            datah = createHisto(sampleColl[-1], var,norm = sampleColl[-1].norm )[0]
        
            #factor obsolete
            if (over == True and fac == True):
                factor = (histos[0].Integral()+histos[0].GetBinContent(0)+histos[0].GetBinContent(var.nbins+1))/(histos[1].Integral()+histos[1].GetBinContent(0)+histos[1].GetBinContent(var.nbins+1))
            elif (over == False and fac == True):
                factor = histos[0].Integral() / histos[1].Integral()

            else:
                factor = 1

            makeStack(datah, histos, var, addFileName = "Pass_cut"+cutStr, rootFile = True, ratioPad = False, over = over ,fact = factor, path = pathDir)
            
            create_datacard_signal(datah.Integral(),histos[0].Integral(),histos[1].Integral(),histos[3].Integral(),histos[2].Integral(), dat_name = var.name+"_Pass")
            
            makeStack(datah0, histos0, var, addFileName = "Fail_cut"+cutStr, rootFile = True, ratioPad = False, over=over, fact = factor, path = pathDir)
            
            create_datacard_control(datah0.Integral(),histos0[0].Integral(),histos0[1].Integral(),histos0[3].Integral(),histos0[2].Integral(),dat_name = var.name+"_Fail")
                
        #Option to copy everything in lxplus and save on the website
        if (saveonwww):
            print("Sinc with the website...")
            os.system("rsync -aP stackPlt/www/ friti@lxplus.cern.ch:/afs/cern.ch/user/f/friti/eos/www/.")
