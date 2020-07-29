#Script for the BDT
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
#from create_datacard import create_datacard_pf
import pickle
from uproot_methods import TLorentzVectorArray
import root_numpy as rp
#from Onlydata_plot import plot_and_save
import ROOT
from variable import variable
from sample import sampleColl
from sample import sample
from cmsstyle import CMS_lumi
from makePlot import makeSinglePlot
from makePlot import makeComparison
from makePlot import makeStack

# fit options
fit = False
feature_importance = False
conf_matr = False


computeFR = False
prob_shapes_plot = False


compCombNorm = False
compSignNorm = False

statErrorsPlots = False
fitPassFail = True
normalPlots = False

compPlot = False
prefitPlot = False


fr = 0.6961
misId_norm=fr/(1-fr)
#fr_BDT = "NAN"

#does not pop up canvases
ROOT.gROOT.SetBatch()
ROOT.gStyle.SetOptStat(0)


#branches to which do the training and test of the bdt
branches=[
    'kpt',
    'Bpt',
    'Bmass',
    'Bcos2D',
    'Blxy_sig',
    'Bmll_raw',
    'Bb_iso03',
    'Bb_iso04',
    'BQ_sq',
    'BDR',
    'Bpt_miss_vec',
#    'mu1_mediumID',
#    'mu2_mediumID',
    'Beta',
    'Bphi',
    'Bsvprob',
    'BE_mu_star',
    'BE_mu_canc',
    'Bpt_var',
    'Bm_miss_sq',
    'mu1pt',
    'mu2pt',
    'kphi',
    'mu1phi',
    'mu2phi',
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
        early_stopping_rounds = 100,
        eval_metric           = 'mlogloss',
        verbose               = True,
        #sample_weight         = train['weight'],
    )
    
    pickle.dump(clf, open("bdtModel/BDT_model21May_medium_onia_newmu.pickle.dat", "wb"))
    print("Model saved.")
    return clf

def inv_mass(pf):
    print("Creating the new inv_mass column...")
#    k_p4= TLorentzVectorArray.from_ptetaphim(pf.kpt,pf.keta,pf.kphi,pf.kmass)
    k_p4= TLorentzVectorArray.from_ptetaphim(pf.kpt,pf.keta,pf.kphi,0.493677)                     
    mu1_p4= TLorentzVectorArray.from_ptetaphim(pf.mu1pt,pf.mu1eta,pf.mu1phi,pf.mu1mass)
    mu2_p4= TLorentzVectorArray.from_ptetaphim(pf.mu2pt,pf.mu2eta,pf.mu2phi,pf.mu2mass)
    pf['inv_mass']= (k_p4+mu1_p4+mu2_p4).mass 
    

def createHisto(sample, var, norm=1):
    his = ROOT.TH1F(sample.histoName,"" ,var.nbins,var.xmin,var.xmax)
    for item in sample.df[var.name]:
            his.Fill(item)
    his.Scale(norm)
    his_int=his.Integral()
    return his, his_int

#passo collezione di p1 p2 p3, che di default è tutti None, ti 6 parametri. E posso scegliere tra pol1 e pol2. quando glieli passo metto a none quelli che non voglio passare
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

        c.SaveAs("png_plots/"+var.name+"_gausFit_"+plt_name+".png")
        print("Saved 'png_plots/"+var.name+"_gausFit_"+plt_name+".png'")

        c.SaveAs("png_plots/"+var.name+"_gausFit_"+plt_name+".pdf")
        print("Saved 'png_plots/"+var.name+"_gausFit_"+plt_name+".pdf'")

        c2 = ROOT.TCanvas()
        fit_gaus.Draw("l")
        c2.SaveAs("png_plots/"+var.name+"_Onlygaus_"+plt_name+".png")
        print("Saved 'png_plots/"+var.name+"_Onlygaus_"+plt_name+".png'")

        f=ROOT.TFile.Open("root_plots/"+var.name+"_gausFit_"+plt_name+".root","RECREATE")
        histo.Write()
        fit_gaus.Write()
        f.Write()
        f.Close()
        print("Saved 'root_plots/"+var.name+"_gausFit_"+plt_name+".root'")
        
        integral = fit_gaus.Integral(var.xmin,var.xmax)
        return integral/((var.xmax-var.xmin)/var.nbins)

    
    

#selections if we want to use the tightId as discriminant for the mis-id background
def selMuonID(df,region):
    if region == 1:
        return df[(df.k_mediumID == 1)]
    if region == 0:
        return df[(df.k_mediumID == 0)]

def seljpsi(df):
    return df[(df.mu1_mediumID == 1) & (df.mu2_mediumID == 1)]
#################################################################################################################
#PREPARE BDT INPUT
#non riapro tutte le volte i file, ma li salvo una volta finito
if(fit):
    print("Computing the BDT test arrays from scratch...")
    print("Opening hd5 files...")
    pf_file=[]
    pf=[]
    files=['hd5_files/OniaX.h5','hd5_files/BcToJpsiTauNu.h5','hd5_files/BcToJpsiMuNu_new.h5']

    for fil in files:
        pf_file.append(pd.read_hdf(fil,'pf'))


    #apply selection on ID
    print("Selection on mediumId to define the misId background in the sideband tightId=0...")
#    print(len(selMuonID(pf_file[1],1))))
#    print(len(seljpsi(selMuonID(pf_file[1],1))))
    pf.append(selMuonID(pf_file[1],1))
    pf.append(selMuonID(pf_file[2],1))
    pf.append(selMuonID(pf_file[0][abs(pf_file[0].k_genpdgId)==13],1)) #OniaX

    #jpsi muons selection
    for p in pf:
        p = seljpsi(p)

    sig_all = pf[0].copy()
    mu_all = pf[1].copy()
    bkg_all = pf[2].copy()
    
    sig_all['target'] = np.ones(sig_all.shape[0]).astype(np.int)
    bkg_all['target'] = np.zeros(bkg_all.shape[0]).astype(np.int)
    mu_all['target'] = np.full((mu_all.shape[0]),2)


    #adapt all at the minimum number of events (inr ealtà dovrei lavorare con i pesi nella bdt, per ora va bene così)
    len_min=min(len(sig_all),len(bkg_all),len(mu_all))
    sig=sig_all[:len_min]
    bkg=bkg_all[:len_min]
    mu=mu_all[:len_min]

    print("Splitting into train, validation and test...")
    data = pd.concat([bkg, sig, mu])
    data['id'] = np.arange(len(data))
    train, test = train_test_split(data, test_size=0.4, random_state=123)
    split_value= 3500 #?
    X_train, X_valid = train[:split_value][branches], train[split_value:][branches]
    y_train, y_valid = train[:split_value]['target'], train[split_value:]['target']
    
    ##############################################################################################################
    #fit
    
    print("Fitting BDT...")
    clf=BDT_fit(X_train, y_train,X_valid,y_valid)
    #plot sulle features importance
        
    
    
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
    tmp_df_bkg=pd.DataFrame(data=clf.predict_proba(bkg_all[len_min:][branches]), columns=['bkg','tau','mu'])
    tmp_df_sig=pd.DataFrame(data=clf.predict_proba(sig_all[len_min:][branches]), columns=['bkg','tau','mu'])

    #Mu non c'è perchè è lui il minimo!
    print("Appending rest of the samples...")
    bkg_all.index = [i for i in range(-len_min,len(bkg_all[len_min:]))]
    tmp_df_bkg = pd.concat([bkg_all[len_min:],tmp_df_bkg],axis=1,sort=False)
    bkg_df = pd.concat([bkg_df,tmp_df_bkg])

    sig_all.index = [i for i in range(-len_min,len(sig_all[len_min:]))]
    tmp_df_sig = pd.concat([sig_all[len_min:],tmp_df_sig],axis=1,sort=False)
    tau_df = pd.concat([tau_df,tmp_df_sig])

    mu_df.index = [i for i in range(0,len(mu_df))]
    tau_df.index = [i for i in range(0,len(tau_df))]
    bkg_df.index = [i for i in range(0,len(bkg_df))]


    mu_df.to_hdf('hd5_files/test/mu_test_onia_newmu.h5','df' ,mode = 'w')
    tau_df.to_hdf('hd5_files/test/tau_test_onia_newmu.h5','df' , mode = 'w')
    bkg_df.to_hdf('hd5_files/test/bkg_test_onia_newmu.h5', 'df' ,mode = 'w')
    print("Saved file 'hd5_files/test/mu_test_onia_newmu.h5'")
    print("Saved file 'hd5_files/test/tau_test_onia_newmu.h5'")
    print("Saved file 'hd5_files/test/bkg_test_onia_newmu.h5'")

    print("Saving complete hd5 files, with BDT results")
    df_list = []
    df_list.append(pd.read_hdf('hd5_files/BcToJpsiMuNu_new.h5','pf' ))
    df_list.append( pd.read_hdf('hd5_files/BcToJpsiTauNu.h5', 'pf'))
    df_list.append( pd.read_hdf('hd5_files/data.h5', 'pf'))
    df_list.append( pd.read_hdf('hd5_files/OniaX.h5', 'pf'))

    df_list[-1] = df_list[-1][abs(df_list[-1].k_genpdgId)==13]

    for i in range(len(df_list)):
        df_list[i] = seljpsi(df_list[i])

    for df,name in zip(df_list,['mu','tau','data','comb']):
        print("Adding new BDT branches...")
        df_pred = clf.predict_proba(df[branches])
        df['bkg'] = [row[0] for row in df_pred]
        df['tau'] = [row[1] for row in df_pred]
        df['mu'] = [row[2] for row in df_pred]

        df.to_hdf('hd5_files/bdt/' + name + '_onia_newmu_bdt.h5','df' ,mode = 'w')
        print("Saved file 'hd5_files/bdt/" + name + "_onia_newmu_bdt.h5'")

    #if you want to continue wit the flow of the code
    mu = df_list[0]
    tau = df_list[1]
    data = df_list[2]
    onia = df_list[3]

else:
    print("Loading BDT dataframes...")
    # there is already a selection on the mediumID of the jpsi muons
    mu = pd.read_hdf('hd5_files/bdt/mu_onia_newmu_bdt.h5','df' )
    tau = pd.read_hdf('hd5_files/bdt/tau_onia_newmu_bdt.h5','df')
    data = pd.read_hdf('hd5_files/bdt/data_onia_newmu_bdt.h5','df' )
    onia = pd.read_hdf('hd5_files/bdt/comb_onia_newmu_bdt.h5','df' )

    print("Loading BDT model...")
    clf = pickle.load(open("bdtModel/BDT_model21May_medium_onia_newmu.pickle.dat", "rb"))
    mu_df = mu.copy()
    mu_df = selMuonID(mu_df,1)
    tau_df = tau.copy()
    tau_df = selMuonID(tau_df,1)
    onia_df = onia.copy()
    onia_df = selMuonID(onia_df,1)
    data_df = data.copy()
    data_df = selMuonID(data_df,1)

# Compute the fake rate with tight ID method
if computeFR:
    print("")
    print("Computing fake rate...")
    data_fr = data.copy()  #already with ID selection on jpsi muons
    inv_mass(data_fr)
    
    var=[]
    var.append(variable("inv_mass","m(3#mu) with mediumId=1","m(3#mu)","[GeV]", 30, 5., 5.45 ))
    var.append(variable("inv_mass","m(3#mu)","m(3#mu)", "[GeV]",30, 5, 5.45 ))
    var.append(variable("inv_mass","m(3#mu) with mediumId=0","m(3#mu)", "[GeV]",60, 3., 8. ))
    var.append(variable("inv_mass","m(3#mu) with mediumId=0","m(3#mu)", "[GeV]",60, 5, 5.45 ))
    sampleColl[-1].df = selMuonID(data_fr,1).copy()

    histo_medium1, t1_int = createHisto(sampleColl[-1],var[0])
    makeSinglePlot(histo_medium1,var[0],sampleColl[-1],"medium1_onia_newmu")
    #medium1_gaus_int=gaus_fit(histo_medium1,var[0],"medium1_onia_newmu",[None,5.28,0.07,None,None,None],pol="pol2")   #without the selection on the jpsi muons
    medium1_gaus_int=gaus_fit(histo_medium1,var[0],"medium1_onia_newmu",[None,5.28,0.07,None,None,None],pol="pol1")

    sampleColl[-1].df = data_fr
    histo_mediumAll, tAll_int = createHisto(sampleColl[-1],var[1])
    makeSinglePlot(histo_mediumAll,var[1],sampleColl[-1],"mediumAll_onia_newmu")
    mediumall_gaus_int=gaus_fit(histo_mediumAll,var[1],"mediumAll_onia_newmu",[None,5.27,0.05,None,None,None])


    fr_mediumId = medium1_gaus_int/mediumall_gaus_int

    sampleColl[-1].df = selMuonID(data_fr,0).copy()
    histo_medium0, t0_int = createHisto(sampleColl[-1],var[2])
    makeSinglePlot(histo_medium0,var[2],sampleColl[-1],"medium0_onia_newmu")


    sampleColl[-1].df = selMuonID(data_fr,0).copy()
    histo_medium0g, t0g_int = createHisto(sampleColl[-1],var[3])
    makeSinglePlot(histo_medium0g,var[3],sampleColl[-1],"medium0g_onia_newmu")
    mediu0g_gaus_int=gaus_fit(histo_medium0g,var[3],"medium0g_onia_newmu",[None,5.27,0.05,None,None,None])


    print("Integral medium 1:",medium1_gaus_int," Integrale all: ", mediumall_gaus_int)
    print("Fake rate computed with mediumId: %s"%fr_mediumId)


else:
    print("")
    print("Using fake rate = ",fr)
    print("MisId norm = ",misId_norm)

if prob_shapes_plot:
    #dfs=[]
    #dfs.append(tau_df)
    #dfs.append(mu_df)
    #dfs.append(bkg_df)

    var=[]
    var.append(variable("tau","prob_{#tau}","prob_{#tau}","", 10, 0., 1 ))
    var.append(variable("mu","prob_{#mu}","prob_{#mu}", "", 10, 0., 1 ))
    var.append(variable("bkg","prob_{comb}","prob_{comb}", "", 10, 0., 1 ))
    
    sampleColl[0].df = mu_df
    sampleColl[1].df = tau_df
    sampleColl[2].df = onia_df
    sampleColl[2].color = ROOT.kBlack
    
    #fare gli histo che mi servono
    for v in var:
        histo = []
        for i in range(len(sampleColl)-2):
            histo.append(createHisto(sampleColl[i],v)[0])
        makeComparison(histo, v, addFile="mediumID_onia_newmu")

        
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
        makeSinglePlot(histo_misid, var, bkg_smpl, "medium_onia_newmu_bkgl01")
        misid_gaus_int=gaus_fit(histo_misid, var, "medium_onia_newmu_bkgl01",[None,5.28,0.05,None,None,None])

        
        #no BDT selection to compute this normalization! We are in the tightId=0 region!
#        histo_misid, misid_int = histo_plot(df,var,"bkgl01","e")
#        misid_gaus_int = gaus_fit(histo_misid,var,"bkgl01",[None,5.27,0.05,None,None,None])

        norm_fact_onia=misid_gaus_int/5.  #5 is the number of processes with K coming from a B+
 
        print("Onia normalization factor: ",norm_fact_onia)

else:
    #norm_fact_onia = 170  #with selection on muonID for jpsi
    norm_fact_onia = 380.41  #with selection on muonID for jpsi
    print("Using normalization factor for combinatorial bkg = ",norm_fact_onia)


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
    #        norm_fact_mu = 18.14777
        #        norm_fact_tau = 0.192634
        norm_fact_tau = 0.1917912973759333
        print("Using normalization factor for signal = ", norm_fact_tau)
        print("Using normalization factor for normalization = ", norm_fact_mu)

if(statErrorsPlots):
        
    bkg_df = selMuonID(data,0).copy()  #jpsi selection

    cut = 0.1
    sampleColl[0].df = mu_df[(mu_df.bkg < cut)]
    sampleColl[1].df = tau_df[(tau_df.bkg < cut)]
    sampleColl[3].df = bkg_df[(bkg_df.bkg < cut)]
    sampleColl[2].df = onia_df[(onia_df.bkg < cut)]
    sampleColl[-1].df = data_df[(data_df.bkg < cut)]
    
    sampleColl[0].norm = norm_fact_mu
    sampleColl[1].norm = norm_fact_tau
    sampleColl[3].norm = misId_norm
    sampleColl[2].norm = norm_fact_onia * 0.2
    
    
    sampleColl[-1].norm = 1
    
    variables=[]
    variables.append(variable("Bpt_miss_vec","p_{T}^{miss}","p_{T}^{miss}","[GeV]",15,0,20))
    
    for var in variables:
            for i in range(len(sampleColl)):
                his = createHisto(sampleColl[i], var, norm = sampleColl[i].norm)[0]
                makeSinglePlot(his, var, sampleColl[i], "_onia_newmu_statSample")
  
if(fitPassFail):
        data_0 = selMuonID(data,0).copy()
        tau_0 = selMuonID(tau,0).copy()
        mu_0 = selMuonID(mu,0).copy()
        onia_0 = selMuonID(onia,0).copy()

        cut = 0.1
        cutStr = "01"
        sampleColl[0].norm = norm_fact_mu #* 2.289
        sampleColl[1].norm = norm_fact_tau #  * 8.6
        #        sampleColl[3].norm = misId_norm
        sampleColl[2].norm = norm_fact_onia # *0.2 * 0.415
        sampleColl[-1].norm = 1

        
        variables=[]
        variables.append(variable("Bpt_miss_vec","p_{T}^{miss}","p_{T}^{miss}","[GeV]",15,0,20))
        variables.append(variable("tau","#tau_{prob}","#tau_{prob}","",10,0,1))
        variables.append(variable("mu","#mu_{prob}","#mu_{prob}","",10,0,1))
        variables.append(variable("bkg","bkg_{prob}","bkg_{prob}","",10,0,0.1))
        
        variables.append(variable("Bpt","p_{T}^{B}","p_{T}^{B}","[GeV]",15,8,50))
        variables.append(variable("kpt","p_{T}^{#mu}","p_{T}^{#mu}","[GeV]",15,2,30))
        variables.append(variable("BDR","#Delta R_{#mu_{1} #mu_{2}}","#Delta R_{#mu_{1} #mu_{2}}","",8,0,1))
        variables.append(variable("Bcos2D","cos2D","cos2D","",8,0.5,1))
        variables.append(variable("Blxy_sig","#sigma_{L_{xy}}","#sigma_{L_{xy}}","[cm]",15,0,20))
        variables.append(variable("Bmll_raw","m_{J/#psi}","m_{J/#psi}","",10,2.5,3.5))
        variables.append(variable("Bb_iso03","B_{iso03}","B_{iso03}","[GeV]",15,0,20))
        variables.append(variable("Bmass","m(3#mu)","m(3#mu)","[GeV]",10,3,7))
        
        
        variables.append(variable("BQ_sq","Q^{2}","Q^{2}","[GeV^{2}]",10,0,20))
        variables.append(variable("Bm_miss_sq","m_{miss}^{2}","m_{miss}^{2}","[GeV^{2}]",10,0,9))
        variables.append(variable("BE_mu_star","E_{#mu}^{*}","E_{#mu}^{*}","[GeV]",12,0,2.2))
        
        '''
        KinVarPlotFunction(data_sel,bkg_sel,tau_sel,mu_sel,onia_sel[(onia_sel.k_tightID==1)  & (onia_sel.k_genpdgId==13)],norm_fact_onia,norm_fact_mu,norm_fact_tau,"Blxy_sig","#sigma_{L_{xy}}","[cm]",rangex=(0,20),bins=15,root_files=True,fake_rate=fake_rate) 
        KinVarPlotFunction(data_sel,bkg_sel,tau_sel,mu_sel,onia_sel[(onia_sel.k_tightID==1)  & (onia_sel.k_genpdgId==13)],norm_fact_onia,norm_fact_mu,norm_fact_tau,"Bmll_raw","m_{J/#psi}","[GeV]",rangex=(2.95,3.25),bins=10,root_files=True,fake_rate=fake_rate) 
        KinVarPlotFunction(data_sel,bkg_sel,tau_sel,mu_sel,onia_sel[(onia_sel.k_tightID==1)  & (onia_sel.k_genpdgId==13)],norm_fact_onia,norm_fact_mu,norm_fact_tau,"Bb_iso03","B_{iso03}","[GeV]",rangex=(0,30),bins=15,root_files=True,fake_rate=fake_rate) 
        KinVarPlotFunction(data_sel,bkg_sel,tau_sel,mu_sel,onia_sel[(onia_sel.k_tightID==1)  & (onia_sel.k_genpdgId==13)],norm_fact_onia,norm_fact_mu,norm_fact_tau,"Bb_iso04","B_{iso04}","[GeV]",rangex=(0,30),bins=15,root_files=True,fake_rate=fake_rate) 
        KinVarPlotFunction(data_sel,bkg_sel,tau_sel,mu_sel,onia_sel[(onia_sel.k_tightID==1)  & (onia_sel.k_genpdgId==13)],norm_fact_onia,norm_fact_mu,norm_fact_tau,"BQ_sq","Q^{2}","[GeV^{2}]",rangex=(0,11),bins=10,root_files=True,fake_rate=fake_rate) 
        KinVarPlotFunction(data_sel,bkg_sel,tau_sel,mu_sel,onia_sel[(onia_sel.k_tightID==1)  & (onia_sel.k_genpdgId==13)],norm_fact_onia,norm_fact_mu,norm_fact_tau,"mu1_mediumID","mediumId_{#mu_{1}}","",rangex=(0,2),bins=2,root_files=True,fake_rate=fake_rate) 
        KinVarPlotFunction(data_sel,bkg_sel,tau_sel,mu_sel,onia_sel[(onia_sel.k_tightID==1)  & (onia_sel.k_genpdgId==13)],norm_fact_onia,norm_fact_mu,norm_fact_tau,"mu2_mediumID","mediumId_{#mu_{2}}","",rangex=(0,2),bins=2,root_files=True,fake_rate=fake_rate) 
        KinVarPlotFunction(data_sel,bkg_sel,tau_sel,mu_sel,onia_sel[(onia_sel.k_tightID==1)  & (onia_sel.k_genpdgId==13)],norm_fact_onia,norm_fact_mu,norm_fact_tau,"Beta","#eta_{B}","",rangex=(-3,3),bins=10,root_files=True,fake_rate=fake_rate) 
        KinVarPlotFunction(data_sel,bkg_sel,tau_sel,mu_sel,onia_sel[(onia_sel.k_tightID==1)  & (onia_sel.k_genpdgId==13)],norm_fact_onia,norm_fact_mu,norm_fact_tau,"Bphi","#phi_{B}","",rangex=(5,5),bins=10,root_files=True,fake_rate=fake_rate) 
        KinVarPlotFunction(data_sel,bkg_sel,tau_sel,mu_sel,onia_sel[(onia_sel.k_tightID==1)  & (onia_sel.k_genpdgId==13)],norm_fact_onia,norm_fact_mu,norm_fact_tau,"Bsvprob","vx_{prob}","",rangex=(0,1),bins=10,root_files=True,fake_rate=fake_rate) 
        KinVarPlotFunction(data_sel,bkg_sel,tau_sel,mu_sel,onia_sel[(onia_sel.k_tightID==1)  & (onia_sel.k_genpdgId==13)],norm_fact_onia,norm_fact_mu,norm_fact_tau,"BE_mu_star","E_{#mu}^{*}","[GeV]",rangex=(0,3),bins=15,root_files=True,fake_rate=fake_rate) 
        KinVarPlotFunction(data_sel,bkg_sel,tau_sel,mu_sel,onia_sel[(onia_sel.k_tightID==1)  & (onia_sel.k_genpdgId==13)],norm_fact_onia,norm_fact_mu,norm_fact_tau,"Bpt_var","p_{T}^{var}","[GeV]",rangex=(-10,30),bins=15,root_files=True,fake_rate=fake_rate) 
        KinVarPlotFunction(data_sel,bkg_sel,tau_sel,mu_sel,onia_sel[(onia_sel.k_tightID==1)  & (onia_sel.k_genpdgId==13)],norm_fact_onia,norm_fact_mu,norm_fact_tau,"Bm_miss_sq","m_{miss}^{2}","[GeV^{2}]",rangex=(-1,9),bins=10,root_files=True,fake_rate=fake_rate) 
        KinVarPlotFunction(data_sel,bkg_sel,tau_sel,mu_sel,onia_sel[(onia_sel.k_tightID==1)  & (onia_sel.k_genpdgId==13)],norm_fact_onia,norm_fact_mu,norm_fact_tau,"mu1pt","p_{T}^{#mu_{1}}","[GeV]",rangex=(3,30),bins=15,root_files=True,fake_rate=fake_rate) 
        KinVarPlotFunction(data_sel,bkg_sel,tau_sel,mu_sel,onia_sel[(onia_sel.k_tightID==1)  & (onia_sel.k_genpdgId==13)],norm_fact_onia,norm_fact_mu,norm_fact_tau,"mu2pt","p_{T}^{#mu_{2}}","[GeV]",rangex=(3,30),bins=15,root_files=True,fake_rate=fake_rate) 
        KinVarPlotFunction(data_sel,bkg_sel,tau_sel,mu_sel,onia_sel[(onia_sel.k_tightID==1)  & (onia_sel.k_genpdgId==13)],norm_fact_onia,norm_fact_mu,norm_fact_tau,"kphi","#phi_{#mu}","",rangex=(-4,4),bins=15,root_files=True,fake_rate=fake_rate) 
        KinVarPlotFunction(data_sel,bkg_sel,tau_sel,mu_sel,onia_sel[(onia_sel.k_tightID==1)  & (onia_sel.k_genpdgId==13)],norm_fact_onia,norm_fact_mu,norm_fact_tau,"mu1phi","#phi_{#mu_{1}}","",rangex=(-4,4),bins=15,root_files=True,fake_rate=fake_rate) 
        KinVarPlotFunction(data_sel,bkg_sel,tau_sel,mu_sel,onia_sel[(onia_sel.k_tightID==1)  & (onia_sel.k_genpdgId==13)],norm_fact_onia,norm_fact_mu,norm_fact_tau,"mu2phi","#phi_{#mu_{2}}","",rangex=(-4,4),bins=15,root_files=True,fake_rate=fake_rate) 
        KinVarPlotFunction(data_sel,bkg_sel,tau_sel,mu_sel,onia_sel[(onia_sel.k_tightID==1)  & (onia_sel.k_genpdgId==13)],norm_fact_onia,norm_fact_mu,norm_fact_tau,"mu1eta","#eta_{#mu_{1}}","",rangex=(-2.5,2.5),bins=15,root_files=True,fake_rate=fake_rate) 
        KinVarPlotFunction(data_sel,bkg_sel,tau_sel,mu_sel,onia_sel[(onia_sel.k_tightID==1)  & (onia_sel.k_genpdgId==13)],norm_fact_onia,norm_fact_mu,norm_fact_tau,"mu2eta","#eta_{#mu_{2}}","",rangex=(-2.5,2.5),bins=15,root_files=True,fake_rate=fake_rate) 
        KinVarPlotFunction(data_sel,bkg_sel,tau_sel,mu_sel,onia_sel[(onia_sel.k_tightID==1)  & (onia_sel.k_genpdgId==13)],norm_fact_onia,norm_fact_mu,norm_fact_tau,"keta","#eta_{#mu}","",rangex=(-2.5,2.5),bins=15,root_files=True,fake_rate=fake_rate) 
        KinVarPlotFunction(data_sel,bkg_sel,tau_sel,mu_sel,onia_sel[(onia_sel.k_tightID==1)  & (onia_sel.k_genpdgId==13)],norm_fact_onia,norm_fact_mu,norm_fact_tau,"tau","#tau_{prob}","",rangex=(0,1),bins=15,root_files=True,fake_rate=fake_rate) 
        KinVarPlotFunction(data_sel,bkg_sel,tau_sel,mu_sel,onia_sel[(onia_sel.k_tightID==1)  & (onia_sel.k_genpdgId==13)],norm_fact_onia,norm_fact_mu,norm_fact_tau,"mu","#mu_{prob}","",rangex=(0,1),bins=15,root_files=True,fake_rate=fake_rate) 

#        variables.append(variable("Bmass","mass","mass","[GeV]",15,3.5,9))
        '''
        for var in variables:

            #Fail region
            sampleColl[0].df = mu_0[(mu_0.bkg<cut)]
            sampleColl[1].df = tau_0[(tau_0.bkg<cut)]
            #        sampleColl[3].df = bkg_df[(bkg_df.bkg < cut)]
            sampleColl[2].df = onia_0[(onia_0.bkg<cut)]
            sampleColl[-1].df = data_0[(data_0.bkg<cut)]

            histos0 = []
            histos0_int = []
            for i in range(len(sampleColl)-1): #data not in this histo list
                if i == 3: #bkg mis_Id
                    continue
                histos0.append(createHisto(sampleColl[i], var, norm = sampleColl[i].norm)[0])

            datah0 = createHisto(sampleColl[-1], var, norm = sampleColl[-1].norm)[0]

            #bkg mis id in the fail region is data - everything else
            bkg_r = ROOT.TH1F("","", var.nbins, var.xmin, var.xmax)
            bkg_r = datah0.Clone()
            bkg_r.SetName("mis_id")

            for i in range(len(histos0)):
                bkg_r.Add(histos0[i],-1)      
            for j in range(bkg_r.GetNbinsX()):
                if(bkg_r.GetBinContent(j) < 0):
                    bkg_r.SetBinContent(j,0.1)
            histos0.append(bkg_r.Clone())
            
            #print("Integral Fail:",datah0.Integral(), histos0[0].Integral(),histos0[1].Integral(),histos0[2].Integral(),histos0[3].Integral())
           # makeStack(datah0, histos0, var, addFileName = "_medium_Fail0_cut"+cutStr, rootFile = True)
           #            create_datacard_5(var.name + "_medium_Fail0_cut"+cutStr+".root", datah0.Integral(),histos0[0].Integral(),histos0[1].Integral(),histos0[3].Integral(),histos0[2].Integral(),dat_name = "Fail")
           
            bkg_df = selMuonID(data,0).copy()  #jpsi selection
            

            #Pass region
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
                histos_int.append(createHisto(sampleColl[i], var, norm = sampleColl[i].norm)[1])
            bkg_r.Scale(misId_norm )#* 0.62 )
            histos.append(bkg_r) #using the subtracted bkg as mis id    
            datah = createHisto(sampleColl[-1], var, norm = sampleColl[-1].norm)[0]
            
            
            factor = histos[0].Integral()/histos[1].Integral()
            print("int",histos[0].Integral(),histos[1].Integral(),factor)
            print("before:",histos[1].Integral(),histos0[1].Integral())
            histos[1].Scale(factor)
            histos0[1].Scale(factor)
            print("after",histos[1].Integral(),histos0[1].Integral())
            
            print("Integrals pass : data",datah.Integral(),"mu:", histos[0].Integral(),"tau:",histos[1].Integral(),"comb:",histos[2].Integral(),"misid:",histos[3].Integral())
            makeStack(datah, histos, var, addFileName = "_medium_onia_newmu_Pass1_cut"+cutStr, rootFile = True, ratioPad = False)
            create_datacard_signal(datah.Integral(),histos[0].Integral(),histos[1].Integral(),histos[3].Integral(),histos[2].Integral(), dat_name = var.name+"_Pass")
           
            makeStack(datah0, histos0, var, addFileName = "_medium_onia_newmu_Fail0_cut"+cutStr, rootFile = True, ratioPad = False)
            create_datacard_control(datah0.Integral(),histos0[0].Integral(),histos0[1].Integral(),histos0[3].Integral(),histos0[2].Integral(),dat_name = var.name+"_Fail")

            #makeSinglePlot(histos0[3],var,sampleColl[2])
            if(normalPlots):
                for i in range(len(sampleColl)-1):
                    histos[i].Scale(1/histos[i].Integral())
                    makeSinglePlot(histos[i], var, sampleColl[i], "_norm_plts_")

            
    # plot normale stack con variabili
if(compPlot):
    variables=[]
    variables.append(variable("Bpt","p_{T}^{B}","p_{T}^{B}","[GeV]",15,8,50))
    variables.append(variable("kpt","p_{T}^{#mu}","p_{T}^{#mu}","[GeV]",15,2,30))
    variables.append(variable("BDR","#Delta R_{#mu_{1} #mu_{2}}","#Delta R_{#mu_{1} #mu_{2}}","",8,0,1))
    variables.append(variable("Bcos2D","cos2D","cos2D","",8,0.5,1))
    variables.append(variable("Blxy_sig","#sigma_{L_{xy}}","#sigma_{L_{xy}}","[cm]",15,0,20))
    variables.append(variable("Bmll_raw","m_{J/#psi}","m_{J/#psi}","",10,2.5,3.5))
    variables.append(variable("Bb_iso03","B_{iso03}","B_{iso03}","[GeV]",15,0,20))
    variables.append(variable("Bmass","m(3#mu)","m(3#mu)","[GeV]",15,3,10))
    
    variables.append(variable("BQ_sq","Q^{2}","Q^{2}","[GeV^{2}]",10,0,11))
    variables.append(variable("Bm_miss_sq","m_{miss}^{2}","m_{miss}^{2}","[GeV^{2}]",10,-1,9))
    variables.append(variable("BE_mu_star","E_{#mu}^{*}","E_{#mu}^{*}","[GeV]",15,0,3))

    sampleColl[0].df = mu_df
    sampleColl[1].df = tau_df
    sampleColl[2].df = onia_df
    sampleColl[2].color = ROOT.kBlack
    
    #fare gli histo che mi servono
    for v in variables:
        histo = []
        for i in range(len(sampleColl)-2):
            histo.append(createHisto(sampleColl[i],v)[0])
        makeComparison(histo, v, addFile="mediumID_onia_newmu")


if(prefitPlot):
        variables=[]
        variables.append(variable("Bpt_miss_vec","p_{T}^{miss}","p_{T}^{miss}","[GeV]",15,0,20))


        sample_coll[0].norm = norm_fact_mu
        sample_coll[1].norm = norm_fact_tau
        sample_coll[2].norm = misId_norm
        sample_coll[3].norm = norm_fact_onia
        sample_coll[-1].norm = 1

         
        for var in variables:
        
            cut = 0.1
            sample_coll[0].df = mu_df[(mu_df.bkg<cut)]
            sample_coll[1].df = tau_df[(tau_df.bkg<cut)]
            sample_coll[3].df = onia[(onia.bkg<cut)]
            sample_coll[2].df = bkg_df[(bkg_df.bkg<cut)]
            sample_coll[-1].df = data[(data.bkg<cut)]
            
            histos = []
            for i in range(len(sample_coll)-1): #data not in this histo list
                histos.append(createHisto(sample_coll[i], var, norm = sample_coll[i].norm)[0])
                datah = createHisto(sample_coll[-1], var, norm = sample_coll[-1].norm)[0]
            makeStack(datah, histos, var, addFileName = "_prefit_", rootFile = True)
        
#        makeComparison([histos[0],histos[2]],var, addFile = "mu_misid")
        #        KinVarPlotFunction(data,bkg_df,tau_df,mu_df,onia[(onia.k_tightID==1)  & (onia.k_genpdgId==13)],norm_fact_onia,norm_fact_mu,norm_fact_tau,"Bpt_miss_vec","p_{T}^{miss}","[GeV]",rangex=(0,20),bins=15,root_files=True,fake_rate=fake_rate)
#        KinVarPlotFunction(data_sel,bkg_sel,tau_sel,mu_sel,onia_sel[(onia_sel.k_tightID==1)  & (onia_sel.k_genpdgId==13)],norm_fact_onia,norm_fact_mu,norm_fact_tau,"Bpt_miss_vec","p_{T}^{miss}","[GeV]",rangex=(0,20),bins=15,root_files=True,fake_rate=fake_rate)
        ''' 
 
        '''
