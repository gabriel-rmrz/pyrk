#Plot of the variables post combine fit
import ROOT
import math

from variable import variable
from makePlot import makeSinglePlot
from makePlot import makeStack
import os

#from create_datacard import create_datacard_5
ROOT.gROOT.SetBatch()

flag = 'v2'
weight = False
saveonwww = True
cutStr = '04'
# file with input histos name from bash
path_combine = "/work/friti/combine/CMSSW_10_2_13/src/HiggsAnalysis/CombinedLimit/Q_sq_2020Dec01_cut04/"
f=ROOT.TFile(path_combine+"fitDiagnostics.root","r")

var = variable("Q_sq", "Q_sq", "Q^{2}", "[GeV]", 15, 3, 8)
#var = variable("pt_miss", "pt_miss", "p_{T}^{miss}", "[GeV]", 15, 0, 20)
#var = variable("tau", "tau", "prob_{#tau}", "", 10, 0, 0.7)
#var = variable("BE_mu_star","E_{#mu}^{*}","E_{#mu}^{*}","[GeV]",12,0.2,2.2)
#var = variable("m_miss_sq", "m_miss_sq", "m_{miss}^2", "[GeV^{2}]",10,0,9)

#var = variable("toy", "toy", "toy", "", 2, 0, 2)
#var = variable("Bmass","m(3#mu)","m(3#mu)","[GeV]",10,3,7)
#var = variable("Bmass","m(3#mu)","m(3#mu)","[GeV]",10,3,13)


if weight:
    scale = 1.333
else:
    scale = 1
 
sample_names = ['mu','tau','chic0','chic1','chic2','hc_mu','jpsi_hc','psi2s_mu','comb','mis_id']


pathDir = "www/" + flag + "/" +"postFit/" + var.name + "_cut" +cutStr + "/" 
if not os.path.exists(pathDir):
    os.makedirs(pathDir)

os.system("cp "+ path_combine + "/datacard.txt " + pathDir + ".")
print(path_combine +var.name+"Pass_cut"+cutStr+".root")
os.system("cp "+ path_combine +var.name+"Pass_cut"+cutStr+".root "+ pathDir + ".")
os.system("cp "+ path_combine +var.name+"Fail_cut"+cutStr+".root "+ pathDir + ".")
os.system("cp "+ path_combine +"fitDiagnostics.root "+ pathDir + ".")



for directoryName in ['ch1','ch2']:
    #for directoryName in ['tight0','tight1']:
    histos = []
    for sname in sample_names:
        print(sname)
        histo = f.Get("shapes_fit_s/" + directoryName + "/" + sname)
        histo_new = ROOT.TH1F(sname,sname, var.nbins, var.xmin, var.xmax)
        for i in range(1,histo.GetNbinsX()+1):
            histo_new.SetBinContent(i,histo.GetBinContent(i))
            histo_new.SetBinError(i,histo.GetBinError(i))
        #histo_new.Scale(scale)
        print("INTEGRAL" + sname,histo_new.Integral())
        histos.append(histo_new)
    
    #DATA
    his_d = f.Get("shapes_fit_s/" + directoryName + "/data")
    histo_d = ROOT.TH1F("data","data", var.nbins, var.xmin, var.xmax)
    for i in range(0,his_d.GetN()):
        histo_d.SetBinContent(i+1,his_d.GetPointY(i))
        #    histo_d.Scale(scale)
    print("data",histo_d.Integral())

    makeStack(histo_d,histos, var, path = pathDir ,fit=True, addFileName = 'post_fit_'+directoryName)

f.Close()



if (saveonwww):
    print("Sinc with the website...")
    os.system("rsync -aP www/"+flag+"/postFit/ friti@lxplus.cern.ch:/afs/cern.ch/user/f/friti/eos/www/"+flag+"/postFit/.")
