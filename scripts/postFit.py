#Plot of the variables post combine fit
import ROOT
import math

from variable import variable
from makePlot import makeSinglePlot
from makePlot import makeStack
from sample import sample_coll
#from create_datacard import create_datacard_5
ROOT.gROOT.SetBatch()

weight = True
asimov = False
PassFail = True
# file with input histos name from bash
#f=ROOT.TFile("/work/friti/new/CMSSW_10_2_15/src/HiggsAnalysis/PassFail/prova/fitDiagnostics.root","r")
f=ROOT.TFile("/work/friti/new/CMSSW_10_2_15/src/HiggsAnalysis/Q_sq_2020Nov02_cut04/fitDiagnostics.root","r")

var = variable("Q_sq", "Q_sq", "Q^{2}", "[GeV]", 12, 1, 10)
#var = variable("pt_miss", "pt_miss", "p_{T}^{miss}", "[GeV]", 15, 0, 20)
#var = variable("tau", "tau", "prob_{#tau}", "", 10, 0, 0.7)
#var = variable("BE_mu_star","E_{#mu}^{*}","E_{#mu}^{*}","[GeV]",12,0.2,2.2)
#var = variable("m_miss_sq", "m_miss_sq", "m_{miss}^2", "[GeV^{2}]",10,0,9)

#var = variable("toy", "toy", "toy", "", 2, 0, 2)
#var = variable("Bmass","m(3#mu)","m(3#mu)","[GeV]",10,3,7)
#var = variable("Bmass","m(3#mu)","m(3#mu)","[GeV]",10,3,13)
if not PassFail:
    his_m = f.Get("shapes_fit_s/control/mc_mu")
    histo_m = ROOT.TH1F("mc_mu","mc_mu", var.nbins, var.xmin, var.xmax)
    for i in range(1,his_m.GetNbinsX()+1):
        histo_m.SetBinContent(i,his_m.GetBinContent(i))
        histo_m.SetBinError(i,his_m.GetBinError(i))
    print("mu",histo_m.Integral())
    makeSinglePlot(histo_m, var, sample_coll[0])


    his_t = f.Get("shapes_fit_s/control/mc_tau")
    histo_t = ROOT.TH1F("mc_tau","mc_tau", var.nbins, var.xmin, var.xmax)
    for i in range(1,his_t.GetNbinsX()+1):
        histo_t.SetBinContent(i,his_t.GetBinContent(i))
        histo_t.SetBinError(i,his_t.GetBinError(i))

    print("tau",histo_t.Integral())
    makeSinglePlot(histo_t, var, sample_coll[1])

    his_x = f.Get("shapes_fit_s/control/mis_id")
    histo_x = ROOT.TH1F("mis_id","mis_id", var.nbins,var.xmin, var.xmax)
    for i in range(1,his_x.GetNbinsX()+1):
        histo_x.SetBinContent(i,his_x.GetBinContent(i))
        histo_x.SetBinError(i,his_x.GetBinError(i))

    print("misid",histo_x.Integral())
    makeSinglePlot(histo_x, var, sample_coll[2])


    his_d = f.Get("shapes_fit_s/control/data")
    histo_d = ROOT.TH1F("data","data_obs", var.nbins, var.xmin, var.xmax)
    for i in range(0,his_d.GetN()):
        histo_d.SetBinContent(i+1,his_d.GetPointY(i))

    print("data",histo_d.Integral())
    makeSinglePlot(histo_d, var, sample_coll[-1])


    his_c=f.Get("shapes_fit_s/control/mc_comb")
    histo_c= ROOT.TH1F("mc_comb","mc_comb", var.nbins, var.xmin, var.xmax)
    for i in range(1,his_c.GetNbinsX()+1):
        histo_c.SetBinContent(i,his_c.GetBinContent(i))
        histo_c.SetBinError(i,his_c.GetBinError(i))

    makeSinglePlot(histo_c, var, sample_coll[3])


    print("comb",histo_c.Integral())
    makeStack(histo_d,[histo_m, histo_t, histo_x, histo_c], var, path = 'fitPlt',fit=False)



    if(asimov):
        fw=ROOT.TFile.Open("root_files/asimov_file_"+var.name+".root","RECREATE")
        histo_m.SetName("mc_mu")
        histo_t.SetName("mc_tau")
        histo_x.SetName("mis_id")
        histo_c.SetName("mc_comb")
        histo_data = ROOT.TH1F("data_obs","data_obs", var.nbins, var.xmin, var.xmax)
        for i in range(1,histo_data.GetNbinsX()+1):
            histo_data.SetBinContent(i,histo_m.GetBinContent(i) + histo_t.GetBinContent(i) + histo_x.GetBinContent(i) + histo_c.GetBinContent(i))
            histo_data.SetBinError(i,math.sqrt(histo_data.GetBinContent(i)))

        fw.cd()
        histo_data.Write()
        histo_c.Write()
        histo_m.Write()
        histo_t.Write()
        histo_x.Write()
        fw.Write()
        
    
        print('Saved file '+'root_files/asimov_file_'+var.name+'.root')
    
        #create datacard
        create_datacard_5('asimov_file_'+var.name+'.root',histo_data.Integral(),histo_m.Integral(),histo_t.Integral(),histo_x.Integral(),histo_c.Integral(),'datacard_asimov_'+var.name)
        print("Saved "+'datacard_asimov_'+var.name+".txt")
        fw.Close()

if weight:
    scale = 1.333
else:
    scale = 1
 
if PassFail:
    for directoryName in ['ch1','ch2']:
    #for directoryName in ['tight0','tight1']:
        his_m = f.Get("shapes_fit_s/" + directoryName + "/mc_mu")
        histo_m = ROOT.TH1F("mc_mu","mc_mu", var.nbins, var.xmin, var.xmax)
        for i in range(1,his_m.GetNbinsX()+1):
            histo_m.SetBinContent(i,his_m.GetBinContent(i))
            histo_m.SetBinError(i,his_m.GetBinError(i))
        histo_m.Scale(scale)
        print("INTEGRAL mu",histo_m.Integral())
        #        makeSinglePlot(histo_m, var, sample_coll[0], saveFile = False)


        his_t = f.Get("shapes_fit_s/" + directoryName + "/mc_tau")
        histo_t = ROOT.TH1F("mc_tau","mc_tau", var.nbins, var.xmin, var.xmax)
        for i in range(1,his_t.GetNbinsX()+1):
            histo_t.SetBinContent(i,his_t.GetBinContent(i))
            histo_t.SetBinError(i,his_t.GetBinError(i))
        histo_t.Scale(scale)
        print("INTEGRAL tau",histo_t.Integral())

        #makeSinglePlot(histo_t, var, sample_coll[1], saveFile = False)

        his_x = f.Get("shapes_fit_s/" + directoryName + "/mis_id")
        histo_x = ROOT.TH1F("mis_id","mis_id", var.nbins,var.xmin, var.xmax)
        for i in range(1,his_x.GetNbinsX()+1):
            histo_x.SetBinContent(i,his_x.GetBinContent(i))
            histo_x.SetBinError(i,his_x.GetBinError(i))
            print(i,histo_x.GetBinContent(i))
        histo_x.Scale(scale)
        
        #makeSinglePlot(histo_x, var, sample_coll[2], saveFile = False)
        print("misid int:",histo_x.Integral())

        his_d = f.Get("shapes_fit_s/" + directoryName + "/data")
        histo_d = ROOT.TH1F("data_obs","data_obs", var.nbins, var.xmin, var.xmax)
        for i in range(0,his_d.GetN()):
            histo_d.SetBinContent(i+1,his_d.GetPointY(i))
        histo_d.Scale(scale)

        print("data",histo_d.Integral())
        #makeSinglePlot(histo_d, var, sample_coll[-1], saveFile = True)


        his_c=f.Get("shapes_fit_s/" + directoryName + "/mc_comb")
        histo_c= ROOT.TH1F("mc_comb","mc_comb", var.nbins, var.xmin, var.xmax)
        for i in range(1,his_c.GetNbinsX()+1):
            histo_c.SetBinContent(i,his_c.GetBinContent(i))
            histo_c.SetBinError(i,his_c.GetBinError(i))
        histo_c.Scale(scale)
        print("comb",histo_c.Integral())

        #makeSinglePlot(histo_c, var, sample_coll[3], saveFile = False)



        makeStack(histo_d,[histo_m, histo_t, histo_x, histo_c], var, path = 'fitPlt',fit=True, addFileName = directoryName)
        '''
        his_m = f.Get("shapes_fit_s/" + directoryName + "/sig")
        histo_m = ROOT.TH1F("mc_mu","mc_mu", var.nbins, var.xmin, var.xmax)
        for i in range(1,his_m.GetNbinsX()+1):
            histo_m.SetBinContent(i,his_m.GetBinContent(i))
            histo_m.SetBinError(i,his_m.GetBinError(i))

        makeSinglePlot(histo_m, var, sample_coll[0], saveFile = False)

        his_x = f.Get("shapes_fit_s/" + directoryName + "/mis_id")
        histo_x = ROOT.TH1F("mis_id","mis_id", var.nbins,var.xmin, var.xmax)
        for i in range(1,his_x.GetNbinsX()+1):
            histo_x.SetBinContent(i,his_x.GetBinContent(i))
            histo_x.SetBinError(i,his_x.GetBinError(i))
            print(i,histo_x.GetBinContent(i))

        makeSinglePlot(histo_x, var, sample_coll[2], saveFile = False)
        print(histo_x.Integral())

        his_d = f.Get("shapes_fit_s/" + directoryName + "/data")
        histo_d = ROOT.TH1F("data_obs","data_obs", var.nbins, var.xmin, var.xmax)
        for i in range(0,his_d.GetN()):
            histo_d.SetBinContent(i+1,his_d.GetPointY(i))

        makeSinglePlot(histo_d, var, sample_coll[-1], saveFile = False)

        makeStack(histo_d,[histo_m,histo_x], var, fit=True, addFileName = directoryName)
        '''

f.Close()
