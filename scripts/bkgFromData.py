import math
from variable import variable

#from utils.pyUtils import MakeRooDataSet as makeds
from utils.cmsstyle import CMS_lumi
from ROOT import gSystem, gROOT, gStyle
from ROOT import TCanvas, TLegend, TH1F, TF1
from root_pandas import read_root

import yaml
import numpy as np
import pandas as pd

dir_plots = 'plots/combinatorial_bkg'
dir_results = 'results/combinatorial_bkg'

def readIntervalsFromFile():
  try:
    with open('%s/sideBandsIntervals.yaml' % (dir_results)) as f:
      intervals=yaml.full_load(f)
      return intervals
  except FileNotFoundError:
    print('sideBandsInterval.yaml not found')

  return None

def createHisto(df, histoName, var, norm=1, over=True):
  histo  = TH1F(histoName, "", var.nbins, var.xmin, var.xmax)
  for varVal in df[var.name]:
    histo.Fill(varVal)

  if(over):
    histo.SetBinContent(1, histo.GetBinContent(0) + histo.GetBinContent(1))
    histo.SetBinError(1, math.sqrt(pow(histo.GetBinError(0), 2) + pow(histo.GetBinError(1), 2)))
    histo.SetBinContent(var.nbins, histo.GetBinContent(var.nbins) + histo.GetBinContent(var.nbins+1))
    histo.SetBinError(var.nbins, math.sqrt(pow(histo.GetBinError(var.nbins), 2) + pow(histo.GetBinError(var.nbins+1), 2)))

  hist_int = histo.Integral()
  return histo, hist_int

def gaus_pol_fit(histo, var,plt_name="", path="", ps=[None,None,None,None,None,None], pol="pol1"):
  c = TCanvas("","",800, 600)
  func= TF1("model", "gaus(0) +"+ pol + "(3)")
  for i,p in enumerate(ps):
    if(p!=None):
      func.SetParameter(i,p)

  for i in range(1, histo.GetNbinsX()+1):
    histo.SetBinError(i, math.sqrt(histo.GetBinContent(i)))
  fit_result = histo.Fit(func, "S")
  histo.Draw()
  CMS_lumi(c, 4, 0)
  c.SetLogy()
  c.SaveAs("%s/dimuon_mass_fit.pdf" % (dir_plots))
  c.SaveAs("%s/dimuon_mass_fit.png" % (dir_plots))
  return fit_result

  #bin_width = (var.xmax - var.xmin)/(1.*nbins)
  #histo.SetTitle(';Events/%2.3f %s ; %s %s; Counts' % (bin_width, var.unit, var.xlabel, var.unit)

def get_int_gaus_pol(fit_result, minVal, maxVal):
  gaus = TF1("gaus","gaus(0)")
  pol = TF1("pol", "pol1(0)")

  gaus.SetParameter(0, fit_result.Parameter(0))
  gaus.SetParameter(1, fit_result.Parameter(1))
  gaus.SetParameter(2, fit_result.Parameter(2))
  
  pol.SetParameter(0, fit_result.Parameter(3))
  pol.SetParameter(1,fit_result.Parameter(4))

  gaus_int = gaus.Integral(minVal, maxVal)
  pol_int = pol.Integral(minVal, maxVal)

  return gaus_int, pol_int

def getSBIntervals(df):
  '''
  Get the interval of the side bands (SB)
  making a fit of the invariant mass of the
  dimuon invariant mass.
  The SB will be defined based on the mean and 
  the sigma of the fit of the resonance
  '''
  var = variable("jpsi_mass", "m(2#mu)", "m(2#mu)", "[GeV]", 40, 2., 4.)
  histoName = "jpsi_mass"
  massHisto, massHisto_int = createHisto(df, histoName, var, 1, True)
  fit_result = gaus_pol_fit(massHisto, var, "Dimuon invariant mass","", [None,3.1,.1,None,None,None], "pol1")
  var_m = fit_result.Parameter(1) #jpsi mass mean
  var_s = fit_result.Parameter(2) #jpsi mass sigma



  ls = 3 # lower number of sigmas from the mean
  us = 7 # upper number of sigmas from the mean

  intervals = {'lsb':{'minVal': var_m - us*var_s , 'maxVal':var_m - ls*var_s}, 'rsb':{'minVal': var_m + ls*var_s, 'maxVal':var_m+ us*var_s}}
  gaus_int, pol_int = get_int_gaus_pol(fit_result, intervals['lsb']['maxVal'], intervals['rsb']['minVal'])
  print('gaus_int = %5.3f, pol_int = %5.3f' % (gaus_int, pol_int))
  with open('%s/sideBandsIntervals.yaml' % (dir_results), 'w') as f:
    yaml.dump(intervals, f)
  return intervals


def plotComparison(df, var_names, categories, cuts):
  for var_name in var_names:
    nStdevs = 2
    nbins =40
    mean = df[var_name].mean(axis = 0, skipna = True)
    sd = df[var_name].std(axis = 0, skipna = True)
    print(mean)
    print(sd)
    rmin = mean - nStdevs * sd
    rmax = mean + nStdevs * sd

    var = variable(var_name, "", var_name, "", nbins, rmin, rmax)
    histos = [createHisto(df.query(cut), cat, var, 1, True) for cat, cut in zip(categories, cuts)]
    c = TCanvas("","",800, 600)
    for i, histo in enumerate(histos):
      histo[0].SetLineColor(i+1)
      histo[0].Scale(1./histo[1])
      histo[0].Draw("same")
    CMS_lumi(c, 4, 0)
    c.SaveAs("%s/comparison_%s.pdf" % (dir_plots, var_name))
    c.SaveAs("%s/comparison_%s.png" % (dir_plots, var_name))
    del histo, histos, c


def main():
  force_fit = True
  inputFile = '/afs/cern.ch/user/g/garamire/work/private/CMSPisa/RJPsiAnalysis/ntuples/2021Jan29/data_UL_2021Jan29.root'
  #inputFile = 'dataframes_local/list_ptmax_UL_flags.root'
  inputTree = 'BTo3Mu'
  data_df = read_root(inputFile, inputTree)
  data_df = data_df[~data_df.event.duplicated(keep=False)]
  sb_intervals = readIntervalsFromFile()
  if (sb_intervals == None or force_fit):
    sb_intervals = getSBIntervals(data_df) # Side bands intervals
  cuts = ["jpsi_mass > %5.3f and jpsi_mass < %5.3f" % (sb_intervals['lsb']['minVal'], sb_intervals['lsb']['maxVal']) ,
          "jpsi_mass > %5.3f and jpsi_mass < %5.3f" % (sb_intervals['rsb']['minVal'], sb_intervals['rsb']['maxVal']) ]
  var_names = ['Q_sq', 'm_miss_sq', 'pt_var']
  categories = ["lsb", "rsb"]
  plotComparison(data_df, var_names, categories, cuts)







if __name__ == '__main__':
  gROOT.SetBatch()
  gStyle.SetOptStat(0)
  main()
