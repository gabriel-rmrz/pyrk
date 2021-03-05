import math
from variable import variable
#from utils.pyUtils import MakeRooDataSet as makeds
from utils.cmsstyle import CMS_lumi
import ROOT
from ROOT import gSystem, gROOT, gStyle
from ROOT import TCanvas, TLegend, TH1F, TF1
from root_pandas import read_root
from uproot_methods import TLorentzVectorArray
from uproot_methods import TLorentzVector
import yaml
import numpy as np
import pandas as pd

dir_plots = 'plots/combinatorial_bkg'
dir_results = 'results/combinatorial_bkg'
def readIntervalsFromFile():
  '''
  Reads the yaml file with the information from the fit.
  '''
  try:
    with open('%s/sideBandsIntervals.yaml' % (dir_results)) as f:
      intervals=yaml.full_load(f)
      return intervals
  except FileNotFoundError:
    print('sideBandsInterval.yaml not found')
  return None

def createHisto(df, histoName, var, norm=1, over=True):
  '''
  Make a histogram form a dataframe
  '''
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
  '''
  Perform a fit with a gaussian for the resonance and 
  a exponential or a polynomial as background
  '''
  c = TCanvas("","",800, 600)
  func= TF1("model", "gaus(0) +"+ pol + "(3)", var.xmin, var.xmax)
  for i,p in enumerate(ps):
    if(p!=None):
      func.SetParameter(i,p)
  for i in range(1, histo.GetNbinsX()+1):
    histo.SetBinError(i, math.sqrt(histo.GetBinContent(i)))
  fit_result = histo.Fit(func, "S")
  histo.Draw()
  CMS_lumi(c, 4, 0)
  #c.SetLogy()
  c.SaveAs("%s/dimuon_mass_fit.pdf" % (dir_plots))
  c.SaveAs("%s/dimuon_mass_fit.png" % (dir_plots))
  return fit_result

def get_int_gaus_pol(fit_result, minVal, maxVal):
  '''
  Get the integral of a gaussian plus a polynomial
  '''
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
  var = variable("jpsi_mass", "m(2#mu)", "m(2#mu)", "[GeV]", 40, 2.97, 3.23)
  histoName = "jpsi_mass"
  massHisto, massHisto_int = createHisto(df, histoName, var, 1, False)
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

def plotComparison(df, var_list, categories, cuts, legends, prefix):
  '''
  Plot a list a variables from a pandas dataframe
  '''
  for var in var_list:
    histos = [createHisto(df.query(cut), cat, var, 1, True) for cat, cut in zip(categories, cuts)]
    c = TCanvas("","",800, 600)
    legend = TLegend(0.67, 0.7, 0.9, 0.86)
    for i, histo in enumerate(histos):
      legend.AddEntry(histo[0], legends[i], "lp")
      histo[0].SetLineColor(i+1)
      histo[0].GetYaxis().SetTitle("a.u.")
      histo[0].GetXaxis().SetTitle("%s [%s]" % (var.xlabel, var.unit))
      histo[0].Scale(1./histo[1])
      histo[0].DrawNormalized("same")
    legend.SetTextFont(43)
    legend.SetTextSize(15)
    legend.SetBorderSize(0)
    legend.SetFillColor(0)
    legend.Draw("SAME")
    CMS_lumi(c, 4, 0)
    c.SaveAs("%s/%s_comparison_%s.pdf" % (dir_plots, prefix, var.name))
    c.SaveAs("%s/%s_comparison_%s.png" % (dir_plots, prefix, var.name))
    del histo, histos, c

def plotComparison2df(df1, df2, var_list, categories, cuts, legends, prefix):
  '''
  Plot a list of variables from 2 pandas dataframes
  '''
  for var in var_list:
    histo1 = createHisto(df1, categories[0], var_list[0], 1, True)
    histo2 = createHisto(df1, categories[1], var_list[1], 1, True) 
    histo3 = createHisto(df2, categories[2], var_list[1], 1, True) 
    c = TCanvas("","",800, 600)
    legend = TLegend(0.17, 0.7, 0.4, 0.86)
    legend.AddEntry(histo1[0], legends[0], "lp")
    histo1[0].SetLineColor(1)
    histo1[0].GetYaxis().SetTitle("a.u.")
    histo1[0].GetXaxis().SetTitle("%s [%s]" % (var_list[0].xlabel, var_list[0].unit))
    histo1[0].DrawNormalized("same")
    legend.AddEntry(histo2[0], legends[1], "lp")
    histo2[0].SetLineColor(2)
    histo2[0].GetYaxis().SetTitle("a.u.")
    histo2[0].GetXaxis().SetTitle("%s [%s]" % (var_list[1].xlabel, var_list[1].unit))
    histo2[0].DrawNormalized("same")
    legend.AddEntry(histo3[0], legends[2], "lp")
    histo3[0].SetLineColor(3)
    histo3[0].GetYaxis().SetTitle("a.u.")
    histo3[0].GetXaxis().SetTitle("%s [%s]" % (var_list[1].xlabel, var_list[1].unit))
    histo3[0].DrawNormalized("same")
    legend.SetTextFont(43)
    legend.SetTextSize(15)
    legend.SetBorderSize(0)
    legend.SetFillColor(0)
    legend.Draw("SAME")
    CMS_lumi(c, 4, 0)
    print('Kolmo: ')
    histo1[0].Scale(1./histo1[0].Integral())
    histo3[0].Scale(1./histo3[0].Integral())
    print(histo3[0].KolmogorovTest(histo1[0], "NX"))
    c.SaveAs("%s/%s_comparison_%s.pdf" % (dir_plots, prefix, var_list[1].name))
    c.SaveAs("%s/%s_comparison_%s.png" % (dir_plots, prefix, var_list[1].name))
    del histo1, histo2, histo3, c

def plotComparisonByCats(df, categories, variables, prefix):
  for var in variables:
    legend = TLegend(0.57, 0.7, 0.8, 0.86)
    c1 = TCanvas("","",800, 600)
    histos = [createHisto(df.query(cat.cuts), cat.name, var, 1, True) for cat in categories]
    for histo, cat in zip(histos, categories):
      legend.AddEntry(histo[0], cat.legend, "lp")
      histo[0].SetLineColor(cat.color)
      histo[0].SetMarkerStyle(cat.marker)
      histo[0].SetMarkerColor(cat.color)
      #histo[0].SetMarkerSize(3)
      histo[0].GetYaxis().SetTitle("a.u.")
      histo[0].GetXaxis().SetTitle("%s [%s]" % (var.xlabel, var.unit))
      histo[0].Scale(1./histo[1])
      histo[0].Draw("same")
    legend.SetTextFont(43)
    legend.SetTextSize(15)
    legend.SetBorderSize(0)
    legend.SetFillColor(0)
    legend.Draw("SAME")
    CMS_lumi(c1, 4, 0)
    c1.SaveAs("%s/%s_comparison_%s.pdf" % (dir_plots, prefix, var.name))
    c1.SaveAs("%s/%s_comparison_%s.png" % (dir_plots, prefix, var.name))
    del c1, legend, histos


class category(object):
  def __init__(self, name, cuts, legend, color, marker):
    self.name = name
    self.cuts = cuts
    self.legend = legend
    self.color = color
    self.marker = marker

def main():
  force_fit = False
  #inputFile = '/afs/cern.ch/user/g/garamire/work/private/CMSPisa/RJPsiAnalysis/ntuples/2021Jan29/data_UL_2021Jan29.root'
  inputFile = 'dataframes_local/data2018_2021Feb12.root'
  inputTree = 'BTo3Mu'
  common_cuts='jpsi_eta < 1.2'
  data_df = read_root(inputFile, inputTree)
  #data_df = data_df[~data_df.event.duplicated(keep=False)]
  data_df = data_df.query(common_cuts).copy()
  data_df['q2_normBykpt2'] = data_df['Q_sq'].divide(data_df['kpt'] * data_df['kpt'])
  data_df['q2_normBydimumass'] = data_df['Q_sq'].divide(data_df['jpsi_mass'] )
  sb_intervals = readIntervalsFromFile()
  if (sb_intervals == None or force_fit):
    sb_intervals = getSBIntervals(data_df) # Side bands intervals
  sb_intervals_psi2s = {'lsb':{'minVal': 3.3 , 'maxVal':3.55}, 'rsb':{'minVal': 3.75, 'maxVal':4}}

  trigger_cuts = "mu1_isDimuon0Trg > 0.5" 
  trigger_cuts2 = "mu1_isJpsiTrk_PsiPrimeTrg > 0.5" 
  cat_jpsi_lsb = category("jpsi_lsb", 
      "%s and jpsi_mass > %5.3f and jpsi_mass < %5.3f" % ( trigger_cuts, sb_intervals['lsb']['minVal'], sb_intervals['lsb']['maxVal']), 
      "J/#psi left side band", 
      1,
      ROOT.kFullCircle)
  cat_jpsi_rsb = category("jpsi_rsb", 
      "%s and jpsi_mass > %5.3f and jpsi_mass < %5.3f" % ( trigger_cuts, sb_intervals['rsb']['minVal'], sb_intervals['rsb']['maxVal']), 
      "J/#psi right side band", 
      2,
      ROOT.kFullSquare)
  cat_psi2s_lsb = category("psi2s_lsb", 
      "%s and jpsi_mass > %5.3f and jpsi_mass < %5.3f" % ( trigger_cuts2, sb_intervals_psi2s['lsb']['minVal'], sb_intervals_psi2s['lsb']['maxVal']), 
      "#psi(2S) left side band", 
      6,
      ROOT.kFullTriangleUp)
  cat_psi2s_res = category("psi2s_res", 
      "%s and jpsi_mass > %5.3f and jpsi_mass < %5.3f" % ( trigger_cuts2, sb_intervals_psi2s['lsb']['maxVal'], sb_intervals_psi2s['rsb']['minVal']), 
      "#psi(2S) resonance region", 
      8,
      ROOT.kFullTriangleDown)
  cat_psi2s_rsb = category("psi2s_rsb", 
      "%s and jpsi_mass > %5.3f and jpsi_mass < %5.3f" % ( trigger_cuts2, sb_intervals_psi2s['rsb']['minVal'], sb_intervals_psi2s['rsb']['maxVal']), 
      "#psi(2S) right side band", 
      9,
      ROOT.kOpenCircle)

  var_list = [
      variable("Q_sq", "Q^{2}", "Q^{2}", "GeV^{2}", 40, 0., 12.),
      variable("m_miss_sq", "m_{miss}^{2}", "m_{miss}^{2}", "GeV^{2}", 40, 0., 10.),
      variable("pt_var", "pt_var", "pt_var", "GeV", 50, 0., 50.),
      variable("E_mu_star", "E^{*}", "E^{*}", "GeV", 40, 0., 4.),
      variable("q2_normBykpt2", "Q^{2}/pt(k)^{2}", "Q^{2}/pt(k)^{2}", "", 40, -0.1, 0.9),
      variable("q2_normBydimumass", "Q^{2}/m(#mu^{+}#mu^{-})", "Q^{2}/m(#mu^{+}#mu^{-})", "", 40, -1, 3.)
      ]
  cats=[cat_jpsi_lsb, cat_jpsi_rsb]
  plotComparisonByCats(data_df, cats, var_list, "jpsi_sidebands")

  cats_psi2s=[cat_psi2s_lsb, cat_psi2s_res, cat_psi2s_rsb]
  plotComparisonByCats(data_df, cats_psi2s, var_list, "psi2s_sidebandsAndResonace")

  cuts = ["jpsi_mass > %5.3f and jpsi_mass < %5.3f" % ( sb_intervals['lsb']['minVal'], sb_intervals['lsb']['maxVal']) ,
          "jpsi_mass > %5.3f and jpsi_mass < %5.3f" % ( sb_intervals['rsb']['minVal'], sb_intervals['rsb']['maxVal']) ]
  legends = [ "#psi(2S) left side band",
             "#psi(2S) right side band"]
  categories = ["lsb", "rsb"]
  prefix = 'psiprimeSB'

  #plotComparison(data_df, var_list, categories, cuts, legends, prefix)
  
  '''
  sb_intervals_psiprime = {'lsb':{'minVal': 3.3 , 'maxVal':3.55}, 'rsb':{'minVal': 3.75, 'maxVal':4}}
  cuts2 = ["jpsi_mass > %5.3f and jpsi_mass < %5.3f" % (sb_intervals_psiprime['lsb']['minVal'],
            sb_intervals_psiprime['lsb']['maxVal']) ,
          "jpsi_mass > %5.3f and jpsi_mass < %5.3f" % (sb_intervals_psiprime['lsb']['maxVal'],
            sb_intervals_psiprime['rsb']['minVal']) ,
          "jpsi_mass > %5.3f and jpsi_mass < %5.3f" % (sb_intervals_psiprime['rsb']['minVal'],
            sb_intervals_psiprime['rsb']['maxVal']) ]
  legends2 = [ "#psi(2S) left side band",
              "#psi(2S) mass region",
             "#psi(2S) right side band"]
  categories2 = ["lsb", "res", "rsb"]
  prefix2 = 'psiprimeAll'
  plotComparison(data_df, var_list, categories2, cuts2, legends2, prefix2)
  '''

  lsb_df = data_df.query(cuts2[0]).copy()
  rsb_df = data_df.query(cuts2[2]).copy()

  mass_mean_lsb = lsb_df['jpsi_mass'].mean()
  mass_mean_rsb = rsb_df['jpsi_mass'].mean()

  BC_MASS_PDG = 6.2756
  jpsi_p4 = TLorentzVectorArray.from_ptetaphim(lsb_df.jpsi_pt, lsb_df.jpsi_eta, lsb_df.jpsi_phi,
      lsb_df.jpsi_mass)
  b_p4_corrected = TLorentzVectorArray.from_ptetaphim(BC_MASS_PDG*lsb_df.Bpt.divide(lsb_df.Bmass), lsb_df.Beta, lsb_df.Bphi, BC_MASS_PDG)
  b_p4 = TLorentzVectorArray.from_ptetaphim(lsb_df.Bpt, lsb_df.Beta, lsb_df.Bphi, lsb_df.Bmass)

  lsb_df['Q_sq_extrapolated'] = (b_p4_corrected - (mass_mean_rsb/mass_mean_lsb)*jpsi_p4).mag2
  var_list2 = [
      variable("Q_sq_extrapolated", "Q^{2}", "Q^{2}", "GeV^{2}", 60, -5., 10.),
      variable("Q_sq", "Q^{2}", "Q^{2}", "GeV^{2}", 60, -5., 10.)
      ]
  cuts3= [ '','','']
  categories3 = ["lsb_corrected", "lsb", "rsb"]
  legends3 = [ "#psi(2S) left side band extrapolated",
               "#psi(2S) left side band",
             "#psi(2S) right side band"]
  prefix3 = 'extrapolated'
  plotComparison2df(lsb_df, rsb_df, var_list2, categories3, cuts3, legends3, prefix3)

if __name__ == '__main__':
  gROOT.SetBatch()
  gStyle.SetOptStat(0)
  main()
