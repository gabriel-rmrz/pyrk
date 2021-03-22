import math
from utils.variable import variable
from utils.category import category
from utils.selections import preselection, pass_id

#from utils.pyUtils import MakeRooDataSet as makeds
from utils.cmsstyle import CMS_lumi
import ROOT
from ROOT import gSystem, gROOT, gStyle
from ROOT import TCanvas, TLegend, TH1F, TF1, TProfile
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

def createProfilePol1Fit(df, var_x, var_y):
  c = TCanvas("","",800, 600)
  prof = TProfile("prof", "", var_x.nbins, var_x.xmin, var_x.xmax, var_y.xmin, var_y.xmax)
  for var_x_val, var_y_val in zip(df[var_x.name],df[var_y.name]) :
    prof.Fill(var_x_val, var_y_val,1)
  
  myPol = TF1("myPol", "[0]+[1]*x", var_x.xmin, var_x.xmax)
  fit_res = prof.Fit("myPol","W").Get()
  print(var_x.xmin)
  print(var_x.xmax)
  print(myPol.GetParameter(0))
  print(myPol.GetParameter(1))
  prof.Draw()
  CMS_lumi(c, 4, 0)
  c.SaveAs("%s/dimuon_mass_%3.1f_%3.1f_vs_q_sq_profile.pdf" % (dir_plots, var_x.xmin, var_x.xmax))
  c.SaveAs("%s/dimuon_mass_%3.1f_%3.1f_vs_q_sq_profile.png" % (dir_plots, var_x.xmin, var_x.xmax))
  return [myPol.GetParameter(0), myPol.GetParameter(1)]


def get_sb_fit_params(categories, vars_x, vars_y):
  '''
  Takes the two categories of the side bands,
  does fit in the two region and extrapolates
  a line that joins the upper limit of the lsb
  and the lower limit of the rsb with line and
  issues the the parameter of the three lines
  '''
  fit_params = []
  for cat, var_x, var_y in zip(categories, vars_x, vars_y):
    fit_params.append(createProfilePol1Fit(cat.df.query(cat.cuts), var_x, var_y))
  
  x1 = vars_x[0].xmax
  x2 = vars_x[1].xmin
  x_r = x1/x2

  b1 = fit_params[0][0]
  m1 = fit_params[0][1]


  b3 = fit_params[1][0]
  m3 = fit_params[1][1]


  b2 = (m1*x1 + b1 - x_r * (m3*x2 + b3))/(1 - x_r)
  m2 = (m1*x1 + b1 -b2)/x1

  return m1 ,b1 , m2, b2, m3, b3

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

def get_int_gaus_pol(fit_result, minVal, maxVal, nBins):
  '''
  Get the integral of a gaussian plus a polynomial
  '''
  print(fit_result.NPar())
  gaus = TF1("gaus","gaus(0)")
  pol = TF1("pol", "pol1(0)")
  gaus.SetParameter(0, fit_result.Parameter(0))
  gaus.SetParameter(1, fit_result.Parameter(1))
  gaus.SetParameter(2, fit_result.Parameter(2))
  pol.SetParameter(0, fit_result.Parameter(3))
  pol.SetParameter(1,fit_result.Parameter(4))
  normVal = (maxVal - minVal)/nBins
  gaus_int = gaus.Integral(minVal, maxVal)/normVal
  pol_int = pol.Integral(minVal, maxVal)/normVal
  return gaus_int, pol_int

def getSBIntervals(cat):
  '''
  Get the interval of the side bands (SB)
  making a fit of the invariant mass of the
  dimuon invariant mass.
  The SB will be defined based on the mean and 
  the sigma of the fit of the resonance
  '''
  #var = variable("jpsi_mass", "m(2#mu)", "m(2#mu)", "[GeV]", 40, 2.97, 3.23)
  var = variable("jpsi_mass", "m(2#mu)", "m(2#mu)", "[GeV]", 40, 3.45, 3.95)
  histoName = "jpsi_mass"
  massHisto, massHisto_int = createHisto(cat.df.query(cat.cuts), histoName, var, 1, False)
  fit_result = gaus_pol_fit(massHisto, var, "Dimuon invariant mass","", [None,3.7,.1,None,None,None], "pol1")
  var_m = fit_result.Parameter(1) #jpsi mass mean
  var_s = fit_result.Parameter(2) #jpsi mass sigma

  ls = 3.0 # lower number of sigmas from the mean
  us = 8.0 # upper number of sigmas from the mean

  intervals = {'lsb':{'minVal': var_m - us*var_s , 'maxVal':var_m - ls*var_s}, 'rsb':{'minVal': var_m + ls*var_s, 'maxVal':var_m+ us*var_s}}
  print(intervals)
  fit_integrals = get_int_gaus_pol(fit_result, intervals['lsb']['maxVal'], intervals['rsb']['minVal'], massHisto.GetNbinsX())
  print('gaus_int = %5.3f, pol_int = %5.3f' % fit_integrals)
  with open('%s/sideBandsIntervals.yaml' % (dir_results), 'w') as f:
    yaml.dump(intervals, f)
  return intervals, fit_integrals

def plotComparisonByCats(categories, variables, prefix, normalize=True):
  for var in variables:
    legend = TLegend(0.52, 0.75, 0.89, 0.86)
    c1 = TCanvas("","",800, 600)
    histos = [createHisto(cat.df.query(cat.cuts), cat.name, var, 1, False) for cat in categories]
    #maxValHisto = 0
    for histo, cat in zip(histos, categories):
      legend.AddEntry(histo[0], cat.legend, "lp")
      if(normalize): 
        histo[0].Scale(1./histo[1])
      else:
        histo[0].Scale(cat.weight)
      #if(maxValHisto < histo[0].GetMaximum()): 
      maxValHisto = histo[0].GetMaximum()
      histo[0].SetLineColor(cat.color)
      histo[0].SetMarkerStyle(cat.marker)
      histo[0].SetMarkerColor(cat.color)
      histo[0].SetMarkerSize(1)
      histo[0].GetYaxis().SetTitle("a.u.")
      histo[0].GetXaxis().SetTitle("%s [%s]" % (var.xlabel, var.unit))
      histo[0].SetMaximum(1.5 * maxValHisto)
      histo[0].Draw("same")
    legend.SetTextFont(43)
    legend.SetTextSize(15)
    legend.SetBorderSize(0)
    legend.SetFillColor(0)
    legend.Draw("SAME")
    CMS_lumi(c1, 4, 0)
    c1.SaveAs("%s/%s_%s.pdf" % (dir_plots, prefix, var.name))
    c1.SaveAs("%s/%s_%s.png" % (dir_plots, prefix, var.name))
    del c1, legend, histos

def main():
  force_fit = True
  inputFile = 'dataframes_local/data2018_2021Feb12.root'
  inputTree = 'BTo3Mu'
  common_cuts=" & ".join([ preselection, pass_id])
  data_df = read_root(inputFile, inputTree)
  #data_df = data_df[~data_df.event.duplicated(keep=False)]
  data_df = data_df.query(common_cuts).copy()
  data_df['q2_normBykpt2'] = data_df['Q_sq'].divide(data_df['kpt'] * data_df['kpt'])
  data_df['q2_normBydimumass'] = data_df['Q_sq'].divide(data_df['jpsi_mass'] )
  #sb_intervals = readIntervalsFromFile()
  #if (sb_intervals == None or force_fit):

  #sb_intervals_psi2s = {'lsb':{'minVal': 3.4 , 'maxVal':3.5}, 'rsb':{'minVal': 3.8, 'maxVal':3.95}}

  trigger_dimuon0 = "mu1_isDimuon0Trg & mu2_isDimuon0Trg & k_isDimuon0Trg" 
  trigger_psiprime = "mu1_isJpsiTrk_PsiPrimeTrg & mu2_isJpsiTrk_PsiPrimeTrg & k_mu1_isJpsiTrk_PsiPrimeTrg" 
  trigger_jpsitrk = "mu1_isJpsiTrkTrg & mu2_isJpsiTrkTrg & k_isJpsiTrkTrg" 
  trigger_nonresonant = "mu1_isJpsiTrk_NonResonantTrg & mu2_isJpsiTrk_NonResonantTrg & k_isJpsiTrk_NonResonantTrg" 
  trigger_psiprime = "mu1_isJpsiTrk_PsiPrimeTrg " 

  trigger_dimuon0_alone = "%s and not ( %s or %s or %s)" % (trigger_dimuon0 , trigger_psiprime, trigger_jpsitrk, trigger_nonresonant)
  trigger_psiprime_alone = "%s and not ( %s or %s or %s)" % (trigger_psiprime, trigger_dimuon0 , trigger_jpsitrk, trigger_nonresonant)

  min_mass_psi2s = 3.35
  max_mass_psi2s = 3.95
  cat_jpsi_jpsitrkTrg = category("jpsi_jpsitrkTrg", 
      data_df,
      "%s and jpsi_mass > 2.95 and jpsi_mass < 3.25" % ( trigger_jpsitrk), 
      "J/#psi mass region + jpsitrkTrg", 
      1,
      ROOT.kFullTriangleUp,
      1.0)
  cat_jpsi_dimuon0Trg = category("jpsi_dimuon0", 
      data_df,
      "%s and jpsi_mass > 2.95 and jpsi_mass < 3.25" % ( trigger_dimuon0), 
      "J/#psi mass region + dimuon0Trg", 
      2,
      ROOT.kFullSquare,
      1.0)
  cat_psi2s = category("psi2s", 
      data_df,
      "%s and jpsi_mass > %5.3f and jpsi_mass < %5.3f" % ( trigger_psiprime_alone, min_mass_psi2s, max_mass_psi2s), 
      "#psi(2S) ", 
      ROOT.kBlue - 4,
      ROOT.kFullCircle,
      1.0)
  
  sb_intervals_psi2s, fit_integrals = getSBIntervals(cat_psi2s)

  cat_psi2s_lsb = category("psi2s_lsb", 
      data_df,
      "%s and jpsi_mass > %5.3f and jpsi_mass < %5.3f" % ( trigger_psiprime_alone, sb_intervals_psi2s['lsb']['minVal'], sb_intervals_psi2s['lsb']['maxVal']), 
      "#psi(2S) left side band", 
      6,
      ROOT.kOpenCircle,
      1.0)
  cat_psi2s_res = category("psi2s_res", 
      data_df,
      "%s and jpsi_mass > %5.3f and jpsi_mass < %5.3f" % ( trigger_psiprime_alone, sb_intervals_psi2s['lsb']['maxVal'], sb_intervals_psi2s['rsb']['minVal']), 
      "#psi(2S) resonance region", 
      8,
      ROOT.kOpenCircle,
      1.0)
  cat_psi2s_rsb = category("psi2s_rsb", 
      data_df,
      "%s and jpsi_mass > %5.3f and jpsi_mass < %5.3f" % ( trigger_psiprime_alone, sb_intervals_psi2s['rsb']['minVal'], sb_intervals_psi2s['rsb']['maxVal']), 
      "#psi(2S) right side band", 
      ROOT.kOrange,
      ROOT.kOpenCircle,
      1.0)

  var_list = [
      variable("Q_sq", "Q^{2}", "Q^{2}", "GeV^{2}", 80, 0., 12.),
      variable("m_miss_sq", "m_{miss}^{2}", "m_{miss}^{2}", "GeV^{2}", 60, 0., 10.),
      variable("pt_var", "pt_var", "pt_var", "GeV", 80, 0., 50.),
      variable("E_mu_star", "E^{*}", "E^{*}", "GeV", 80, 0., 2.),
      #variable("q2_normBykpt2", "Q^{2}/pt(k)^{2}", "Q^{2}/pt(k)^{2}", "", 40, -0.1, 0.9),
      #variable("q2_normBydimumass", "Q^{2}/m(#mu^{+}#mu^{-})", "Q^{2}/m(#mu^{+}#mu^{-})", "", 40, -1, 3.)
      ]

###############
###### bkg from sidebands BEGIN #####
  psi2s_lsb_df = data_df.query(cat_psi2s_lsb.cuts).copy()
  psi2s_rsb_df = data_df.query(cat_psi2s_rsb.cuts).copy()
  psi2s_res_df = data_df.query(cat_psi2s_res.cuts)

  dimuon_mass_mean_lsb = psi2s_lsb_df['jpsi_mass'].mean()
  dimuon_mass_mean_rsb = psi2s_rsb_df['jpsi_mass'].mean()
  dimuon_mass_mean_res = psi2s_res_df['jpsi_mass'].mean()

  lsb_middle = (sb_intervals_psi2s['lsb']['maxVal'] + sb_intervals_psi2s['lsb']['minVal']) /2.0
  rsb_middle = (sb_intervals_psi2s['rsb']['maxVal'] + sb_intervals_psi2s['rsb']['minVal']) /2.0
  res_middle = (sb_intervals_psi2s['rsb']['minVal'] + sb_intervals_psi2s['lsb']['maxVal']) /2.0

  BC_MASS_PDG = 6.2756
  dimuon_lsb_p4 = TLorentzVectorArray.from_ptetaphim(psi2s_lsb_df.jpsi_pt, psi2s_lsb_df.jpsi_eta, psi2s_lsb_df.jpsi_phi, psi2s_lsb_df.jpsi_mass)
  k_lsb_p4 = TLorentzVectorArray.from_ptetaphim(psi2s_lsb_df.kpt, psi2s_lsb_df.keta, psi2s_lsb_df.kphi, psi2s_lsb_df.kmass)
  b_p4_lsb = k_lsb_p4 + (res_middle/lsb_middle)*dimuon_lsb_p4
  b_p4_lsb_corrected = TLorentzVectorArray.from_ptetaphim(BC_MASS_PDG*np.divide(b_p4_lsb.pt, b_p4_lsb.mass), b_p4_lsb.eta, b_p4_lsb.phi, BC_MASS_PDG)
  #b_p4_lsb_corrected = TLorentzVectorArray.from_ptetaphim(BC_MASS_PDG*psi2s_lsb_df.Bpt.divide(psi2s_lsb_df.Bmass), psi2s_lsb_df.Beta, psi2s_lsb_df.Bphi, BC_MASS_PDG)

  psi2s_res_extrapolated_from_lsb_df = pd.DataFrame()
  psi2s_res_extrapolated_from_lsb_mass_scaled_df = pd.DataFrame()
  #psi2s_res_extrapolated_from_lsb_mass_scaled_df['Q_sq'] = (b_p4_lsb_corrected - (dimuon_mass_mean_res/dimuon_mass_mean_lsb)*dimuon_lsb_p4).mag2
  psi2s_res_extrapolated_from_lsb_mass_scaled_df['Q_sq'] = (b_p4_lsb_corrected - (res_middle/lsb_middle)*dimuon_lsb_p4).mag2
  psi2s_res_extrapolated_from_lsb_df['Q_sq'] = (b_p4_lsb_corrected - dimuon_lsb_p4).mag2

  dimuon_rsb_p4 = TLorentzVectorArray.from_ptetaphim(psi2s_rsb_df.jpsi_pt, psi2s_rsb_df.jpsi_eta, psi2s_rsb_df.jpsi_phi, psi2s_rsb_df.jpsi_mass)
  k_rsb_p4 = TLorentzVectorArray.from_ptetaphim(psi2s_rsb_df.kpt, psi2s_rsb_df.keta, psi2s_rsb_df.kphi, psi2s_rsb_df.kmass)
  b_p4_rsb = k_rsb_p4 + (res_middle/rsb_middle)*dimuon_rsb_p4
  b_p4_rsb_corrected = TLorentzVectorArray.from_ptetaphim(BC_MASS_PDG*np.divide(b_p4_rsb.pt, b_p4_rsb.mass), b_p4_rsb.eta, b_p4_rsb.phi, BC_MASS_PDG)
  #b_p4_rsb_corrected = TLorentzVectorArray.from_ptetaphim(BC_MASS_PDG*psi2s_rsb_df.Bpt.divide(psi2s_rsb_df.Bmass), psi2s_rsb_df.Beta, psi2s_rsb_df.Bphi, BC_MASS_PDG)

  psi2s_res_extrapolated_from_rsb_df = pd.DataFrame()
  psi2s_res_extrapolated_from_rsb_mass_scaled_df = pd.DataFrame()
  psi2s_res_extrapolated_from_rsb_df['Q_sq'] = (b_p4_rsb_corrected - dimuon_rsb_p4).mag2
  #psi2s_res_extrapolated_from_rsb_mass_scaled_df['Q_sq'] = (b_p4_rsb_corrected - (dimuon_mass_mean_res/dimuon_mass_mean_rsb)*dimuon_rsb_p4).mag2
  psi2s_res_extrapolated_from_rsb_mass_scaled_df['Q_sq'] = (b_p4_rsb_corrected - (res_middle/rsb_middle)*dimuon_rsb_p4).mag2
  
  psi2s_res_extrapolated_from_sb_df = pd.concat([psi2s_res_extrapolated_from_rsb_mass_scaled_df, psi2s_res_extrapolated_from_lsb_mass_scaled_df])

  var_list3 = [variable("Q_sq", "Q^{2}", "Q^{2}", "GeV^{2}", 60, 2., 8.)]

  cats_psi2s=[cat_psi2s_lsb, cat_psi2s_res, cat_psi2s_rsb]
  n_lsb = (cat_psi2s_lsb.df.query(cat_psi2s_lsb.cuts)).shape[0]
  n_rsb = (cat_psi2s_rsb.df.query(cat_psi2s_rsb.cuts)).shape[0]
  n_res = (cat_psi2s_res.df.query(cat_psi2s_res.cuts)).shape[0]


  weight = fit_integrals[1]/(fit_integrals[0] + fit_integrals[1]) * n_res / (n_lsb + n_rsb)
  

  cat_psi2s_res_from_sbs = category("cat_psi2s_res_from_sbs",
      data_df,
      " or ".join([cat_psi2s_lsb.cuts, cat_psi2s_rsb.cuts]),
      "#psi(2S) bkg from sbs",
      ROOT.kOrange + 3,
      ROOT.kFullStar,
      weight)
      #fit_integrals[1]/(n_lsb + n_rsb))
  
  cat_psi2s_res_extrapolated_from_sb = category("psi2s_res_extrapolated_from_sb", 
      psi2s_res_extrapolated_from_sb_df,
      "Q_sq > -100",
      "#psi(2S) sbs extrapolated", 
      ROOT.kRed + 3,
      ROOT.kFullDiamond,
      weight)

  cat_psi2s_res_extrapolated_from_lsb = category("psi2s_res_extrapolated_from_lsb", 
      psi2s_res_extrapolated_from_lsb_df,
      "Q_sq > -100",
      "#psi(2S) lsb extrapolated", 
      ROOT.kAzure - 3,
      ROOT.kFullDiamond,
      weight)
  cat_psi2s_res_extrapolated_from_rsb = category("psi2s_res_extrapolated_from_rsb", 
      psi2s_res_extrapolated_from_rsb_df,
      "Q_sq > -100",
      "#psi(2S) rsb extrapolated", 
      ROOT.kPink - 3,
      ROOT.kFullDiamond,
      weight)

  cat_psi2s_res_extrapolated_from_lsb_mass_scaled = category("psi2s_res_extrapolated_from_lsb_mass_scaled", 
      psi2s_res_extrapolated_from_lsb_mass_scaled_df,
      "Q_sq > -100",
      "#psi(2S) lsb extrapolated (mass_scaled)", 
      ROOT.kAzure - 3,
      ROOT.kFullDiamond,
      weight)
  cat_psi2s_res_extrapolated_from_rsb_mass_scaled = category("psi2s_res_extrapolated_from_rsb_mass_scaled", 
      psi2s_res_extrapolated_from_rsb_mass_scaled_df,
      "Q_sq > -100",
      "#psi(2S) rsb extrapolated (mass_scaled)", 
      ROOT.kPink - 3,
      ROOT.kFullDiamond,
      weight)
  cats_psi2s_fromSB = [cat_psi2s_res, cat_psi2s_res_from_sbs, cat_psi2s_res_extrapolated_from_lsb, cat_psi2s_res_extrapolated_from_rsb]
  plotComparisonByCats(cats_psi2s_fromSB, var_list3, "psi2s_fromSBS", False)
  plotComparisonByCats(cats_psi2s_fromSB, var_list3, "psi2s_fromSBS_normalized", True)
  exit()


  cats_psi2s_mass_scaled_fromSB = [cat_psi2s_res, cat_psi2s_res_extrapolated_from_sb, cat_psi2s_res_extrapolated_from_lsb_mass_scaled, cat_psi2s_res_extrapolated_from_rsb_mass_scaled]
  plotComparisonByCats(cats_psi2s_mass_scaled_fromSB, var_list3, "psi2s_fromSBS_mass_scaled", False)
  plotComparisonByCats(cats_psi2s_mass_scaled_fromSB, var_list3, "psi2s_fromSBS_mass_scaled_normalized", True)

  var_psi2s_Q_sq = variable("Q_sq", "Q^{2}", "Q^{2}", "GeV^{2}", 80, 0., 12.)
  var_psi2s_lsb_mass =  variable("jpsi_mass", "m(#mu#mu)", "m(#mu#mu)", "GeV", 40, sb_intervals_psi2s['lsb']['minVal'], sb_intervals_psi2s['lsb']['maxVal'])
  var_psi2s_rsb_mass =  variable("jpsi_mass", "m(#mu#mu)", "m(#mu#mu)", "GeV", 40, sb_intervals_psi2s['rsb']['minVal'], sb_intervals_psi2s['rsb']['maxVal'])
  
  cats_jpsi2s_sb = [cat_psi2s_lsb,cat_psi2s_rsb]
  vars_x = [var_psi2s_lsb_mass,var_psi2s_rsb_mass]
  vars_y = [var_psi2s_Q_sq, var_psi2s_Q_sq]
  m1, b1, m2, b2, m3, b3= get_sb_fit_params(cats_jpsi2s_sb, vars_x, vars_y)

  dimuon_mass_lsb_weight_nominator = b2 + (res_middle/rsb_middle)*m2 * psi2s_lsb_df.jpsi_mass
  dimuon_mass_rsb_weight_nominator = b2 + m2 * (res_middle/rsb_middle)*psi2s_rsb_df.jpsi_mass
  dimuon_mass_lsb_weight_denominator = b1 + (res_middle/rsb_middle)*psi2s_lsb_df.jpsi_mass * m1
  dimuon_mass_rsb_weight_denominator = b3 + (res_middle/rsb_middle)*psi2s_rsb_df.jpsi_mass * m3
  dimuon_mass_lsb_weight = dimuon_mass_lsb_weight_nominator.divide(dimuon_mass_lsb_weight_denominator)
  print(dimuon_mass_lsb_weight)
  dimuon_mass_rsb_weight = dimuon_mass_rsb_weight_nominator.divide(dimuon_mass_rsb_weight_denominator)
  dimuon_lsb_rw_p4 = TLorentzVectorArray.from_ptetaphim(dimuon_mass_lsb_weight.multiply(psi2s_lsb_df.jpsi_pt), psi2s_lsb_df.jpsi_eta, psi2s_lsb_df.jpsi_phi, dimuon_mass_lsb_weight.multiply(psi2s_lsb_df.jpsi_mass))
  dimuon_rsb_rw_p4 = TLorentzVectorArray.from_ptetaphim(dimuon_mass_rsb_weight.multiply(psi2s_rsb_df.jpsi_pt), psi2s_rsb_df.jpsi_eta, psi2s_rsb_df.jpsi_phi, dimuon_mass_rsb_weight.multiply(psi2s_rsb_df.jpsi_mass))
  psi2s_res_extrapolated_from_lsb_mass_rw_df = pd.DataFrame()
  psi2s_res_extrapolated_from_rsb_mass_rw_df = pd.DataFrame()
  psi2s_res_extrapolated_from_lsb_mass_rw_df['Q_sq'] = (b_p4_lsb_corrected - dimuon_lsb_rw_p4).mag2
  psi2s_res_extrapolated_from_rsb_mass_rw_df['Q_sq'] = (b_p4_rsb_corrected - dimuon_rsb_rw_p4).mag2

  psi2s_res_extrapolated_from_sb_rw_df = pd.concat([psi2s_res_extrapolated_from_lsb_mass_rw_df, psi2s_res_extrapolated_from_rsb_mass_rw_df])

  cat_psi2s_res_extrapolated_from_lsb_mass_rw = category("psi2s_res_extrapolated_from_lsb_mass_rw", 
      psi2s_res_extrapolated_from_lsb_mass_rw_df,
      "Q_sq > -100",
      "#psi(2S) lsb reweighted", 
      ROOT.kAzure - 3,
      ROOT.kFullDiamond,
      weight)
  cat_psi2s_res_extrapolated_from_rsb_mass_rw= category("psi2s_res_extrapolated_from_rsb_mass_rw", 
      psi2s_res_extrapolated_from_rsb_mass_rw_df,
      "Q_sq > -100",
      "#psi(2S) rsb reweited", 
      ROOT.kPink - 3,
      ROOT.kFullDiamond,
      weight)
  cat_psi2s_res_extrapolated_from_rw_sb = category("psi2s_res_extrapolated_from_rw_sb", 
      psi2s_res_extrapolated_from_sb_rw_df,
      "Q_sq > -100",
      "#psi(2S) sbs reweited", 
      ROOT.kRed + 3,
      ROOT.kFullDiamond,
      weight)
  cats_psi2s_mass_scaled_fromSB_rw = [cat_psi2s_res, cat_psi2s_res_extrapolated_from_rw_sb, cat_psi2s_res_extrapolated_from_lsb_mass_rw, cat_psi2s_res_extrapolated_from_rsb_mass_rw]
  plotComparisonByCats(cats_psi2s_mass_scaled_fromSB_rw, var_list3, "psi2s_fromSBS_reweighted", False)
  #plotComparisonByCats(cats_psi2s_mass_scaled_fromSB_rw, var_list3, "psi2s_fromSBS_reweighted_normalized", True)

###### bkg from sidebands END ####

  plotComparisonByCats(cats_psi2s, var_list, "psi2s_sidebandsAndResonace", True)

  cats_jpsi_psi2s=[cat_psi2s_lsb, cat_psi2s_res, cat_psi2s_rsb, cat_jpsi_jpsitrkTrg]
  plotComparisonByCats(cats_jpsi_psi2s, var_list, "jpsi_psi2s", True)

  cats_jpsi_trgs=[cat_jpsi_jpsitrkTrg, cat_jpsi_dimuon0Trg]
  plotComparisonByCats(cats_jpsi_trgs, var_list, "jpsi_trgs", True)

  var_psi2s_mass =  [variable("jpsi_mass", "m(#mu#mu)", "m(#mu#mu)", "GeV", 40, 3.3, 4.)]
  cats_psi2s_all=[cat_psi2s, cat_psi2s_lsb, cat_psi2s_res, cat_psi2s_rsb]
  plotComparisonByCats(cats_psi2s_all, var_psi2s_mass, "psi2s_sideBandsAndResonace", False) 

if __name__ == '__main__':
  gROOT.SetBatch()
  gStyle.SetOptStat(0)
  main()
