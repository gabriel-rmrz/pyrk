import ROOT

class sample(object):
    def __init__(self, name, histoName, legendName, color):#, unit,nbins, xmin, xmax):
        self.name = name
        self.histoName = histoName
        self.legendName = legendName
        self.color=color
    
#ROOT.gStyle.SetPalette(ROOT.kRainBow)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)

mu = sample("mu", "mu", "B_{c} #rightarrow J/#psi + #mu  #nu", ROOT.kBlue-3)

tau = sample("tau", "tau", "B_{c} #rightarrow J/#psi + #tau_{#mu} #nu", ROOT.kGreen+2)
chic0 = sample("chic0", "chic0", "B_{c} #rightarrow #Chi_{c}^{0} + #mu_{#mu} #nu",ROOT.kRed-10)
chic1 = sample("chic1", "chic1", "B_{c} #rightarrow #Chi_{c}^{1} + #mu_{#mu} #nu", ROOT.kOrange+7)
chic2 = sample("chic2", "chic2", "B_{c} #rightarrow #Chi_{c}^{2} + #mu_{#mu} #nu", ROOT.kCyan+1)
hc_mu = sample("hc_mu", "hc_mu", "B_{c} #rightarrow h_{c} + #mu_{#mu} #nu", ROOT.kSpring-5)
jpsi_3pi = sample("jpsi_3pi", "jpsi_3pi", "B_{c} #rightarrow J/#psi + 3 #pi", ROOT.kPink)
jpsi_hc = sample("jpsi_hc", "jpsi_hc", "B_{c} #rightarrow J/#psi + D^{s}", ROOT.kViolet)
psi2s_mu = sample("psi2s_mu", "psi2s_mu", "B_{c} #rightarrow #psi(2S) + #mu_{#mu} #nu", ROOT.kGray)
psi2s_tau = sample("psi2s_tau", "psi2s_tau", "B_{c} #rightarrow #psi(2S) + #tau_{#mu} #nu", ROOT.kAzure)
misid = sample("mis_id", "mis_id", "misId bkg", ROOT.kYellow)
comb = sample("comb", "comb", "comb J/#psi + #mu", ROOT.kRed+1)
data = sample("data", "data_obs", "data", ROOT.kBlack)

#missing jpsi_3pi!
#missing psi2s tau
#sample_collection = [mu, tau,chic0,chic1,chic2,hc_mu,jpsi_hc,psi2s_mu,psi2s_tau, comb,  data,misid]
sample_collection = [mu, tau,chic0,chic1,chic2,hc_mu,jpsi_hc,psi2s_mu, comb,  data,misid]

sample_dic = {}

for collection in (sample_collection):
    sample_dic[collection.name] = collection
