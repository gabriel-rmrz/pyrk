import ROOT

class sample(object):
    def __init__(self, name, histoName, legendName, color):#, unit,nbins, xmin, xmax):
        self.name = name
        self.histoName = histoName
        self.legendName = legendName
        self.color=color

    
    
data = sample("data", "data_obs", "data", ROOT.kBlack)
mu = sample("mu", "mc_mu", "B_{c} #rightarrow J/#psi(#mu#mu) + #mu  #nu", ROOT.kRed) 
tau = sample("tau", "mc_tau", "B_{c} #rightarrow J/#psi(#mu#mu) + #tau_{#mu} #nu", ROOT.kGreen) 
misid = sample("mis_id", "mis_id", "misId bkg", ROOT.kBlue)
comb = sample("comb", "mc_comb", "combinatorial J/#psi(#mu#mu) + #mu", ROOT.kYellow) 


#inserisci altri campi, come il file hd5 e le selezioni diverse con i richiami ai dataframe salvati che si possono usare

sample_coll = [mu, tau,  misid, comb, data]
sampleColl = [mu, tau, comb, misid, data]
