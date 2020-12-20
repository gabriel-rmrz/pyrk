#macro that makes a stack plot with CMS requirements
import ROOT
from cmsstyle import CMS_lumi
#from sample import sample_collection
from sample import sample_dic
import os
import math
#in input:
#histo data, array of histos of mc, class var,  

def makeComparison(mch_tmp, var,path, addFile = ""):
    order = ["mu", "tau", "mis_id", "comb"]
    his = [histo for i,name in enumerate(order) for histo in mch_tmp if histo.GetName() == name ]
    sm = [i for i,name in enumerate(order) for histo in mch_tmp if histo.GetName() == name ]
    names = [name for i,name in enumerate(order) for histo in mch_tmp if histo.GetName() == name ]

    print(his)
    if(len(his) != len(mch_tmp)):
        print("ERROR: Problem in histos names!")
        return 0
    
    maxy = []
    for i,histo,sname in zip(sm,his,names):
        print(i,histo)
        histo.SetFillColor(0)
        histo.SetLineColor(sample_dic[sname].color)
        histo.GetXaxis().SetRangeUser(var.xmin,var.xmax)
        histo.Scale(1/histo.Integral()) #drawnormalized
        maxy.append(histo.GetMaximum())

    maximum = max(maxy)
    his[0].SetMaximum(maximum*1.35)
    his[0].SetMinimum(0.)

    #legend
    legend = ROOT.TLegend(0.57,0.70,0.90,0.86)
    for i,histo,sname in zip(sm,his,names):
        legend.AddEntry(histo, sample_dic[sname].legendName, "f")
    legend.SetTextFont(43)
    legend.SetTextSize(15)
    legend.SetBorderSize(0)
    legend.SetFillColor(0)

    ########
    #canvas#
    ########
    c = ROOT.TCanvas("","",700, 700) 
    c.cd() 

    c.SetTicks(True)
    c.SetBottomMargin(2)
    c.SetLeftMargin(0.15)
    c.SetRightMargin(0.15)
    c.Draw()
    c.cd()

    his[0].Draw("hist")
    for i in range(1, len(his)):
        his[i].Draw("histSAME")
    legend.Draw("SAME")

    his[0].GetYaxis().SetTitle("Events / bin")
    his[0].GetYaxis().SetLabelOffset(0.01)
    his[0].GetYaxis().SetTitleOffset(2)
    his[0].GetYaxis().SetLabelFont(43)
    his[0].GetYaxis().SetLabelSize(15)
    his[0].GetYaxis().SetTitleFont(43)
    his[0].GetYaxis().SetTitleSize(18)

    his[0].GetXaxis().SetLabelOffset(0.01)
    his[0].GetXaxis().SetTitleOffset(1.6)
    his[0].GetXaxis().SetLabelFont(43)
    his[0].GetXaxis().SetLabelSize(15)
    his[0].GetXaxis().SetTitleFont(43)
    his[0].GetXaxis().SetTitleSize(18)
    his[0].GetXaxis().SetTitle(var.xlabel + " " + var.unit)
    
    CMS_lumi(c, 4, 1)
    c.RedrawAxis()

    folder = path
    c.SaveAs(folder +"BDT_output/"+ var.name + "_comparison_"   + addFile + ".png")
    c.SaveAs(folder + "BDT_output/"+ var.name + "_comparison_"  + addFile + ".pdf")
    


def makeSinglePlot(his,var,smpl, addFile="", saveFile = True, retur = False):
    ROOT.gROOT.SetBatch()
    ROOT.gStyle.SetOptStat(0)
    his.SetMarkerStyle(20)
    his.SetMarkerSize(0.9)
    his.SetLineColor(ROOT.kBlack)
    #legend
    legend = ROOT.TLegend(0.51,0.70,0.78,0.86)
    legend.AddEntry(his, smpl.legendName, "lp")
    legend.SetTextFont(43)
    legend.SetTextSize(15)
    legend.SetBorderSize(0)
    legend.SetFillColor(0)

    ########
    #canvas#
    ########
    c = ROOT.TCanvas("","",700, 700)

    c.cd() 
    c.SetTicks(True)
    c.SetBottomMargin(2)
    c.SetLeftMargin(0.15)
    c.SetRightMargin(0.15)
    c.Draw()
    c.cd()

    his.Draw("e")
    legend.Draw("SAME")
    his.SetTitle("")
    his.GetYaxis().SetTitle("Events / bin")
    his.GetYaxis().SetLabelOffset(0.01)
    his.GetYaxis().SetTitleOffset(2)
    his.GetYaxis().SetLabelFont(43)
    his.GetYaxis().SetLabelSize(15)
    his.GetYaxis().SetTitleFont(43)
    his.GetYaxis().SetTitleSize(18)
    his.SetMaximum(his.GetMaximum()*1.5)

    his.GetXaxis().SetLabelOffset(0.01)
    his.GetXaxis().SetTitleOffset(1.6)
    his.GetXaxis().SetLabelFont(43)
    his.GetXaxis().SetLabelSize(15)
    his.GetXaxis().SetTitleFont(43)
    his.GetXaxis().SetTitleSize(18)
    his.GetXaxis().SetTitle(var.xlabel + " " + var.unit)
    CMS_lumi(c, 4, 1)
    c.RedrawAxis()
    
    if saveFile:
        folder = "pngPlots/"
        c.SaveAs(folder + var.name + "_" + smpl.name + addFile + ".png")
        c.SaveAs(folder + var.name + "_" + smpl.name + addFile + ".pdf")
        print("Saved '"+folder + var.name + "_" + smpl.name + addFile + ".png'")
    if (retur == True):
        return his

def makeStack_old(datah, mch_tmp, var, path , fit=False,  addFileName = "", rootFile = False, ratioPad = True, over = False, fact = 1):
    os.system('set LD_PRELOAD=libtcmalloc.so')
    ROOT.gROOT.SetBatch()
    order = ["mc_mu", "mc_tau", "mis_id", "mc_comb"]
    his = [histo for i,name in enumerate(order) for histo in mch_tmp if histo.GetName() == name ]
    sm = [i for i,name in enumerate(order) for histo in mch_tmp if histo.GetName() == name ]
    if(len(his) != len(mch_tmp)):
        print("ERROR: Problem in histos names!")
        return 0

    stack= ROOT.THStack("","")
    for i,histo,name in zip(sm,his,order):
        if name == 'mc_tau':  #to multiply the over and underflow for the factor to have tau = mu
            f = fact
        else:
            f = 1
        if over:
            histo.SetBinContent(1, histo.GetBinContent(0) + histo.GetBinContent(1))
            histo.SetBinError(1, math.sqrt(pow(histo.GetBinError(0),2) + pow(histo.GetBinError(1),2)))
            histo.SetBinContent(var.nbins, histo.GetBinContent(var.nbins) + histo.GetBinContent(var.nbins+1))
            histo.SetBinError(var.nbins, math.sqrt(pow(histo.GetBinError(var.nbins),2) + pow(histo.GetBinError(var.nbins+1),2)))

        histo.Scale(f)
        histo.SetFillColor(sample_coll[i].color)
        histo.SetLineColor(sample_coll[i].color)
        histo.GetXaxis().SetRangeUser(var.xmin,var.xmax)
        stack.Add(histo)

    #stack
    datah.SetMarkerStyle(20)
    datah.SetMarkerSize(0.9)
    datah.SetBinContent(1, datah.GetBinContent(0) + datah.GetBinContent(1))
    datah.SetBinError(1, math.sqrt(pow(datah.GetBinError(0),2) + pow(datah.GetBinError(1),2)))
    datah.SetBinContent(var.nbins, datah.GetBinContent(var.nbins) + datah.GetBinContent(var.nbins+1))
    datah.SetBinError(var.nbins, math.sqrt(pow(datah.GetBinError(var.nbins),2) + pow(datah.GetBinError(var.nbins+1),2)))
    
    maximum = max(datah.GetMaximum(), stack.GetMaximum())
    stack.SetMaximum(maximum*1.35)
    stack.SetMinimum(0.)

    #legend
    legend = ROOT.TLegend(0.45,0.70,0.78,0.86)
    legend.AddEntry(datah, sample_coll[-1].legendName, "lp")
    for i,histo in zip(sm,his):
        legend.AddEntry(histo, sample_coll[i].legendName, "f")
    legend.SetTextFont(43)
    legend.SetTextSize(15)
    legend.SetBorderSize(0)
    legend.SetFillColor(0)

    ########
    #canvas#
    ########

    if not ratioPad:
        main_pad = ROOT.TCanvas("","",700, 700) 
        main_pad.SetBottomMargin(2)

    if ratioPad:
        c = ROOT.TCanvas("","",700, 700)                                                                              
        c.cd()
        ############
        # MAIN PAD #
        ############
        main_pad = ROOT.TPad("","",0.  , 0.25, 1. , 1.  )
        main_pad.SetBottomMargin(0.02)
    
    main_pad.SetTicks(True)
    main_pad.SetLeftMargin(0.15)
    main_pad.SetRightMargin(0.15)
    main_pad.Draw()
    main_pad.cd()
        
    stack.Draw("hist")
    legend.Draw("SAME")

    stack.GetYaxis().SetTitle("Events / bin")
    stack.GetYaxis().SetLabelOffset(0.01)
    stack.GetYaxis().SetTitleOffset(2)
    stack.GetYaxis().SetLabelFont(43)
    stack.GetYaxis().SetLabelSize(15)
    stack.GetYaxis().SetTitleFont(43)
    stack.GetYaxis().SetTitleSize(18)

    if not ratioPad:
        stack.GetXaxis().SetLabelOffset(0.01)
    if ratioPad:
        stack.GetXaxis().SetLabelOffset(0.8)

    stack.GetXaxis().SetTitleOffset(1.6)
    stack.GetXaxis().SetLabelFont(43)
    stack.GetXaxis().SetLabelSize(15)
    stack.GetXaxis().SetTitleFont(43)
    stack.GetXaxis().SetTitleSize(18)
    stack.GetXaxis().SetTitle(var.xlabel + " " + var.unit)
    stack.SetMinimum(0.1)
    #ROOT.gPad.SetLogy()

    # Stat. Errors
    h_err = stack.GetStack().Last().Clone("h_err")
    h_err.SetLineWidth(100)
    h_err.SetFillStyle(3154)
    h_err.SetMarkerSize(0)
    h_err.SetFillColor(ROOT.kGray+2)
    h_err.Draw("e2same0")
    if fit:
        legend.AddEntry(h_err, "Stat. + Syst. Unc.", "f")
    else:
        legend.AddEntry(h_err, "Stat. Unc.", "f")

    datah.Draw("E same")


    CMS_lumi(main_pad, 4, 1)
    main_pad.RedrawAxis()

    if ratioPad:
    
        #############
        # RATIO PAD #
        #############
        c.cd() 
        ratio_pad = ROOT.TPad("","",0.  , 0.  , 1. , 0.25)
        
        ratio_pad.SetLeftMargin(0.15)
        ratio_pad.SetRightMargin(0.15)
        ratio_pad.SetTopMargin(0.05)   
        ratio_pad.SetBottomMargin(0.5)
        ratio_pad.SetGridy()
        c.cd()
        ratio_pad.Draw()
        ratio_pad.cd()

        
        #f1 = ROOT.TLine(0, 1., 1, 1.)
        
        f1 = ROOT.TLine(var.xmin, 1., var.xmax, 1.)
        f1.SetLineColor(ROOT.kBlack)
        f1.SetLineStyle(ROOT.kDashed)
        f1.Draw("same")
        
        hratio = stack.GetStack().Last()
        ratio = datah.Clone("ratio")
        ratio.SetLineColor(ROOT.kBlack)
        ratio.SetMaximum(1.2)
        ratio.SetMinimum(0.8)
        ratio.Sumw2()
        ratio.SetStats(0)
    
        ratio.Divide(hratio)
        ratio.SetMarkerStyle(20)
        ratio.SetMarkerSize(0.9)
        ratio.Draw("epx0e0")
        ratio.SetTitle("")


        ratio.GetYaxis().SetTitle("Data / MC")
        ratio.GetYaxis().SetNdivisions(503)
        ratio.GetYaxis().SetLabelFont(43)
        ratio.GetYaxis().SetTitleFont(43)
        ratio.GetYaxis().SetTitleOffset(1.5)
        ratio.GetYaxis().SetLabelSize(15)
        ratio.GetYaxis().SetTitleSize(18)
        ratio.GetYaxis().SetRangeUser(0.8, 1.2)
        ratio.GetYaxis().SetLabelOffset(0.01)
        
        ratio.GetXaxis().SetLabelFont(43)
        ratio.GetXaxis().SetTitleFont(43)
        ratio.GetXaxis().SetTitleOffset(4)
        ratio.GetXaxis().SetLabelSize(15)
        ratio.GetXaxis().SetTitleSize(18)
        ratio.GetXaxis().SetRangeUser(var.xmin, var.xmax)
        ratio.GetXaxis().SetTitle(var.xlabel + " " + var.unit)
        ratio.GetXaxis().SetLabelOffset(0.04)
        
        h_bkg_err = hratio.Clone("h_err")
        h_bkg_err.Reset()
        h_bkg_err.Sumw2()
    
        h_bkg_err.SetLineWidth(100)
        
        for i in range(1,hratio.GetNbinsX()+1):
            h_bkg_err.SetBinContent(i,1)
            if(hratio.GetBinContent(i)):
                h_bkg_err.SetBinError(i, (hratio.GetBinError(i)/hratio.GetBinContent(i)))
            else:
                h_bkg_err.SetBinError(i, 10^(-99))
        h_bkg_err.SetLineWidth(100)
            
        h_bkg_err.SetMarkerSize(0)
        h_bkg_err.SetFillColor(ROOT.kGray+1)
        h_bkg_err.SetFillStyle(3154)
        h_bkg_err.Draw("e20same")
        
        ratio.Draw("epx0e0same")
        
    #ROOT.TGaxis.SetMaxDigits(3)
    if ratioPad:
        c.Update()
    if not ratioPad:    
        main_pad.Update()

    if (rootFile == True):
        f = ROOT.TFile.Open("rootFiles/"+ var.name + addFileName +  ".root", "RECREATE")
        f.cd()
        for h in his:
            h.Write()
        datah.SetName("data_obs")
        
        datah.Write()
        f.Close()
        print("Saved rootFiles/"+ var.name + addFileName +  ".root")
    if fit:
        folder = "fitPlt/"
    else:
        folder = path + "/"

    if not ratioPad:
        main_pad.SaveAs(folder + var.name  +"_"+ addFileName+ ".png")
        main_pad.SaveAs(folder + var.name  + "_"+ addFileName+".pdf")
        
        
    if ratioPad:
        c.SaveAs(folder + var.name + "_"+ addFileName+ ".png")
        c.SaveAs(folder + var.name + "_"+ addFileName+ ".pdf")


    stack.Delete()
    legend.Delete()
    h_err.Delete()
    main_pad.Delete()
    if ratioPad:
        ratio.Delete()
        f1.Delete()
        h_bkg_err.Delete()
        ratio_pad.Delete()
        c.Delete()

def makeStack(datah, mch_tmp, var, path , fit=False,  addFileName = "", rootFile = False, ratioPad = True, over = False,rootPath = "rootFiles"):
    os.system('set LD_PRELOAD=libtcmalloc.so')
    ROOT.gROOT.SetBatch()
    #order = ['mu','tau','chic0','chic1','chic2','hc_mu','jpsi_hc','psi2s_mu','psi2s_tau','comb','mis_id']
    order = ['mu','tau','chic0','chic1','chic2','hc_mu','jpsi_hc','psi2s_mu','comb','mis_id']
    #    order = ["mc_mu", "mc_tau", "mis_id", "mc_comb"]
    his = [histo for i,name in enumerate(order) for histo in mch_tmp if histo.GetName() == name ]
    sm = [i for i,name in enumerate(order) for histo in mch_tmp if histo.GetName() == name ]
    if(len(his) != len(mch_tmp)):
        print("ERROR: Problem in histos names!")
        return 0
   
    stack= ROOT.THStack("","")
    for i,histo,sname in zip(sm,his,order):
        '''
        if over:
            histo.SetBinContent(1, histo.GetBinContent(0) + histo.GetBinContent(1))
            histo.SetBinError(1, math.sqrt(pow(histo.GetBinError(0),2) + pow(histo.GetBinError(1),2)))
            histo.SetBinContent(var.nbins, histo.GetBinContent(var.nbins) + histo.GetBinContent(var.nbins+1))
            histo.SetBinError(var.nbins, math.sqrt(pow(histo.GetBinError(var.nbins),2) + pow(histo.GetBinError(var.nbins+1),2)))
        '''
        #histo.Scale(f)
        histo.SetFillColor(sample_dic[sname].color)
        histo.SetLineColor(sample_dic[sname].color)
        histo.GetXaxis().SetRangeUser(var.xmin,var.xmax)
        stack.Add(histo)
    '''
    if over: #for data!
        datah.SetBinContent(1, datah.GetBinContent(0) + datah.GetBinContent(1))
        datah.SetBinError(1, math.sqrt(pow(datah.GetBinError(0),2) + pow(datah.GetBinError(1),2)))
        datah.SetBinContent(var.nbins, datah.GetBinContent(var.nbins) + datah.GetBinContent(var.nbins+1))
        datah.SetBinError(var.nbins, math.sqrt(pow(datah.GetBinError(var.nbins),2) + pow(datah.GetBinError(var.nbins+1),2)))
    '''

    #stack
    datah.SetMarkerStyle(20)
    datah.SetMarkerSize(0.9)
    '''
    datah.SetBinContent(1, datah.GetBinContent(0) + datah.GetBinContent(1))
    datah.SetBinError(1, math.sqrt(pow(datah.GetBinError(0),2) + pow(datah.GetBinError(1),2)))
    datah.SetBinContent(var.nbins, datah.GetBinContent(var.nbins) + datah.GetBinContent(var.nbins+1))
    datah.SetBinError(var.nbins, math.sqrt(pow(datah.GetBinError(var.nbins),2) + pow(datah.GetBinError(var.nbins+1),2)))
    '''
    maximum = max(datah.GetMaximum(), stack.GetMaximum())
    stack.SetMaximum(maximum*1.5)
    stack.SetMinimum(0.)
    print("INTEGRALE DATA",datah.Integral())
    #legend
    legend = ROOT.TLegend(0.18,0.65,0.85,0.81)
    legend.SetNColumns(3);
    legend.AddEntry(datah, sample_dic["data"].legendName, "lp")
    for histo,sname in zip(his,order):
        legend.AddEntry(histo, sample_dic[sname].legendName, "f")
    legend.SetTextFont(43)
    legend.SetTextSize(15)
    legend.SetBorderSize(0)
    legend.SetFillColor(0)

    ########
    #canvas#
    ########

    if not ratioPad:
        main_pad = ROOT.TCanvas("","",700, 700) 
        main_pad.SetBottomMargin(2)

    if ratioPad:
        c = ROOT.TCanvas("","",700, 700)                                                                              
        c.cd()
        ############
        # MAIN PAD #
        ############
        main_pad = ROOT.TPad("","",0.  , 0.25, 1. , 1.  )
        main_pad.SetBottomMargin(0.02)
    
    main_pad.SetTicks(True)
    main_pad.SetLeftMargin(0.15)
    main_pad.SetRightMargin(0.15)
    main_pad.Draw()
    main_pad.cd()
        
    stack.Draw("hist")
    legend.Draw("SAME")

    stack.GetYaxis().SetTitle("Events / bin")
    stack.GetYaxis().SetLabelOffset(0.01)
    stack.GetYaxis().SetTitleOffset(2)
    stack.GetYaxis().SetLabelFont(43)
    stack.GetYaxis().SetLabelSize(15)
    stack.GetYaxis().SetTitleFont(43)
    stack.GetYaxis().SetTitleSize(18)

    if not ratioPad:
        stack.GetXaxis().SetLabelOffset(0.01)
    if ratioPad:
        stack.GetXaxis().SetLabelOffset(0.8)

    stack.GetXaxis().SetTitleOffset(1.6)
    stack.GetXaxis().SetLabelFont(43)
    stack.GetXaxis().SetLabelSize(15)
    stack.GetXaxis().SetTitleFont(43)
    stack.GetXaxis().SetTitleSize(18)
    stack.GetXaxis().SetTitle(var.xlabel + " " + var.unit)
    stack.SetMinimum(0.1)
    #ROOT.gPad.SetLogy()

    # Stat. Errors
    h_err = stack.GetStack().Last().Clone("h_err")
    h_err.SetLineWidth(100)
    h_err.SetFillStyle(3154)
    h_err.SetMarkerSize(0)
    h_err.SetFillColor(ROOT.kGray+2)
    h_err.Draw("e2same0")
    if fit:
        legend.AddEntry(h_err, "Stat. + Syst. Unc.", "f")
    else:
        legend.AddEntry(h_err, "Stat. Unc.", "f")

    datah.Draw("E same")


    CMS_lumi(main_pad, 4, 1)
    main_pad.RedrawAxis()

    if ratioPad:
    
        #############
        # RATIO PAD #
        #############
        c.cd() 
        ratio_pad = ROOT.TPad("","",0.  , 0.  , 1. , 0.25)
        
        ratio_pad.SetLeftMargin(0.15)
        ratio_pad.SetRightMargin(0.15)
        ratio_pad.SetTopMargin(0.05)   
        ratio_pad.SetBottomMargin(0.5)
        ratio_pad.SetGridy()
        c.cd()
        ratio_pad.Draw()
        ratio_pad.cd()

        
        #f1 = ROOT.TLine(0, 1., 1, 1.)
        
        f1 = ROOT.TLine(var.xmin, 1., var.xmax, 1.)
        f1.SetLineColor(ROOT.kBlack)
        f1.SetLineStyle(ROOT.kDashed)
        f1.Draw("same")
        
        hratio = stack.GetStack().Last()
        ratio = datah.Clone("ratio")
        ratio.SetLineColor(ROOT.kBlack)
        ratio.SetMaximum(1.2)
        ratio.SetMinimum(0.8)
        ratio.Sumw2()
        ratio.SetStats(0)
    
        ratio.Divide(hratio)
        ratio.SetMarkerStyle(20)
        ratio.SetMarkerSize(0.9)
        ratio.Draw("epx0e0")
        ratio.SetTitle("")


        ratio.GetYaxis().SetTitle("Data / MC")
        ratio.GetYaxis().SetNdivisions(503)
        ratio.GetYaxis().SetLabelFont(43)
        ratio.GetYaxis().SetTitleFont(43)
        ratio.GetYaxis().SetTitleOffset(1.5)
        ratio.GetYaxis().SetLabelSize(15)
        ratio.GetYaxis().SetTitleSize(18)
        ratio.GetYaxis().SetRangeUser(0.8, 1.2)
        ratio.GetYaxis().SetLabelOffset(0.01)
        
        ratio.GetXaxis().SetLabelFont(43)
        ratio.GetXaxis().SetTitleFont(43)
        ratio.GetXaxis().SetTitleOffset(4)
        ratio.GetXaxis().SetLabelSize(15)
        ratio.GetXaxis().SetTitleSize(18)
        ratio.GetXaxis().SetRangeUser(var.xmin, var.xmax)
        ratio.GetXaxis().SetTitle(var.xlabel + " " + var.unit)
        ratio.GetXaxis().SetLabelOffset(0.04)
        
        h_bkg_err = hratio.Clone("h_err")
        h_bkg_err.Reset()
        h_bkg_err.Sumw2()
    
        h_bkg_err.SetLineWidth(100)
        
        for i in range(1,hratio.GetNbinsX()+1):
            h_bkg_err.SetBinContent(i,1)
            if(hratio.GetBinContent(i)):
                h_bkg_err.SetBinError(i, (hratio.GetBinError(i)/hratio.GetBinContent(i)))
            else:
                h_bkg_err.SetBinError(i, 10^(-99))
        h_bkg_err.SetLineWidth(100)
            
        h_bkg_err.SetMarkerSize(0)
        h_bkg_err.SetFillColor(ROOT.kGray+1)
        h_bkg_err.SetFillStyle(3154)
        h_bkg_err.Draw("e20same")
        
        ratio.Draw("epx0e0same")
        
    #ROOT.TGaxis.SetMaxDigits(3)
    if ratioPad:
        c.Update()
    if not ratioPad:    
        main_pad.Update()

    if (rootFile == True):
        f = ROOT.TFile.Open(rootPath+"/"+ var.name + addFileName +  ".root", "RECREATE")
        f.cd()
        for h in his:
            h.Write()
        datah.SetName("data_obs")
        
        datah.Write()
        f.Close()
        print("Saved "+rootPath+"/"+ var.name + addFileName +  ".root")
    if fit:
        folder = path + "/"
    else:
        folder = path + "/"

    if not ratioPad:
        main_pad.SaveAs(folder + var.name  +"_"+ addFileName+ ".png")
        main_pad.SaveAs(folder + var.name  + "_"+ addFileName+".pdf")
        
        
    if ratioPad:
        c.SaveAs(folder + var.name + "_"+ addFileName+ ".png")
        c.SaveAs(folder + var.name + "_"+ addFileName+ ".pdf")


    stack.Delete()
    legend.Delete()
    h_err.Delete()
    main_pad.Delete()
    if ratioPad:
        ratio.Delete()
        f1.Delete()
        h_bkg_err.Delete()
        ratio_pad.Delete()
        c.Delete()

