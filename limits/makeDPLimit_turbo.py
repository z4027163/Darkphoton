import ROOT
import os,sys
from ROOT import TGraphAsymmErrors
from ROOT import TGraphErrors
from ROOT import TColor
#from ROOT import TGraph
from array import array
from ROOT import *
from operator import truediv
import random
import math
from glob import glob
import re 
import sys
from math import sqrt

year = sys.argv[1]

fff = open("eps2.txt", "w")
#fff2 = open("eps2_obs.txt", "w")


limit1=array('d')
limiteps2=array('d')
limit190=array('d')
limiteps290=array('d')
limit195up=array('d')
limit195down=array('d')
limit168up=array('d')
limit168down=array('d')
limit1Observed=array('d')
limit2=array('d')
limit2eps2=array('d')
limit290=array('d')
limit2eps290=array('d')
limit295up=array('d')
limit295down=array('d')
limit268up=array('d')
limit268down=array('d')
limit2Observed=array('d')
mass1=array('d')
mass2=array('d')
masserr1=array('d')
masserr2=array('d')

if year == "2018":
        lumi = 61.3
 #       lumi = 6.6
 #       lumi = 1.1


if year == "2017":
        lumi = 35.3

if year == "bothYears":
        lumi = 96.6

if year == "CL90":
        lumi = 96.6

#lumi_project = 96.6
lumi_project = lumi

#ACCEPTANCE
m_acc = array('d',[1.1, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5,8.5,9.5,15])
#NNLO
a_acc = array('d',[0.000457,0.001013,0.003603,0.01086,0.01441,0.01903,0.026,0.0369,0.07352,0.06841,0.1424])
#NLO
#a_acc = array('d',[0.000138,0.0002246,0.0005281,0.001168,0.002222,0.003816,0.006513,0.01196,0.03348,0.07715,0.1781])
#pythia
#a_acc = array('d',[0.00013,0.000205,0.0004193,0.001192,0.0028,0.005065,0.00891,0.01516,0.03266,0.0629,0.188])
accgraph = TGraph(len(m_acc),m_acc,a_acc);

#THEO CROSS SECTION FOR EPS=0.02
m= array('d',[1.1,2.0,  3.0,  4.0,  5.0,  6.0,  7.0,  8.0,   9.0,   10.0,  12.5,  20.0])
#xSec= array('d',[8541, 7514, 3323, 2055, 1422, 1043, 793.6, 621.1, 484.3, 292.8, 98.95])#[pb], model-dependent
#xSec= array('d',[26900,21800, 10990, 3894, 2166, 1421, 1013, 760.5, 589.6, 457.4, 275.1, 92.97])
xSec= array('d',[33500,27400,13100,5181,3126,1960,1346,975,738,559,327,112])

xsecgraph = TGraph(len(m),m,xSec);
eps2scale = 1.
base_eps = 0.02 #epsilon for which the cross sections are computed

a = (base_eps**2)*eps2scale/sqrt(lumi_project/lumi) # for lumi projection (6.6->100)

files = glob("combine_output/"+year+"/higgsCombineasympMassIndex_*.AsymptoticLimits.mH*.root")

d_m = {}
for fname in files:
        m = float(re.search(r"mH(\d+\.?\d+).root", fname).group(1))
        d = int(re.search(r"Index_(\d+).Asymptotic", fname).group(1))
        d_m[d] = [m, fname]


d_m = sorted(d_m.items())
print d_m

for d,m_fname in d_m:
    m, fname = m_fname
    #if (m == 1.716 or m ==7.632):
    #        continue
    #file90=glob.glob("higgsCombineIterV9_CL90_ForPress_2018_"+str(d)+".AsymptoticLimits.mH*.root")

    acc = accgraph.Eval(m,0,"S")
    xsec = xsecgraph.Eval(m,0,"S")


    if (m < 3):
            f=ROOT.TFile.Open(fname)
            tree=f.Get("limit")
            tree.GetEntry(2)
            limit1.append(tree.limit*a/(acc*xsec))
            a_exp=tree.limit*a/(acc*xsec)
             
            tree.GetEntry(0)
            limit195up.append(abs(tree.limit*a/(acc*xsec)-limit1[-1]))
            tree=f.Get("limit")
            tree.GetEntry(4)
            limit195down.append(abs(tree.limit*a/(acc*xsec)-limit1[-1]))
            
            
            tree.GetEntry(1)
            limit168up.append(abs(tree.limit*a/(acc*xsec)-limit1[-1]))
            tree=f.Get("limit")
            tree.GetEntry(3)
            limit168down.append(abs(tree.limit*a/(acc*xsec)-limit1[-1]))

            tree.GetEntry(5)
            limit1Observed.append(abs(tree.limit*a/(acc*xsec)))
            a_obs=tree.limit*a/(acc*xsec)            

            tree=f.Get("limit")
            mass1.append(m)
            masserr1.append(0.)
#            fff.write("{0} {1}\n".format(m, math.sqrt(tree.limit*a)))
            fff.write("{0} {1}\n".format(m, a_exp))
#            fff.write("{0} {1}\n".format(m, a_obs))
            
    if (m > 3):
            f=ROOT.TFile.Open(fname)
            tree=f.Get("limit")
            tree.GetEntry(2)
            limit2.append(tree.limit*a/(acc*xsec))
            a_exp=tree.limit*a/(acc*xsec)

            tree.GetEntry(0)
            limit295up.append(abs(tree.limit*a/(acc*xsec)-limit2[-1]))
            tree=f.Get("limit")
            tree.GetEntry(4)
            limit295down.append(abs(tree.limit*a/(acc*xsec)-limit2[-1]))
            
            
            tree.GetEntry(1)
            limit268up.append(abs(tree.limit*a/(acc*xsec)-limit2[-1]))
            tree=f.Get("limit")
            tree.GetEntry(3)
            limit268down.append(abs(tree.limit*a/(acc*xsec)-limit2[-1]))
            
            tree.GetEntry(5)
            limit2Observed.append(abs(tree.limit*a/(acc*xsec)))
            a_obs=tree.limit*a/(acc*xsec)           
 
            tree=f.Get("limit")
            mass2.append(m)
            masserr2.append(0.)
#            fff.write("{0} {1}\n".format(m, math.sqrt(tree.limit*a)))
            fff.write("{0} {1}\n".format(m, a_exp))
#            fff.write("{0} {1}\n".format(m, a_obs))   

c1=ROOT.TCanvas("c1","c1",700,500)
#c1.SetGrid()
c1.SetLogy()
c1.SetLogx()

mg=ROOT.TMultiGraph()
mgeps=ROOT.TMultiGraph()
graph_limit1=ROOT.TGraph(len(mass1),mass1,limit1)
graph_limit1.SetTitle("graph_limit1")
#graph_limit1.SetMarkerSize(1)
#graph_limit1.SetMarkerStyle(20)
#graph_limit1.SetMarkerColor(kBlack)
#graph_limit1.SetLineWidth(2)
#graph_limit1.SetLineStyle(7)
graph_limit1.SetLineColor(kRed+3)
graph_limit1.SetLineWidth(2)
graph_limit1.SetLineStyle(2)
graph_limit1.GetYaxis().SetRangeUser(0,100)
graph_limit1.GetXaxis().SetRangeUser(10,70)
graph_limit1.GetXaxis().SetMoreLogLabels()
graph_limit1.GetYaxis().SetTitle("#sigma(pp#rightarrow A)#times BR(A#rightarrow #mu#mu)[pb]")
graph_limit1.GetYaxis().SetTitleSize(2)
graph_limit1.GetXaxis().SetTitle("Dark Photon Mass [GeV]")

graph_limit2=ROOT.TGraph(len(mass2),mass2,limit2)
graph_limit2.SetTitle("graph_limit2")
#graph_limit2.SetMarkerSize(1)
#graph_limit2.SetMarkerStyle(20)
#graph_limit2.SetMarkerColor(kBlack)
#graph_limit2.SetLineWidth(2)
#graph_limit2.SetLineStyle(7)
graph_limit2.SetLineColor(kRed+3)
graph_limit2.SetLineWidth(2)
graph_limit2.SetLineStyle(2)

#graph_limit=ROOT.TGraph(len(mass),mass,limitObserved)
#graph_limit.Draw("same")
graph_limit195up=ROOT.TGraphAsymmErrors(len(mass1),mass1,limit1,masserr1,masserr1,limit195up,limit195down)
graph_limit195up.SetTitle("graph_limit195up")
graph_limit195up.SetFillColor(ROOT.TColor.GetColor(252,241,15))

graph_limit295up=ROOT.TGraphAsymmErrors(len(mass2),mass2,limit2,masserr2,masserr2,limit295up,limit295down)
graph_limit295up.SetTitle("graph_limit295up")
graph_limit295up.SetFillColor(ROOT.TColor.GetColor(252,241,15))


graph_limit168up=ROOT.TGraphAsymmErrors(len(mass1),mass1,limit1,masserr1,masserr1,limit168up,limit168down)
graph_limit168up.SetTitle("graph_limit168up")
graph_limit168up.SetMarkerColor(kGreen)
graph_limit168up.SetFillColor(kGreen);

graph_limitObs1=ROOT.TGraph(len(mass1),mass1,limit1Observed)
graph_limitObs1.SetTitle("graph_limitObs1")
graph_limitObs1.SetMarkerSize(1)
graph_limitObs1.SetMarkerStyle(1)
graph_limitObs1.SetMarkerColor(kBlack)
graph_limitObs1.SetLineWidth(2)
graph_limitObs1.SetLineStyle(1)

graph_limitObs2=ROOT.TGraph(len(mass2),mass2,limit2Observed)
graph_limitObs2.SetTitle("graph_limitObs2")
graph_limitObs2.SetMarkerSize(1)
graph_limitObs2.SetMarkerStyle(1)
graph_limitObs2.SetMarkerColor(kBlack)
graph_limitObs2.SetLineWidth(2)
graph_limitObs2.SetLineStyle(1)

graph_limit268up=ROOT.TGraphAsymmErrors(len(mass2),mass2,limit2,masserr2,masserr2,limit268up,limit268down)

graph_limit268up.SetTitle("graph_limit268up")
graph_limit268up.SetFillColor(kGreen);
graph_limit268up.SetMarkerColor(kGreen)

mg.Add(graph_limit195up,"3")
mg.Add(graph_limit168up,"3")
mg.Add(graph_limit1,"pl")
mg.Add(graph_limitObs1,"l")
mg.Add(graph_limit295up,"3")
mg.Add(graph_limit268up,"3")
mg.Add(graph_limit2,"pl")
mg.Add(graph_limitObs2,"l")

mg.Draw("APC")
mg.GetXaxis().SetRangeUser(1.2,9.)
mg.GetYaxis().SetRangeUser(1e-7,1e-2)
#mg.GetYaxis().SetTitle("xSec*BR [pb]")
#mg.GetXaxis().SetTitle("Dark Photon Mass [GeV]")
mg.GetYaxis().SetTitle("#epsilon^{2}")
mg.GetYaxis().SetTitleOffset(0.9)
mg.GetYaxis().SetTitleSize(0.05)
mg.GetXaxis().SetTitle("Dark Photon Mass [GeV]")
mg.GetXaxis().SetTitleSize(0.05)
mg.GetXaxis().SetMoreLogLabels()
c1.Update()
legend=ROOT.TLegend(0.5,0.6,0.8,0.9)
cmsTag=ROOT.TLatex(0.13,0.917,"#scale[1.1]{CMS}")
cmsTag.SetNDC()
cmsTag.SetTextAlign(11)
cmsTag.Draw()
cmsTag2=ROOT.TLatex(0.215,0.917,"#scale[0.825]{#bf{#it{Preliminary}}}")
cmsTag2.SetNDC()
cmsTag2.SetTextAlign(11)
#cmsTag.SetTextFont(61)
cmsTag2.Draw()
cmsTag3=ROOT.TLatex(0.90,0.917,"#scale[0.9]{#bf{"+str(lumi)+" fb^{-1} (13 TeV)}}")
cmsTag3.SetNDC()
cmsTag3.SetTextAlign(31)
#cmsTag.SetTextFont(61)
cmsTag3.Draw()
leg=ROOT.TLegend(0.65, 0.65,0.85, 0.85)  
leg.SetBorderSize( 0 )
leg.SetFillStyle( 1001 )
leg.SetFillColor(kWhite) 
leg.AddEntry( graph_limitObs1 , "Observed",  "LP" )
leg.AddEntry( graph_limit1 , "Expected",  "LP" )
leg.AddEntry( graph_limit168up, "#pm 1#sigma",  "F" ) 
leg.AddEntry( graph_limit195up, "#pm 2#sigma",  "F" ) 
leg.Draw("same")
c1.SaveAs("limit"+year+"DarkPhoton_eps2_turbo.root")
c1.SaveAs("limit"+year+"DarkPhoton_eps2_turbo.pdf")
c1.SaveAs("limit"+year+"DarkPhoton_eps2_turbo.png")
c2=ROOT.TCanvas("c2","c2",700,500)
c2.SetLogy()
cmsTag.Draw()
cmsTag2.Draw()
cmsTag3.Draw()
mgeps.Draw("APC")
leg2=ROOT.TLegend(0.65, 0.65,0.8, 0.85)  
leg2.SetBorderSize( 0 )
leg2.SetFillStyle( 1001 )
leg2.SetFillColor(kWhite) 

leg2.Draw("same")
c2.SaveAs("thep.root")

