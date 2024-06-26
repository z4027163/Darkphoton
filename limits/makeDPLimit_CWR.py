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

fff.write("mass  -2              -1              0              +1              +2              Obs\n")


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

#xsec=1.#pb


#some counters
#m=0.2
#t=0
#make the loop

#a=10#1./0.1143
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

print(lumi)
#ACCEPTANCE
acc_file = TFile.Open("acceptances_dy.root")
acc_teff = acc_file.Get("cmsacc")
nbins_acc = acc_teff.GetPassedHistogram().GetNbinsX()
acceptances = array('d')
m_acceptances = array('d')
for j in range(nbins_acc):
	acceptances.append(acc_teff.GetEfficiency(j+1))
	m_acceptances.append(acc_teff.GetPassedHistogram().GetBinCenter(j+1))
accgraph = TGraph(nbins_acc,m_acceptances,acceptances);

#THEO CROSS SECTION FOR EPS=0.02
eps2scale = 1.

a = eps2scale/sqrt(lumi_project/lumi) # for lumi projection (6.6->100)

#dark photon
#files = glob("dual/reso0p2/combine_output/"+year+"/higgsCombineasympMassIndex_*.AsymptoticLimits.mH*.root")
#2hdm
#files = glob("ptcut/IPsig/combine_output/"+year+"/higgsCombineasympMassIndex_*.AsymptoticLimits.mH*.root")
files = glob("combine_output/"+year+"/higgsCombineasympMassIndex_*.AsymptoticLimits.mH*.root")

d_m = {}
for fname in files:
        print(fname)
        m = float(re.search(r"mH(\d+\.?\d+).root", fname).group(1))
        d = int(re.search(r"Index_(\d+).Asymptotic", fname).group(1))
	d_m[d] = [m, fname]


d_m = sorted(d_m.items())
print d_m

for d,m_fname in d_m:
	m, fname = m_fname
        ##if (m == 1.716 or m ==7.632):
        ##        continue
	#file90=glob.glob("higgsCombineIterV9_CL90_ForPress_2018_"+str(d)+".AsymptoticLimits.mH*.root")

	acc = accgraph.Eval(m)
        if (m < 3):
                f=ROOT.TFile.Open(fname)
                tree=f.Get("limit")
                tree.GetEntry(2)
                limit1.append(tree.limit*a)
                print("Mass: " + str(m) + "         xSec "+ str(tree.limit*a))
                a2=tree.limit*a  
               
                tree.GetEntry(0)
                limit195up.append(abs(tree.limit*a-limit1[-1]))
                a0=tree.limit*a

                tree=f.Get("limit")
                tree.GetEntry(4)
                limit195down.append(abs(tree.limit*a-limit1[-1]))
                a4=tree.limit*a

                
                tree.GetEntry(1)
                limit168up.append(abs(tree.limit*a-limit1[-1]))
                a1=tree.limit*a
 
                tree=f.Get("limit")
                tree.GetEntry(3)
                limit168down.append(abs(tree.limit*a-limit1[-1]))
                a3=tree.limit*a

                tree.GetEntry(5)
		limit1Observed.append(abs(tree.limit*a))
                tree=f.Get("limit")
                #print("Mass: " + str(m) + "         xSec "+ str(tree.limit*a))



                mass1.append(m)
                masserr1.append(0.)
                fff.write("{0} {1} {2} {3} {4} {5} {6}\n".format(m, a0,a1,a2,a3,a4,tree.limit*a))
                #exp
                #fff.write("{0} {1}\n".format(m, a2))
                #obs
                #fff.write("{0} {1}\n".format(m, tree.limit*a))

        if (m > 3):
                f=ROOT.TFile.Open(fname)
                tree=f.Get("limit")
                tree.GetEntry(2)
                limit2.append(tree.limit*a)
                a2=tree.limit*a 
                print("Mass: " + str(m) + "         xSec "+ str(tree.limit*a))
               
                tree.GetEntry(0)
                limit295up.append(abs(tree.limit*a-limit2[-1]))
                a0=tree.limit*a

                tree=f.Get("limit")
                tree.GetEntry(4)
                limit295down.append(abs(tree.limit*a-limit2[-1]))
                a4=tree.limit*a                
                
                tree.GetEntry(1)
                limit268up.append(abs(tree.limit*a-limit2[-1]))
                a1=tree.limit*a

                tree=f.Get("limit")
                tree.GetEntry(3)
                limit268down.append(abs(tree.limit*a-limit2[-1]))
                a3=tree.limit*a                

                tree.GetEntry(5)
                limit2Observed.append(abs(tree.limit*a))
                #limit2Observed.append(abs(tree.limit*a-limit2[-1]))
                tree=f.Get("limit")
		
                mass2.append(m)
                masserr2.append(0.)
                fff.write("{0} {1} {2} {3} {4} {5} {6}\n".format(m, a0,a1,a2,a3,a4,tree.limit*a)) 
                #exp
                #fff.write("{0} {1}\n".format(m, a2))
                #obs
                #fff.write("{0} {1}\n".format(m, tree.limit*a))

print limit1

c1=ROOT.TCanvas("c1","c1",700,500)
#c1.SetGrid()
c1.SetLogy()
c1.SetLogx()
c1.SetBottomMargin(0.13);

mg=ROOT.TMultiGraph()
mgeps=ROOT.TMultiGraph()
graph_limit1=ROOT.TGraph(len(mass1),mass1,limit1)
graph_limit1.SetTitle("graph_limit1")
#graph_limit1.SetMarkerSize(1)
#graph_limit1.SetMarkerStyle(20)
#graph_limit1.SetMarkerColor(kBlack)
graph_limit1.SetLineColor(kRed+3)
graph_limit1.SetLineWidth(2)
graph_limit1.SetLineStyle(2)
graph_limit1.GetYaxis().SetRangeUser(0,100)
graph_limit1.GetXaxis().SetRangeUser(10,70)
graph_limit1.GetXaxis().SetMoreLogLabels()
#graph_limit1.GetYaxis().SetTitle("#sigma(pp#rightarrow A)#times BR(A#rightarrow #mu#mu)[pb]")
#briefing
graph_limit1.GetYaxis().SetTitle("#sigma(pp#rightarrow A)#times B(A#rightarrow #mu#mu)[pb]")

graph_limit1.GetYaxis().SetTitleSize(2)
graph_limit1.GetXaxis().SetTitle("Dark Photon Mass [GeV]")

graph_limit2=ROOT.TGraph(len(mass2),mass2,limit2)
graph_limit2.SetTitle("graph_limit2")
#graph_limit2.SetMarkerSize(1)
#graph_limit2.SetMarkerStyle(20)
#graph_limit2.SetMarkerColor(kBlack)
graph_limit2.SetLineWidth(2)
graph_limit2.SetLineStyle(2)
graph_limit2.SetLineColor(kRed+3)

#graph_limit=ROOT.TGraph(len(mass),mass,limitObserved)
#graph_limit.Draw("same")
graph_limit195up=ROOT.TGraphAsymmErrors(len(mass1),mass1,limit1,masserr1,masserr1,limit195up,limit195down)
graph_limit195up.SetTitle("graph_limit195up")
graph_limit195up.SetFillColor(ROOT.TColor.GetColor(252,241,15))

graph_limit295up=ROOT.TGraphAsymmErrors(len(mass2),mass2,limit2,masserr2,masserr2,limit295up,limit295down)
graph_limit295up.SetTitle("graph_limit295up")
graph_limit295up.SetFillColor(ROOT.TColor.GetColor(252,241,15))

graph_limit168up=ROOT.TGraphAsymmErrors(len(mass1),mass1,limit1,masserr1,masserr1,limit168up,limit168down)
graph_limit168up.SetTitle("graph_limit68up")
graph_limit168up.SetFillColor(kGreen);
graph_limit168up.SetMarkerColor(kGreen)

graph_limit268up=ROOT.TGraphAsymmErrors(len(mass2),mass2,limit2,masserr2,masserr2,limit268up,limit268down)
graph_limit268up.SetTitle("graph_limit68up")
graph_limit268up.SetFillColor(kGreen);
graph_limit268up.SetMarkerColor(kGreen)

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

mg.Add(graph_limit195up,"3")
mg.Add(graph_limit168up,"3")
mg.Add(graph_limit1,"l")
mg.Add(graph_limitObs1,"l")
mg.Add(graph_limit295up,"3")
mg.Add(graph_limit268up,"3")
mg.Add(graph_limit2,"l")
mg.Add(graph_limitObs2,"l")
#mg.Add(graph_limit,"pl")

mg.Draw("APCE5")
mg.GetXaxis().SetRangeUser(1.2,9.)
#mg.GetYaxis().SetRangeUser(0.01,10)  #dark photon
mg.GetYaxis().SetRangeUser(0.001,4)  #2hdm
#mg.GetYaxis().SetTitle("xSec*BR [pb]")
#mg.GetXaxis().SetTitle("Dark Photon Mass [GeV]")
#mg.GetYaxis().SetTitle("#sigma(pp#rightarrow X)#times BR(X#rightarrow #mu#mu) #times Acc.[pb]")
#briefing
mg.GetYaxis().SetTitle("#sigma(pp#rightarrow X)#times B(X#rightarrow #mu#mu) #times Acc.[pb]")

mg.GetYaxis().SetTitleOffset(1.0)
mg.GetYaxis().SetTitleSize(0.05)
mg.GetXaxis().SetTitle("m_{#mu#mu} [GeV]")
mg.GetXaxis().SetMoreLogLabels()
mg.GetXaxis().SetTitleSize(0.05)
mg.GetXaxis().SetTitleOffset(1.2)

mg.GetYaxis().SetLabelSize(0.05);
#mg.GetYaxis().SetLabelFont(43);
#mg.GetYaxis().SetLabelOffset(0.0);
mg.GetXaxis().SetLabelSize(0.05);
#mg.GetXaxis().SetLabelFont(43);
mg.GetXaxis().SetLabelOffset(0.0);

c1.Update()
c1.SetBottomMargin(0.13);
legend=ROOT.TLegend(0.5,0.6,0.8,0.9)
cmsTag=ROOT.TLatex(0.13,0.83,"#scale[1.1]{CMS}")
cmsTag.SetNDC()
cmsTag.SetTextAlign(11)
cmsTag.Draw()
cmsTag2=ROOT.TLatex(0.215,0.917,"#scale[1.0]{#bf{#it{Preliminary}}}")
cmsTag2.SetNDC()
cmsTag2.SetTextAlign(11)
#cmsTag.SetTextFont(61)
#cmsTag2.Draw()
cmsTag3=ROOT.TLatex(0.90,0.917,"#scale[0.9]{#bf{"+str(lumi)+" fb^{-1} (13 TeV)}}")
cmsTag3.SetNDC()
cmsTag3.SetTextAlign(31)
#cmsTag.SetTextFont(61)
cmsTag3.Draw()
#leg=ROOT.TLegend(0.6, 0.65,0.8, 0.85)  
leg=ROOT.TLegend(0.55, 0.60,0.8, 0.85)
leg.SetBorderSize( 0 )
leg.SetFillStyle( 1001 )
leg.SetFillColor(kWhite) 
leg.SetTextSize(0.05)
leg.AddEntry( graph_limitObs1, "Observed",  "L" )
leg.AddEntry( graph_limit1 , "Median expected",  "L" )
leg.AddEntry( graph_limit168up, "68% expected",  "F" ) 
leg.AddEntry( graph_limit195up, "95% expected",  "F" )

leg.Draw("same")
c1.SaveAs("limit"+year+"DarkPhoton.root")
c1.SaveAs("limit"+year+"DarkPhoton.pdf")
c1.SaveAs("limit"+year+"DarkPhoton.png")
c2=ROOT.TCanvas("c2","c2",700,500)
c2.SetLogy()




cmsTag.Draw()
cmsTag2.Draw()
cmsTag3.Draw()
mgeps.Draw("APC")
leg2=ROOT.TLegend(0.65, 0.65,0.87, 0.85)  
leg2.SetBorderSize( 0 )
leg2.SetFillStyle( 1001 )
leg2.SetFillColor(kWhite) 
leg2.Draw("same")
c2.SaveAs("thep.root")


graph_limit1.SaveAs("modelIndepedantLimit"+year+"Lower.root")
graph_limit2.SaveAs("modelIndepedantLimit"+year+"Upper.root")
