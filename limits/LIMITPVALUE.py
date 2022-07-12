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

#files = glob("combine_output/"+year+"/higgsCombineasympMassIndex_*.AsymptoticLimits.mH*.root")
files = glob("combine_output/"+year+"/higgsCombineasympMassIndex_*.Significance.mH*.root")

d_m = {}
for fname in files:
        m = float(re.search(r"mH(\d+\.?\d+).root", fname).group(1))
        d = int(re.search(r"Index_(\d+).Significance", fname).group(1))
	d_m[d] = [m, fname]


d_m = sorted(d_m.items())
print d_m

for d,m_fname in d_m:
	m, fname = m_fname
        ##if (m == 1.716 or m ==7.632):
        ##        continue
	#file90=glob.glob("higgsCombineIterV9_CL90_ForPress_2018_"+str(d)+".AsymptoticLimits.mH*.root")
        print "mass=",m
	acc = accgraph.Eval(m)

        if (m < 3):
                f=ROOT.TFile.Open(fname)
                tree=f.Get("limit")
                tree.GetEntry(0)
                limit1.append(tree.limit*a)


                print("Mass: " + str(m) + "         pVal: "+ str(tree.limit*a))
                mass1.append(m)
                masserr1.append(0.)
                fff.write("{0} {1} {2}\n".format(m, tree.limit*a, math.sqrt(tree.limit*a)))
        if (m > 3):
                f=ROOT.TFile.Open(fname)
                tree=f.Get("limit")
                tree.GetEntry(0)
                limit2.append(tree.limit)

                print("Mass: " + str(m) + "         pVal: "+ str(tree.limit*a))


                mass2.append(m)
                masserr2.append(0.)
                fff.write("{0} {1} {2}\n".format(m, tree.limit*a, math.sqrt(tree.limit*a)))


print limit1

c1=ROOT.TCanvas("c1","c1",700,500)
#c1.SetGrid()
c1.SetLogy()
c1.SetLogx()


mg=ROOT.TMultiGraph()
mgeps=ROOT.TMultiGraph()
graph_limit1=ROOT.TGraph(len(mass1),mass1,limit1)
graph_limit1.SetTitle("graph_limit1")
graph_limit1.SetMarkerSize(0)
graph_limit1.SetMarkerStyle(20)
graph_limit1.SetMarkerColor(kBlack)
graph_limit1.SetLineWidth(2)
graph_limit1.SetLineStyle(1)
graph_limit1.GetYaxis().SetRangeUser(0,100)
graph_limit1.GetXaxis().SetRangeUser(10,70)
graph_limit1.GetXaxis().SetMoreLogLabels()
graph_limit1.GetYaxis().SetTitle("P-Value")
graph_limit1.GetYaxis().SetTitleSize(2)
graph_limit1.GetXaxis().SetTitle("Dark Photon Mass [GeV]")

graph_limit2=ROOT.TGraph(len(mass2),mass2,limit2)
graph_limit2.SetTitle("graph_limit2")
graph_limit2.SetMarkerSize(0)
graph_limit2.SetMarkerStyle(20)
graph_limit2.SetMarkerColor(kBlack)
graph_limit2.SetLineWidth(2)
graph_limit2.SetLineStyle(1)


mg.Add(graph_limit1,"pc")
mg.Add(graph_limit2,"pc")
#mg.Add(graph_limit,"pl")

mg.Draw("APCE5")
mg.GetXaxis().SetRangeUser(1.2,9.)
#mg.GetYaxis().SetRangeUser(0.01,5)
#mg.GetYaxis().SetTitle("xSec*BR [pb]")
#mg.GetXaxis().SetTitle("Dark Photon Mass [GeV]")
mg.GetYaxis().SetTitle("P-Value")
mg.GetYaxis().SetTitleOffset(0.9)
mg.GetYaxis().SetTitleSize(0.05)
mg.GetXaxis().SetTitle("Dimuon Mass [GeV]")
mg.GetXaxis().SetMoreLogLabels()
mg.GetXaxis().SetTitleSize(0.05)
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
leg=ROOT.TLegend(0.6, 0.35,0.8, 0.55)  
leg.SetBorderSize( 0 )
leg.SetFillStyle( 1001 )
leg.SetFillColor(kWhite) 
#leg.AddEntry( obse , "Observed",  "LP" )
leg.AddEntry( graph_limit1 , "Local P-Value",  "PL" )
leg.Draw("same")
c1.SaveAs("pval"+year+"DarkPhoton.root")
c1.SaveAs("pval"+year+"DarkPhoton.pdf")
c1.SaveAs("pval"+year+"DarkPhoton.png")
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
