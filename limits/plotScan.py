from ROOT import *

fEnv = TFile.Open("higgsCombineEnvelope.MultiDimFit.mH1.803.root","READ")
tEnv = fEnv.Get("limit")
tEnv.GetEntry(0)
offset=2*(tEnv.deltaNLL+tEnv.nll+tEnv.nll0)
print offset
n=tEnv.Draw("r:2*(deltaNLL+nll+nll0)-"+str(offset),"quantileExpected>0 && deltaNLL<100","goff");
#n=tEnv.Draw("r:(2*(deltaNLL+nll0)-"+str(offset)+"+5485286.)","quantileExpected>0 && deltaNLL<10","goff");
gEnv = TGraph(n,tEnv.GetV1(),tEnv.GetV2())

f = {}
t = {}
g = {}

for i in range(4):
 f[i] = TFile.Open("higgsCombinefixed_pdf_"+str(i)+".MultiDimFit.mH1.803.root","READ")
 t[i] = f[i].Get("limit")

 n = t[i].Draw("r:2*(deltaNLL+nll+nll0)-"+str(offset),"quantileExpected>0 && deltaNLL<100","goff"); 
 #n = t[i].Draw("r:(2*(deltaNLL+nll0)-"+str(offset)+"+5485286.)","quantileExpected>0 && deltaNLL<10","goff"); 
 g[i] =TGraph(n,t[i].GetV1(),t[i].GetV2())
 g[i].SetLineColor(i+2)
 g[i].SetMarkerColor(i+2)
 g[i].SetMarkerStyle(20)
 g[i].SetMarkerSize(1.2)
 g[i].SetLineWidth(2)
 
c = TCanvas("c","c",800,800)
c.cd()
gEnv.SetLineColor(1)
gEnv.SetMarkerColor(1)
gEnv.SetMarkerStyle(20)
gEnv.SetMarkerSize(1.2)
gEnv.SetLineWidth(2)
gEnv.GetYaxis().SetTitle("Inverse Log-Likelihood")
#gEnv.GetYaxis().SetRangeUser(-10,50)


gEnv.Draw("ALP")

leg = TLegend(0.7,0.7,0.9,0.9)
leg.AddEntry(gEnv,"pdf Envelope","lp")
for i in range(4):
  g[i].Draw("LPSame")
  leg.AddEntry(g[i],"pdf_index="+str(i),"lp")
gEnv.Draw("LPSame")
leg.Draw("same")

c.SaveAs("test.png")
