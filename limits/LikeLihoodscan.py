from ROOT import *

fEnv2017 = TFile.Open("higgsCombineEnvelope2017.MultiDimFit.mH1.584.root","READ")
fEnv2018 = TFile.Open("higgsCombineEnvelope2018.MultiDimFit.mH1.584.root","READ")
fEnvBoth = TFile.Open("higgsCombineEnvelopeBoth.MultiDimFit.mH1.584.root","READ")

tEnv2017 = fEnv2017.Get("limit")
tEnv2018 = fEnv2018.Get("limit")
tEnvBoth = fEnvBoth.Get("limit")

tEnv2017.GetEntry(0)
tEnv2018.GetEntry(0)
tEnvBoth.GetEntry(0)

offset2017=2*(tEnv2017.deltaNLL+tEnv2017.nll+tEnv2017.nll0)
offset2018=2*(tEnv2018.deltaNLL+tEnv2018.nll+tEnv2018.nll0)
offsetBoth=2*(tEnvBoth.deltaNLL+tEnvBoth.nll+tEnvBoth.nll0)

print offset2017
n2017=tEnv2017.Draw("r:2*(deltaNLL+nll+nll0)-"+str(offset2017),"quantileExpected>0 && deltaNLL<100","goff");
n2018=tEnv2018.Draw("r:2*(deltaNLL+nll+nll0)-"+str(offset2018),"quantileExpected>0 && deltaNLL<100","goff");
nBoth=tEnvBoth.Draw("r:2*(deltaNLL+nll+nll0)-"+str(offsetBoth),"quantileExpected>0 && deltaNLL<100","goff");

#n=tEnv.Draw("r:(2*(deltaNLL+nll0)-"+str(offset)+"+5485286.)","quantileExpected>0 && deltaNLL<10","goff");
gEnv2017 = TGraph(n2017,tEnv2017.GetV1(),tEnv2017.GetV2())
gEnv2018 = TGraph(n2018,tEnv2018.GetV1(),tEnv2018.GetV2())
gEnvBoth = TGraph(nBoth,tEnvBoth.GetV1(),tEnvBoth.GetV2())


c = TCanvas("c","c",800,800)
c.cd()
gEnv2017.SetLineColor(4)
gEnv2017.SetMarkerColor(4)
gEnv2017.SetMarkerStyle(20)
gEnv2017.SetMarkerSize(1.2)
gEnv2017.SetLineWidth(2)
gEnv2017.SetTitle("1.584 GeV Excess Compatibility")
gEnv2017.GetYaxis().SetTitle("Inverse Log-Likelihood")

gEnv2018.SetLineColor(2)
gEnv2018.SetMarkerColor(2)
gEnv2018.SetMarkerStyle(20)
gEnv2018.SetMarkerSize(1.2)
gEnv2018.SetLineWidth(2)

gEnvBoth.SetLineColor(1)
gEnvBoth.SetMarkerColor(1)
gEnvBoth.SetMarkerStyle(20)
gEnvBoth.SetMarkerSize(1.2)
gEnvBoth.SetLineWidth(3)


gEnv2017.Draw("ALP")
gEnv2018.Draw("LP")
gEnvBoth.Draw("LP")

leg = TLegend(0.7,0.7,0.9,0.9)
leg.AddEntry(gEnv2017,"2017","lp")
leg.AddEntry(gEnv2018,"2018","lp")
leg.AddEntry(gEnvBoth,"2017+2018","lp")


leg.Draw("same")


c.SaveAs("test1_584.png")
