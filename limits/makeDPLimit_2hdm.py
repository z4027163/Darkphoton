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

fff_obs = open("limit_2hdm_obs.txt", "w")
fff_exp = open("limit_2hdm_exp.txt", "w")


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

#for theoretical
limit1Observed_up=array('d')
limit1Observed_dow=array('d')
limit2Observed_up=array('d')
limit2Observed_dow=array('d')

limit1Observed_up_rel=array('d')
limit1Observed_dow_rel=array('d')
limit2Observed_up_rel=array('d')
limit2Observed_dow_rel=array('d')

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


m= array('d',[1.176, 1.187, 1.199, 1.211, 1.223, 1.235, 1.248, 1.26, 1.273, 1.286, 1.299, 1.312, 1.325, 1.338, 1.351, 1.365, 1.378, 1.392, 1.406, 1.42, 1.434, 1.449, 1.463,
1.478, 1.493, 1.508, 1.523, 1.538, 1.553, 1.569, 1.584, 1.6, 1.616, 1.632, 1.649, 1.665, 1.682, 1.699, 1.716, 1.733, 1.75, 1.768, 1.785, 1.803, 1.821, 1.839,
1.858, 1.876, 1.895, 1.914, 1.933, 1.953, 1.972, 1.992, 2.012, 2.032, 2.052, 2.073, 2.094, 2.114, 2.136, 2.157, 2.179, 2.2, 2.222, 2.245, 2.267, 2.29, 2.313,
2.336, 2.359, 2.383, 2.406, 2.43, 2.455, 2.479, 2.504, 2.529, 2.554, 2.58, 2.5801, 4.2009, 4.201, 4.243, 4.286, 4.328, 4.372, 4.415, 4.46, 4.504, 4.549, 4.595, 4.641, 4.687,
4.734, 4.781, 4.829, 4.877, 4.926, 4.975, 5.025, 5.075, 5.126, 5.177, 5.229, 5.282, 5.334, 5.388, 5.442, 5.496, 5.551, 5.606, 5.663, 5.719, 5.776, 5.834, 5.892,
5.951, 6.011, 6.071, 6.132, 6.193, 6.255, 6.318, 6.381, 6.445, 6.509, 6.574, 6.64, 6.706, 6.773, 6.841, 6.909, 6.978, 7.048, 7.119, 7.19, 7.262, 7.334, 7.408,
7.482, 7.557, 7.632, 7.709, 7.786, 7.864])
xSec= array('d',[3383220.0, 3006680.0, 2658290.0, 2368340.0, 2129990.0, 1936410.0, 1769280.0, 1647010.0, 1541050.0, 1454200.0, 1377770.0, 1304270.0, 1232400.0, 1162790.0, 1096050.0, 1028110.0, 969336.0, 911430.0, 859849.0, 814678.0, 775188.0, 738297.0, 708101.0, 679381.0, 653485.0, 629478.0, 606876.0, 585602.0, 565597.0, 545588.0, 528013.0, 510459.0, 494065.0, 478757.0, 463602.0, 450305.0, 437120.0, 424821.0, 413321.0, 402532.0, 392366.0, 382183.0, 373024.0, 363710.0, 354718.0, 346042.0, 337231.0, 329215.0, 321108.0, 313370.0, 306006.0, 298662.0, 292077.0, 285564.0, 279483.0, 273816.0, 268528.0, 263343.0, 258495.0, 254153.0, 249642.0, 245554.0, 241454.0, 237671.0, 233809.0, 229887.0, 226262.0, 222618.0, 219140.0, 215845.0, 212747.0, 209744.0, 207102.0, 204610.0, 202320.0, 200434.0, 198817.0, 197556.0, 196616.0, 195931.0, 195928.0, 88277.4, 88272.1, 86159.4, 84170.2, 82386.3, 80674.9, 79146.7, 77688.1, 76389.0, 75177.9, 74049.7, 73018.9, 72071.9, 71176.9, 70340.9, 69532.9, 68755.8, 67978.5, 67201.0, 66390.9, 65556.7, 64683.5, 63790.8, 62863.8, 61904.7, 60953.1, 59957.2, 58957.1, 57956.4, 56940.4, 55931.3, 54896.9, 53896.2, 52897.5, 51906.3, 50944.8, 50001.8, 49084.1, 48210.9, 47367.5, 46567.0, 45795.1, 45051.9, 44347.9, 43670.6, 43029.3, 42412.4, 41819.1, 41256.6, 40714.5, 40191.5, 39693.2, 39209.8, 38739.6, 38280.2, 37835.5, 37396.3, 36965.7, 36528.7, 36093.7, 35651.5, 35204.1, 34735.6, 34253.8, 33748.2])
Acc= array('d',[0.0031132290489434052, 0.0031950108551093218, 0.0032842080310042064, 0.0033733850261104172, 0.0034625418404279525, 0.003551678473956814, 0.0036482203868203056, 0.003737314977039429, 0.0038338113429840365, 0.00393028402453083, 0.0040267330216798, 0.004123158334430952, 0.004219559962784284, 0.004315937906739798, 0.004412292166297493, 0.0045160318046898774, 0.004612336873575177, 0.004716023537397291, 0.004819682732923708, 0.004923314460154432, 0.005026918719089458, 0.00513789280083457, 0.005241440161157086, 0.005352353280103077, 0.005463234866566765, 0.005574084920548148, 0.0056849034420472285, 0.005795690431064008, 0.005906445887598482, 0.006024550285432539, 0.0061352405748369195, 0.006253275461065544, 0.006371274470336523, 0.006489237602649857, 0.006614534120238175, 0.006732423256326371, 0.006857641152925477, 0.006982818547802869, 0.007107955440958545, 0.00723305183239251, 0.007358107722104755, 0.007490475695148383, 0.007615448198962985, 0.0077477278810561614, 0.007879962156374824, 0.008012151024918965, 0.008151634458748763, 0.008283729991145287, 0.008423114903485933, 0.008562449223710526, 0.008701732951819068, 0.00884829274615698, 0.008987472627290468, 0.009133923109023062, 0.009280317533009338, 0.009426655899249296, 0.009572938207742937, 0.009726474299511786, 0.009879948587615319, 0.010026057117095375, 0.010186711752826446, 0.010340000629934036, 0.010500522689095059, 0.01065368501584034, 0.010814074498431326, 0.01098168188733988, 0.011141932627008732, 0.011309394966498689, 0.01147678316961914, 0.011644097236370084, 0.011811337166751529, 0.011985769357031597, 0.012152857791353616, 0.01232713189878907, 0.012508581588943619, 0.012682690886604899, 0.012863968899911352, 0.01304515932298918, 0.013226262155838392, 0.013414516186287051, 0.013415240057361914, 0.024963651972764275, 0.024964353127838894, 0.02525871435754353, 0.025559828074126927, 0.025853688988445686, 0.026161278126496677, 0.026461614462283033, 0.026775642444122144, 0.027082417623696613, 0.027395884327350732, 0.02771602363559511, 0.028035866398361475, 0.02835541261564982, 0.028681599204324296, 0.029007476214094732, 0.02933996724811054, 0.029672135389507553, 0.030010890647431704, 0.030349309418733574, 0.030694287838266894, 0.031038915896885726, 0.031390075574862834, 0.0317408707373445, 0.03209816892973384, 0.03246194828950977, 0.03281848129366187, 0.033188326012491075, 0.033557762070349634, 0.03392678946723754, 0.034302230622072084, 0.034677247840200105, 0.03506545470726263, 0.03544640751076122, 0.035833711707963664, 0.036227343334183426, 0.03662050351475667, 0.03701995861085435, 0.03742568381692375, 0.03783090450327629, 0.038242361664952755, 0.038653297349444095, 0.039070435314034105, 0.03949374949187073, 0.03991650743671953, 0.040345406278434776, 0.04077373108882773, 0.04120816091912949, 0.04164866844118872, 0.04208856549439056, 0.042534503241237295, 0.042986453512711416, 0.04343775575663811, 0.04389503240592468, 0.044358254450687425, 0.04482739246060951, 0.04529582400278364, 0.04577013170911722, 0.046243712907058546, 0.046729692078412, 0.04721490381921836, 0.04770588941879717, 0.04819608670631838, 0.04869853579667821, 0.04920015397109325, 0.04970743949456272])

BR = array('d',[0.02418, 0.02418, 0.02418, 0.02545, 0.02651, 0.02651, 0.02651, 0.02651, 0.02651, 0.02651, 0.02651, 0.02738, 0.02738, 0.02738, 0.02738, 0.02738, 0.02738, 0.02738, 0.02804, 0.02804, 0.02804, 0.02804, 0.02851, 0.02851, 0.02851, 0.02878, 0.02878, 0.02878, 0.02878, 0.02878, 0.02878, 0.02887, 0.02887, 0.02887, 0.02879, 0.02879, 0.02879, 0.02879, 0.02854, 0.02854, 0.02854, 0.02854, 0.02854, 0.02814, 0.02814, 0.02814, 0.02758, 0.02758, 0.02758, 0.02688, 0.02688, 0.02688, 0.02688, 0.02688, 0.02605, 0.02605, 0.02605, 0.02508, 0.02508, 0.02399, 0.02399, 0.02399, 0.02399, 0.02276, 0.02141, 0.02141, 0.02141, 0.02141, 0.01991, 0.01991, 0.01933, 0.01873, 0.01811, 0.01746, 0.0168, 0.01611, 0.01539, 0.01465, 0.01389, 0.01309, 0.01309, 0.00395, 0.00395, 0.00402, 0.00401, 0.00397, 0.00381, 0.00374, 0.00361, 0.00347, 0.0034, 0.00328, 0.0031, 0.00304, 0.00289, 0.00284, 0.0028, 0.00265, 0.00252, 0.00252, 0.00242, 0.00242, 0.00232, 0.00232, 0.00225, 0.00225, 0.00218, 0.00218, 0.00213, 0.00213, 0.00208, 0.00204, 0.00204, 0.002, 0.002, 0.00197, 0.00197, 0.00195, 0.00191, 0.00191, 0.00189, 0.00189, 0.00187, 0.00186, 0.00186, 0.00185, 0.00184, 0.00184, 0.00183, 0.00182, 0.00182, 0.00182, 0.00181, 0.00181, 0.0018, 0.0018, 0.0018, 0.00179, 0.00179, 0.00179, 0.00179, 0.00179, 0.00178, 0.00178, 0.00178, 0.00178])
#xSec_dow= array('d',[14.1,17.52,22.11,26.53,32.19,45.46,47.92,41,35.83,32.54,32.5,35.79,32.88,32.54,5.125])
#xSec_up=array('d',[225.9,194.7,177.3,155.4,140.9,105.4,70.21,51.88,42.92,37.83,36,50.96,63.33,42.5,5.017])

fidxsec = array('d')
for i in range(len(m)):
    fidxsec.append(xSec[i]*Acc[i]*BR[i])


m_unc = array('d',[1.05,1.176,  1.3,  1.4,  1.5,  1.6,  1.7,  1.8,  1.9  ,2.0,  2.5, 2.6 ,4.0,  4.2, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0])
xsec_unc = array('d',[1.02,0.89,0.76, 0.67, 0.58, 0.52, 0.46, 0.42, 0.38, 0.32, 0.17, 0.14, 0.10,0.10,0.10,0.10,0.10,0.10,0.10,0.10])

xsecgraph = TGraph(len(m),m,fidxsec);

uncgraph = TGraph(len(m_unc),m_unc,xsec_unc)

eps2scale = 1.
base_eps = 0.02 #epsilon for which the cross sections are computed

a = 1 # for lumi projection (6.6->100)

files = glob("FINAL/boosted/combine_output/"+year+"/higgsCombineasympMassIndex_*.AsymptoticLimits.mH*.root")
#files = glob("combine_output/"+year+"/higgsCombineasympMassIndex_*.AsymptoticLimits.mH*.root")

d_m = {}
for fname in files:
        m = float(re.search(r"mH(.*?)\.root", fname).group(1))
        d = int(re.search(r"Index_(\d+).Asymptotic", fname).group(1))
        d_m[d] = [m, fname]


d_m = sorted(d_m.items())
print d_m

for d,m_fname in d_m:
    m, fname = m_fname
    #if (m == 1.716 or m ==7.632):
    #        continue
    #file90=glob.glob("higgsCombineIterV9_CL90_ForPress_2018_"+str(d)+".AsymptoticLimits.mH*.root")
    xsec = xsecgraph.Eval(m,0,"S")

    unc_up=uncgraph.Eval(m,0,"S")
    unc_dow=uncgraph.Eval(m,0,"S")
   
    if (m>3):
        unc_up=0.10
        unc_dow=0.10
    #add 0.3 to acceptance uncertainty
    unc_up = sqrt(unc_up*unc_up+0.3*0.3)
    unc_dow = sqrt(unc_dow*unc_dow+0.3*0.3)


    if (m < 3 and m > 1.17):
            f=ROOT.TFile.Open(fname)
            tree=f.Get("limit")
            tree.GetEntry(2)
            limit1.append(sqrt(tree.limit*a/(xsec)))
            a_exp=sqrt(tree.limit*a/(xsec))
             
            tree.GetEntry(0)
            limit195up.append(abs(sqrt(tree.limit*a/(xsec))-limit1[-1]))
            tree=f.Get("limit")
            tree.GetEntry(4)
            limit195down.append(abs(sqrt(tree.limit*a/(xsec))-limit1[-1]))
            
            
            tree.GetEntry(1)
            limit168up.append(abs(sqrt(tree.limit*a/(xsec))-limit1[-1]))
            tree=f.Get("limit")
            tree.GetEntry(3)
            limit168down.append(abs(sqrt(tree.limit*a/(xsec))-limit1[-1]))

            tree.GetEntry(5)
            limit1Observed.append(sqrt(abs(tree.limit*a/(xsec))))
            a_obs=sqrt(tree.limit*a/(xsec))            
           
            
            limit1Observed_up_rel.append(sqrt(abs(tree.limit*a/(xsec)))*(sqrt(1+unc_up)-1))
            limit1Observed_dow_rel.append(sqrt(abs(tree.limit*a/(xsec)))*(1-sqrt(1-unc_dow)))

            tree=f.Get("limit")
            mass1.append(m)
            masserr1.append(0.)
            #fff.write("{0} {1}\n".format(m, math.sqrt(tree.limit*a)))
            fff_exp.write("{0} {1}\n".format(m, a_exp))
            fff_obs.write("{0} {1}\n".format(m, a_obs))
            
    if (m > 3):
            f=ROOT.TFile.Open(fname)
            tree=f.Get("limit")
            tree.GetEntry(2)
            limit2.append(sqrt(tree.limit*a/(xsec)))
            a_exp=sqrt(tree.limit*a/(xsec))

            tree.GetEntry(0)
            limit295up.append(abs(sqrt(tree.limit*a/(xsec))-limit2[-1]))
            tree=f.Get("limit")
            tree.GetEntry(4)
            limit295down.append(abs(sqrt(tree.limit*a/(xsec))-limit2[-1]))
            
            
            tree.GetEntry(1)
            limit268up.append(abs(sqrt(tree.limit*a/(xsec))-limit2[-1]))
            tree=f.Get("limit")
            tree.GetEntry(3)
            limit268down.append(abs(sqrt(tree.limit*a/(xsec))-limit2[-1]))
            
            tree.GetEntry(5)
            limit2Observed.append(abs(sqrt(tree.limit*a/(xsec))))
            a_obs=sqrt(tree.limit*a/(xsec))           


            limit2Observed_up_rel.append(sqrt(abs(tree.limit*a/(xsec)))*(sqrt(1+unc_up)-1))
            limit2Observed_dow_rel.append(sqrt(abs(tree.limit*a/(xsec)))*(1-sqrt(1-unc_dow)))

 
            tree=f.Get("limit")
            mass2.append(m)
            masserr2.append(0.)
            fff_exp.write("{0} {1}\n".format(m, a_exp))
            fff_obs.write("{0} {1}\n".format(m, a_obs))


c1=ROOT.TCanvas("c1","c1",700,500)
#c1.SetGrid()
c1.SetLogy()
#qier
c1.SetLogx()

mg=ROOT.TMultiGraph()
mgeps=ROOT.TMultiGraph()
graph_limit1=ROOT.TGraph(len(mass1),mass1,limit1)
graph_limit1.SetTitle("graph_limit1")
graph_limit1.SetMarkerSize(0)
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
graph_limit1.GetXaxis().SetTitle("m(a) [GeV]")

graph_limit2=ROOT.TGraph(len(mass2),mass2,limit2)
graph_limit2.SetTitle("graph_limit2")
graph_limit2.SetMarkerSize(0)
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
graph_limit195up.SetLineColor(ROOT.TColor.GetColor(252,241,15))

graph_limit295up=ROOT.TGraphAsymmErrors(len(mass2),mass2,limit2,masserr2,masserr2,limit295up,limit295down)
graph_limit295up.SetTitle("graph_limit295up")
graph_limit295up.SetFillColor(ROOT.TColor.GetColor(252,241,15))


graph_limit168up=ROOT.TGraphAsymmErrors(len(mass1),mass1,limit1,masserr1,masserr1,limit168up,limit168down)
graph_limit168up.SetTitle("graph_limit168up")
graph_limit168up.SetMarkerColor(kGreen)
graph_limit168up.SetFillColor(kGreen);
graph_limit168up.SetLineColor(kGreen);

graph_limit1_unc=ROOT.TGraphAsymmErrors(len(mass1),mass1,limit1Observed,masserr1,masserr1,limit1Observed_up_rel,limit1Observed_dow_rel)
graph_limit1_unc.SetTitle("graph_limit1_unc")
graph_limit1_unc.SetFillColorAlpha(15, 0.7);
#graph_limit1_unc.SetLineWidth(-802);
graph_limit1_unc.SetLineColorAlpha(15,0.7)

graph_limit2_unc=ROOT.TGraphAsymmErrors(len(mass2),mass2,limit2Observed,masserr2,masserr2,limit2Observed_up_rel,limit2Observed_dow_rel)
graph_limit2_unc.SetTitle("graph_limit2_unc")
#graph_limit2_unc.SetFillColor(ROOT.TColor.GetColor(15))
graph_limit2_unc.SetFillColorAlpha(15, 0.7);
#graph_limit1_unc.SetLineColor(15)


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

#theory up
col=16
tem_t=1;

graph_limit268up=ROOT.TGraphAsymmErrors(len(mass2),mass2,limit2,masserr2,masserr2,limit268up,limit268down)

graph_limit268up.SetTitle("graph_limit268up")
graph_limit268up.SetFillColor(kGreen);
graph_limit268up.SetMarkerColor(kGreen)

mg.Add(graph_limit195up,"3")
mg.Add(graph_limit168up,"3")
mg.Add(graph_limit1_unc,"3")
mg.Add(graph_limit1,"l")
mg.Add(graph_limitObs1,"l")
mg.Add(graph_limit295up,"3")
mg.Add(graph_limit268up,"3")
mg.Add(graph_limit2_unc,"3")
mg.Add(graph_limit2,"l")
mg.Add(graph_limitObs2,"l")


mg.Draw("APC")
mg.GetXaxis().SetRangeUser(1.16,9.)
mg.GetYaxis().SetRangeUser(0.002,10.99)
#mg.GetYaxis().SetTitle("xSec*BR [pb]")
#mg.GetXaxis().SetTitle("Dark Photon Mass [GeV]")
mg.GetYaxis().SetTitle("sin(#it{#theta}_{H})")
mg.GetYaxis().SetTitleOffset(1.0)
mg.GetYaxis().SetTitleSize(0.05)
mg.GetXaxis().SetTitle("#it{m_{a}} [GeV]")
mg.GetXaxis().SetTitleSize(0.05)
mg.GetXaxis().SetTitleOffset(1.2)
mg.GetXaxis().SetMoreLogLabels()

mg.GetYaxis().SetLabelSize(0.05);
mg.GetXaxis().SetLabelSize(0.05);
mg.GetXaxis().SetLabelOffset(0.0);

c1.SetBottomMargin(0.13);
c1.Update()
legend=ROOT.TLegend(0.5,0.6,0.8,0.9)
cmsTag=ROOT.TLatex(0.13,0.83,"#scale[1.2]{CMS}")
cmsTag.SetNDC()
cmsTag.SetTextAlign(11)
cmsTag.Draw()
#cmsTag2=ROOT.TLatex(0.215,0.917,"#scale[1.0]{#bf{#it{Preliminary}}}")
cmsTag2=ROOT.TLatex(0.23,0.83,"#scale[1.0]{#bf{#it{Supplementary}}}")
cmsTag2.SetNDC()
cmsTag2.SetTextAlign(11)
#cmsTag.SetTextFont(61)
cmsTag2.Draw()
cmsTag3=ROOT.TLatex(0.90,0.917,"#scale[0.9]{#bf{"+str(lumi)+" fb^{-1} (13 TeV)}}")
cmsTag3.SetNDC()
cmsTag3.SetTextAlign(31)
#cmsTag.SetTextFont(61)
cmsTag3.Draw()
leg=ROOT.TLegend(0.55, 0.54,0.80, 0.85)  
leg.SetBorderSize( 0 )
leg.SetFillStyle( 1001 )
leg.SetFillColor(kWhite) 
leg.SetTextSize(0.05)
leg.AddEntry( graph_limitObs1 , "Observed",  "LP" )
leg.AddEntry( graph_limit1 , "Median expected",  "L" )
leg.AddEntry( graph_limit168up, "68% expected",  "F" )
leg.AddEntry( graph_limit195up, "95% expected",  "F" )
leg.AddEntry( graph_limit1_unc, "theoretical #pm 1#sigma",  "F" )


#qier
tail=""
leg.Draw("same")
c1.SetBottomMargin(0.13);
c1.SaveAs("limit"+year+"2HDM.root")
c1.SaveAs("limit"+year+"2HDM"+tail+".pdf")
c1.SaveAs("limit"+year+"2HDM"+tail+".png")
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

