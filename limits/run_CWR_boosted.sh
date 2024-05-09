root -l -b -q MultiMakeCardsAndWS_ptcut_CWR.cpp
root -l -b -q makeDual_CWR_ptcutkk.cpp
root -l -b -q makeDual_CWR_ptcutkk_offpeak.cpp
root -l -b -q makeDual_CWR_ptcutkpi.cpp
root -l -b -q makeDual_CWR_ptcutkpi_offpeak.cpp
root -l -b -q makeDual2_CWR_ptcut.cpp

root -l -b -q MultiMakeCardsAndWS_ptcut_edgeCWR.cpp

mv output_edge/* output/

mkdir -p output_dual/dual2card
mv output_dual/*dual2.txt output_dual/dual2card

mkdir -p combine_output/dual
python bothYearsLimitProcessing_dual_other.py

python bothYearsLimitProcessing.py 

mv combine_output/dual/* combine_output/bothYears


#pvalue
python check_pvalue.py
python check_pvalue_dual.py
cp combine_output/test/* combine_output/pval/

#CL90
python bothYearsLimitProcessing_CL90.py
python bothYearsLimitProcessing_CL90_dual.py
mv combine_output/dual_CL90/* combine_output/CL90


# separate year
python pyDPLimitsProcessing.py 2017
python pyDPLimitsProcessing.py 2018
python pyDPLimitsProcessing_dual_other.py 2017
cp combine_output/dual_2017/* combine_output/2017
python pyDPLimitsProcessing_dual_other.py 2018
cp combine_output/dual_2018/* combine_output/2018


