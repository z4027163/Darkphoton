root -l -b -q MultiMakeCardsAndWS_CWR.cpp
root -l -b -q makeDual_CWR_1.cpp
root -l -b -q makeDual_CWR_3.cpp
root -l -b -q makeDual_CWR_kk.cpp
root -l -b -q makeDual_CWR_bothD.cpp
root -l -b -q makeDual_CWR_kpi.cpp
root -l -b -q makeDual2_CWR_bothD.cpp

mkdir -p output_dual/dual2card
mv output_dual/*dual2.txt output_dual/dual2card

#bothyear limit
mkdir -p combine_output/dual
python bothYearsLimitProcessing_dual_other.py
python bothYearsLimitProcessing.py

cp combine_output/dual/* combine_output/bothYears


#CL90
python bothYearsLimitProcessing_CL90.py
python bothYearsLimitProcessing_CL90_dual.py
mv combine_output/dual_CL90/* combine_output/CL90

#pval
python check_pvalue.py
python check_pvalue_dual.py
cp combine_output/test/* combine_output/pval/


#separate year
python pyDPLimitsProcessing.py 2017
python pyDPLimitsProcessing_dual_other.py 2017
mv combine_output/dual_2017/higgsCombineasympMassIndex* combine_output/2017
python pyDPLimitsProcessing.py 2018
python pyDPLimitsProcessing_dual_other.py 2018
mv combine_output/dual_2018/higgsCombineasympMassIndex* combine_output/2018


