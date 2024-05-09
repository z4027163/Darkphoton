root -l -b -q makeDual_1.cpp
root -l -b -q makeDual_2.cpp
root -l -b -q makeDual_3.cpp
root -l -b -q makeDual_kk.cpp
root -l -b -q makeDual_bothD.cpp
root -l -b -q makeDual_kpi.cpp
root -l -b -q makeDual2_bothD.cpp

python pyDPLimitsProcessing_dual_other.py 2017
mv higgsCombineasympMassIndex* combine_output/2017
python pyDPLimitsProcessing_dual_other.py 2018
mv higgsCombineasympMassIndex* combine_output/2018
python bothYearsLimitProcessing_dual_other.py
mv higgsCombineasympMassIndex* combine_output/bothYears
python bothYearsLimitProcessing_CL90_dual.py
mv higgsCombineasympMassIndex* combine_output/CL90
