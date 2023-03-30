root -l -b -q MultiMakeCardsAndWS.cpp

python pyDPLimitsProcessing.py 2018

python makeDPLimit.py 2018

python bothYearsLimitProcessing.py

root -l makeDual_1.cpp   //make 207 SR

root -l makeDual_2.cpp   //make 207 high ip region

root -l makeDual_3.cpp   //make 215,216 SR

root -l makeDual_kk.cpp  // make 202-206 SR

root -l makeDual_bothD.cpp //make 208-214 SR

root -l makeDual_kpi.cpp   //make 217-222 SR

root -l makeDual2_bothD.cpp  //make 202-222 CR

python pyDPLimitsProcessing_dual_other.py year

mv higgs* file to combine_output

python bothYearsLimitProcessing_dual_other.py

mv combine file to combine_output


