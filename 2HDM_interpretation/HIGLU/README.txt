# HIGLU and HDECAY
The total cross section for Higgs boson production via the dominant gluon fusion mechanism including NLO [two-loop] QCD corrections can be calculated numerically with the program HIGLU.
The QCD corrections are included for arbitrary Higgs and quark masses and increase the cross section at the LHC by up to a factor of 2.
The source code HIGLU provides the evaluation of the production of the Standard Model [SM] Higgs boson as well as the neutral Higgs bosons of the minimal supersymmetric extension [MSSM].
The program HDECAY determines the decay widths and branching ratios of the Higgs bosons within the SM and the MSSM, including the dominant higher-order corrections. 

1) To setup HIGLU:
   ```
   mkdir HIGLU; cd HIGLU
   wget http://tiger.web.psi.ch/higlu/higlu.tar.gz
   tar -xzvf higlu.tar.gz
   cd higlu
   make
   ```
Note: Download and install LHAPDF-6.2.3 from here:

https://lhapdf.hepforge.org/downloads?
https://lhapdf.hepforge.org/lhapdf5/install

Then edit the makefile and include the correct path.

2) To run HIGLU, settings can be changed in the higlu.in configuration file. Then: 
   ```
   make
   ./run
   ```
   
3) Documentation:

HIGLU
http://tiger.web.psi.ch/higlu/
https://arxiv.org/abs/hep-ph/9610350

HDECAY
https://arxiv.org/pdf/hep-ph/9704448.pdf
https://arxiv.org/pdf/1810.00768.pdf








