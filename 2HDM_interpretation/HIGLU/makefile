# Makro-Definitions:

LIBS =

FFLAGS =

LIBS =

HDECAY = hdecay.o haber.o feynhiggs.o hsqsq.o susylha.o hgaga.o dmb.o elw.o hgg.o

OBJS = higlu.o Cteq61Pdf.o $(HDECAY)

LIBS = /home/elisa/work/DarkPhotonSearch_Scouting/2HDM_Limits/Grids/lib/libLHAPDF.so

#FFLAGS= -fno-emulate-complex -fno-automatic -ffixed-line-length-none -ffast-math -march=pentiumpro -malign-double -Wall -fno-silent

#FC=f77

#FFLAGS= -fno-fixed-form -fno-automatic -ffixed-line-length-none -ffast-math -march=pentiumpro -malign-double -Wall

FC=gfortran

#FFLAGS= -pc 64 -g77libs

#FC=pgf77

# Commands:
.f.o:
	$(FC) -c $(FFLAGS) $*.f

higlu: $(OBJS)
	$(FC) $(FFLAGS) $(OBJS) $(LIBS) -o run

higlu0: $(OBJS)
	$(FC) $(FFLAGS) $(OBJS) $(LIBS) -o run0

clean:
	rm -f $(OBJS)
