##############################################################################
################################# makefile ###################################
##############################################################################
#									     #
#   makefile of QPPnltMP						     #
#									     #
#   Input:  $(CC)     = compiler command				     #
#           $(SW)     = compiler options                                     #
#           $(OPTH)   = the OPT****.h files to include			     #
#           $(OPTINC) = the -I$(OPT****.h directory)			     #
#   	    $(NDOH)   = the .h files to include for NDOSolver/FiOracle	     #
#           $(NDOINC) = the -I$(NDOSolver/FiOracle directories)		     #
#           $(MPSH)   = the .h files to include for MPSolver (only)	     #
#           $(BDLINC) = the -I$(Bundle/MPSolver directory)		     #
#           $(QMPDIR) = the directory where the QPPnltMP source is           #
#									     #
#   Output: $(QMPOBJ) = the final object(s) / library			     #
#           $(QMPLIB) = external libreries + -L<libdirs>		     #
#           $(QMPH)   = the .h files to include				     #
#           $(QMPINC) = the -I$(QPPnltMP directory)			     #
#									     #
#                                VERSION 2.00				     #
#                	        09 - 05 - 2012				     #
#									     #
#			      Antonio Frangioni				     #
#   			   Operations Research Group			     #
#			  Dipartimento di Informatica			     #
#   			     Universita' di Pisa			     #
#									     #
##############################################################################

# clean - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

clean::
	rm -f $(QMPDIR)*.o $(QMPDIR)*~

# macroes to be exported- - - - - - - - - - - - - - - - - - - - - - - - - - -

# comment out BMinQuad.o to build the unconstrained version
QMPOBJ = $(QMPDIR)QPPnltMP.o $(QMPDIR)MinQuad.o\
#	 $(QMPDIR)BMinQuad.o
QMPLIB =
# comment out BMinQuad.h to build the unconstrained version
QMPH   = $(QMPDIR)QPPnltMP.h $(QMPDIR)MinQuad.h\
#	 $(QMPDIR)BMinQuad.h
QMPINC = -I$(QMPDIR)

# dependencies: every .o from its .C + every recursively included .h- - - - -

$(QMPDIR)MinQuad.o: $(QMPDIR)MinQuad.C $(QMPDIR)MinQuad.h $(OPTH)
	$(CC) -c $*.C -o $@ $(OPTINC) $(SW)

# comment the following three lines to build the unconstrained version
#$(QMPDIR)BMinQuad.o: $(QMPDIR)BMinQuad.C $(QMPDIR)BMinQuad.h\
#	$(QMPDIR)MinQuad.h $(OPTH)
#	$(CC) -c $*.C -o $@ $(OPTINC) $(SW)

$(QMPDIR)QPPnltMP.o: $(QMPDIR)QPPnltMP.C $(QMPH) $(MPSH) $(NDOH) $(OPTH)
	$(CC) -c $*.C -o $@ $(NDOINC) $(BDLINC) $(OPTINC) $(SW)

########################### End of makefile ##################################
