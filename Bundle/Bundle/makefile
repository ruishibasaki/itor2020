##############################################################################
################################ makefile ####################################
##############################################################################
#									     #
#   makefile of Bundle							     #
#									     #
#   Input:  $(CC)     = compiler command				     #
#           $(SW)     = compiler options                                     #
#           $(OPTH)   = the OPT****.h files to include			     #
#           $(OPTINC) = the -I$(OPT****.h directory)			     #
#   	    $(NDOH)   = the .h files to include for NDOSolver/FiOracle	     #
#           $(FIOH)   = the .h files to include for FiOracle (only)	     #
#           $(NDOINC) = the -I$(NDOSolver/FiOracle directories)		     #
#           $(BDLDIR) = the directory where the Bundle/MPSolver source is    #
#									     #
#   Output: $(BDLOBJ) = the final object(s) / library			     #
#           $(BDLLIB) = external libreries + -L<libdirs>		     #
#           $(BDLH)   = the .h files to include for Bundle/MPSolver	     #
#           $(MPSH)   = the .h files to include for MPSolver (only)	     #
#           $(BDLINC) = the -I$(Bundle/MPSolver directory)		     #
#									     #
#                                VERSION 3.00				     #
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
	rm -f $(BDLDIR)*.o $(BDLDIR)*~

# macroes to be exported- - - - - - - - - - - - - - - - - - - - - - - - - - -

BDLLIB =
BDLOBJ = $(BDLDIR)Bundle.o
BDLINC = -I$(BDLDIR)
BDLH   = $(BDLDIR)Bundle.h $(BDLDIR)MPSolver.h
MPSH   = $(BDLDIR)MPSolver.h

# dependencies: every .o from its .C + every recursively included .h- - - - -

$(BDLOBJ): $(BDLDIR)Bundle.C $(BDLH) $(NDOH) $(OPTH)
	$(CC) -c $*.C -o $@ $(NDOINC) $(OPTINC) $(SW)

########################## End of makefile ###################################
