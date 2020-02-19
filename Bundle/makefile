##############################################################################
################################ makefile ####################################
##############################################################################
#									     #
#   makefile of TestNDO							     #
#									     #
#                                VERSION 2.00				     #
#                	        20 - 04 - 2005				     #
#									     #
# 		               Implementation by:			     #
#									     #
#			       Antonio Frangioni			     #
#   			   Operations Research Group			     #
#			  Dipartimento di Informatica			     #
#   			     Universita' di Pisa			     #
#									     #
##############################################################################

# module name
NAME = BdlKnap	

# basic directory
DIR = ./

# debug switches
#SW = -g -Wall
# production switches
SW = -O3

# compiler
CC = g++

# default target- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

default: $(NAME)

# clean - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

clean::
	rm -f $(DIR)*.o $(DIR)*~ $(NAME)

# define & include the necessary modules- - - - - - - - - - - - - - - - - - -
# if a module is not used in the current configuration, just comment out the
# corresponding include line
# each module outputs some macros to be used here:
# *OBJ is the final object(s) / library
# *LIB is the external libreries + -L<libdirs>
# *H   is the list of all include files
# *INC is the -I<include directories>

# OPT****.h
OPTH	= $(DIR)OPTUtils.h $(DIR)OPTvect.h $(DIR)OPTtypes.h
OPTINC	= -I$(DIR)

# NDOSlver
NDODIR = ./NDOSlver/
include $(NDODIR)makefile

# Bundle
BDLDIR = ./Bundle/
include $(BDLDIR)makefile


# QPPnltMP
QMPDIR = ./QPPnltMP/
include $(QMPDIR)makefile


#SYSTEM = x86-64_osx
#LIBFORMAT = static_pic
#CONCERTDIR = /Applications/CPLEX_Studio128/concert
#CPLEXDIR = 	/Applications/CPLEX_Studio128/cplex
#CPLEXLIBDIR = $(CPLEXDIR)/lib/$(SYSTEM)/$(LIBFORMAT)
#CONCERTLIBDIR = $(CONCERTDIR)/lib/$(SYSTEM)/$(LIBFORMAT)
#CPLEXFLAGS =   -L$(CPLEXLIBDIR) -lilocplex -lcplex -L$(CONCERTLIBDIR) -lconcert -lm -ldl -lpthread

SYSTEM = x86-64_linux
LIBFORMAT = static_pic
CONCERTDIR = /opt/ibm128/ILOG/CPLEX_Studio128/concert
CPLEXDIR = 	/opt/ibm128/ILOG/CPLEX_Studio128/cplex
CPLEXLIBDIR = $(CPLEXDIR)/lib/$(SYSTEM)/$(LIBFORMAT)
CONCERTLIBDIR = $(CONCERTDIR)/lib/$(SYSTEM)/$(LIBFORMAT)
CPLEXFLAGS =   -L$(CPLEXLIBDIR) -lilocplex -lcplex -L$(CONCERTLIBDIR) -lconcert -lm -ldl -lpthread

# main module (linking phase) - - - - - - - - - - - - - - - - - - - - - - - -

# object files
MOBJ =	$(BDLOBJ)\
	$(OMPOBJ)\
	$(QMPOBJ)\
	UtilsMethods.o \
	MCND_y.o \
	TestFi.o 

# libreries
MLIB =	$(BDLLIB)\
	$(OMPLIB)\
	$(QMPLIB) 

$(NAME): $(MOBJ) $(DIR)Main.o
	$(CC) -o $(NAME) $(DIR)Main.o $(MOBJ) $(MLIB) $(CPLEXFLAGS) $(SW)
	#rm -f $(DIR)*.o $(DIR)*~
# dependencies: every .o from its .C + every recursively included .h- - - - -

# include dirs
MINC =	$(NDOINC)\
	$(BDLINC)\
	$(QMPINC)\
	 -I$(DIR)

# includes
MH =	$(NDOH)\
	$(QMPH)\
	$(BDLH)\
	$(OPTH) \
	MCND_y.hpp \
	structures.hpp \
	UtilsMethods.hpp \
	TestFi.h 
	
$(TSFDIR)TestFi.o: $(TSFDIR)TestFi.C $(TSFH) $(FIOH) $(OPTH)
	$(CC) -c $*.C -o $@  $(MINC)  $(SW)

$(DIR)Main.o: $(DIR)Main.C $(MH)
	$(CC) -c $*.C -o $@ $(MINC) -DIL_STD -I$(CPLEXDIR)/include -I$(CONCERTDIR)/include $(SW)
	
$(DIR)UtilsMethods.o: $(DIR)UtilsMethods.cpp $(DIR)UtilsMethods.hpp
	$(CC) -c $(DIR)UtilsMethods.cpp -o $@ -DIL_STD -I$(CPLEXDIR)/include -I$(CONCERTDIR)/include $(SW)	


$(DIR)MCND_y.o: $(DIR)MCND_y.cpp $(DIR)MCND_y.hpp 
	$(CC) -c $(DIR)MCND_y.cpp -o $@  $(MINC)  -DIL_STD -I$(CPLEXDIR)/include -I$(CONCERTDIR)/include $(SW)	


############################ End of makefile #################################
