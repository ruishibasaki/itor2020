# compiler
CC = g++

# module name
NAME = mcnd_knap

# basic directory
DIR = ./

# debug switches
#SW = -g -Wall -fPIC -fexceptions -DNDEBUG -DIL_STD
# production switches
SW = -O3 

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

VOLIB = 
VOLOBJ = $(VOLDIR)VolVolume.o
VOLH   = $(VOLDIR)VolVolume.hpp
VOLINC = -I$(VOLDIR) 


# main module (linking phase) - - - - - - - - - - - - - - - - - - - - - - - -

# object files

$(DIR)VolVolume.o: $(DIR)VolVolume.cpp $(DIR)VolVolume.hpp
	$(CC) -c $(DIR)VolVolume.cpp -o $@ $(SW)
	
$(DIR)UtilsMethods.o: $(DIR)UtilsMethods.cpp $(DIR)UtilsMethods.hpp
	$(CC) -c $(DIR)UtilsMethods.cpp -o $@ -I$(DIR) -DIL_STD -I$(CPLEXDIR)/include -I$(CONCERTDIR)/include $(SW)	
		
$(DIR)ylp.o: $(DIR)MCND_y.cpp $(DIR)MCND_y.hpp $(DIR)UtilsMethods.hpp $(DIR)mcnd.hpp
	$(CC) -c $(DIR)MCND_y.cpp -o $@ -I$(DIR) -DIL_STD -I$(CPLEXDIR)/include -I$(CONCERTDIR)/include $(SW)	

$(DIR)mcnd.o: $(DIR)mcnd.cpp  $(DIR)mcnd.hpp $(DIR)structures.hpp 
	$(CC) -c $(DIR)mcnd.cpp -o $@ -I$(DIR) $(SW)
	
$(DIR)main.o: $(DIR)main.cpp  $(DIR)mcnd.hpp $(DIR)structures.hpp 
	$(CC) -c $(DIR)main.cpp -o $@ -I$(DIR)  -DIL_STD -I$(CPLEXDIR)/include -I$(CONCERTDIR)/include $(SW)
	

$(NAME): $(DIR)mcnd.o $(DIR)ylp.o  $(DIR)UtilsMethods.o $(DIR)VolVolume.o $(DIR)main.o
	$(CC) -o $(NAME) $(DIR)main.o $(DIR)ylp.o $(DIR)UtilsMethods.o $(DIR)mcnd.o $(DIR)VolVolume.o $(CPLEXFLAGS) $(SW)

	#rm -f $(DIR)*.o $(DIR)*~
