# configuration -- change this to install in different directories,
#                  or if "make" doesn't work.
# prefix = /usr/local

#########################
#       TARGET          #
#########################

TARGET= reliability
SUF90=f90
SUF77=f

.SUFFIXES: .f90 .f .o .mod

GSL_prefix = /usr/local

#########################
#      COMPILATION      #
#########################

FAD	= mpif90
F90	= mpif90
F77	= mpif77


CC = gcc
FC = ifort

FFLAGS  = -r8 -openmp -O4  #-debug extended -ftrapuv

CFLAGS = -fPIC -Wall -O0 -g -I$(GSL_prefix)/include
LFLAGS =  -L$(GSL_prefix)/lib -lgsl -lgslcblas -lm
LIBS = -ldl -lstdc++

#export:
#    ar rvs pcestimate.a *.o

SRCS = dimrel.o mpi.o main.o randomroutines.o

OBJS =  ${SRCS:.$(SUF)=.o}

all:  $(TARGET)

$(TARGET): $(OBJS) 
	$(F90) $(FFLAGS) -o $(TARGET) $(OBJS)  $(LFLAGS) -Wl,-rpath=.
	@echo " ----------- ${TARGET} created ----------- "

######################################
####### Compilation
######################################
%.o : %.mod

.$(SUF90).o:
	$(F90) $(FFLAGS) -c $<
.$(SUF77).o:
	$(F77) $(FFLAGS) -c $<

##################################
# Clean
##################################

clean:
	rm -f *.o *~ *.so pc /output/output*
