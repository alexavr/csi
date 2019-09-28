################################################
# Change this part
NETCDF=/opt/netcdf4-serial
FC = ifort
################################################

################################################
# Do not edit part
# (unless you know what you are doing)

# Program Name
PROG = csi.exe

# Scr folder
VPATH = ./src/

# Compiler Flags
FFLAGS = -I$(NETCDF)/include 
FLINK = -O3 -L$(NETCDF)/lib -lnetcdf -lnetcdff
# LINKER = $(FC) -traceback -o 
LINKER = $(FC) -o 

# Object files
OBJS = gaussian_filter.o module_globals.o module_io.o module_csi.o csi.o

LPT: $(PROG)

# Create the LPT
$(PROG): $(OBJS)
	@echo "--------------------------------------"
	@echo "Creating the executable for the CSI"
	@echo "--------------------------------------"
	$(LINKER) $(PROG) $(OBJS) $(FLINK)
	mv *.o *.mod $(VPATH)

%.o: %.f90
	@echo "--------------------------------------"
	@echo "Compiling $<"
	@echo "--------------------------------------"
	$(FC) -c $(FFLAGS) $<
	
# Clean up everything
clean:
	@echo "--------------------------------------"
	@echo "Cleaning everything up in the CSI"
	@echo "--------------------------------------"
	rm -f $(VPATH)*.o $(VPATH)*.mod $(VPATH)*.exe $(PROG)

gaussian_filter.o   : gaussian_filter.f90
module_globals.o    : module_globals.f90
module_io.o         : module_io.f90 module_globals.o
module_csi.o        : module_globals.o gaussian_filter.o
csi.o               : csi.f90 module_io.o module_csi.o
