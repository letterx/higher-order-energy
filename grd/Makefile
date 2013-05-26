# Petter Strandmark 2011

# You might have to change these to g++-4.5 or something similar if you have several versions installed
COMPILER = g++ -c
LINKER = g++

# Temporary directory to store object files; needs to exist
OBJDIR = ./obj/submodular/

# Location of Clp if not installed where gcc can find it
CLPLIBDIR = ~/src/coin-Clp/lib/
CLPINCLUDEDIR = ~/src/coin-Clp/include/

# Names of executables
TARGET  = bin/submodular

##################################################################################

SUBMOBJ = $(OBJDIR)PseudoBoolean.o $(OBJDIR)PseudoBoolean_create_g.o $(OBJDIR)PseudoBoolean_heuristic.o $(OBJDIR)PseudoBoolean_minimize.o $(OBJDIR)PseudoBoolean_minimize_lp.o $(OBJDIR)PseudoBoolean_minimize_reduction.o $(OBJDIR)PseudoBoolean_reduce.o $(OBJDIR)PseudoBoolean_complete.o $(OBJDIR)PseudoBoolean_branchandbound.o $(OBJDIR)GeneratorPseudoBoolean.o $(OBJDIR)VertexPacking.o
QPBOOBJ = ../qpbo/QPBO.o 
MAXFLOWOBJ = $(OBJDIR)graph.o $(OBJDIR)maxflow.o
LIBDIR = $(OBJDIR)
INCLUDE = -I source/library -I thirdparty/Petter -I thirdparty/HOCR -I ../qpbo -I thirdparty/maxflow-v3.01.src -I thirdparty/ibfs
OPTIONS = -std=c++0x -O3




all : $(TARGET)

$(TARGET): $(OBJDIR)main_program.o $(OBJDIR)submodular_tests.o  $(OBJDIR)libpetter.a $(QPBOOBJ) $(MAXFLOWOBJ)
	mkdir -p bin
	LD_RUN_PATH=$$LD_RUN_PATH:$(CLPLIBDIR) \
	$(LINKER) -L $(LIBDIR) -L $(CLPLIBDIR) $(OBJDIR)main_program.o  $(OBJDIR)submodular_tests.o  $(QPBOOBJ) $(MAXFLOWOBJ) -o $(TARGET) -lpetter -lClp -lCoinUtils -llapack

clean:
	rm -f $(TARGET)
	rm -f $(OBJDIR)*.o
	rm -f $(OBJDIR)*.la
	rm -f $(OBJDIR)*.a

test: $(TARGET)
	./$(TARGET)

# Main program; compile with -Wall -pedantic
$(OBJDIR)main_program.o: source/main_program.cpp source/library/PseudoBoolean.h source/library/Posiform.h
	mkdir -p $(OBJDIR)
	$(COMPILER) $(OPTIONS) -Wall -pedantic $(INCLUDE) source/main_program.cpp -o $(OBJDIR)main_program.o
$(OBJDIR)submodular_tests.o: source/submodular_tests.cpp source/library/PseudoBoolean.h
	mkdir -p $(OBJDIR)
	$(COMPILER) $(OPTIONS) -Wall -pedantic $(INCLUDE) source/submodular_tests.cpp -o $(OBJDIR)submodular_tests.o
#$(OBJDIR)posiform_tests.o: source/posiform_tests.cpp source/library/PseudoBoolean.h source/library/Posiform.h source/library/VertexPacking.h
#	$(COMPILER) $(OPTIONS) -Wall -Werror -pedantic-errors $(INCLUDE) source/posiform_tests.cpp -o $(OBJDIR)posiform_tests.o
# Main program; -pedantic-errors not possible since it includes QPBO.h
$(OBJDIR)image_denoising.o: source/image_denoising.cpp
	mkdir -p $(OBJDIR)
	$(COMPILER) $(OPTIONS) -Wall -Werror $(INCLUDE) source/image_denoising.cpp -o $(OBJDIR)image_denoising.o 

$(LIBDIR)libpetter.a : $(SUBMOBJ) $(OBJDIR)Petter-Color.o
	mkdir -p $(LIBDIR)
	ar rcs $(LIBDIR)libpetter.a $(OBJDIR)Petter-Color.o $(SUBMOBJ)

# Library
$(OBJDIR)Petter-Color.o: thirdparty/Petter/Petter-Color.cc thirdparty/Petter/Petter-Color.h
	mkdir -p $(OBJDIR)
	$(COMPILER) $(OPTIONS) -Wall -pedantic $(INCLUDE) thirdparty/Petter/Petter-Color.cc -o $(OBJDIR)Petter-Color.o	

# Library
$(OBJDIR)%.o: source/library/%.cpp source/library/PseudoBoolean.h
	mkdir -p $(OBJDIR)
	$(COMPILER) $(OPTIONS) -Wall $(INCLUDE) -I $(CLPINCLUDEDIR)  $< -o $@

../qpbo/QPBO.o: ../qpbo/QPBO.cpp
	$(COMPILER) $(OPTIONS) $(INCLUDE) $< -o $@

# Will generate warnings in g++ 4.5
$(OBJDIR)%.o: thirdparty/maxflow-v3.01.src/%.cpp
	mkdir -p $(OBJDIR)
	$(COMPILER) $(OPTIONS) $(INCLUDE) $< -o $@


