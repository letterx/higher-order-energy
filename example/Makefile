HIGHER_ORDER_DIR = ../include
QPBO_DIR = ../qpbo
CLPLIBDIR = ~/src/coin-Clp/lib/

CXX ?= g++
DEFS =
INCLUDES = -I$(QPBO_DIR) -I$(HIGHER_ORDER_DIR) -I../hocr -I../grd/source/library/
OPT ?= -O3
CXX_FLAGS = $(OPT) -Wall --std=c++11 $(INCLUDES) $(DEFS)
LD_FLAGS = -L../grd/obj/submodular/ -L $(CLPLIBDIR)
LIBS = -lboost_program_options -lpetter -lClp -lCoinUtils -llapack

EXAMPLE_SRCS = higher-order-example.cpp \
			   image.cpp \
			   foe-cliques.cpp
EXAMPLE_OBJS = $(EXAMPLE_SRCS:.cpp=.o)

QPBO_SRCS = $(QPBO_DIR)/QPBO.cpp
QPBO_OBJS = $(QPBO_SRCS:.cpp=.o)

MAXFLOW_OBJS = ../grd/obj/submodular/graph.o ../grd/obj/submodular/maxflow.o

SRCS = $(EXAMPLE_SRCS) $(QPBO_SRCS) foe-9.cpp test-pairwise.cpp

LIBPETTER = ../grd/obj/submodular/libpetter.a

TARGET = higher-order-example

.PHONY: all
all: $(TARGET) foe-9 test-pairwise

$(TARGET): $(EXAMPLE_OBJS) $(QPBO_OBJS) $(LIBPETTER)
	LD_RUN_PATH=$$LD_RUN_PATH:$(CLPLIBDIR) \
	$(CXX) $(CXX_FLAGS) $(LD_FLAGS) -o $@ $(EXAMPLE_OBJS) $(QPBO_OBJS) $(MAXFLOW_OBJS) $(LIBS)

test-pairwise: test-pairwise.o $(QPBO_OBJS) $(LIBPETTER)
	LD_RUN_PATH=$$LD_RUN_PATH:$(CLPLIBDIR) \
	$(CXX) $(CXX_FLAGS) $(LD_FLAGS) -o $@ test-pairwise.o $(QPBO_OBJS) $(MAXFLOW_OBJS) $(LIBS)

foe-9: foe-9.o image.o foe-cliques.o $(QPBO_OBJS) $(LIBPETTER)
	LD_RUN_PATH=$$LD_RUN_PATH:$(CLPLIBDIR) \
	$(CXX) $(CXX_FLAGS) $(LD_FLAGS) -o $@ foe-9.o image.o foe-cliques.o $(QPBO_OBJS) $(MAXFLOW_OBJS) $(LIBS)

$(LIBPETTER): $(QPBO_OBJS)
	cd ../grd; make

%.o: %.cpp
	$(CXX) $(CXX_FLAGS) -MMD -o $@ -c $<
	@cp $*.d $*.P; \
	 sed -e 's/#.*//' -e 's/^[^:]*: *//' -e 's/ *\\$$//' \
    	 -e '/^$$/ d' -e 's/$$/ :/' < $*.d >> $*.P; \
	 rm $*.d

-include $(SRCS:.cpp=.P)

.PHONY: clean
clean: 
	rm -f $(EXAMPLE_OBJS)
	rm -f $(QPBO_OBJS)
	rm -f *.o
	rm -f *.P
	rm -f *.d
	cd ../grd; make clean

.PHONY: distclean
distclean: clean
	rm -f $(TARGET)
