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

MAIN_SRCS = random-graph.cpp
MAIN_OBJS = $(MAIN_SRCS:.cpp=.o)

QPBO_SRCS = $(QPBO_DIR)/QPBO.cpp
QPBO_OBJS = $(QPBO_SRCS:.cpp=.o)

MAXFLOW_OBJS = ../grd/obj/submodular/maxflow.o

SRCS = $(MAIN_SRCS) $(QPBO_SRCS)

LIBPETTER = ../grd/obj/submodular/libpetter.a

TARGET = random-graph

.PHONY: all
all: $(TARGET) 

$(TARGET): $(MAIN_OBJS) $(QPBO_OBJS) $(LIBPETTER)
	LD_RUN_PATH=$$LD_RUN_PATH:$(CLPLIBDIR) \
	$(CXX) $(CXX_FLAGS) $(LD_FLAGS) -o $@ $(MAIN_OBJS) $(QPBO_OBJS) $(MAXFLOW_OBJS) $(LIBS)

$(LIBPETTER): $(QPBO_OBJS)
	cd ../grd; $(MAKE)

%.o: %.cpp
	$(CXX) $(CXX_FLAGS) -MMD -o $@ -c $<
	@cp $*.d $*.P; \
	 sed -e 's/#.*//' -e 's/^[^:]*: *//' -e 's/ *\\$$//' \
    	 -e '/^$$/ d' -e 's/$$/ :/' < $*.d >> $*.P; \
	 rm $*.d

-include $(SRCS:.cpp=.P)

.PHONY: clean
clean: 
	rm -f $(OBJS)
	rm -f *.o
	rm -f *.P
	rm -f *.d
	cd ../grd; make clean

.PHONY: distclean
distclean: clean
	rm -f $(TARGET)
