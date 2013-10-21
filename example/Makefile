HIGHER_ORDER_DIR = ../include
FUSION_DIR = ../fusion-move
QPBO_DIR = ../qpbo

CXX ?= g++
DEFS =
OPT ?= -O3
INCLUDES = -I$(QPBO_DIR) -I$(HIGHER_ORDER_DIR) -I$(FUSION_DIR)
CXX_FLAGS = --std=c++11 $(OPT) -Wall $(INCLUDES) $(DEFS)
LD_FLAGS = 
LIBS = -lopencv_core -lopencv_imgproc -lopencv_highgui

EXAMPLE_SRCS = higher-order-example.cpp \
			   foe-cliques.cpp
EXAMPLE_OBJS = $(EXAMPLE_SRCS:.cpp=.o)

QPBO_SRCS = $(QPBO_DIR)/QPBO.cpp
QPBO_OBJS = $(QPBO_SRCS:.cpp=.o)

TARGET = higher-order-example

.PHONY: all
all: $(TARGET)

$(TARGET): $(EXAMPLE_OBJS) $(QPBO_OBJS)
	$(CXX) $(CXX_FLAGS) $(LD_FLAGS) -o $@ $(EXAMPLE_OBJS) $(QPBO_OBJS) $(LIBS)

%.o: %.cpp
	$(CXX) $(CXX_FLAGS) -o $@ -c $<

.PHONY: test
test: $(TARGET)
	./$(TARGET) medium-test.pgm medium-out.pgm

.PHONY: clean
clean: 
	rm -f $(EXAMPLE_OBJS)
	rm -f $(QPBO_OBJS)
	rm -f *.o
	rm -f *.P
	rm -f *.d

.PHONY: distclean
distclean: clean
	rm -f $(TARGET)
