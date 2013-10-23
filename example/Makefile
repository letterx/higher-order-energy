HIGHER_ORDER_DIR = ../include
QPBO_DIR = ../qpbo

CXX ?= g++
DEFS =
INCLUDES = -I$(QPBO_DIR) -I$(HIGHER_ORDER_DIR)
CXX_FLAGS = -std=c++11 -O3 -Wall $(INCLUDES) $(DEFS)
LD_FLAGS = 
LIBS = 

EXAMPLE_SRCS = higher-order-example.cpp \
			   image.cpp \
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
	./$(TARGET) medium-test.pgm medium-test-out.pgm

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
