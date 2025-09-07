CXX = mpicxx
CXXFLAGS = -O3 -Wall -std=c++11

HYPRE_DIR = $(PWD)/hypre/src/hypre

CXXFLAGS += -I$(HYPRE_DIR)/include
LDFLAGS  += -L$(HYPRE_DIR)/lib -lHYPRE

EXEC = diffusion

OBJS = methodes.o main.o

all: $(EXEC)

$(EXEC): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $(EXEC) $(OBJS) $(LDFLAGS)

methodes.o: methodes.cc methodes.hh
	$(CXX) $(CXXFLAGS) -c methodes.cc

main.o: main.cc methodes.hh
	$(CXX) $(CXXFLAGS) -c main.cc

clean:
	rm -f *.o $(EXEC)
