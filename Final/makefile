EXE = myplacer
CXX = g++
GCC = gcc

SRC = $(wildcard *.cpp)
DEP = $(wildcard *.h)
OBJ = flute_new.o

CFLAGS = -std=c++0x -Wall -O3
LFLAGS = -fopenmp -static -L./libs -lui-timer
INC = -I./includes

$(EXE): $(SRC) $(DEP) $(OBJ)	
	/bin/rm -f $(EXE) 
	$(CXX) $(CFLAGS) $(SRC) $(OBJ) $(INC) $(LFLAGS) -o $(EXE)

flute_new.o: Flute/flute.h Flute/flute_new.cpp
	/bin/rm -f flute_new.o
	$(CXX) $(CFLAGS) Flute/flute_new.cpp -c 

.PHONY: test
test: $(EXE) 
	./$(EXE) ICCAD15.parm simple/simple.iccad2015 1.0 ||true

.PHONY: benchmark
benchmark: $(EXE)
	./$(EXE) ICCAD15.parm superblue18/superblue18.iccad2015 1.0 ||true

.PHONY: clean
clean:
	/bin/rm -f $(EXE) *.gch $(OBJ)

.PHONY: backup
backup:
	@mkdir -p script-backup/`date +%m%d`
	@tar -cvzf iccad2015_evaluate_solution_`date +%m%d`.tar.gz $(SRC) $(DEP) $(OBJ) $(EXE) makefile *.parm simple Flute *.dat includes libs
	@mv iccad2015_evaluate_solution_`date +%m%d`.tar.gz script-backup/`date +%m%d`
