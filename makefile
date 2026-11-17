FC = gfortran
FFLAGS = -O2

TARGET = ode_test
OBJS = precision_mod.o problem_mod.o solvers_mod.o main.o

all: $(TARGET)

$(TARGET): $(OBJS)
	$(FC) $(FFLAGS) -o $@ $(OBJS)

precision_mod.o: precision_mod.f90
	$(FC) $(FFLAGS) -c $<

problem_mod.o: problem_mod.f90 precision_mod.o
	$(FC) $(FFLAGS) -c $<

solvers_mod.o: solvers_mod.f90 problem_mod.o precision_mod.o
	$(FC) $(FFLAGS) -c $<

main.o: main.f90 solvers_mod.o problem_mod.o precision_mod.o
	$(FC) $(FFLAGS) -c $<

run: $(TARGET)
	./$(TARGET)
	python3 check_ode.py

check:
	python3 check_ode.py

clean:
	rm -f $(OBJS) $(TARGET) *.mod

cleanall: clean
	rm -f rk.dat ae.dat ai.dat

.PHONY: all run check clean cleanall
