CXX = g++
OBJECTS = main.o RK4CoupledStep.o
Python = python3

Function_and_Derivative_Plot: foo Plotter.py data.dat zeros.dat
	${Python} Plotter.py

zeros.dat: ZeroDetector.py data.dat
	${Python} ZeroDetector.py

data.dat: foo
	./foo > data.dat



foo: ${OBJECTS}
	${CXX} ${OBJECTS} -o foo



clean:
	rm data.dat
	rm  *.o *.pdf
	rm   *.pdf
	rm  foo
