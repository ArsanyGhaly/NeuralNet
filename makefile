.PHONY: NeuralNetDriver NeuralNetEnsembleDriver

DATAFILE1=test1.txt
DATAFILE2=test2.txt

CC 	= g++
CFLAGS 	= -g -O0 -std=c++11

AA	= NeuralNetDriver.exe
BB	= NeuralNetArchitectureDriver.exe


default: NeuralNetDriver


NeuralNetDriver:	Neurons.o  NeuralNet.o NeuralNetDriver.o 
		$(CC) $(CFLAGS) -o $(AA) Neurons.o NeuralNetDriver.o NeuralNet.o
		./NeuralNetDriver $(DATAFILE1)

NeuralNetArchitectureDriver:	Neurons.o NeuralNet.o NeuralNetArchitecture.o NeuralNetArchitectureDriver.o
		$(CC) $(CFLAGS) -o $(BB) Neurons.o NeuralNet.o NeuralNetArchitecture.o NeuralNetArchitectureDriver.o
		./NeuralNetArchitectureDriver $(DATAFILE2)


clean:
	$(RM) *.o *.exe

Neurons.o:	Neurons.cpp Neurons.h
			$(CC) $(CFLAGS) -c Neurons.cpp

NeuralNet.o: NeuralNet.cpp NeuralNet.h Neurons.h
			$(CC) $(CFLAGS) -c NeuralNet.cpp

NeuralNetArchitecture.o: NeuralNetArchitecture.cpp NeuralNet.h Neurons.h
			$(CC) $(CFLAGS) -c NeuralNetArchitecture.cpp

NeuralNetDriver.o:	NeuralNetDriver.cpp NeuralNet.h
			$(CC) $(CFLAGS) -c NeuralNetDriver.cpp

NeuralNetArchitectureDriver.o: NeuralNetArchitectureDriver.cpp NeuralNet.h Neurons.h NeuralNetArchitecture.h
			$(CC) $(CFLAGS) -c NeuralNetArchitectureDriver.cpp
