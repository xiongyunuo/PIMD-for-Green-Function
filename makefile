CPPX = g++
CPP_FLAG = -O2 -std=c++11

prog : main.cpp
	$(CPPX) -o prog main.cpp XRandom.cpp XDiffEq.cpp $(CPP_FLAG)

.PHONY : clean
clean :
	rm prog