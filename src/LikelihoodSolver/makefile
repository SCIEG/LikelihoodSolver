
CXX := g++
CXXFLAGS := -g -O2 -fPIC -c
INCLUDES := -I ../../lib/boost_1_53_0

all: CachingSolver.o NoSuspectFourUnknownsLikelihoodSolver.o \
            NoSuspectOneUnknownLikelihoodSolver.o \
            NoSuspectThreeUnknownsLikelihoodSolver.o \
            NoSuspectTwoUnknownsLikelihoodSolver.o \
            OneSuspectNoUnknownsLikelihoodSolver.o \
            OneSuspectOneUnknownLikelihoodSolver.o \
            OneSuspectThreeUnknownsLikelihoodSolver.o \
            OneSuspectTwoUnknownsLikelihoodSolver.o \
            UnknownsSolver.o LikelihoodSolver.o

CachingSolver.o: CachingSolver.cpp
	$(CXX) $(INCLUDES) $(CXXFLAGS) CachingSolver.cpp

NoSuspectFourUnknownsLikelihoodSolver.o: NoSuspectFourUnknownsLikelihoodSolver.cpp
	$(CXX) $(INCLUDES) $(CXXFLAGS) NoSuspectFourUnknownsLikelihoodSolver.cpp

NoSuspectOneUnknownLikelihoodSolver.o: NoSuspectOneUnknownLikelihoodSolver.cpp
	$(CXX) $(INCLUDES) $(CXXFLAGS) NoSuspectOneUnknownLikelihoodSolver.cpp

NoSuspectThreeUnknownsLikelihoodSolver.o: NoSuspectThreeUnknownsLikelihoodSolver.cpp
	$(CXX) $(INCLUDES) $(CXXFLAGS) NoSuspectThreeUnknownsLikelihoodSolver.cpp

NoSuspectTwoUnknownsLikelihoodSolver.o: NoSuspectTwoUnknownsLikelihoodSolver.cpp
	$(CXX) $(INCLUDES) $(CXXFLAGS) NoSuspectTwoUnknownsLikelihoodSolver.cpp

OneSuspectNoUnknownsLikelihoodSolver.o: OneSuspectNoUnknownsLikelihoodSolver.cpp
	$(CXX) $(INCLUDES) $(CXXFLAGS) OneSuspectNoUnknownsLikelihoodSolver.cpp

OneSuspectOneUnknownLikelihoodSolver.o: OneSuspectOneUnknownLikelihoodSolver.cpp
	$(CXX) $(INCLUDES) $(CXXFLAGS) OneSuspectOneUnknownLikelihoodSolver.cpp

OneSuspectThreeUnknownsLikelihoodSolver.o: OneSuspectThreeUnknownsLikelihoodSolver.cpp
	$(CXX) $(INCLUDES) $(CXXFLAGS) OneSuspectThreeUnknownsLikelihoodSolver.cpp

OneSuspectTwoUnknownsLikelihoodSolver.o: OneSuspectTwoUnknownsLikelihoodSolver.cpp
	$(CXX) $(INCLUDES) $(CXXFLAGS) OneSuspectTwoUnknownsLikelihoodSolver.cpp

UnknownsSolver.o: UnknownsSolver.cpp
	$(CXX) $(INCLUDES) $(CXXFLAGS) UnknownsSolver.cpp

LikelihoodSolver.o: LikelihoodSolver.cpp
	$(CXX) $(INCLUDES) $(CXXFLAGS) LikelihoodSolver.cpp

clean:
	rm -f *.o
