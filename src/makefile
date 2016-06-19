
CXX := g++
CXXFLAGS := -g -O2 -fPIC
INCLUDES := -I ../lib/boost_1_53_0

all: copy-labr copy-sampler

copy-labr: labr
	cp labr ../target/labr

labr: deps lrmain.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) lrmain.cpp Configuration.o \
		LikelihoodSolver/LikelihoodSolver.o \
		LikelihoodSolver/CachingSolver.o \
		LikelihoodSolver/NoSuspectFourUnknownsLikelihoodSolver.o \
		LikelihoodSolver/NoSuspectOneUnknownLikelihoodSolver.o \
		LikelihoodSolver/NoSuspectThreeUnknownsLikelihoodSolver.o \
		LikelihoodSolver/NoSuspectTwoUnknownsLikelihoodSolver.o \
		LikelihoodSolver/OneSuspectNoUnknownsLikelihoodSolver.o \
		LikelihoodSolver/OneSuspectOneUnknownLikelihoodSolver.o \
		LikelihoodSolver/OneSuspectThreeUnknownsLikelihoodSolver.o \
		LikelihoodSolver/OneSuspectTwoUnknownsLikelihoodSolver.o \
		LikelihoodSolver/UnknownsSolver.o \
		utils/LikelihoodUtil.o utils/FileReaderUtil.o utils/ProbabilityUtil.o \
		utils/DebugUtil.o utils/InputParserUtil.o utils/StringUtil.o\
		-o labr

copy-sampler: sampler
	cp sampler ../target/sampler

sampler: deps sampler-main.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) sampler-main.cpp Configuration.o \
		LikelihoodSolver/LikelihoodSolver.o LikelihoodSolver/UnknownsSolver.o \
		utils/FileReaderUtil.o utils/StringUtil.o utils/InputParserUtil.o \
		utils/ProbabilityUtil.o \
		-o sampler

deps: Configuration.o LikelihoodSolver utils

Configuration.o: Configuration.cpp
	$(CXX) $(CXXFLAGS) -c Configuration.cpp

LikelihoodSolver: force
	$(MAKE) -C LikelihoodSolver

utils: force
	$(MAKE) -C utils

force:

clean:
	$(MAKE) -C LikelihoodSolver clean
	$(MAKE) -C utils clean
	rm -f *.o labr ../target/labr*

