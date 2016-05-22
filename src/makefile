
all: copy-labr copy-sampler

copy-labr: labr
	cp labr ../target/labr

labr: LikelihoodSolver utils Configuration.o lrmain.cpp
	g++ -I ../lib/boost_1_53_0 -g lrmain.cpp Configuration.o LikelihoodSolver/LikelihoodSolver.o \
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

sampler: Configuration.o utils sampler-main.cpp
	g++ -I ../lib/boost_1_53_0 -g sampler-main.cpp \
		utils/FileReaderUtil.o utils/StringUtil.o utils/InputParserUtil.o \
		-o sampler

Configuration.o: Configuration.cpp
	g++ -g -c Configuration.cpp

LikelihoodSolver: force
	$(MAKE) -C LikelihoodSolver

utils: force
	$(MAKE) -C utils

force:

clean:
	$(MAKE) -C LikelihoodSolver clean
	$(MAKE) -C utils clean
	rm -f *.o labr ../target/labr*

