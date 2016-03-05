
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
		utils/DebugUtil.o \
		-o labr
	cp labr ../target/labr

Configuration.o: Configuration.cpp
	g++ -g -c Configuration.cpp

LikelihoodSolver: force
	$(MAKE) -C LikelihoodSolver

utils: force
	$(MAKE) -C utils

force:
	true

clean:
	rm -rf *.o */*.o labr ../target/labr*

