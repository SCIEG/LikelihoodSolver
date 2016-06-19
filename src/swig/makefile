

all: libLabRetriever.so

libLabRetriever.so : deps lr-swig_wrap.cxx
	g++ -O2 -c -fPIC lr-swig.cpp lr-swig_wrap.cxx
	g++ -I ../../lib/boost_1_53_0 -shared lr-swig.o lr-swig_wrap.o \
		../LikelihoodSolver/LikelihoodSolver.o \
		../LikelihoodSolver/UnknownsSolver.o \
		../utils/*.o ../*.o -o libLabRetriever.so

clean:
	rm -f *~ *.o *.so

clean-all:
	rm -f *~ StringVector.cs StringVectorVector.cs LabRetriever.cs LabRetrieverPINVOKE.cs *.o \
		*_wrap.cxx *.so

deps:
	cd .. && make deps

swig: lr-swig.i
	swig -csharp -c++ lr-swig.i

