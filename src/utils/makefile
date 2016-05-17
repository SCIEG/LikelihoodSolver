all: DebugUtil.o FileReaderUtil.o LikelihoodUtil.o ProbabilityUtil.o InputParserUtil.o

DebugUtil.o: DebugUtil.cpp DebugUtil.h
	g++ -I ../../lib/boost_1_53_0 -g -c DebugUtil.cpp

FileReaderUtil.o: FileReaderUtil.cpp FileReaderUtil.h
	g++ -I ../../lib/boost_1_53_0 -g -c FileReaderUtil.cpp

LikelihoodUtil.o: LikelihoodUtil.cpp LikelihoodUtil.h
	g++ -I ../../lib/boost_1_53_0 -g -c LikelihoodUtil.cpp

ProbabilityUtil.o: ProbabilityUtil.cpp ProbabilityUtil.h
	g++ -I ../../lib/boost_1_53_0 -g -c ProbabilityUtil.cpp

InputParseUtil.o: InputParserUtil.cpp InputParserUtil.h
	g++ -I ../../lib/boost_1_53_0 -g -c InputParserUtil.cpp

StringUtil.o: StringUtil.cpp StringUtil.h
	g++ -I ../../lib/boost_1_53_0 -g -c StringUtil.cpp

clean:
	rm -f *.o
