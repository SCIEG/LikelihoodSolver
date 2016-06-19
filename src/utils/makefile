
CXX := g++
CXXFLAGS := -g -O2 -fPIC -c
INCLUDES := -I ../../lib/boost_1_53_0

all: DebugUtil.o FileReaderUtil.o LikelihoodUtil.o ProbabilityUtil.o InputParserUtil.o \
     StringUtil.o

DebugUtil.o: DebugUtil.cpp DebugUtil.h
	$(CXX) $(INCLUDES) $(CXXFLAGS) DebugUtil.cpp

FileReaderUtil.o: FileReaderUtil.cpp FileReaderUtil.h
	$(CXX) $(INCLUDES) $(CXXFLAGS) FileReaderUtil.cpp

LikelihoodUtil.o: LikelihoodUtil.cpp LikelihoodUtil.h
	$(CXX) $(INCLUDES) $(CXXFLAGS) LikelihoodUtil.cpp

ProbabilityUtil.o: ProbabilityUtil.cpp ProbabilityUtil.h
	$(CXX) $(INCLUDES) $(CXXFLAGS) ProbabilityUtil.cpp

InputParseUtil.o: InputParserUtil.cpp InputParserUtil.h
	$(CXX) $(INCLUDES) $(CXXFLAGS) InputParserUtil.cpp

StringUtil.o: StringUtil.cpp StringUtil.h
	$(CXX) $(INCLUDES) $(CXXFLAGS) StringUtil.cpp

clean:
	rm -f *.o
