SAT=../modules/sat-engine
DFA=../modules/dfa-sat
LOCAL=/opt/homebrew/include

CC=g++
CCFLAGS=-I$(LOCAL) -I$(SAT)/include -I$(DFA)/src -Wall -O3 -std=c++14

strips:
	$(CC) $(CCFLAGS) -o strips strips.cc

clean:
	rm -f strips *~

strips:	strips.cc strips.h utils.h
strips:	$(SAT)/include/sat/encoder/sat.h $(SAT)/include/sat/encoder/theory.h
strips:	$(DFA)/src/dfa.h

