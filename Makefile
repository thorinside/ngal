CC=g++
CCFLAGS=-Wall -O3 -s
VER=1.9.1
MESSAGE="SGAL Simple Genetic Algorithm Library 1.9.0"
PROGS=gal gal2 gal3 opt
HEADERS=ACG.h RNG.h chrom.h pool.h operator.h engine.h eval.h pRandom.h ga.h
SOURCES=ACG.cc RNG.cc pool.cc operator.cc engine.cc pRandom.cc opt.cc gal.cc gal2.cc gal3.cc
OBJS=ACG.o RNG.o pool.o operator.o engine.o pRandom.o
LIBS=-lm

all: $(PROGS)

gal: $(OBJS) gal.o
	$(CC) $(CCFLAGS) -o $@ $@.o $(OBJS) $(LIBS)

gal2: $(OBJS) gal2.o
	$(CC) $(CCFLAGS) -o $@ $@.o $(OBJS) $(LIBS)

gal3: $(OBJS) gal3.o
	$(CC) $(CCFLAGS) -o $@ $@.o $(OBJS) $(LIBS)

opt: $(OBJS) opt.o
	$(CC) $(CCFLAGS) -o $@ $@.o $(OBJS) $(LIBS)

$(HEADERS):
	co $@

.cc.o:
	$(CC) $(CCFLAGS) -c $<

depend dep: co
	$(CC) -M *.cc > .depend

checkout co:
	co $(HEADERS)
	co $(SOURCES)
	co copyright

editall:
	co -l $(HEADERS) $(SOURCES) copyright Makefile

saveall:
	ci -f$(VER) $(HEADERS) $(SOURCES) copyright Makefile

test:
	echo $(REVISION)

clean:
	rm -f $(PROGS) *.o

dist: clean co
	mkdir ngal-$(VER)
	cp Makefile ngal-$(VER)/
	cp $(HEADERS) ngal-$(VER)/
	cp $(SOURCES) ngal-$(VER)/
	cp copyright ngal-$(VER)/
	chmod 0755 ngal-$(VER)
	chmod 0644 ngal-$(VER)/*
	tar zcvf ngal-$(VER).tar.gz ngal-$(VER)
	rm -rf ngal-$(VER)

# This dummy stuff makes sure dependencies are included
# To keep everything up to date.

dummy:
ifeq (.depend,$(wildcard .depend))
include .depend
endif
