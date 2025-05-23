#############################################################################
# Makefile for building: prank
#############################################################################

####### Compiler, tools and options

CC            = gcc
CXX           = g++
DEFINES       = -DNDEBUG
CFLAGS        = -m64 -pipe -O3 $(DEFINES)
CXXFLAGS      = -m64 -pipe -O3 $(DEFINES)
INCPATH       = -I. -I/usr/include
LINK          = g++
LFLAGS        = -m64
LIBS          = $(SUBLIBS)    
AR            = ar cqs
RANLIB        = 
TAR           = tar -cf
COMPRESS      = gzip -9f
COPY          = cp -f
SED           = sed
COPY_FILE     = $(COPY)
COPY_DIR      = $(COPY) -r
STRIP         = strip
INSTALL_FILE  = install -m 644 -p
INSTALL_DIR   = $(COPY_DIR)
INSTALL_PROGRAM = install -m 755 -p
DEL_FILE      = rm -f
SYMLINK       = ln -f -s
DEL_DIR       = rmdir
MOVE          = mv -f
CHK_DIR_EXISTS= test -d
MKDIR         = mkdir -p
POD2MAN       = pod2man --utf8 --center "The Probabilistic Alignment Kit" -release="v.121211"

####### Output directory

OBJECTS_DIR   = ./

####### Files

SOURCES       = writefile.cpp \
		treenode.cpp \
		translatesequences.cpp \
		terminalsequence.cpp \
		terminalnode.cpp \
		site.cpp \
		sequence.cpp \
		readnewick.cpp \
		readfile.cpp \
		readalignment.cpp \
		pwsite.cpp \
		pwhirschberg.cpp \
		progressivealignment.cpp \
		prank.cpp \
		postprobability.cpp \
		phylomatchscore.cpp \
		node.cpp \
		intmatrix.cpp \
		hmmodel.cpp \
		hirschberg.cpp \
		guidetree.cpp \
		fullprobability.cpp \
		flmatrix.cpp \
		eigen.cpp \
		dbmatrix.cpp \
		characterprobability.cpp \
		boolmatrix.cpp \
		ancestralsequence.cpp \
		ancestralnode.cpp \
		check_version.cpp \
		exonerate_reads.cpp \
		mafft_alignment.cpp \
		bppancestors.cpp \
		raxmlancestors.cpp \
		fasttree_tree.cpp \
		raxmlrebl.cpp 
OBJECTS       = writefile.o \
		treenode.o \
		translatesequences.o \
		terminalsequence.o \
		terminalnode.o \
		site.o \
		sequence.o \
		readnewick.o \
		readfile.o \
		readalignment.o \
		pwsite.o \
		pwhirschberg.o \
		progressivealignment.o \
		prank.o \
		postprobability.o \
		phylomatchscore.o \
		node.o \
		intmatrix.o \
		hmmodel.o \
		hirschberg.o \
		guidetree.o \
		fullprobability.o \
		flmatrix.o \
		eigen.o \
		dbmatrix.o \
		characterprobability.o \
		boolmatrix.o \
		ancestralsequence.o \
		ancestralnode.o \
		check_version.o \
		exonerate_reads.o \
		mafft_alignment.o \
		bppancestors.o \
		raxmlancestors.o \
		fasttree_tree.o \
		raxmlrebl.o
DESTDIR       = 
TARGET        = prank

MANPAGES      = prank.1

first: all
####### Implicit rules

.SUFFIXES: .o .c .cpp .cc .cxx .C

.cpp.o:
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o "$@" "$<"

.cc.o:
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o "$@" "$<"

.cxx.o:
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o "$@" "$<"

.C.o:
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o "$@" "$<"

.c.o:
	$(CC) -c $(CFLAGS) $(INCPATH) -o "$@" "$<"

####### Build rules


all: Makefile $(TARGET) $(MANPAGES)

$(TARGET):  $(OBJECTS)  
	$(LINK) $(LFLAGS) -o $(TARGET) $(OBJECTS) $(OBJCOMP) $(LIBS)


clean:compiler_clean 
	-$(DEL_FILE) $(OBJECTS)
	-$(DEL_FILE) *~ core *.core


####### Sub-libraries

distclean: clean
	-$(DEL_FILE) $(TARGET) 
	-$(DEL_FILE) Makefile


check: first

compiler_clean: 

####### Compile

writefile.o: writefile.cpp writefile.h \
		ancestralnode.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o writefile.o writefile.cpp

treenode.o: treenode.cpp config.h \
		hmmodel.h \
		dbmatrix.h \
		flmatrix.h \
		intmatrix.h \
		ancestralnode.h \
		writefile.h \
		treenode.h \
		sequence.h \
		site.h \
		boolmatrix.h \
		hirschberg.h \
		phylomatchscore.h \
		ancestralsequence.h \
		terminalsequence.h \
		fullprobability.h \
		postprobability.h \
		characterprobability.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o treenode.o treenode.cpp

translatesequences.o: translatesequences.cpp translatesequences.h \
		config.h \
		hmmodel.h \
		dbmatrix.h \
		flmatrix.h \
		intmatrix.h \
		ancestralnode.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o translatesequences.o translatesequences.cpp

terminalsequence.o: terminalsequence.cpp terminalsequence.h \
		sequence.h \
		site.h \
		intmatrix.h \
		flmatrix.h \
		dbmatrix.h \
		boolmatrix.h \
		config.h \
		hmmodel.h \
		ancestralnode.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o terminalsequence.o terminalsequence.cpp

terminalnode.o: terminalnode.cpp terminalnode.h \
		treenode.h \
		sequence.h \
		site.h \
		intmatrix.h \
		flmatrix.h \
		dbmatrix.h \
		boolmatrix.h \
		terminalsequence.h \
		config.h \
		hmmodel.h \
		ancestralnode.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o terminalnode.o terminalnode.cpp

site.o: site.cpp site.h \
		intmatrix.h \
		flmatrix.h \
		dbmatrix.h \
		boolmatrix.h \
		hmmodel.h \
		ancestralnode.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o site.o site.cpp

sequence.o: sequence.cpp sequence.h \
		site.h \
		intmatrix.h \
		flmatrix.h \
		dbmatrix.h \
		boolmatrix.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o sequence.o sequence.cpp

readnewick.o: readnewick.cpp readnewick.h \
		treenode.h \
		sequence.h \
		site.h \
		intmatrix.h \
		flmatrix.h \
		dbmatrix.h \
		boolmatrix.h \
		ancestralnode.h \
		node.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o readnewick.o readnewick.cpp

readfile.o: readfile.cpp readfile.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o readfile.o readfile.cpp

readalignment.o: readalignment.cpp readalignment.h \
		sequence.h \
		site.h \
		intmatrix.h \
		flmatrix.h \
		dbmatrix.h \
		boolmatrix.h \
		treenode.h \
		phylomatchscore.h \
		ancestralsequence.h \
		terminalsequence.h \
		config.h \
		hmmodel.h \
		ancestralnode.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o readalignment.o readalignment.cpp

pwsite.o: pwsite.cpp pwsite.h \
		flmatrix.h \
		intmatrix.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o pwsite.o pwsite.cpp

pwhirschberg.o: pwhirschberg.cpp config.h \
		hmmodel.h \
		dbmatrix.h \
		flmatrix.h \
		intmatrix.h \
		ancestralnode.h \
		pwhirschberg.h \
		pwsite.h \
		exonerate_reads.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o pwhirschberg.o pwhirschberg.cpp

progressivealignment.o: progressivealignment.cpp readnewick.h \
		treenode.h \
		sequence.h \
		site.h \
		intmatrix.h \
		flmatrix.h \
		dbmatrix.h \
		boolmatrix.h \
		readfile.h \
		writefile.h \
		ancestralnode.h \
		guidetree.h \
		progressivealignment.h \
		config.h \
		hmmodel.h \
		translatesequences.h \
		node.h \
		mafft_alignment.h \
		exonerate_reads.h \
		bppancestors.h \
		raxmlancestors.h \
		fasttree_tree.h \
		readalignment.h \
		phylomatchscore.h \
		ancestralsequence.h \
		terminalsequence.h \
		hirschberg.h \
		raxmlrebl.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o progressivealignment.o progressivealignment.cpp

prank.o: prank.cpp progressivealignment.h \
		config.h \
		hmmodel.h \
		dbmatrix.h \
		flmatrix.h \
		intmatrix.h \
		ancestralnode.h \
		treenode.h \
		sequence.h \
		site.h \
		boolmatrix.h \
		readfile.h \
		writefile.h \
		readnewick.h \
		guidetree.h \
		translatesequences.h \
		node.h \
		mafft_alignment.h \
		exonerate_reads.h \
		bppancestors.h \
		raxmlancestors.h \
		fasttree_tree.h \
		readalignment.h \
		phylomatchscore.h \
		ancestralsequence.h \
		terminalsequence.h \
		hirschberg.h \
		check_version.h \
		prank.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o prank.o prank.cpp

postprobability.o: postprobability.cpp config.h \
		hmmodel.h \
		dbmatrix.h \
		flmatrix.h \
		intmatrix.h \
		ancestralnode.h \
		postprobability.h \
		sequence.h \
		site.h \
		boolmatrix.h \
		phylomatchscore.h \
		ancestralsequence.h \
		terminalsequence.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o postprobability.o postprobability.cpp

phylomatchscore.o: phylomatchscore.cpp config.h \
		hmmodel.h \
		dbmatrix.h \
		flmatrix.h \
		intmatrix.h \
		ancestralnode.h \
		phylomatchscore.h \
		sequence.h \
		site.h \
		boolmatrix.h \
		ancestralsequence.h \
		terminalsequence.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o phylomatchscore.o phylomatchscore.cpp

node.o: node.cpp node.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o node.o node.cpp

intmatrix.o: intmatrix.cpp intmatrix.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o intmatrix.o intmatrix.cpp

hmmodel.o: hmmodel.cpp hmmodel.h \
		dbmatrix.h \
		flmatrix.h \
		intmatrix.h \
		ancestralnode.h \
		eigen.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o hmmodel.o hmmodel.cpp

hirschberg.o: hirschberg.cpp config.h \
		hmmodel.h \
		dbmatrix.h \
		flmatrix.h \
		intmatrix.h \
		ancestralnode.h \
		exonerate_reads.h \
		hirschberg.h \
		sequence.h \
		site.h \
		boolmatrix.h \
		phylomatchscore.h \
		ancestralsequence.h \
		terminalsequence.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o hirschberg.o hirschberg.cpp

guidetree.o: guidetree.cpp guidetree.h \
		flmatrix.h \
		intmatrix.h \
		fasttree_tree.h \
		config.h \
		hmmodel.h \
		dbmatrix.h \
		ancestralnode.h \
		pwhirschberg.h \
		pwsite.h \
		translatesequences.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o guidetree.o guidetree.cpp

fullprobability.o: fullprobability.cpp config.h \
		hmmodel.h \
		dbmatrix.h \
		flmatrix.h \
		intmatrix.h \
		ancestralnode.h \
		fullprobability.h \
		sequence.h \
		site.h \
		boolmatrix.h \
		phylomatchscore.h \
		ancestralsequence.h \
		terminalsequence.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o fullprobability.o fullprobability.cpp

flmatrix.o: flmatrix.cpp flmatrix.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o flmatrix.o flmatrix.cpp

eigen.o: eigen.cpp eigen.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o eigen.o eigen.cpp

dbmatrix.o: dbmatrix.cpp dbmatrix.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o dbmatrix.o dbmatrix.cpp

characterprobability.o: characterprobability.cpp config.h \
		hmmodel.h \
		dbmatrix.h \
		flmatrix.h \
		intmatrix.h \
		ancestralnode.h \
		characterprobability.h \
		sequence.h \
		site.h \
		boolmatrix.h \
		ancestralsequence.h \
		terminalsequence.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o characterprobability.o characterprobability.cpp

boolmatrix.o: boolmatrix.cpp boolmatrix.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o boolmatrix.o boolmatrix.cpp

ancestralsequence.o: ancestralsequence.cpp ancestralsequence.h \
		sequence.h \
		site.h \
		intmatrix.h \
		flmatrix.h \
		dbmatrix.h \
		boolmatrix.h \
		config.h \
		hmmodel.h \
		ancestralnode.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o ancestralsequence.o ancestralsequence.cpp

ancestralnode.o: ancestralnode.cpp config.h \
		hmmodel.h \
		dbmatrix.h \
		flmatrix.h \
		intmatrix.h \
		ancestralnode.h \
		writefile.h \
		treenode.h \
		sequence.h \
		site.h \
		boolmatrix.h \
		hirschberg.h \
		phylomatchscore.h \
		ancestralsequence.h \
		terminalsequence.h \
		fullprobability.h \
		postprobability.h \
		characterprobability.h \
		terminalnode.h \
		readalignment.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o ancestralnode.o ancestralnode.cpp

check_version.o: check_version.cpp check_version.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o check_version.o check_version.cpp

exonerate_reads.o: exonerate_reads.cpp exonerate_reads.h \
		config.h \
		hmmodel.h \
		dbmatrix.h \
		flmatrix.h \
		intmatrix.h \
		ancestralnode.h \
		translatesequences.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o exonerate_reads.o exonerate_reads.cpp

mafft_alignment.o: mafft_alignment.cpp mafft_alignment.h \
		config.h \
		hmmodel.h \
		dbmatrix.h \
		flmatrix.h \
		intmatrix.h \
		ancestralnode.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o mafft_alignment.o mafft_alignment.cpp

bppancestors.o: bppancestors.cpp bppancestors.h \
		ancestralnode.h \
		config.h \
		hmmodel.h \
		dbmatrix.h \
		flmatrix.h \
		intmatrix.h \
		readfile.h \
		readnewick.h \
		treenode.h \
		sequence.h \
		site.h \
		boolmatrix.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o bppancestors.o bppancestors.cpp

raxmlancestors.o: raxmlancestors.cpp raxmlancestors.h \
		ancestralnode.h \
		config.h \
		hmmodel.h \
		dbmatrix.h \
		flmatrix.h \
		intmatrix.h \
		readfile.h \
		readnewick.h \
		treenode.h \
		sequence.h \
		site.h \
		boolmatrix.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o raxmlancestors.o raxmlancestors.cpp

fasttree_tree.o: fasttree_tree.cpp fasttree_tree.h \
		config.h \
		hmmodel.h \
		dbmatrix.h \
		flmatrix.h \
		intmatrix.h \
		ancestralnode.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o fasttree_tree.o fasttree_tree.cpp

raxmlrebl.o: raxmlrebl.cpp raxmlrebl.h \
		ancestralnode.h \
		config.h \
		hmmodel.h \
		dbmatrix.h \
		flmatrix.h \
		intmatrix.h \
		node.h \
		readfile.h \
		readnewick.h \
		treenode.h \
		sequence.h \
		site.h \
		boolmatrix.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o raxmlrebl.o raxmlrebl.cpp

####### Manpages

%.1: %.1.pod
	$(POD2MAN) --section=1 --name="$(shell basename $@ .1 | tr a-z A-Z)" $< >$@

doc: $(MANPAGES)


####### Install

install:   FORCE

uninstall:   FORCE

FORCE:

