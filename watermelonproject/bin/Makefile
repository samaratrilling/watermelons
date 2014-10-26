JFLAGS = -g
JC = /cygdrive/c/Program\ Files/Java/jdk1.8.0_20/bin/javac.exe
JVM = /cygdrive/c/Program\ Files/Java/jdk1.8.0_20/jre/bin/java.exe

.SUFFIXES: .java .class

all:
	$(JC) $(JFLAGS) watermelon/sim/*.java
#	$(JC) $(JFLAGS) watermelon/dumb/*.java
	$(JC) $(JFLAGS) watermelon/group6/*.java

run:
#	$(JVM) watermelon.sim.watermelon map0 40 50 true dumb 0.01
	$(JVM) watermelon.sim.watermelon map0 40 50 true group6 0.01

.PHONY: clean
clean:
	$(RM) watermelon/sim/*.class	\
	      watermelon/dumb/*.class	\
	      watermelon/group6/*.class
