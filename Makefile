include Make.inc
#CPP = clang++

TESTS  = testChooseTransition testParticleFilter testTrajParticle
TESTS += testReadTrajectory

#SRCS  = $(wildcard *.c)
#SRCS += $(wildcard *.cpp)

MODELS_SRC = $(wildcard models/*.cpp)
MODELS = $(MODELS_SRC:.cpp=.o)

OBJS  = MCTraj.o Trajectory.o ParticleFilter.o TrajParticle.o 
OBJS += StateTransition.o TransitionType.o TrajParticleFilter.o
OBJS += pfLik.o Model.o EpiState.o BranchState.o TreeNode.o
OBJS += $(MODELS)

CFLAGS += -fopenmp -Imodels -I. -I.. -isystem ./include

clean:
	rm -rf $(OBJS)

distclean: clean
	rm -f mctraj.a
	rm -f $(TESTS)

tests: all

lib: mctraj.a

mctraj.a: $(OBJS)
	ar rcs mctraj.a $(OBJS)

cdream/dream.a: cdream/*.cpp
	cd cdream; make dream.a

##############################################################################

models/SIS.o: models/SIS.h models/SIS.cpp
	$(CPP) $(CFLAGS) -c -o models/SIS.o models/SIS.cpp

models/SIR.o: models/SIR.h models/SIR.cpp
	$(CPP) $(CFLAGS) -c -o models/SIR.o models/SIR.cpp

models/SEIS.o: models/SEIS.h models/SEIS.cpp
	$(CPP) $(CFLAGS) -c -o models/SEIS.o models/SEIS.cpp

##############################################################################

testChooseTransition: testChooseTransition.cpp MCTraj.o ParticleFilter.o
	$(CPP) $(CFLAGS) -o testChooseTransition testChooseTransition.cpp mctraj.a -lgsl 

testParticleFilter: testParticleFilter.cpp ParticleFilter.o
	$(CPP) $(CFLAGS) -o testParticleFilter testParticleFilter.cpp mctraj.a -lgsl 

testTrajParticle: testTrajParticle.cpp mctraj.a Tree.o
	$(CPP) $(CFLAGS) -o testTrajParticle testTrajParticle.cpp mctraj.a \
		Tree.o -lgsl $(LDFLAGS)

testReadTrajectory: testReadTrajectory.cpp
	$(CPP) $(CFLAGS) -o testReadTrajectory testReadTrajectory.cpp mctraj.a -lgsl 

calcPfLik: calcPfLik.cpp mctraj.a Tree.o
	$(CPP) $(CFLAGS) -o $@ $^ -lgsl $(LDFLAGS)

simulateTrajectory: simulateTrajectory.cpp mctraj.a Tree.o
	$(CPP) $(CFLAGS) -o $@ $^ -lgsl $(LDFLAGS)

metrop_pf: metrop_pf.cpp mctraj.a Tree.o
	$(CPP) $(CFLAGS) -o $@ $^ -lgsl $(LDFLAGS)

pfSPSA: pfSPSA.cpp mctraj.a Tree.o
	$(CPP) $(CFLAGS) -o $@ $^ -lgsl $(LDFLAGS)

pfDREAM: pfDREAM.cpp mctraj.a cdream/dream.a Tree.o
	$(CPP) $(CFLAGS) -o $@ $^ -lgsl $(LDFLAGS)

pfMH: pfMH.cpp mctraj.a cdream/dream.a Tree.o
	$(CPP) $(CFLAGS) -o $@ $^ -lgsl $(LDFLAGS)

##############################################################################

testSingleTraj: testSingleTraj.cpp mctraj.a $(EXPOMV) ../Tree.o ../expoTreeSIR.o ../expoTree.o
	$(CPP) $(CFLAGS) -o testSingleTraj testSingleTraj.cpp mctraj.a \
		../expoTreeSIR.o ../expoTree.o $(EXPOMV) ../Tree.o -lgsl $(LDFLAGS)

clockPfLik: clockPfLik.cpp mctraj.a $(EXPOMV) ../Tree.o ../expoTreeSIR.o ../expoTree.o
	$(CPP) $(CFLAGS) -o clockPfLik clockPfLik.cpp mctraj.a \
		../expoTreeSIR.o ../expoTree.o $(EXPOMV) ../Tree.o -lgsl $(LDFLAGS)

##############################################################################

testRng: testRng.cpp
	$(CPP) $(CFLAGS) -o $@ $< -lgsl $(LDFLAGS)

testBranchState: testBranchState.cpp BranchState.o
	$(CPP) $(CFLAGS) -o $@ $< BranchState.o -lgsl $(LDFLAGS)

##############################################################################

.c.o: $<
	$(CC) $(CFLAGS) -c $<

.f.o: $<
	$(FC) $(FFLAGS) -c $<

.cpp.o: $<
	$(CPP) $(CPPFLAGS) -c $<

