include Make.inc

TESTS  = testChooseTransition testParticleFilter testTrajParticle
TESTS += testReadTrajectory

#SRCS  = $(wildcard *.c)
#SRCS += $(wildcard *.cpp)

MODELS_SRC = $(wildcard models/*.cpp)
MODELS = $(MODELS_SRC:.cpp=.o)

OBJS  = MCTraj.o Trajectory.o ParticleFilter.o TrajParticle.o 
OBJS += StateTransition.o TransitionType.o
OBJS += TrajParticleFilter.o HistoryFilter.o StaticFilter.o
OBJS += pfLik.o Model.o EpiState.o BranchState.o Tree.o TreeNode.o
OBJS += $(MODELS)

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

cpso/libpso.a: cpso/*.cpp
	cd cpso; make libpso.a

cpso/libpso_mpi.a: cpso/*.cpp
	cd cpso; make libpso_mpi.a

##############################################################################

models/SIS.o: models/SIS.h models/SIS.cpp
	$(CPP) $(CPPFLAGS) -c -o models/SIS.o models/SIS.cpp

models/SIR.o: models/SIR.h models/SIR.cpp
	$(CPP) $(CPPFLAGS) -c -o models/SIR.o models/SIR.cpp

models/SEIS.o: models/SEIS.h models/SEIS.cpp
	$(CPP) $(CPPFLAGS) -c -o models/SEIS.o models/SEIS.cpp

##############################################################################

testChooseTransition: testChooseTransition.cpp MCTraj.o ParticleFilter.o mctraj.a
	$(CPP) $(CPPFLAGS) -o $@ $^ -lgsl 

testParticleFilter: testParticleFilter.cpp ParticleFilter.o
	$(CPP) $(CPPFLAGS) -o testParticleFilter testParticleFilter.cpp mctraj.a -lgsl 

testTrajParticle: testTrajParticle.cpp mctraj.a Tree.o
	$(CPP) $(CPPFLAGS) -o testTrajParticle testTrajParticle.cpp mctraj.a \
		Tree.o -lgsl $(LDFLAGS)

testReadTrajectory: testReadTrajectory.cpp
	$(CPP) $(CPPFLAGS) -o testReadTrajectory testReadTrajectory.cpp mctraj.a -lgsl 

calcPfLik: calcPfLik.cpp mctraj.a Tree.o
	$(CPP) $(CPPFLAGS) -o $@ $^ -lgsl $(LDFLAGS)

simulateTrajectory: simulateTrajectory.cpp mctraj.a Tree.o
	$(CPP) $(CPPFLAGS) -o $@ $^ -lgsl $(LDFLAGS)

metrop_pf: metrop_pf.cpp mctraj.a Tree.o
	$(CPP) $(CPPFLAGS) -o $@ $^ -lgsl $(LDFLAGS)

pfSPSA: pfSPSA.cpp mctraj.a Tree.o
	$(CPP) $(CPPFLAGS) -o $@ $^ -lgsl $(LDFLAGS)

pfDREAM: pfDREAM.cpp mctraj.a cdream/dream.a Tree.o
	$(CPP) $(CPPFLAGS) -o $@ $^ -lgsl $(LDFLAGS)

pfMH: pfMH.cpp mctraj.a cdream/dream.a Tree.o
	$(CPP) $(CPPFLAGS) -o $@ $^ -lgsl $(LDFLAGS)

pfPSO: pfPSO.cpp mctraj.a cpso/libpso.a Tree.o
	$(CPP) $(CPPFLAGS) -o $@ $^ -lgsl $(LDFLAGS)

pfPSO_mpi: pfPSO.cpp mctraj.a cpso/libpso_mpi.a Tree.o
	$(MPICPP) $(MPICFLAGS) -o $@ $^ -lgsl $(MPILDFLAGS)

##############################################################################

testSingleTraj: testSingleTraj.cpp mctraj.a $(EXPOMV) ../Tree.o ../expoTreeSIR.o ../expoTree.o
	$(CPP) $(CPPFLAGS) -o testSingleTraj testSingleTraj.cpp mctraj.a \
		../expoTreeSIR.o ../expoTree.o $(EXPOMV) ../Tree.o -lgsl $(LDFLAGS)

clockPfLik: clockPfLik.cpp mctraj.a $(EXPOMV) ../Tree.o ../expoTreeSIR.o ../expoTree.o
	$(CPP) $(CPPFLAGS) -o clockPfLik clockPfLik.cpp mctraj.a \
		../expoTreeSIR.o ../expoTree.o $(EXPOMV) ../Tree.o -lgsl $(LDFLAGS)

##############################################################################

testRng: testRng.cpp
	$(CPP) $(CPPFLAGS) -o $@ $< -lgsl $(LDFLAGS)

testBranchState: testBranchState.cpp BranchState.o
	$(CPP) $(CPPFLAGS) -o $@ $< BranchState.o -lgsl $(LDFLAGS)

##############################################################################

.c.o: $<
	$(CC) $(CFLAGS) -c $<

.f.o: $<
	$(FC) $(FFLAGS) -c $<

.cpp.o: $<
	$(CPP) $(CPPFLAGS) -c $<

