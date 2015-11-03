include Make.inc

TESTS  = testChooseTransition testParticleFilter testTrajParticle
TESTS += testReadTrajectory

MODELS_SRC = $(wildcard models/*.cpp)
MODELS = $(MODELS_SRC:.cpp=.o)

OBJS  = MCTraj.o Trajectory.o ParticleFilter.o TrajParticle.o 
OBJS += StateTransition.o TransitionType.o
OBJS += TrajParticleFilter.o HistoryFilter.o StaticFilter.o
OBJS += pfLik.o Model.o EpiState.o BranchState.o Tree.o TreeNode.o
OBJS += Parameters.o JSON.o $(MODELS)

.PHONY: all debug

all: calcPfLik

DEBUG_FLAGS := -g -DDEBUG

clean: clean_build clean_debug

clean_build:
	rm -f $(OBJS:%.o=build/%.o)

clean_debug:
	rm -f $(OBJS:%.o=debug/%.o)

distclean: clean
	rm -f build/libmctraj.a
	rm -f $(TESTS)

tests: all

lib: build/libmctraj.a

build/libmctraj.a: $(OBJS:%.o=build/%.o)
	ar rcs $@ $(OBJS:%.o=build/%.o)

debug/libmctraj.a: $(OBJS:%.o=debug/%.o)
	ar rcs $@ $(OBJS:%.o=debug/%.o)

cdream/dream.a: cdream/*.cpp
	cd cdream; make dream.a

cpso/libpso.a: cpso/*.cpp
	cd cpso; make libpso.a

cpso/libpso_mpi.a: cpso/*.cpp
	cd cpso; make libpso_mpi.a

##############################################################################

.c.o: $<
	$(CC) $(CFLAGS) -c $<

.f.o: $<
	$(FC) $(FFLAGS) -c $<

.cpp.o: $<
	$(CPP) $(CPPFLAGS) -c $<

models/%.o: models/%.cpp models/%.h 
	$(CPP) $(CPPFLAGS) -c -o $@ $<

build/%.o: %.cpp
	$(CPP) $(CPPFLAGS) -o $@ -c $<

debug/%.o: %.cpp
	$(CPP) $(DEBUG_FLAGS) $(CPPFLAGS) -o $@ -c $<

build/%: %.cpp build/libmctraj.a
	$(CPP) $(CPPFLAGS) -o $@ $^ -lgsl $(LDFLAGS)

debug/%: %.cpp debug/libmctraj.a
	$(CPP) $(DEBUG_FLAGS) $(CPPFLAGS) -o $@ $^ -lgsl $(LDFLAGS)

##############################################################################

testChooseTransition: testChooseTransition.cpp MCTraj.o ParticleFilter.o libmctraj.a
	$(CPP) $(CPPFLAGS) -o $@ $^ -lgsl 

testParticleFilter: testParticleFilter.cpp ParticleFilter.o
	$(CPP) $(CPPFLAGS) -o testParticleFilter testParticleFilter.cpp libmctraj.a -lgsl 

testTrajParticle: testTrajParticle.cpp libmctraj.a Tree.o
	$(CPP) $(CPPFLAGS) -o testTrajParticle testTrajParticle.cpp libmctraj.a \
		Tree.o -lgsl $(LDFLAGS)

testReadTrajectory: testReadTrajectory.cpp
	$(CPP) $(CPPFLAGS) -o testReadTrajectory testReadTrajectory.cpp libmctraj.a -lgsl 

simulateTrajectory: simulateTrajectory.cpp libmctraj.a
	$(CPP) $(CPPFLAGS) -o $@ $^ -lgsl $(LDFLAGS)

metrop_pf: metrop_pf.cpp libmctraj.a Tree.o
	$(CPP) $(CPPFLAGS) -o $@ $^ -lgsl $(LDFLAGS)

build/pfSPSA: pfSPSA.cpp libmctraj.a Tree.o
	$(CPP) $(CPPFLAGS) -o $@ $^ -lgsl $(LDFLAGS)

build/pfMH: pfMH.cpp build/libmctraj.a cdream/dream.a Tree.o
	$(CPP) $(CPPFLAGS) -o $@ $^ -lgsl $(LDFLAGS)

# DREAM

build/pfDREAM: pfDREAM.cpp build/libmctraj.a cdream/dream.a Tree.o
	$(CPP) $(CPPFLAGS) -o $@ $^ -lgsl $(LDFLAGS)

debug/pfDREAM: pfDREAM.cpp debug/libmctraj.a cdream/dream.a Tree.o
	$(CPP) $(CPPFLAGS) -o $@ $^ -lgsl $(LDFLAGS)

# PSO

build/pfPSO: pfPSO.cpp build/libmctraj.a cpso/libpso.a Tree.o
	$(CPP) $(CPPFLAGS) -o $@ $^ -lgsl $(LDFLAGS)

debug/pfPSO: pfPSO.cpp debug/libmctraj.a cpso/libpso.a Tree.o
	$(CPP) $(CPPFLAGS) -o $@ $^ -lgsl $(LDFLAGS)

pfPSO_mpi: pfPSO.cpp libmctraj.a cpso/libpso_mpi.a Tree.o
	$(MPICPP) $(MPICFLAGS) -o $@ $^ -lgsl $(MPILDFLAGS)

##############################################################################

testSingleTraj: testSingleTraj.cpp libmctraj.a $(EXPOMV) ../Tree.o ../expoTreeSIR.o ../expoTree.o
	$(CPP) $(CPPFLAGS) -o testSingleTraj testSingleTraj.cpp libmctraj.a \
		../expoTreeSIR.o ../expoTree.o $(EXPOMV) ../Tree.o -lgsl $(LDFLAGS)

clockPfLik: clockPfLik.cpp libmctraj.a $(EXPOMV) ../Tree.o ../expoTreeSIR.o ../expoTree.o
	$(CPP) $(CPPFLAGS) -o clockPfLik clockPfLik.cpp libmctraj.a \
		../expoTreeSIR.o ../expoTree.o $(EXPOMV) ../Tree.o -lgsl $(LDFLAGS)

##############################################################################

testRng: testRng.cpp
	$(CPP) $(CPPFLAGS) -o $@ $< -lgsl $(LDFLAGS)

testBranchState: testBranchState.cpp BranchState.o
	$(CPP) $(CPPFLAGS) -o $@ $< BranchState.o -lgsl $(LDFLAGS)

##############################################################################


