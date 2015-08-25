include ../Make.inc

CPPFLAGS += -isystem ../include

../build/mctraj.a:
	cd ..
	make build/mctraj.a

build/%: %.cpp ../build/mctraj.a
	$(CPP) $(CPPFLAGS) -o $@ $^ -lgsl $(LDFLAGS)


