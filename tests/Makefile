include ../Make.inc

CPPFLAGS += -isystem ../include

../debug/libmctraj.a:
	cd ..
	make debug/libmctraj.a

build/%: %.cpp ../build/libmctraj.a
	cd ..; make build/libmctraj.a
	$(CPP) $(CPPFLAGS) -o $@ $^ -lgsl $(LDFLAGS)

debug/%: %.cpp
	cd ..; make debug/libmctraj.a
	$(CPP) -DDEBUG $(CPPFLAGS) -o $@ $^ -lgsl $(LDFLAGS)

