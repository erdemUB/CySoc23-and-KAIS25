include makefile.in

TARGET     = nucleus

SRCS	= $(CFILES:%=%.c) $(CXXFILES:%=%.cpp)

nucleus: main.o ktruss.o hierarchy.o util.o bucket.o graph.o
	$(LD) $(LDFLAGS) -o ../$@ main.o ktruss.o hierarchy.o util.o bucket.o graph.o  $(LIBS)
	/bin/rm *.o
	
clean: 
	/bin/rm -f ../$(TARGET) *.o 
