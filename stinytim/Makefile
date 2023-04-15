LFLAGS_NORMAL = "LFLAGS = -lm"
CFLAGSSPARC = "CFLAGS = -O"
CFLAGSHP = "CFLAGS = -Ae -O"
CFLAGSLINUX = "CFLAGS = -O3 -Wcast-align -march=i686 -ffloat-store"
CFLAGSOTHER = "CFLAGS =  -O"
CFLAGSTHREADEDSOLARIS = "CFLAGS = -mt -fsingle -O -DTT_THREADED -D_REENTRANT -D_POSIX_C_SOURCE=199506L"
LFLAGS_THREADEDSOLARIS = "LFLAGS = -lm -lthread -lpthread" 
CFLAGSTHREADEDLINUX = "CFLAGS = -O3 -D_REENTRANT -DTT_THREADED -ffloat-store"
LFLAGS_THREADEDLINUX = "LFLAGS = -lm -pthread" 
 
nothing: 
	@echo
	@echo "To compile Tiny Tim/SIRTF, use one of the following :"
	@echo
	@echo "          make sparcstation"
	@echo "          make linux"
	@echo "          make hp"
	@echo
	@echo "If you are compiling on something else, then"
	@echo "you must edit the Makefile according to the"
	@echo "manual and enter"
	@echo
	@echo "          make other"
	@echo

sparcstation:
	@make -f make.st1 $(CFLAGSSPARC) $(LFLAGS_NORMAL)
	@make -f make.st2 $(CFLAGSSPARC) $(LFLAGS_NORMAL)
	@echo

linux: 
	@make -f make.st1 $(CFLAGSLINUX) $(LFLAGS_NORMAL)
	@make -f make.st2 $(CFLAGSLINUX) $(LFLAGS_NORMAL)
	@echo

hp: 
	@make -f make.st1 $(CFLAGSHP) $(LFLAGS_NORMAL)
	@make -f make.st2 $(CFLAGSHP) $(LFLAGS_NORMAL)
	@echo

other:
	@make -f make.st1 $(CFLAGSOTHER) $(LFLAGS_NORMAL)
	@make -f make.st2 $(CFLAGSOTHER) $(LFLAGS_NORMAL)
	@echo

threadedsolaris: 
	@make -f make.st1 $(CFLAGSTHREADEDSOLARIS) $(LFLAGS_THREADEDSOLARIS)
	@make -f make.st2 $(CFLAGSTHREADEDSOLARIS) $(LFLAGS_THREADEDSOLARIS)
	@echo

threadedlinux: 
	@make -f make.st1 $(CFLAGSTHREADEDLINUX) $(LFLAGS_THREADEDLINUX)
	@make -f make.st2 $(CFLAGSTHREADEDLINUX) $(LFLAGS_THREADEDLINUX)
	@echo

clean:
	rm -f stiny1 stiny2 *.o

