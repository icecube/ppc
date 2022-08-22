ifdef opcl
# opts	=	-I$(opcl)/include -L$(opcl)/lib/x86_64 -lOpenCL
opts	=	-I$(opcl)/include -L$(opcl)/lib64 -lOpenCL
else
opts	=	-framework OpenCL
endif

gcpp	=	$(CXX) ppc.cxx -O2 --fast-math $(opts) -Wall -Wextra

mlib	=	-fPIC -DXLIB -c -o ppc.o && $(CC) -shared \
		-fPIC -Wl,-soname,xppc ppc.o -o libxppc.so

warn	=	2>&1 | grep -v assuming

all:
	@echo "	make obj: compile the object code"
	@echo "	make lib: compile the xppc library"
	@echo "	make ocl: compile the ppc executable"

obj:
	$(gcpp) -o ppc.o -DXLIB -c

ocl:
	$(gcpp) -o ppc

lib:
	$(gcpp) $(mlib)

clean:
	rm ppc.o ppc libxppc.so || true
