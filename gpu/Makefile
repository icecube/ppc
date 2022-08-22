ifdef arch
archf	=	-arch=sm_$(arch)
else
archf	=
endif

nvcc	=	nvcc ppc.cu -Xptxas=-v $(archf) --maxrregcount=64 \
		-O2 --use_fast_math --compiler-options=-O2,--fast-math

gcpp	=	$(CXX) -x c++ ppc.cu -O2 --fast-math

mlib	=	-fPIC -DXLIB -c -o ppc.o && $(CC) -shared \
		-fPIC -Wl,-soname,xppc ppc.o -o libxppc.so

warn	=	2>&1 | grep -v assuming

all:
	@echo "	make go:   compile the ppc object code for GPU"
	@echo "	make co:   compile the ppc object code for CPU"
	@echo "	make gpu:  compile the ppc  executable for GPU"
	@echo "	make cpu:  compile the ppc  executable for CPU"
	@echo "	make glib: compile the libxppc library for GPU"
	@echo "	make clib: compile the libxppc library for CPU"

go:
	$(nvcc) -o ppc.o -DXLIB -c $(warn)

gpu:
	$(nvcc) -o ppc $(warn)

glib:
	$(nvcc),$(mlib)

link:
	ln -s ppc.cu ppc.cxx || true
	ln -s pro.cu pro.cxx || true

co:
	$(gcpp) -o ppc.o -DXLIB -c

cpu:
	$(gcpp) -o ppc

clib:
	$(gcpp) $(mlib)

clib_osx:	link
	$(gcpp) -fPIC -DXLIB -dynamic -dynamiclib -o libxppc.dylib

clean:
	rm ppc.o ppc libxppc.so ppc.cxx pro.cxx || true
