gcpp	=	nvc++ ppc.cxx -O2

mlib	=	-fPIC -DXLIB -c -o ppc.o && $(CC) -shared \
		-fPIC -Wl,-soname,xppc ppc.o -o libxppc.so

all:
	@echo "	make go:   compile the ppc object code for GPU"
	@echo "	make co:   compile the ppc object code for CPU"
	@echo "	make gpu:  compile the ppc  executable for GPU"
	@echo "	make cpu:  compile the ppc  executable for CPU"
	@echo "	make glib: compile the libxppc library for GPU"
	@echo "	make clib: compile the libxppc library for CPU"

go:
	$(gcpp) -o ppc.o -DXLIB -c -stdpar=gpu

gpu:
	$(gcpp) -o ppc -stdpar=gpu

glib:
	$(gcpp) -stdpar=gpu $(mlib)

co:
	$(gcpp) -o ppc.o -DXLIB -c -stdpar=multicore -DXCPU

cpu:
	$(gcpp) -o ppc -stdpar=multicore -DXCPU

clib:
	$(gcpp) -stdpar=multicore -DXCPU $(mlib)

clean:
	rm ppc.o ppc libxppc.so || true
