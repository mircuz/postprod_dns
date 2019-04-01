CCX = mpicxx
CFLAGS = /usr/local/fri/complex.o /usr/local/fri/fft.o /usr/local/fri/rbmat.o 


all: translation compilation

translation:
	/usr/local/fri/fri /usr/local/fri/cpl postprod.cpl

compilation: 
	$(CC) -O3 -fPIC -c initialization.c
	$(CC) -O3 -fPIC -c postprod_support.c
	$(CC) -O3 -fPIC -c postprod_oper.c
	$(CC) -O3 -fPIC -c postprod.c
	$(CC) -O3 -fPIC -c fft_support.c
	$(CC) -O3 -fPIC -c convol_trasp.c
	
	$(CCX) $(CFLAGS) -o postprod_exe postprod.o initialization.o postprod_oper.o postprod_support.o fft_support.o convol_trasp.o /home/mirco/Scrivania/fftmpi-1Oct18/src/libfft3dmpi.a -L/home/mirco/Scrivania/hdf5/lib /home/mirco/Scrivania/hdf5/lib/libhdf5_hl.a /home/mirco/Scrivania/hdf5/lib/libhdf5.a -lz -ldl -lm -Wl,-rpath -Wl,/home/mirco/Scrivania/hdf5/lib -I/usr/local/include -pthread -Wl,-rpath -Wl,/usr/local/lib -Wl,--enable-new-dtags -L/usr/local/lib -lmpi
	make remove_useless
		#--> Executable ready <--
		#--> run as mpiexec -n "#procs" postprod_exe <--

remove_useless:
	mv fft_support.c fft_support.cost
	rm *.c
	mv fft_support.cost fft_support.c
	rm *.o
	rm *.d

clean: 
	rm postprod_exe
	rm results/Processed_Field.*
