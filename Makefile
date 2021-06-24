
objs=main.o init_rho.o init_x.o save_timestep.o mk_matrix_A.o mk_matrix_B.o pressure.o timestep.o lhs.o
cc=gfortran
bin=shock
libs= -lm 
cflags= -g

all: $(objs)
	$(cc) $? -o $(bin) $(libs) $(cflags)

%.o: %.f90
	$(cc) -c -o $@ $< $(cflags)

clean:
	rm *.o  $(bin) 
