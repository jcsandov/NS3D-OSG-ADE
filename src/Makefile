FC      = mpif90
FFLAGS  = -c -O3 -g
LDFLAGS = -O3
SRCS   := precision.F90      \
          global.F90         \
		  global_app.F90     \
		  global_mpi.F90     \
		  global_param.F90   \
		  global_osg.F90     \
                  checksum.F90       \
		  wf_mpi.F90         \
		  main.F90           \
		  init.F90           \
		  mg_driver.F90      \
		  mg_nlevel.F90      \
		  mg_inject.F90      \
		  mg_metrics.F90     \
		  solver_daf.F90     \
		  bcond_fm.F90       \
		  mg_brhs.F90        \
		  nonreflect_ibc.F90 \
		  des_eddy.F90       \
		  rsmooth.F90        \
		  sgtsv.f            \
		  xerbla.f           \
		# -llapack -lblas

OBJS   := $(SRCS:.F90=.o)

channel: $(OBJS)
	$(FC) $(LDFLAGS) -o $@ $(OBJS) 

%.o : %.F90
	$(FC) $(FFLAGS) -o $@ $<

clean:
	@- /bin/rm -fv *.o *.mod channel

.PHONY: clean


###### 
 
# 
# mpif90 -O3 -g precision.F90 
# 
# mpif90 -O3 -g global.F90  
# 
# mpif90 -O3 -g global_app.F90  
# 
# mpif90 -O3 -g global_mpi.F90  
# 
# mpif90 -O3 -g global_param.F90  
# 
# mpif90 -O3 -g global_osg.F90
# 
# mpif90 -O3 -g global_turb.F90  
# 
# mpif90 -O3 -g checksum.F90  
# 
# mpif90 -O3 -g wf_mpi.F90  
# 
# mpif90 -O3 -g main.F90 
# 
# mpif90 -O3 -g init.F90 
# 
# mpif90 -O3 -g mg_driver.F90 
# 
# mpif90 -O3 -g mg_nlevel.F90 
# 
# mpif90 -O3 -g mg_inject.F90  
# 
# mpif90 -O3 -g mg_metrics.F90  
# 
# mpif90 -O3 -g solver_daf.F90 
# 
# mpif90 -O3 -g bcond_fm.F90 
# 
# mpif90 -O3 -g mg_brhs.F90 
# 
# mpif90 -O3 -g nonreflect_ibc.F90 
# 
# mpif90 -O3 -g des_eddy.F90 
# 		
# mpif90 -O3 -g rsmooth.F90 
# 
# mpif90 -O3 -g sgtsv.f
# 
# mpif90 -O3 -g xerbla.f
# 
# mpif90 -O3 -g -o channel precision.o global.o global_param.o global_app.o \
#                 global_mpi.o global_osg.o checksum.o wf_mpi.o  global_turb.o \
# 	main.o \
# 	init.o \
# 	mg_nlevel.o \
# 	mg_driver.o \
# 	mg_inject.o \
# 	mg_metrics.o \
# 	mg_brhs.o \
# 	nonreflect_ibc.o \
# 	solver_daf.o \
# 	bcond_fm.o \
# 	des_eddy.o \
#         sgtsv.o \
#         xerbla.o \
# 	rsmooth.o
# 
# 
