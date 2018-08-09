OBJS = functions.o Basis_1orb_SF.o Model_1orb_SF.o Spin_dynamics_VNE_1orbHubbard_engine.o Space_Time_FourierTransform.o main.o
#DEBUG = -g3
OPTFLAG = -O3
CPPFLAGS = -std=c++11
CC = g++ $(OPTFLAG) $(CPPFLAGS)
CFLAGS = -c $(DEBUG) #-DWITH_COMPLEX
LFLAGS = $(DEBUG)
MKL_LIB += -ldl -lpthread -lm
MKL_LIB2 = -llapack -lblas #/opt/intel/mkl/lib/libmkl_core.a  /opt/intel/mkl/lib/libmkl_intel_lp64.a  /opt/intel/mkl/lib/libmkl_sequential.a
#/opt/intel/mkl/lib/libmkl_intel_thread.a #/opt/intel/mkl/lib/libmkl_gnu_thread.a  #/opt/intel/mkl/lib/libmkl_sequential.a
MKL_include = #-I/opt/intel/mkl/include
OPENMP = #/opt/intel/compilers_and_libraries_2016.3.170/mac/compiler/lib/
LIBS_1 =  #-L$(OPENMP)
LIBS_1 += #-fopenmp #-liomp5
 

all :	$(OBJS) 
	$(CC) $(LFLAGS) $(OBJS) -o dynamics $(MKL_include) $(MKL_LIB2) $(MKL_LIB) $(LIBS_1)
	cp dynamics NI_Skw

functions.o : functions.cpp
	$(CC) $(CFLAGS) functions.cpp $(MKL_include) $(MKL_LIB) $(LIBS_1)

Basis_1orb_SF.o : src/Model_1orbHubbard/Basis_1orb_SF.cpp
	$(CC) $(CFLAGS) src/Model_1orbHubbard/Basis_1orb_SF.cpp $(MKL_include) $(MKL_LIB) $(LIBS_1)

Model_1orb_SF.o : tensor_type.h functions.h src/Model_1orbHubbard/Model_1orb_SF.h src/Model_1orbHubbard/Model_1orb_SF.cpp
	$(CC) $(CFLAGS) src/Model_1orbHubbard/Model_1orb_SF.cpp $(MKL_include) $(MKL_LIB) $(LIBS_1)

#Spin_dynamics_ED_engine.o : tensor_type.h functions.h src/Model_1orbHubbard/Basis_1orb_SF.h src/Model_1orbHubbard/Model_1orb_SF.h Spin_dynamics_ED_engine.h Spin_dynamics_ED_engine.cpp
#	$(CC) $(CFLAGS) Spin_dynamics_ED_engine.cpp $(MKL_include) $(MKL_LIB) $(LIBS_1)

Spin_dynamics_VNE_1orbHubbard_engine.o : src/Model_1orbHubbard/Spin_dynamics_VNE_1orbHubbard_engine.cpp
	$(CC) $(CFLAGS) src/Model_1orbHubbard/Spin_dynamics_VNE_1orbHubbard_engine.cpp $(MKL_include) $(MKL_LIB) $(LIBS_1)

#Spin_dynamics_VNE_3orbPnictides_engine.o : src/Model_3orbPnictides/Spin_dynamics_VNE_3orbPnictides_engine.cpp
#	$(CC) $(CFLAGS) src/Model_3orbPnictides/Spin_dynamics_VNE_3orbPnictides_engine.cpp $(MKL_include) $(MKL_LIB) $(LIBS_1)

Space_Time_FourierTransform.o : Space_Time_FourierTransform.cpp
	$(CC) $(CFLAGS) Space_Time_FourierTransform.cpp $(MKL_include) $(MKL_LIB) $(LIBS_1)

main.o : main.cpp
	$(CC) $(CFLAGS) main.cpp $(MKL_include) $(MKL_LIB) $(LIBS_1)

clean:
	rm *.o dynamics NI_Skw ST_Fourier

NI_Skw : dynamics
	cp dynamics NI_Skw

ST_Fourier : dynamics
	cp dynamics ST_Fourier
