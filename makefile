OBJS = functions.o Basis_1orb_SF.o Model_1orb_SF.o Spin_dynamics_VNE_1orbHubbard_engine.o Coordinates.o ParametersEngine.o MFParams.o Hamiltonian.o Observables.o Spin_dynamics_VNE_3orbPnictides_engine.o Space_Time_FourierTransform.o Spin_dynamics_VNE_1orb_engine_MCMF.o Coordinates_MCMF.o ParametersEngine_MCMF.o Hamiltonian_MCMF.o Observables_MCMF.o MFParams_MCMF.o Spin_dynamics_VNE_MultiOrbSF.o Coordinates_MultiOrbSF.o ParametersEngine_MultiOrbSF.o Hamiltonian_MultiOrbSF.o Observables_MultiOrbSF.o MFParams_MultiOrbSF.o Spin_dynamics_Rotor.o ParametersEngine_Rotor.o ParametersEngine_GenericHamilSpinsO3.o Spin_dynamics_VNE_GenericHamilSpinsO3.o main.o
DEBUG = #-g3
OPTFLAG = -O3
CPPFLAGS = -std=c++11
CC = g++ $(OPTFLAG) $(CPPFLAGS)
CFLAGS = -c $(DEBUG) #-DWITH_COMPLEX
LFLAGS = $(DEBUG)
MKL_LIB += -llapack -lblas #-ldl -lpthread -lm
MKL_LIB2 = #-llapack -lblas #/opt/intel/mkl/lib/libmkl_core.a  /opt/intel/mkl/lib/libmkl_intel_lp64.a  /opt/intel/mkl/lib/libmkl_sequential.a
MKL_LIB2 = #/opt/intel/mkl/lib/libmkl_core.a  /opt/intel/mkl/lib/libmkl_intel_lp64.a  /opt/intel/mkl/lib/libmkl_sequential.a
#/opt/intel/mkl/lib/libmkl_intel_thread.a #/opt/intel/mkl/lib/libmkl_gnu_thread.a  #/opt/intel/mkl/lib/libmkl_sequential.a
MKL_include = #-I/opt/intel/mkl/include
#MKL_include = /opt/intel/mkl/lib/libmkl_core.a  /opt/intel/mkl/lib/libmkl_intel_lp64.a  /opt/intel/mkl/lib/libmkl_sequential.a
OPENMP = #/opt/intel/compilers_and_libraries_2016.3.170/mac/compiler/lib/
LIBS_1 =  #-L$(OPENMP)
LIBS_1 += #-fopenmp #-liomp5  #-lgomp #-liomp5
#LIBS_1 += -qopenmp -liomp5  #-lgomp #-liomp5

all :	$(OBJS) 
	$(CC) $(LIBS_1) $(LFLAGS) $(OBJS) -o dynamics $(MKL_include) $(MKL_LIB2) $(MKL_LIB)
	cp dynamics NI_Skw
	cp dynamics ST_Fourier
	cp dynamics Calculate_Fqw
	cp dynamics Calculate_2pFqw
	cp dynamics Calculate_Sqw
	cp dynamics convolution_w
	cp dynamics sum_rule

functions.o : functions.cpp
	$(CC) $(LIBS_1) $(CFLAGS) functions.cpp $(MKL_include) $(MKL_LIB) 

Basis_1orb_SF.o : src/Model_1orbHubbard_old/Basis_1orb_SF.cpp
	$(CC) $(LIBS_1) $(CFLAGS) src/Model_1orbHubbard_old/Basis_1orb_SF.cpp $(MKL_include) $(MKL_LIB)

Model_1orb_SF.o : tensor_type.h functions.h src/Model_1orbHubbard_old/Model_1orb_SF.h src/Model_1orbHubbard_old/Model_1orb_SF.cpp
	$(CC) $(LIBS_1) $(CFLAGS) src/Model_1orbHubbard_old/Model_1orb_SF.cpp $(MKL_include) $(MKL_LIB) 

#Spin_dynamics_ED_engine.o : tensor_type.h functions.h src/Model_1orbHubbard_old/Basis_1orb_SF.h src/Model_1orbHubbard_old/Model_1orb_SF.h Spin_dynamics_ED_engine.h Spin_dynamics_ED_engine.cpp
#	$(CC) $(CFLAGS) Spin_dynamics_ED_engine.cpp $(MKL_include) $(MKL_LIB) $(LIBS_1)

Spin_dynamics_VNE_1orbHubbard_engine.o : src/Model_1orbHubbard_old/Spin_dynamics_VNE_1orbHubbard_engine.cpp
	$(CC) $(LIBS_1) $(CFLAGS) src/Model_1orbHubbard_old/Spin_dynamics_VNE_1orbHubbard_engine.cpp $(MKL_include) $(MKL_LIB)

Space_Time_FourierTransform.o : src/Model_1orbHubbard_old/Space_Time_FourierTransform.cpp
	$(CC) $(LIBS_1) $(CFLAGS) src/Model_1orbHubbard_old/Space_Time_FourierTransform.cpp $(MKL_include) $(MKL_LIB)


Spin_dynamics_VNE_3orbPnictides_engine.o : src/Model_3orbPnictides/Spin_dynamics_VNE_3orbPnictides_engine.cpp
	$(CC) $(LIBS_1) $(CFLAGS) src/Model_3orbPnictides/Spin_dynamics_VNE_3orbPnictides_engine.cpp $(MKL_include) $(MKL_LIB)

Coordinates.o : src/Model_3orbPnictides/Coordinates.cpp
	$(CC) $(LIBS_1) $(CFLAGS) src/Model_3orbPnictides/Coordinates.cpp $(MKL_include) $(MKL_LIB)

ParametersEngine.o : src/Model_3orbPnictides/ParametersEngine.cpp
	$(CC) $(LIBS_1) $(CFLAGS) src/Model_3orbPnictides/ParametersEngine.cpp $(MKL_include) $(MKL_LIB)

Hamiltonian.o : src/Model_3orbPnictides/Hamiltonian.cpp
	$(CC) $(LIBS_1) $(CFLAGS) src/Model_3orbPnictides/Hamiltonian.cpp $(MKL_include) $(MKL_LIB)

Observables.o : src/Model_3orbPnictides/Observables.cpp
	$(CC) $(LIBS_1) $(CFLAGS) src/Model_3orbPnictides/Observables.cpp $(MKL_include) $(MKL_LIB)

MFParams.o : src/Model_3orbPnictides/MFParams.cpp
	$(CC) $(LIBS_1) $(CFLAGS) src/Model_3orbPnictides/MFParams.cpp $(MKL_include) $(MKL_LIB)


Spin_dynamics_VNE_1orb_engine_MCMF.o : src/MCMF_1orbHubbard/Spin_dynamics_VNE_1orb_engine_MCMF.cpp
	$(CC) $(LIBS_1) $(CFLAGS) src/MCMF_1orbHubbard/Spin_dynamics_VNE_1orb_engine_MCMF.cpp $(MKL_include) $(MKL_LIB)

Coordinates_MCMF.o : src/MCMF_1orbHubbard/Coordinates_MCMF.cpp
	$(CC) $(LIBS_1) $(CFLAGS) src/MCMF_1orbHubbard/Coordinates_MCMF.cpp $(MKL_include) $(MKL_LIB)

ParametersEngine_MCMF.o : src/MCMF_1orbHubbard/ParametersEngine_MCMF.cpp
	$(CC) $(LIBS_1) $(CFLAGS) src/MCMF_1orbHubbard/ParametersEngine_MCMF.cpp $(MKL_include) $(MKL_LIB)

Hamiltonian_MCMF.o : src/MCMF_1orbHubbard/Hamiltonian_MCMF.cpp
	$(CC) $(LIBS_1) $(CFLAGS) src/MCMF_1orbHubbard/Hamiltonian_MCMF.cpp $(MKL_include) $(MKL_LIB)

Observables_MCMF.o : src/MCMF_1orbHubbard/Observables_MCMF.cpp
	$(CC) $(LIBS_1) $(CFLAGS) src/MCMF_1orbHubbard/Observables_MCMF.cpp $(MKL_include) $(MKL_LIB)

MFParams_MCMF.o : src/MCMF_1orbHubbard/MFParams_MCMF.cpp
	$(CC) $(LIBS_1) $(CFLAGS) src/MCMF_1orbHubbard/MFParams_MCMF.cpp $(MKL_include) $(MKL_LIB)


Spin_dynamics_VNE_MultiOrbSF.o : src/MultiOrb_SpinFermion/Spin_dynamics_VNE_MultiOrbSF.cpp
	$(CC) $(LIBS_1) $(CFLAGS) src/MultiOrb_SpinFermion/Spin_dynamics_VNE_MultiOrbSF.cpp $(MKL_include) $(MKL_LIB)

Coordinates_MultiOrbSF.o : src/MultiOrb_SpinFermion/Coordinates_MultiOrbSF.cpp
	$(CC) $(LIBS_1) $(CFLAGS) src/MultiOrb_SpinFermion/Coordinates_MultiOrbSF.cpp $(MKL_include) $(MKL_LIB)

ParametersEngine_MultiOrbSF.o : src/MultiOrb_SpinFermion/ParametersEngine_MultiOrbSF.cpp
	$(CC) $(LIBS_1) $(CFLAGS) src/MultiOrb_SpinFermion/ParametersEngine_MultiOrbSF.cpp $(MKL_include) $(MKL_LIB)

Hamiltonian_MultiOrbSF.o : src/MultiOrb_SpinFermion/Hamiltonian_MultiOrbSF.cpp
	$(CC) $(LIBS_1) $(CFLAGS) src/MultiOrb_SpinFermion/Hamiltonian_MultiOrbSF.cpp $(MKL_include) $(MKL_LIB)

Observables_MultiOrbSF.o : src/MultiOrb_SpinFermion/Observables_MultiOrbSF.cpp
	$(CC) $(LIBS_1) $(CFLAGS) src/MultiOrb_SpinFermion/Observables_MultiOrbSF.cpp $(MKL_include) $(MKL_LIB)

MFParams_MultiOrbSF.o : src/MultiOrb_SpinFermion/MFParams_MultiOrbSF.cpp
	$(CC) $(LIBS_1) $(CFLAGS) src/MultiOrb_SpinFermion/MFParams_MultiOrbSF.cpp $(MKL_include) $(MKL_LIB)

Spin_dynamics_Rotor.o : src/RotorDynamics/Spin_dynamics_Rotor.cpp
	$(CC) $(LIBS_1) $(CFLAGS) src/RotorDynamics/Spin_dynamics_Rotor.cpp $(MKL_include) $(MKL_LIB)

ParametersEngine_Rotor.o : src/RotorDynamics/ParametersEngine_Rotor.cpp
	$(CC) $(LIBS_1) $(CFLAGS) src/RotorDynamics/ParametersEngine_Rotor.cpp $(MKL_include) $(MKL_LIB)

ParametersEngine_GenericHamilSpinsO3.o : src/GenericHamilSpinsO3/ParametersEngine_GenericHamilSpinsO3.cpp
	$(CC) $(LIBS_1) $(CFLAGS) src/GenericHamilSpinsO3/ParametersEngine_GenericHamilSpinsO3.cpp $(MKL_INCLUDE) $(MKL_LIB)

Spin_dynamics_VNE_GenericHamilSpinsO3.o : src/GenericHamilSpinsO3/Spin_dynamics_VNE_GenericHamilSpinsO3.cpp
	$(CC) $(LIBS_1) $(CFLAGS) src/GenericHamilSpinsO3/Spin_dynamics_VNE_GenericHamilSpinsO3.cpp $(MKL_INCLUDE) $(MKL_LIB)


main.o : main.cpp
	$(CC) $(LIBS_1) $(CFLAGS) main.cpp $(MKL_include) $(MKL_LIB)

clean:
	rm *.o dynamics NI_Skw ST_Fourier Calculate_Fqw Calculate_Sqw Calculate_2pFqw

sum_rule : dynamics
	cp dynamics sum_rule

NI_Skw : dynamics
	cp dynamics NI_Skw

ST_Fourier : dynamics
	cp dynamics ST_Fourier

Calculate_Fqw: dynamics
	cp dynamics Calculate_Fqw

Calculate_2pFqw: dynamics
	cp dynamics Calculate_2pFqw

Calculate_Sqw: dynamics
	cp dynamics Calculate_Sqw

convolution_w: dynamics
	cp dynamics convolution_w
