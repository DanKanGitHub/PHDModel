include ./Makefile.export.Trilinos

CC=g++
CFLAGS= -Wall -g -I. ${Trilinos_INCLUDE_DIRS} ${Trilinos_TPL_INCLUDE_DIRS} ${Trilinos_MPI_INCLUDE_DIRS} ${Trilinos_MPI_INCLUDE_DIRS}
LINKS = ${Trilinos_LIBRARY_DIRS} ${Trilinos_TPL_LIBRARY_DIRS} ${Trilinos_LIBRARIES} ${Trilinos_TPL_LIBRARIES} ${Trilinos_MPI_LIBRARIES}

fem2d: main.o VelBioMesh2d.o PreMesh2d.o SparseAssembly.o QuadPtWt.o InitialVel.o InitialStr.o ElementNeigh.o Conformation.o \
	DepartureFoot.o Shape2d.o FeetSearch.o functions2.o NatBoundary2d.o PreEssenBoundary2d.o \
	VelEssenBoundary2d.o InitMatZero.o InitialBio.o InitMatZeroBio.o BioSparseAssembly.o BioDiffSparseAssembly.o ProcNodePartitionClass.o BioWeightFunc.o \
	RetardationDividedByRelaxation.o PolymericViscosityFunctionofBioVolFrac.o WriteVelPreData.o WriteBioData.o WriteStrData.o VectorNorm.o InitialGuess.o
	${CC} -o $@ $^  $(LINKS)

main.o: main.cpp
	${CC} main.cpp -c ${CFLAGS}

VelBioMesh2d.o:  VelBioMesh2d.cpp VelBioMesh2d.h
	${CC} VelBioMesh2d.cpp -c ${CFLAGS}

PreMesh2d.o: PreMesh2d.cpp PreMesh2d.h
	${CC} PreMesh2d.cpp -c ${CFLAGS}

SparseAssembly.o: SparseAssembly.cpp SparseAssembly.h Shape2d.h DepartureFoot.h NatBoundary2d.h PreEssenBoundary2d.h VelEssenBoundary2d.h BioWeightFunc.h RetardationDividedByRelaxation.h
	${CC} SparseAssembly.cpp -c ${CFLAGS}

BioSparseAssembly.o: BioSparseAssembly.cpp BioSparseAssembly.h Shape2d.h
	${CC} BioSparseAssembly.cpp -c ${CFLAGS}

BioDiffSparseAssembly.o: BioDiffSparseAssembly.cpp BioDiffSparseAssembly.h Shape2d.h
	${CC} BioDiffSparseAssembly.cpp -c ${CFLAGS}

DepartureFoot.o: DepartureFoot.cpp DepartureFoot.h Shape2d.h FeetSearch.h
	${CC} DepartureFoot.cpp -c ${CFLAGS}
	
QuadPtWt.o: QuadPtWt.cpp QuadPtWt.h
	${CC} QuadPtWt.cpp -c ${CFLAGS} 

InitialVel.o: InitialVel.cpp InitialVel.h
	${CC} InitialVel.cpp -c ${CFLAGS}

InitialStr.o: InitialStr.cpp InitialStr.h BioWeightFunc.h RetardationDividedByRelaxation.h
	${CC} InitialStr.cpp -c ${CFLAGS}

ElementNeigh.o: ElementNeigh.cpp ElementNeigh.h
	${CC} ElementNeigh.cpp -c ${CFLAGS}

Conformation.o: Conformation.cpp Conformation.h Shape2d.h DepartureFoot.h
	${CC} Conformation.cpp -c ${CFLAGS}

Shape2d.o: Shape2d.cpp Shape2d.h
	${CC} Shape2d.cpp -c ${CFLAGS}

FeetSearch.o: FeetSearch.cpp FeetSearch.h
	${CC} FeetSearch.cpp -c ${CFLAGS}

InitMatZero.o: InitMatZero.cpp InitMatZero.h
	${CC} InitMatZero.cpp -c ${CFLAGS}

NatBoundary2d.o: NatBoundary2d.cpp NatBoundary2d.h functions2.h
	${CC} NatBoundary2d.cpp -c ${CFLAGS}

VelEssenBoundary2d.o: VelEssenBoundary2d.cpp VelEssenBoundary2d.h functions2.h
	${CC} VelEssenBoundary2d.cpp -c ${CFLAGS}

PreEssenBoundary2d.o: PreEssenBoundary2d.cpp PreEssenBoundary2d.h functions2.h
	${CC} PreEssenBoundary2d.cpp -c ${CFLAGS}

functions2.o: functions2.cpp functions2.h
	${CC} functions2.cpp -c ${CFLAGS}

InitialBio.o: InitialBio.cpp InitialBio.h
	${CC} InitialBio.cpp -c ${CFLAGS}

InitMatZeroBio.o: InitMatZeroBio.cpp InitMatZeroBio.h
	${CC} InitMatZeroBio.cpp -c ${CFLAGS}

ProcNodePartitionClass.o: ProcNodePartitionClass.cpp ProcNodePartitionClass.h
	${CC} ProcNodePartitionClass.cpp -c ${CFLAGS}

BioWeightFunc.o: BioWeightFunc.cpp BioWeightFunc.h
	${CC} BioWeightFunc.cpp -c ${CFLAGS}

RetardationDividedByRelaxation.o: RetardationDividedByRelaxation.cpp RetardationDividedByRelaxation.h
	${CC} RetardationDividedByRelaxation.cpp -c ${CFLAGS}

PolymericViscosityFunctionofBioVolFrac.o: PolymericViscosityFunctionofBioVolFrac.cpp PolymericViscosityFunctionofBioVolFrac.h
	${CC} PolymericViscosityFunctionofBioVolFrac.cpp -c ${CFLAGS}

WriteBioData.o: WriteBioData.cpp WriteBioData.h
	${CC} WriteBioData.cpp -c ${CFLAGS}

WriteVelPreData.o: WriteVelPreData.cpp WriteVelPreData.h
	${CC} WriteVelPreData.cpp -c ${CFLAGS}

WriteStrData.o: WriteStrData.cpp WriteStrData.h
	${CC} WriteStrData.cpp -c ${CFLAGS}

VectorNorm.o: VectorNorm.cpp VectorNorm.h
	${CC} VectorNorm.cpp -c ${CFLAGS}

InitialGuess.o: InitialGuess.cpp InitialGuess.h
	${CC} InitialGuess.cpp -c ${CFLAGS}

clean:
	rm -rf *.o
	rm -rf fem2d