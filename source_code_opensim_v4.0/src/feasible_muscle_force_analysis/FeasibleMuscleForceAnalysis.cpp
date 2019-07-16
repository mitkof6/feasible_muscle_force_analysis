#include "FeasibleMuscleForceAnalysis.h"
#include <cassert>
#include <OpenSim/Simulation/Model/Model.h>
#include <OpenSim/Simulation/Model/Muscle.h>
extern "C" {
#include "lrslib.h"
};

using namespace std;
using namespace OpenSim;
using namespace SimTK;

/******************************************************************************/

Matrix calculateMomentArm(const State& s, const Model& model) {
    auto& coordinateSet = model.getCoordinateSet();
    auto& muscleSet = model.getMuscles();
    Matrix R(muscleSet.getSize(), coordinateSet.getSize());
    for (int i = 0; i < muscleSet.getSize(); i++) {
        for (int j = 0; j < coordinateSet.getSize(); j++) {
            R[i, j] = muscleSet[i].computeMomentArm(s, coordinateSet[j]);
        }
    }
    return R;
}

Vector calcMaxForce(const State& s, const Model& model) {
    vector<double> optimalForce;
    auto& muscleSet = model.getMuscles();
    Vector fMax(muscleSet.getSize(), 0.0);
    for (int i = 0; i < muscleSet.getSize(); i++) {
        fMax[i] = muscleSet[i].getMaxIsometricForce();
    }
    return fMax;
}

Matrix calculateNullSpace(const Matrix& A, double atol=1e-13, double rtol=0) {
    Vector s;
    Matrix U, VT;
    FactorSVD svd(A);
    svd.getSingularValuesAndVectors(s, U, VT);
    auto tol = max(atol, rtol * s[0]);
    int nz = 0;
    for (auto v : s) {
        if (v >= tol) nz++;
    }
    Matrix N = VT(nz, 0, VT.nrow() - nz, VT.ncol()).transpose();
    return N;
}

void testNullSpace() {
    Real a[24] = {0.03197292, 0.87063039, 0.13821516, 0.13098836, 0.75546764, 0.42618434,
                  0.87926417, 0.09365083, 0.10514934, 0.57219935, 0.59429898, 0.27111457,
                  0.48186018, 0.73522382, 0.05211834, 0.11899367, 0.08551981, 0.50789553,
                  0.74345284, 0.39507762, 0.64844182, 0.92448094, 0.81439658, 0.66590865};
    Real n[12] = {0.27668194, -0.19129641,
                  -0.21901866, -0.39516696,
                  0.49314911, -0.41430876,
                  -0.75577206, -0.0449566 ,
                  0.19822795,  0.10174348,
                  0.14763525,  0.78944549};
    Matrix A(4, 6, a);
    Matrix Nr(6, 2, n);
    auto N = calculateNullSpace(A);

    // cout << A << endl;
    // cout << N << endl;
    // cout << Nr << endl;

    // check if result is correct
    for (int i = 0; i < Nr.nrow(); i++) {
        assert((N[i] - Nr[i]).norm() < 0.0001);
    }

    // check if this is a null space matrix
    auto AN  = A * N;
    for (int i = 0; i < AN.nrow(); i++) {
        assert(AN[i].norm() < 0.0001);
    }
}

void constructLinearMuscleInequality(const Matrix& NR,
                                     const Vector& fmPar,
                                     const Vector& fmMax,
                                     Matrix& Z,
                                     Vector& b) {
    Z = Matrix(2 * NR.nrow(), NR.ncol());
    Z(0, 0, NR.nrow(), NR.ncol()) = -1.0 * NR;
    Z(NR.nrow(), 0, NR.nrow(), NR.ncol()) = NR;
    b = Vector(2 * fmPar.size());
    b(0, fmPar.size()) = fmPar;
    b(fmPar.size(), fmPar.size()) = fmMax - fmPar;

    // tests
    // assert((b(0, m) - fmPar).norm() < 0.0001);
    // assert((b(m, m) - (fmMax - fmPar)).norm() < 0.0001);
    // for (int i = 0; i < Z.nrow() / 2; i++) {
    //  assert((Z[i] + NR[i]).norm() < 0.0001);
    // }
    // for (int i = Z.nrow() / 2; i < Z.nrow(); i++) {
    //  assert((Z[i] - NR[i]).norm() < 0.0001);
    // }
}

void vertexEnumeration(const Matrix& A, const Vector& b) {
    // variables
    lrs_dic* P;
    lrs_dat* Q;
    lrs_mp_vector output;
    lrs_mp_matrix Lin;   
    long i;
    long col;

    // initialize
    if (!lrs_init((char*) "lrs")) {
	throw OpenSim::Exception("could not initialize lrs library");
    }
        
    // initialize Q data
    Q = lrs_alloc_dat("lrs_data");
    if (Q == NULL) {
	throw OpenSim::Exception("could not allocate Q data");
    }
    Q->m = A.nrow();
    Q->n = A.ncol() + 1;

    // initialize P dictionary
    P = lrs_alloc_dic(Q); 
    if (P == NULL) {
	throw OpenSim::Exception("could not allocate P dictionary");
    }

    // populate P and Q from A and b
    lrs_mp_vector num, den;
    num = lrs_alloc_mp_vector(Q->n);
    den = lrs_alloc_mp_vector(Q->n);
    for (int i = 1; i <= Q->m; i++){
	// TODO populate num, den
	lrs_set_row_mp(P, Q, i, num, den, GE);
    }

    // pivot to a starting dictionary
    output = lrs_alloc_mp_vector(Q->n);
    if (!lrs_getfirstbasis(&P, Q, &Lin, false)) {
	// print non redundant coefficients
        for (col = 0L; col < Q->nredundcol; col++) 
            lrs_printoutput(Q, Lin[col]);

        do {
            for (col = 0; col <= P->d; col++)
                if (lrs_getsolution(P, Q, output, col))
                    lrs_printoutput(Q, output);
        } while (lrs_getnextbasis(&P, Q, FALSE));

        lrs_printtotals(P, Q);

	// free space
        lrs_clear_mp_vector(output, Q->n);
	lrs_clear_mp_vector(num, Q->n);
        lrs_clear_mp_vector(den, Q->n);
        lrs_free_dic(P, Q);
        lrs_free_dat(Q);
	lrs_close((char*) "lrs");
    } else {
	// free space
	lrs_clear_mp_vector(output, Q->n);
	lrs_clear_mp_vector(num, Q->n);
        lrs_clear_mp_vector(den, Q->n);
        lrs_free_dic(P, Q);
        lrs_free_dat(Q);
	lrs_close((char*) "lrs");
	throw OpenSim::Exception("empty set");
    }
}

/******************************************************************************/

FeasibleMuscleForceAnalysis::FeasibleMuscleForceAnalysis() : Analysis() {
    setNull();
}

FeasibleMuscleForceAnalysis::FeasibleMuscleForceAnalysis(const std::string& aFileName)
    : Analysis(aFileName, false) {
    setNull();
    updateFromXMLDocument();
}

int FeasibleMuscleForceAnalysis::begin(const State& s) {
    if (!proceed()) return 0;

    soStorage = new Storage(get_so_results());

    return record(s);
}

int FeasibleMuscleForceAnalysis::step(const State& s, int stepNumber) {
    if (!proceed(stepNumber)) return 0;
    return record(s);
}

int FeasibleMuscleForceAnalysis::end(const SimTK::State & s) {
    if (!proceed()) return 0;
    return record(s);
}

int FeasibleMuscleForceAnalysis::printResults(const string& aBaseName,
                                              const string& aDir,
                                              double aDT,
                                              const string& aExtension) {
    return 0;
}

void FeasibleMuscleForceAnalysis::setNull() {
    setName("FeasibleMuscleForceAnalysis");
    constructProperty_id_results("");
    constructProperty_so_results("");
}

int FeasibleMuscleForceAnalysis::record(const State& s) {
    _model->realizePosition(s);
    // get quantities
    auto m = _model->getMuscles().getSize();
    auto fmMax = calcMaxForce(s, *_model);
    auto R = calculateMomentArm(s, *_model);
    auto NR = calculateNullSpace(R);
    Vector fmPar(m);
    soStorage->getDataAtTime(s.getTime(), m, fmPar);

    // formulate inequlaties
    Matrix Z;
    Vector b;
    constructLinearMuscleInequality(NR, fmPar, fmMax, Z, b);



    exit(0);

    // testNullSpace();
    // exit(0);

    return 0;
}
