/*
 * TODO:
 *     1. Check if coordinates and muscles are disabled (DONE)
 *     2. Use nonlinear muscle model
 *     3. Use generalized forces instead of so results
 *     4. Use largest polytope instead of the smallest
 **/
#include "FeasibleMuscleForceAnalysis.h"
#include <cassert>
#include <ctime>
#include <limits>
#include <OpenSim/Simulation/Model/Model.h>
#include <OpenSim/Simulation/Model/Muscle.h>
#ifdef __cplusplus
extern "C" {
#endif
// #undef __cplusplus		     // required to work with gmp
#include "lrslib.h"	
#ifdef __cplusplus
}
#endif

using namespace std;
using namespace OpenSim;
using namespace SimTK;

/******************************************************************************/

// start measuring time
#define START_TIME()							      \
    clock_t begin = clock();

// end measuring time
#define END_TIME()					            	      \
    clock_t end = clock();  					  	      \
    double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;               \
    cout << "Elapsed time: " << elapsed_secs << endl;

// used as macro in order to insert __FILE__ and __LINE__
#define THROW_EXCEPTION(msg) throw FileLineException(msg, __FILE__, __LINE__)

/**
 * \brief An exception that prints the file and line number along with the
 * message.
 */
class FileLineException : public std::exception {
    std::string msg;
public:
    FileLineException(const std::string &arg, const char* file, int line)
        : std::exception() {
        std::ostringstream o;
        o << file << ":" << line << ": " << arg;
        msg = o.str();
    }
    ~FileLineException() throw() {};
    const char* what() const throw() { return msg.c_str(); }
};

/**
 * Calculates the muscle moment arm matrix R[n, m] (n DoFs, m muscles) given the
 * current state. The calculations assume that all the coordinates and muscles
 * are enabled.
 */
Matrix calculateMomentArm(const State& s, const Model& model) {
    const auto& coordinateSet = model.getCoordinateSet();
    const auto& muscleSet = model.getMuscles();
    Matrix R(coordinateSet.getSize(), muscleSet.getSize());
    for (int i = 0; i < coordinateSet.getSize(); i++) {
	for (int j = 0; j < muscleSet.getSize(); j++) {
            R[i][j] = muscleSet[j].computeMomentArm(s, coordinateSet[i]);
        }
    }
    return R;
}

/**
 * Get the maximum isometric force of the muscles
 */
Vector getMaxIsometricForce(const Model& model) {
    const auto& muscleSet = model.getMuscles();
    Vector fMax(muscleSet.getSize());
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
    return VT(nz, 0, VT.nrow() - nz, VT.ncol()).transpose();
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

    // test
    // assert((b(0, m) - fmPar).norm() < 0.0001);
    // assert((b(m, m) - (fmMax - fmPar)).norm() < 0.0001);
    // for (int i = 0; i < Z.nrow() / 2; i++) {
    //  assert((Z[i] + NR[i]).norm() < 0.0001);
    // }
    // for (int i = Z.nrow() / 2; i < Z.nrow(); i++) {
    //  assert((Z[i] - NR[i]).norm() < 0.0001);
    // }
}

long gcd(long a, long b) {
    if (a == 0)
        return b;
    else if (b == 0)
        return a;

    if (a < b)
        return gcd(a, b % a);
    else
        return gcd(b, a % b);
}

void fraction(const double& input, long& numerator, long& denominator) {
    double integral = std::floor(input);
    double frac = input - integral;

    const long precision = 1000000000; // This is the accuracy.

    long gcd_ = gcd(round(frac * precision), precision);

    denominator = precision / gcd_;
    numerator = round(frac * precision) / gcd_;

    // std::cout << integral << " + ";
    // std::cout << numerator << " / " << denominator << std::endl;
}

Matrix vertexEnumeration(const Matrix& A, const Vector& b) {
    // initialize
    if (!lrs_init((char*) "lrs")) {
	THROW_EXCEPTION("could not initialize lrs library");
    }

    // initialize Q data
    lrs_dat* Q = lrs_alloc_dat("lrs_data");
    if (Q == NULL) {
	THROW_EXCEPTION("could not allocate lrs Q data");
    }
    Q->m = A.nrow();
    Q->n = A.ncol() + 1;

    // initialize P dictionary
    lrs_dic* P = lrs_alloc_dic(Q);
    if (P == NULL) {
	THROW_EXCEPTION("could not allocate lrs P dictionary");
    }

    // populate P and Q from A and b
    lrs_mp_vector num, den;
    num = lrs_alloc_mp_vector(Q->n);
    den = lrs_alloc_mp_vector(Q->n);
    for (int i = 1; i <= Q->m; i++){
        for (int j = 0; j < Q->n; j++) {
            if (j > 0) {
		double value = -A[i - 1][j - 1];
		// cout << value << endl;
		if (std::abs(value) < 0.001) value = 0.0;
       /*         mpq_t op;
                mpq_init(op);
                mpq_set_d(op, value);
                mpq_canonicalize(op);
		itomp(mpz_get_si(mpq_numref(op)), num[j]);
                itomp(mpz_get_si(mpq_denref(op)), den[j]);
                mpq_clear(op);*/
        long nnum, dden;
        fraction(value, nnum, dden);
        *num[j] = nnum;
        *den[j] = dden;
            } else {
		double value = b[i - 1];
		// cout << value << endl;
		if (std::abs(value) < 0.001) value = 0.0;
                /*mpq_t op;
                mpq_init(op);
                mpq_set_d(op, value);
                mpq_canonicalize(op);
		itomp(mpz_get_si(mpq_numref(op)), num[j]);
                itomp(mpz_get_si(mpq_denref(op)), den[j]);
                mpq_clear(op);*/
        long nnum, dden;
        fraction(value, nnum, dden);
        *num[j] = nnum;
        *den[j] = dden;
	    }
        }
	 lrs_printoutput (Q, num);
	 lrs_printoutput (Q, den);
        lrs_set_row_mp(P, Q, i, num, den, GE);
    }
    lrs_clear_mp_vector(num, Q->n);
    lrs_clear_mp_vector(den, Q->n);

    // try pivot to a starting dictionary
    lrs_mp_matrix Lin;
    if (!lrs_getfirstbasis(&P, Q, &Lin, FALSE)) {
        // free space
	if (Q->nredundcol > 0) lrs_clear_mp_matrix(Lin, Q->nredundcol, Q->n);
        lrs_free_dic(P, Q);
        lrs_free_dat(Q);
        lrs_close((char*) "lrs");
	return Matrix(0, 0);
        // throw OpenSim::Exception("empty set: infeasible problem");
    }

    // print non redundant coefficients
    // for (int col = 0L; col < Q->nredundcol; col++) {
    //     lrs_printoutput(Q, Lin[col]);
    // }

    // obtain the vertices of the polytope
    lrs_mp_vector output = lrs_alloc_mp_vector(Q->n);
    vector<double> result;
    int rows = 0;
    do {
	// get solution into output
        for (int col = 0; col <= P->d; col++) {
	    // if (lrs_getsolution (P, Q, output, col)) {
	    // 	 lrs_printoutput (Q, output);
	    // }
	    lrs_getsolution(P, Q, output, col);
	}
	// transform as double
	for (int i = 1; i < Q->n; i++) {
	    double res;
	    rattodouble(output[i], output[0], &res);
	    result.push_back(res);
	}
	rows++;
    } while (lrs_getnextbasis(&P, Q, FALSE));
    lrs_printtotals(P, Q); 

    // free space
    lrs_clear_mp_vector(output, Q->n);
    if (Q->nredundcol > 0) lrs_clear_mp_matrix(Lin, Q->nredundcol, Q->n); 
    lrs_free_dic(P, Q);
    lrs_free_dat(Q);
    lrs_close((char*) "lrs");
    
    return Matrix(rows, P->d, &result[0]);
}

void testVertexEnumeration() {
    // define a cube
    Real a[18] = { 0,  0, -1,
                   0, -1,  0,
                   1,  0,  0,
                   -1,  0,  0,
                   0,  1,  0,
                   0,  0,  1};
    Real s[24] = {-0.5,  0.5,  0.5,
                  0.5, -0.5,  0.5,
                  0.5, -0.5, -0.5,
                  0.5,  0.5,  0.5,
                  0.5,  0.5, -0.5,
                  -0.5, -0.5,  0.5,
                  -0.5, -0.5, -0.5,
                  -0.5,  0.5, -0.5};
    Matrix A(6, 3, a);
    Vector b(6, 0.5);
    Matrix SR(8, 3, s);
    Matrix S = vertexEnumeration(A, b);
    cout << S<< endl;
    for (int i = 0; i < S.nrow(); i++) {
	assert((S[i] - SR[i]).norm() < 0.0001);
    }
}

/******************************************************************************/

FeasibleMuscleForceAnalysis::FeasibleMuscleForceAnalysis() : Analysis() {
    setNull();
}

FeasibleMuscleForceAnalysis::FeasibleMuscleForceAnalysis(const std::string& fileName)
    : Analysis(fileName, false) {
    setNull();
    updateFromXMLDocument();
}

int FeasibleMuscleForceAnalysis::begin(const State& s) {
    if (!proceed()) return 0;
    testVertexEnumeration();
    exit(0);
    // check if muscles are disabled
    const auto& muscleSet = _model->getMuscles();
    for (int i = 0; i < muscleSet.getSize(); i++) {
	if (!muscleSet[i].appliesForce(s)) {
	    THROW_EXCEPTION("disabled muscles are not supported");
	}
    }
    
    // setup internal variables
    fmMax = getMaxIsometricForce(*_model);
    feasibleMuscleForces.clear();
    smallestPolytope = numeric_limits<int>::max();
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

int FeasibleMuscleForceAnalysis::printResults(const string& baseName,
                                              const string& dir,
                                              double dt,
                                              const string& extension) {
    // create column labels
    Array<string> labels;
    labels.append("time");
    const auto& muscleSet = _model->getMuscles();
    for (int i = 0; i < muscleSet.getSize(); i++) {
	labels.append(muscleSet[i].getName());
    }

    // storage description
    string description;
    description = "\nThis file contains the feasible muscle forces.\n";
    description += "\nUnits are S.I. units (seconds, meters, Newtons, ...)\n\n";

    for (int i = 0; i < smallestPolytope; i++) {
	Storage storage;
	storage.setName(to_string(i));
	storage.setColumnLabels(labels);
	storage.setDescription(description);
	for (int j = 0; j < feasibleMuscleForces.size(); j++) {
	    storage.append(feasibleMuscleForces[j].first, feasibleMuscleForces[j].second[i]);
	}
	Storage::printResult(&storage, baseName + "_" + storage.getName(), dir, dt, extension);
    }
    
    return 0;
}

void FeasibleMuscleForceAnalysis::setNull() {
    setName("FeasibleMuscleForceAnalysis");
    constructProperty_id_results("");
    constructProperty_so_results("");
}

int FeasibleMuscleForceAnalysis::record(const State& s) {
    cout << endl << "Analyze at time: " << s.getTime() << endl;
    
    // get necessary quantities
    _model->realizePosition(s);
    auto m = _model->getMuscles().getSize();
    auto R = calculateMomentArm(s, *_model);
    auto NR = calculateNullSpace(R);
    Vector fmPar(m);
    soStorage->getDataAtTime(s.getTime(), m, fmPar);
  
    // formulate feasible inequalities
    Matrix Z;
    Vector b;
    constructLinearMuscleInequality(NR, fmPar, fmMax, Z, b);

    // find the vertices of the feasible polytope
    START_TIME();
    auto fm0 = vertexEnumeration(Z, b);
    END_TIME();

    // check if vertex enumeration returned valid solution
    if (fm0.nrow() == 0 || fm0.ncol() == 0) {
	cout << "*** warning: invalid vertex enumeration solution at: "
	     << s.getTime() << " ***" << endl;
	return 0;
    }

    // calculate the feasible muscle forces
    vector<Vector> fm;
    for (int i = 0; i < fm0.nrow(); i++) {
	fm.push_back(fmPar + NR * ~fm0[i]);
    }
    feasibleMuscleForces.push_back({s.getTime(), fm});

    // keep track of the smallest polytope
    if (fm.size() < smallestPolytope) {
	smallestPolytope = fm.size();
    }
    
    return 0;
}
