/*
 * TODO:
 *     1. Check if coordinates and muscles are disabled
 *     2. Use nonlinear muscle model
 *     3. Use generalized forces instead of so results
 **/
#include "FeasibleMuscleForceAnalysis.h"
#include <cassert>
#include <ctime>
#include <limits>
#include <OpenSim/Simulation/Model/Model.h>
#include <OpenSim/Simulation/Model/Muscle.h>
extern "C" {
#undef __cplusplus		// required to work with gmp
#include "lrslib.h"		// must be enclosed in extern "C"
};

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


Matrix calcMomentArm(const State& s, const Model& model) {
    const auto& coordinateSet = model.getCoordinateSet();
    const auto& muscleSet = model.getMuscles();
    Matrix R(coordinateSet.getSize(), muscleSet.getSize(), 0.0);
    for (int i = 0; i < coordinateSet.getSize(); i++) {
	for (int j = 0; j < muscleSet.getSize(); j++) {
            R[i][j] = muscleSet[j].computeMomentArm(s, coordinateSet[i]);
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

Matrix calcNullSpace(const Matrix& A, double atol=1e-13, double rtol=0) {
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
    auto N = calcNullSpace(A);

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

Matrix vertexEnumeration(const Matrix& A, const Vector& b) {
    // initialize
    if (!lrs_init((char*) "lrs")) {
        throw OpenSim::Exception("could not initialize lrs library");
    }

    // initialize Q data
    lrs_dat* Q = lrs_alloc_dat("lrs_data");
    if (Q == NULL) {
        throw OpenSim::Exception("could not allocate Q data");
    }
    Q->m = A.nrow();
    Q->n = A.ncol() + 1;

    // initialize P dictionary
    lrs_dic* P = lrs_alloc_dic(Q);
    if (P == NULL) {
        throw OpenSim::Exception("could not allocate P dictionary");
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
                mpq_t op;
                mpq_init(op);
                mpq_set_d(op, value);
                mpq_canonicalize(op);
		itomp(mpz_get_si(mpq_numref(op)), num[j]);
                itomp(mpz_get_si(mpq_denref(op)), den[j]);
                mpq_clear(op);
            } else {
		double value = b[i - 1];
		// cout << value << endl;
		if (std::abs(value) < 0.001) value = 0.0;
                mpq_t op;
                mpq_init(op);
                mpq_set_d(op, value);
                mpq_canonicalize(op);
		itomp(mpz_get_si(mpq_numref(op)), num[j]);
                itomp(mpz_get_si(mpq_denref(op)), den[j]);
                mpq_clear(op);
	    }
        }
	// lrs_printoutput (Q, num);
	// lrs_printoutput (Q, den);
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
        throw OpenSim::Exception("empty set: infeasible problem");
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

FeasibleMuscleForceAnalysis::FeasibleMuscleForceAnalysis(const std::string& aFileName)
    : Analysis(aFileName, false) {
    setNull();
    updateFromXMLDocument();
}

int FeasibleMuscleForceAnalysis::begin(const State& s) {
    if (!proceed()) return 0;

    fs.clear();
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

int FeasibleMuscleForceAnalysis::printResults(const string& aBaseName,
                                              const string& aDir,
                                              double aDT,
                                              const string& aExtension) {
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
	for (int j = 0; j < fs.size(); j++) {
	    storage.append(fs[j].first, fs[j].second[i]);
	}
	Storage::printResult(&storage, aBaseName + "_" + storage.getName(),
			     aDir, aDT, aExtension);
    }
    
    return 0;
}

void FeasibleMuscleForceAnalysis::setNull() {
    setName("FeasibleMuscleForceAnalysis");
    constructProperty_id_results("");
    constructProperty_so_results("");
}

int FeasibleMuscleForceAnalysis::record(const State& s) {
    // if (s.getTime() < 0.05) return 0;
    cout << endl << "Analyze at time: " << s.getTime() << endl;
    
    // get necessary quantities
    _model->realizePosition(s);
    auto m = _model->getMuscles().getSize();
    auto fmMax = calcMaxForce(s, *_model);
    auto R = calcMomentArm(s, *_model);
    auto NR = calcNullSpace(R);
    Vector fmPar(m);
    soStorage->getDataAtTime(s.getTime(), m, fmPar);

    // double n[] = {
    // 		  -0.46146, -0.0238014, 0.368391, 0.0103206, -0.0358229, 0.0337891, -0.211974, -0.236186, 0.171261, 0.115441, 0.385872, 0.126023,
    // 		  0.0501787, 0.879929, 0.0433098, -0.00133284, 0.00462632, -0.00436367, 0.0273752, 0.0305021, -0.0221173, -0.0149085, -0.0498331, -0.0162751,
    // 		  0.0302064, 0.194628, -0.725281, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    // 		  -0.459017, -0.0415017, -0.448695, -0.288797, 0.0828552, -0.329365, 0.0598763, -0.131832, -0.00476989, 0.0886512, 0.279753, 0.262193,
    // 		  -0.133447, 0.283149, -0.0807475, 0.298292, -0.130292, 0.370154, -0.347349, -0.191602, 0.237714, 0.0698152, 0.249674, -0.0865959,
    // 		  -0.116413, 0.302133, 0.332922, -0.138576, 0.217269, -0.321435, 0.187709, -0.17668, -0.299387, -0.0112508, -0.109504, 0.267366,
    // 		  0.0686642, -0.092554, -0.100786, 0.343647, 0.199077, 0.104364, -0.102035, -0.465032, -0.301166, -0.0165876, -0.163447, 0.191762,
    // 		  -0.106521, 0.0574287, 0.0602131, -0.261443, -0.230864, -0.00862777, -0.0443883, 0.330641, 0.394781, -0.388529, -0.103044, 0.219663,
    // 		  -0.0444961, -0.0470447, -0.0541152, 0.114782, -0.032149, 0.122772, -0.18565, -0.184082, 0.105189, -0.503853, -0.337065, 0.518146,
    // 		  -0.183061, -0.00114376, -0.00793102, 0.618069, -0.0311461, -0.294781, 0.200284, 0.274326, 0.0701886, -0.0385367, 0.0353692, 0.0548506,
    // 		  0.00510827, 3.19163e-05, 0.000221313, -0.0362911, 0.833371, 0.124004, -0.0976747, 0.174014, 0.261375, -0.0929982, 0.0545053, -0.0503825,
    // 		  -0.172368, -0.00107695, -0.00746779, -0.289328, 0.128696, 0.635141, 0.258017, 0.0668163, -0.182348, 0.0524662, -0.0190498, 0.0977353,
    // 		  -0.163062, -0.00101881, -0.0070646, 0.211395, -0.0934018, 0.263622, 0.759881, -0.124947, 0.165746, -0.0712031, 0.053258, 0.0297087,
    // 		  -0.283448, -0.00177098, -0.0122803, 0.300433, 0.181252, 0.0829554, -0.118897, 0.536361, -0.252098, 0.0692146, -0.0212613, 0.151859,
    // 		  0.178358, 0.00111438, 0.00772729, 0.0686845, 0.256404, -0.178452, 0.175232, -0.22899, 0.576215, 0.16451, -0.109354, 0.0124899,
    // 		  -0.0113545, -7.09424e-05, -0.000491927, -0.02127, -0.0931609, 0.0683861, -0.0565118, 0.0943308, 0.14778, 0.564561, -0.33224, 0.297817,
    // 		  0.38126, 0.0023821, 0.0165179, 0.0359162, 0.0437746, -0.00717825, 0.0768864, 0.0339584, -0.11302, -0.296245, 0.585048, 0.186212,
    // 		  0.419579, 0.00262152, 0.0181781, 0.0155006, -0.0610769, 0.0731817, 0.0201245, 0.150767, 0.0473816, 0.334952, 0.26942, 0.565735};
    // double a[] ={0.0921391, 5.40254, 0.0654329, 83.7385, 26.6831, 0.117335, 7.79152, 2.41927, 61.5669, 284.912, 0.139795, 225.933, 0.0587067, 0.32315, 822.66, 0.135622, 50.1176, 5.44083};
    // NR = Matrix(18, 12, n);
    // fmPar = Vector(18, a);
    
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
	cout << "Invalid vertex enumeration solution at: " << s.getTime() << endl;
	return 0;
    }

    // calculate the feasible muscle forces
    vector<Vector> fm;
    for (int i = 0; i < fm0.nrow(); i++) {
	if (abs(s.getTime() - 0.199999) < 0.000001) {
	    cout << NR << endl;
	    cout << fmPar << endl;
	    //cout << fm0 << endl;
	}
	fm.push_back(fmPar + NR * ~fm0[i]);
    }
    fs.push_back({s.getTime(), fm});

    // keep track of the smallest polytope
    if (fm.size() < smallestPolytope) {
	smallestPolytope = fm.size();
    }
    
    // testVertexEnumeration();
    // testNullSpace();
    
    return 0;
}
