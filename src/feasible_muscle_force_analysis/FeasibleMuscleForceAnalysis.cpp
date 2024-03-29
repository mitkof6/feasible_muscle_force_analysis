/*
 * Things to improve:
 *     1. Check if coordinates and muscles are disabled (DONE)
 *     2. Use nonlinear muscle model (DONE: not solving dynamics)
 *     3. Use generalized forces instead of so results (DONE)
 *     4. Use largest polytope instead of the smallest (DONE)
 *     5. Sample the feasible polytope to generate internal solutions (DONE)
 **/
#include "FeasibleMuscleForceAnalysis.h"
#include <cassert>
#include <ctime>
#include <limits>
#include <regex>
#include <OpenSim/Simulation/Model/Model.h>
#include <OpenSim/Simulation/Model/Muscle.h>
extern "C" {
#undef __cplusplus
#include "lrslib/lrslib.h"
    FILE* lrs_ofp = stdout;
    FILE* lrs_ifp = stdin;
}

using namespace std;
using namespace OpenSim;
using namespace SimTK;

/******************************************************************************/

namespace {

    // start measuring time
#define START_TIME()				\
    clock_t begin = clock();

    // end measuring time
#define END_TIME()						\
    clock_t end = clock();                                      \
    double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC; \
    cout << "Elapsed time: " << elapsed_secs << endl;

    // used as macro in order to insert __FILE__ and __LINE__
#define THROW_EXCEPTION(msg) throw FileLineException(msg, __FILE__, __LINE__)

    /**
     * \brief An exception that prints the file and line number along with the
     * message.
     */
    class FileLineException : public exception {
        string msg;
    public:
        FileLineException(const string &arg, const char* file, int line)
            : exception() {
            ostringstream o;
            o << file << ":" << line << ": " << arg;
            msg = o.str();
        }
        ~FileLineException() throw() {};
        const char* what() const throw() { return msg.c_str(); }
    };

    Matrix calculateMomentArm(const State& s, const Model& model,
                              const vector<int>& activeCoordinateIndices,
                              const vector<int>& activeMuscleIndices) {
        const auto& coordinateSet = model.getCoordinateSet();
        const auto& muscleSet = model.getMuscles();
        Matrix R(activeCoordinateIndices.size(), activeMuscleIndices.size());
        for (int i = 0; i < activeCoordinateIndices.size(); i++) {
            for (int j = 0; j < activeMuscleIndices.size(); j++) {
                R[i][j] = muscleSet[activeMuscleIndices[j]].computeMomentArm(
                    s, coordinateSet[activeCoordinateIndices[i]]);
            }
        }
        return R;
    }

    Vector getMuscleRelatedValues(const Model& model,
                                  const vector<int>& activeMuscleIndices,
                                  function<double(const Muscle&)> operation) {
        const auto& muscleSet = model.getMuscles();
        Vector vec(activeMuscleIndices.size());
        for (int i = 0; i < activeMuscleIndices.size(); i++) {
            vec[i] = operation(muscleSet[activeMuscleIndices[i]]);
        }
        return vec;
    }

    Matrix pInv(const Matrix& A) {
        Matrix AInv;
        FactorSVD svd(A);
        svd.inverse(AInv);
        return AInv;
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
            SimTK_ASSERT_ALWAYS((N[i] - Nr[i]).norm() < 0.0001, "");
        }

        // check if this is a null space matrix
        auto AN  = A * N;
        for (int i = 0; i < AN.nrow(); i++) {
            SimTK_ASSERT_ALWAYS(AN[i].norm() < 0.0001, "");
        }
    }

    /**
     * Feasible inequality for a linear muscle model
     *
     *     fm = fmMax am
     *
     *     Z = [-NR; NR]
     *     b = [fmPar; fmMax - fmPar]
     *
     * For detailed derivations:
     * <a href="https://doi.org/10.1371/journal.pone.0209171">[Publication]</a>
     **/
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

        // // validation
        // int m = fmPar.size();
        // SimTK_ASSERT_ALWAYS((b(0, m) - fmPar).norm() < 0.0001, "");
        // SimTK_ASSERT_ALWAYS((b(m, m) - (fmMax - fmPar)).norm() < 0.0001, "");
        // for (int i = 0; i < m; i++) {
        //      SimTK_ASSERT_ALWAYS((Z[i] + Z[i + m]).norm() < 0.0001, "");
        // }
    }

    /**
     * Feasible inequality for a nonlinear muscle model
     *
     *     fm = fmMax (am fmL fmV + fmPE)
     *
     *     Z = [-NR; NR]
     *     b = [fmPar - fmMax o cosAlpha o fmPE; fmMax o (fmL o fmV + fmPE) o cosAlpha - fmPar]
     *
     * For detailed derivations:
     * <a href="https://doi.org/10.1371/journal.pone.0209171">[Publication]</a>
     **/
    void constructNonlinearMuscleInequality(const Matrix& NR,
                                            const Vector& fmPar,
                                            const Vector& fmMax,
                                            const Vector& fmL,
                                            const Vector& fmV,
                                            const Vector& fmPE,
                                            const Vector& cosAlpha,
                                            Matrix& Z,
                                            Vector& b) {
        Z = Matrix(2 * NR.nrow(), NR.ncol());
        Z(0, 0, NR.nrow(), NR.ncol()) = -1.0 * NR;
        Z(NR.nrow(), 0, NR.nrow(), NR.ncol()) = NR;
        b = Vector(2 * fmPar.size());
        b(0, fmPar.size()) = fmPar -
            fmMax.elementwiseMultiply(fmPE).elementwiseMultiply(cosAlpha);
        b(fmPar.size(), fmPar.size()) =
            fmMax.elementwiseMultiply(
                fmL.elementwiseMultiply(fmV) + fmPE).elementwiseMultiply(cosAlpha) - fmPar;
    }

    /**
     * Find the V-representations from the H-representation (vertex
     * enumeration).  In our problem Z fm0 <= b is convex and bounded, therefore
     * through vertex enumeration we can find the extreme vertices of a
     * polytope.
     **/
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
                double value;
                if (j > 0) {
                    value = -A[i - 1][j - 1];
                } else {
                    value = b[i - 1];
                }
                if (abs(value) < 0.001) value = 0.0;
                mpq_t op;
                mpq_init(op);
                mpq_set_d(op, value);
                mpq_canonicalize(op);
                itomp(mpz_get_si(mpq_numref(op)), num[j]);
                itomp(mpz_get_si(mpq_denref(op)), den[j]);
                // cout << value << endl;
                // cout << mpz_get_si(mpq_numref(op)) << " / "
                //      << mpz_get_si(mpq_denref(op)) <<  endl;
            }
            //lrs_printoutput (Q, num);
            //lrs_printoutput (Q, den);
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
                /*           if (lrs_getsolution (P, Q, output, col)) {
                             lrs_printoutput (Q, output);
                             }*/
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
                      0.5,  0.5,  0.5,
                      0.5, -0.5,  0.5,
                      0.5, -0.5, -0.5,
                      0.5,  0.5, -0.5,
                      -0.5, -0.5,  0.5,
                      -0.5, -0.5, -0.5,
                      -0.5,  0.5, -0.5};
        Matrix A(6, 3, a);
        Vector b(6, 0.5);
        Matrix SR(8, 3, s);
        Matrix S = vertexEnumeration(A, b);
        // cout << S<< endl;
        for (int i = 0; i < S.nrow(); i++) {
            SimTK_ASSERT_ALWAYS((S[i] - SR[i]).norm() < 0.0001,
                                "testVertexEnumeration failed");
        }
    }
};

/******************************************************************************/

FeasibleMuscleForceAnalysis::FeasibleMuscleForceAnalysis() : Analysis() {
    constructProperties();
}

FeasibleMuscleForceAnalysis::FeasibleMuscleForceAnalysis(const string& fileName)
    : Analysis(fileName, false) {
    constructProperties();
    updateFromXMLDocument();
}

int FeasibleMuscleForceAnalysis::begin(const State& s) {
    if (!proceed()) return 0;

    // initialize active model coordinates
    const auto& coordinateSet = _model->getCoordinateSet();
    auto excludedCoordinatesRegex = regex(get_excluded_coordiantes());
    cout << "Active coordinates:" << endl;
    for (int i = 0; i < coordinateSet.getSize(); i++) {
        if (!regex_match(coordinateSet[i].getName(), excludedCoordinatesRegex)) {
            activeCoordinateIndices.push_back(i);
            cout << coordinateSet[i].getName() << endl;
        }
    }

    // check if muscles are disabled and append active muscle coordinates
    const auto& muscleSet = _model->getMuscles();
    auto excludedMusclesRegex = regex(get_excluded_muscles());
    cout << "Active muscles:" << endl;
    for (int i = 0; i < muscleSet.getSize(); i++) {
        if (!muscleSet[i].appliesForce(s)) {
            THROW_EXCEPTION("disabled muscles are not supported");
        }
        if (!regex_match(muscleSet[i].getName(), excludedMusclesRegex)) {
            activeMuscleIndices.push_back(i);
            cout << muscleSet[i].getName() << endl;
        }
    }

    // setup internal variables
    fmMax = getMuscleRelatedValues(*_model, activeMuscleIndices,
				   [](const Muscle& m)->double {
				       return m.getMaxIsometricForce();
				   });
    feasibleMuscleForces.clear();
    largestPolytope = 0;
    idStorage = Storage(get_id_results());

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

    // create feasible force storages: create M files where M is the dimension
    // of the large polytope and each file containing a feasible muscle force
    // solution for each instance of the analysis.
    for (int i = 0; i < largestPolytope; i++) {
        Storage storage;
        storage.setName(to_string(i));
        storage.setColumnLabels(labels);
        storage.setDescription(description);
        for (int j = 0; j < feasibleMuscleForces.size(); j++) {
            Vector temp;
            if (i < feasibleMuscleForces[j].second.size()) {
                // if within the bounds get the current vector i
                temp = feasibleMuscleForces[j].second[i];
            } else {
                // else use the last solution that satisfies the solution at t
                temp = feasibleMuscleForces[j].second.back();
            }
            // fill inactive muscles with 0
            Vector fm(muscleSet.getSize(), 0.0);
            for (int k = 0; k < activeMuscleIndices.size(); k++) {
                fm[activeMuscleIndices[k]] = temp[k];
            }
            storage.append(feasibleMuscleForces[j].first, fm);
        }
        Storage::printResult(&storage, baseName + "_" + storage.getName(),
                             dir, dt, extension);
    }

    return 0;
}

void FeasibleMuscleForceAnalysis::constructProperties() {
    setName("FeasibleMuscleForceAnalysis");
    constructProperty_id_results("");
    constructProperty_use_linear_muscle_model(true);
    constructProperty_excluded_coordiantes("");
    constructProperty_excluded_muscles("");
    constructProperty_convex_sampling_depth(0);
}

int FeasibleMuscleForceAnalysis::record(const State& s) {
    cout << endl << "Analyze at time: " << s.getTime() << endl;

    // get necessary quantities
    auto m = _model->getMuscles().getSize();
    auto n = s.getNQ();
    auto R = calculateMomentArm(s, *_model,
                                activeCoordinateIndices,
                                activeMuscleIndices);
    auto NR = calculateNullSpace(R);
    auto RInv = pInv(R);

    // fmPar = R + tau
    Vector tauTemp(n);
    idStorage.getDataAtTime(s.getTime(), n, tauTemp);
    Vector tau(activeCoordinateIndices.size());
    for (int i = 0; i < activeCoordinateIndices.size(); i++) {
        tau[i] = tauTemp[activeCoordinateIndices[i]];
    }
    auto fmPar = RInv * tau;

    // formulate feasible inequalities
    Matrix Z;
    Vector b;
    if (get_use_linear_muscle_model()) {
        constructLinearMuscleInequality(NR, fmPar, fmMax, Z, b);
    } else{
        auto fmL = getMuscleRelatedValues(*_model, activeMuscleIndices,
                                          [&](const Muscle& m)->double {
                                              return m.getActiveForceLengthMultiplier(s);
                                          });
        auto fmV = getMuscleRelatedValues(*_model, activeMuscleIndices,
					  [&](const Muscle& m)->double {
					      return m.getForceVelocityMultiplier(s);
					  });
        auto fmPE = getMuscleRelatedValues(*_model, activeMuscleIndices,
					   [&](const Muscle& m)->double {
					       return m.getPassiveForceMultiplier(s);
					   });
        auto cosAlpha = getMuscleRelatedValues(*_model, activeMuscleIndices,
					       [&](const Muscle& m)->double {
						   return m.getCosPennationAngle(s);
					       });
        constructNonlinearMuscleInequality(NR, fmPar, fmMax, fmL, fmV, fmPE,
                                           cosAlpha, Z, b);
    }

    // find the vertices of the feasible polytope
    // START_TIME();
    auto fm0 = vertexEnumeration(Z, b);
    // END_TIME();

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

    // since the feasible space is a convex set we can find additional solution
    // in the form z = a * x_i + (1 - a) x_j
    for (int i = 0; i < get_convex_sampling_depth(); i++) {
        int n = fm.size();
        for (int j = 0; j < n; j++) {
            for (int k = 0; k < n; k++) {
                if (i == j) continue;
                double a = 0.5;  // midpoint interpolation
                auto x1 = fm[j];
                auto x2 = fm[k];
                fm.push_back(a * x1 + (1 - a) * x2);
            }
        }
    }
    feasibleMuscleForces.push_back({s.getTime(), fm});

    // keep track of the largest polytope
    if (fm.size() > largestPolytope) {
        largestPolytope = fm.size();
    }

    return 0;
}
