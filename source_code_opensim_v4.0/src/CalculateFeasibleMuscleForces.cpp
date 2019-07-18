/**
 * @file CalculateFeasibleMuscleForces.h
 *
 * \brief This executable initiates the FeasibleMuscleForceAnalysis, which can
 * also be used as an OpenSim plugin.
 *
 * @author Dimitar Stanev <stanev@ece.upatras.gr>
 *
 * @see <a href="https://simtk.org/projects/redundancy">[SimTK Project]</a>, <a
 * href="https://doi.org/10.1371/journal.pone.0209171">[Publication]</a>
 */
#include <iostream>
#include <OpenSim/OpenSim.h>
#include "feasible_muscle_force_analysis/FeasibleMuscleForceAnalysis.h"
#include "Settings.h"

using namespace OpenSim;
using namespace std;

void run() {
    auto subjectDir = DATA_DIR + "/gait1018/";
    auto resultsDir = subjectDir + "results/";
    auto outputDir = resultsDir + "feasible_muscle_forces/";
    auto modelFile = subjectDir + "subject01_scaled.osim";
    auto ikFile = resultsDir + "subject01_walk1_ik.mot";
    auto idFile = resultsDir + "subject01_inverse_dynamics.sto";
    auto soFile = resultsDir + "subject01_StaticOptimization_force.sto";

    Model model(modelFile);
    Storage motion(ikFile);

    FeasibleMuscleForceAnalysis* feasibleForces = new FeasibleMuscleForceAnalysis();
    feasibleForces->setModel(model);
    feasibleForces->set_so_results(soFile);
    feasibleForces->set_id_results(idFile);
    model.addAnalysis(feasibleForces);

    AnalyzeTool analysis;
    analysis.setName("FeasibleMuscleForces");
    analysis.setModel(model);
    analysis.setInitialTime(motion.getFirstTime());
    analysis.setFinalTime(motion.getLastTime());
    analysis.setLowpassCutoffFrequency(6);
    analysis.setCoordinatesFileName(ikFile);
    analysis.setLoadModelAndInput(true);
    analysis.setResultsDir(outputDir);
    analysis.run();
}

int main(int argc, char *argv[]) {
    try {
        run();
        getchar();
    }
    catch (exception &e) {
        cout << e.what() << endl;
        getchar();
        return -1;
    }
    return 0;
}
