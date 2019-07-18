/**
 * @file FeasibleMuscleForceAnalysis.h
 *
 * \brief This analysis calculates the feasible muscle forces that satisfy both
 * the motion under study and the physiological muscle constraints.
 *
 * @author Dimitar Stanev <stanev@ece.upatras.gr>
 *
 * @see <a href="https://simtk.org/projects/redundancy">[SimTK Project]</a>, <a
 * href="https://doi.org/10.1371/journal.pone.0209171">[Publication]</a>
 */
#ifndef FEASIBLE_MUSCLE_FORCE_ANALYSIS_H
#define FEASIBLE_MUSCLE_FORCE_ANALYSIS_H

#include "internal/FeasibleMuscleForceAnalysisExports.h"
#include <OpenSim/Simulation/Model/Analysis.h>

namespace OpenSim {
    /**
     * \brief TODO
     */
    class FeasibleMuscleForceAnalysis_API FeasibleMuscleForceAnalysis : public Analysis {
        OpenSim_DECLARE_CONCRETE_OBJECT(FeasibleMuscleForceAnalysis, Analysis);
    public:
        OpenSim_DECLARE_PROPERTY(id_results, std::string,
                                 "Generalized forces output storage from inverse dynamics tool");
        OpenSim_DECLARE_PROPERTY(so_results, std::string,
                                 "Muscle forces output storage from static optimization tool");
    public:
        FeasibleMuscleForceAnalysis();
        FeasibleMuscleForceAnalysis(const std::string& fileName);
        int begin(const SimTK::State& s) override;
        int step(const SimTK::State& s, int step) override;
        int end(const SimTK::State& s) override;
        int printResults(const std::string& baseName,
                         const std::string& dir = "",
                         double dt= -1.0,
                         const std::string& extension = ".sto") override;
    private:
        void setNull();
        int record(const SimTK::State& s);
    private:
        SimTK::ReferencePtr<OpenSim::Storage> idStorage, soStorage;
	SimTK::Vector fmMax;
	// [{t0, feasible muscle forces}, ..., {tf, feasible muscle forces}]
	std::vector<std::pair<double, std::vector<SimTK::Vector>>> feasibleMuscleForces;
	int smallestPolytope;
    };
}

#endif
