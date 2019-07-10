/**
 * @file FeasibleMuscleForceAnalysis.h
 *
 * \brief TODO
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
		FeasibleMuscleForceAnalysis(const std::string& aFileName);
		int begin(const SimTK::State& s) override;
		int step(const SimTK::State& s, int stepNumber) override;
		int end(const SimTK::State& s) override;
		int printResults(const std::string& aBaseName,
					     const std::string& aDir = "",
						 double aDT = -1.0,
						 const std::string& aExtension = ".sto") override;
	private:
		void setNull();
		int record(const SimTK::State& s);
	private:
		SimTK::ReferencePtr<OpenSim::Storage> idStorage, soStorage;
    };
}

#endif
