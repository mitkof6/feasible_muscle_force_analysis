#include "FeasibleMuscleForceAnalysis.h"
#include <OpenSim/Simulation/Model/Model.h>
#include <OpenSim/Simulation/Model/Muscle.h>

using namespace std;
using namespace OpenSim;
using namespace SimTK;

Matrix calculateMomentArmMatrix(const State& s, const Model& model) {
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
	const string& aDir, double aDT, const string& aExtension) {
	return 0;
}

void FeasibleMuscleForceAnalysis::setNull() {
	setName("FeasibleMuscleForceAnalysis");
	constructProperty_id_results("");
	constructProperty_so_results("");
}

int FeasibleMuscleForceAnalysis::record(const State& s) {
	_model->realizePosition(s);

	auto fMax = calcMaxForce(s, *_model);
	auto R = calculateMomentArmMatrix(s, *_model);

	//cout << s.getTime() << endl;
	//cout << fMax << endl;
	//cout << R << endl;

	return 0;
}