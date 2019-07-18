/**
 * @file CalculateFeasibleJointReactionLoads.h
 *
 * \brief This executable is called after CalculateFeasibleMuscleForces to
 * calculate the feasible joint reaction loads for a particular movement.
 *
 * @author Dimitar Stanev <stanev@ece.upatras.gr>
 *
 * @see <a href="https://simtk.org/projects/redundancy">[SimTK Project]</a>, <a
 * href="https://doi.org/10.1371/journal.pone.0209171">[Publication]</a>
 */
#include <iostream>
#include <OpenSim/OpenSim.h>
#include "Settings.h"
#include <dirent.h>

using namespace OpenSim;
using namespace std;

/**
 * \brief Performs joint reaction analysis.
 */
class JointReactionAnalysis {
public:
    struct Parameters {
        std::string subjectName;
        OpenSim::Model model;
        std::string inverseKinematicsMotion;
        std::string groundReactionXMLTemplate;
        std::string groundReactionForces;
        std::string reserveActuators;
        std::string muscleForces;
        /** specify joint of interest list or ALL keyword */
        std::vector<std::string> jointNames;
        /** specify express in frames (e.g. child or parent) */
        std::vector<std::string> applyOnBody;
        /** specify express in frames (e.g. ground, child or parent) */
        std::vector<std::string> expressInFrames;
        std::string resultsDir;
    };
    JointReactionAnalysis(Parameters& parameters) : parameters(parameters) {
    }
    void run() {
        string tempExternalLoadsXML = parameters.resultsDir +
            parameters.subjectName + "_external_loads.xml";
        // prepare external loads
        ExternalLoads externalLoads(parameters.groundReactionXMLTemplate, false);
        // update template
        externalLoads.setExternalLoadsModelKinematicsFileName(
            parameters.inverseKinematicsMotion);
        externalLoads.setDataFileName(parameters.groundReactionForces);
        externalLoads.setLowpassCutoffFrequencyForLoadKinematics(6);
        externalLoads.print(tempExternalLoadsXML);

        // // add reserve actuators
        // ForceSet forceSet(parameters.reserveActuators);
        // for (int i = 0; i < forceSet.getSize(); i++) {
        //     parameters.model.updForceSet().append(forceSet.get(i));
        // }

        // add static optimization
        Storage motion(parameters.inverseKinematicsMotion);
        auto jointReaction = new JointReaction();
        jointReaction->setName("JointReaction");
        jointReaction->setStartTime(motion.getFirstTime());
        jointReaction->setEndTime(motion.getLastTime());
        jointReaction->setForcesFileName(parameters.muscleForces);
        Array<string> jointNames, expressInFrames, applyOnBody;
        for (int i = 0; i < parameters.jointNames.size(); i++) {
            jointNames.append(parameters.jointNames[i]);
            expressInFrames.append(parameters.expressInFrames[i]);
            applyOnBody.append(parameters.applyOnBody[i]);
        }
        jointReaction->setJointNames(jointNames);
        jointReaction->setInFrame(expressInFrames);
        jointReaction->setOnBody(applyOnBody);
        parameters.model.addAnalysis(jointReaction);

        AnalyzeTool tool;
        tool.setName(parameters.subjectName);
        tool.setModel(parameters.model);
        tool.setInitialTime(motion.getFirstTime());
        tool.setFinalTime(motion.getLastTime());
        tool.setCoordinatesFileName(parameters.inverseKinematicsMotion);
        // let createExternalLoads to handle the creation of the external loads
        tool.setExternalLoadsFileName(tempExternalLoadsXML);
        // _loadModelAndInput be set so that .mot file is loaded correctly
        tool.setLoadModelAndInput(true);
        tool.setLowpassCutoffFrequency(6);
        tool.setResultsDir(parameters.resultsDir);
        tool.run();
    }
private:
    Parameters parameters;
};

void run(int argc, char *argv[]) {
    cout << "OpenSim JointReaction analysis probably has a memory leak meaning "
         << "that this batch analysis will overflow the RAM if the feasible "
         << "muscle force set is large. To avoid this problem either use the "
         << "python scripts or run this executable providing the previous "
         << "iteration as a command line argument." << endl;

    // read command line arguments
    int previousIteration;
    int maxIterations;
    if(argc == 3){
        previousIteration = atoi(argv[1]);
        maxIterations = atoi(argv[2]);
	assert(previousIteration >= 0 && maxIterations > 0);
    }
    else{
	auto message = string("this executable accepts 2 integer arguments: ") +
	    "previousIteration maxIterations";
        throw OpenSim::Exception(message);
    }

    auto subjectDir = DATA_DIR + "/gait1018/";
    auto resultsDir = subjectDir + "results/";
    auto feasibleForcesDir = resultsDir + "feasible_muscle_forces/";
    auto outputDir = resultsDir + "joint_reaction_analyses/";
    auto modelFile = subjectDir + "subject01_scaled.osim";
    auto groundReactionXMLTemplate = subjectDir + "subject01_walk1_grf.xml";
    auto groundReactionForces = subjectDir + "subject01_walk1_grf.mot";
    auto reserveActuators = subjectDir + "subject01_actuators.xml";
    auto ikFile = resultsDir + "subject01_walk1_ik.mot";

    // find all feasible muscle forces and perform joint reaction analysis
    DIR *dir;
    struct dirent *ent;
    int currentIteration = 0;
    int skip = previousIteration;
    if ((dir = opendir(feasibleForcesDir.c_str())) != NULL) {
        while ((ent = readdir(dir)) != NULL) {
	    auto file = string(ent->d_name);
	    if (file.find(".sto") == string::npos) continue; // not a sto file
	    if (skip > 0) {
		skip--;
		continue;
	    }
            // perform JRA
            Model model(modelFile);
            JointReactionAnalysis::Parameters jraParameters
		{
		 to_string(currentIteration + previousIteration) + "_subject01",
		 model,
		 ikFile,
		 groundReactionXMLTemplate,
		 groundReactionForces,
		 reserveActuators,
		 feasibleForcesDir + file,
		 {"ALL"},
		 {"child"},
		 {"ground"},
		 outputDir
		};
            JointReactionAnalysis jra(jraParameters);
            jra.run();
	    
	    // update counters
	    currentIteration++;
	    if (currentIteration == maxIterations) {
		closedir(dir);
		auto message = string("restart analysis with ") +
		    "previousIteration = " +
		    to_string(previousIteration + maxIterations);
		throw OpenSim::Exception(message);
	    }
        }
    }
    closedir(dir);
}

int main(int argc, char *argv[]) {
    try {
        run(argc, argv);
    }
    catch (exception &e) {
        cout << e.what() << endl;
        return -1;
    }
    return 0;
}
