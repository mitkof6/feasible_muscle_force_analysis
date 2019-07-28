# \brief Useful OpenSim utilities.
#
# @author Dimitar Stanev (stanev@ece.upatras.gr)
import os
import opensim


def construct_ik_task_set(model, marker_data, task_set):
    """Construct OpenSim Inverse Kinematics task set.

    In older versions of OpenSim (e.g. 3.3) IK will not execute when there are
    virtual markers that do not exist in the marker data.

    """
    virtual_markers = model.getMarkerSet()
    marker_names = marker_data.getMarkerNames()
    for i in range(0, marker_names.getSize()):
        marker_name = marker_names.get(i)
        exists = False
        for j in range(0, virtual_markers.getSize()):
            if virtual_markers.get(j).getName() == marker_name:
                task = opensim.IKMarkerTask()
                task.setName(marker_name)
                task.setApply(True)
                task.setWeight(1)
                task_set.adoptAndAppend(task)
                exists = True
                break

        if not exists:
            task = opensim.IKMarkerTask()
            task.setName(marker_name)
            task.setApply(False)
            task.setWeight(1)
            task_set.adoptAndAppend(task)


def perform_ik(model_file, trc_file, results_dir):
    """Performs Inverse Kinematics using OpenSim.

    Parameters
    ----------
    model_file: str
        OpenSim model (.osim)
    trc_file: str
        the experimentally measured marker trajectories (.trc)
    results_dir: str
        directory to store the results

    """
    model = opensim.Model(model_file)
    # model.set_assembly_accuracy(1e-3)
    model.initSystem()
    marker_data = opensim.MarkerData(trc_file)
    name = os.path.basename(trc_file)[:-4]
    ik_tool = opensim.InverseKinematicsTool()
    task_set = ik_tool.getIKTaskSet()
    construct_ik_task_set(model, marker_data, task_set)
    ik_tool.setName(name)
    ik_tool.setModel(model)
    ik_tool.setStartTime(marker_data.getStartFrameTime())
    ik_tool.setEndTime(marker_data.getLastFrameTime())
    ik_tool.setMarkerDataFileName(trc_file)
    ik_tool.setResultsDir(results_dir)
    ik_file = results_dir + name + '_ik.mot'
    ik_tool.setOutputMotionFileName(ik_file)
    ik_tool.run()
    return ik_file


def perform_so(model_file, ik_file, grf_file, grf_xml, reserve_actuators,
               results_dir):
    """Performs Static Optimization using OpenSim.

    Parameters
    ----------
    model_file: str
        OpenSim model (.osim)
    ik_file: str
        kinematics calculated from Inverse Kinematics
    grf_file: str
        the ground reaction forces
    grf_xml: str
        xml description containing how to apply the GRF forces
    reserve_actuators: str
        path to the reserve actuator .xml file
    results_dir: str
        directory to store the results
    """
    # model
    model = opensim.Model(model_file)

    # prepare external forces xml file
    name = os.path.basename(grf_file)[:-8]
    external_loads = opensim.ExternalLoads(model, grf_xml)
    external_loads.setExternalLoadsModelKinematicsFileName(ik_file)
    external_loads.setDataFileName(grf_file)
    external_loads.setLowpassCutoffFrequencyForLoadKinematics(6)
    external_loads.printToXML(results_dir + name + '.xml')

    # add reserve actuators
    force_set = opensim.ForceSet(model, reserve_actuators)
    force_set.setMemoryOwner(False)  # model will be the owner
    for i in range(0, force_set.getSize()):
        model.updForceSet().append(force_set.get(i))

    # construct static optimization
    motion = opensim.Storage(ik_file)
    static_optimization = opensim.StaticOptimization()
    static_optimization.setStartTime(motion.getFirstTime())
    static_optimization.setEndTime(motion.getLastTime())
    static_optimization.setUseModelForceSet(True)
    static_optimization.setUseMusclePhysiology(True)
    static_optimization.setActivationExponent(2)
    static_optimization.setConvergenceCriterion(0.0001)
    static_optimization.setMaxIterations(100)
    # model.addAnalysis(static_optimization)
    model.updAnalysisSet().adoptAndAppend(static_optimization)

    # analysis
    analysis = opensim.AnalyzeTool(model)
    analysis.setName(name)
    analysis.setModel(model)
    analysis.setSolveForEquilibrium(True)
    analysis.setInitialTime(motion.getFirstTime())
    analysis.setFinalTime(motion.getLastTime())
    analysis.setLowpassCutoffFrequency(6)
    analysis.setCoordinatesFileName(ik_file)
    analysis.setExternalLoadsFileName(results_dir + name + '.xml')
    analysis.setLoadModelAndInput(True)
    analysis.setResultsDir(results_dir)
    analysis.run()
    so_force_file = results_dir + name + '_StaticOptimization_force.sto'
    so_activations_file = results_dir + name + \
        '_StaticOptimization_activation.sto'
    return (so_force_file, so_activations_file)


def perform_jra(model_file, ik_file, grf_file, grf_xml, reserve_actuators,
                muscle_forces_file, results_dir, prefix='',
                joint_names=['All'],
                apply_on_bodies=['parent'],
                express_in_frame=['ground']):
    """Performs Static Optimization using OpenSim.

    Parameters
    ----------
    model_file: str
        OpenSim model (.osim)
    ik_file: str
        kinematics calculated from Inverse Kinematics
    grf_file: str
        the ground reaction forces
    grf_xml: str
        xml description containing how to apply the GRF forces
    reserve_actuators: str
        path to the reserve actuator .xml file
    muscle_forces_file: str
        path to the file containing the muscle forces from SO
    results_dir: str
        directory to store the results
    prefix: str
        prefix of the resultant joint reaction loads
    joint_names: list
        joint names of interest
    apply_on_bodies: list
        apply on child or parent
    express_in_frame: list
        frame to express results
    """
    assert(len(joint_names) == len(apply_on_bodies) == len(express_in_frame))

    # model
    model = opensim.Model(model_file)

    # prepare external forces xml file
    name = os.path.basename(grf_file)[:-8]
    external_loads = opensim.ExternalLoads(grf_xml, False)
    external_loads.setExternalLoadsModelKinematicsFileName(ik_file)
    external_loads.setDataFileName(grf_file)
    external_loads.setLowpassCutoffFrequencyForLoadKinematics(6)
    external_loads.printToXML(results_dir + name + '.xml')

    # # TODO this may not be needed
    # # add reserve actuators (must not be appended when performing JRA)
    # force_set = opensim.ForceSet(model, reserve_actuators)
    # force_set.setMemoryOwner(False)  # model will be the owner
    # for i in range(0, force_set.getSize()):
    #     model.updForceSet().append(force_set.get(i))
    #     # model.addForce(force_set.get(i))

    # construct joint reaction analysis
    motion = opensim.Storage(ik_file)
    joint_reaction = opensim.JointReaction(model)
    joint_reaction.setName('JointReaction')
    joint_reaction.setStartTime(motion.getFirstTime())
    joint_reaction.setEndTime(motion.getLastTime())
    joint_reaction.setForcesFileName(muscle_forces_file)
    joint_names_arr = opensim.ArrayStr()
    apply_on_bodies_arr = opensim.ArrayStr()
    express_in_frame_arr = opensim.ArrayStr()
    for j, b, f in zip(joint_names, apply_on_bodies, express_in_frame):
        joint_names_arr.append(j)
        apply_on_bodies_arr.append(b)
        express_in_frame_arr.append(f)

    joint_reaction.setJointNames(joint_names_arr)
    joint_reaction.setOnBody(apply_on_bodies_arr)
    joint_reaction.setInFrame(express_in_frame_arr)
    # model.addAnalysis(joint_reaction)
    model.updAnalysisSet().adoptAndAppend(joint_reaction)
    model.initSystem()

    # analysis
    analysis = opensim.AnalyzeTool(model)
    analysis.setName(prefix + name)
    analysis.setModel(model)
    analysis.setSolveForEquilibrium(True)
    analysis.setModelFilename(model_file)
    analysis.setInitialTime(motion.getFirstTime())
    analysis.setFinalTime(motion.getLastTime())
    analysis.setLowpassCutoffFrequency(6)
    analysis.setCoordinatesFileName(ik_file)
    analysis.setExternalLoadsFileName(results_dir + name + '.xml')
    analysis.setLoadModelAndInput(True)
    analysis.setResultsDir(results_dir)
    analysis.run()
    jra_file = results_dir + name + '_JointReaction_ReactionLoads.sto'
    return jra_file
