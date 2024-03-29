# This script processes the joint reaction loads and computes the
# minimum and maximum bounds for a given joint of interest. The
# subject directory, the folder containing the joint reaction results
# form perform_joint_reaction_batch.py, the joint of interest and the
# mass of the subject must be provided.
#
# @author Dimitar Stanev (stanev@ece.upatras.gr)
import os
import numpy as np
from tqdm import tqdm
import matplotlib.pyplot as plt
plt.rcParams['font.size'] = 9

###############################################################################
# parameters

subject_dir = os.getcwd() + '/../data/gait1018/'
os_jra_file = subject_dir + \
    'results/subject01_JointReaction_ReactionLoads.sto'
jra_results_dir = subject_dir + 'results/joint_reaction_analyses/'
figures_dir = subject_dir + 'results/fig/'

collect = True
joint = 'knee_r'
joints = 3
mass = 72.6  # kg
g = 9.8  # m/s^2
body_weight = mass * g

if not os.path.isfile(os_jra_file):
    raise RuntimeError('required files do not exist')

if not (os.path.isdir(jra_results_dir) and
        os.path.isdir(figures_dir)):
    raise RuntimeError('required folders do not exist')


###############################################################################
# utilities

def readMotionFile(filename):
    """Reads OpenSim .sto files.

    Parameters
    ----------
    filename: str
        absolute path to the .sto file

    Returns
    -------
    header: list of str
        the header of the .sto
    labels: list of str
        the labels of the columns
    data: list of lists
        an array of the data

    """

    if not os.path.exists(filename):
        print('file do not exists')

    file_id = open(filename, 'r')

    # read header
    next_line = file_id.readline()
    header = [next_line]
    nc = 0
    nr = 0
    while 'endheader' not in next_line:
        if 'datacolumns' in next_line:
            nc = int(next_line[next_line.index(' ') + 1:len(next_line)])
        elif 'datarows' in next_line:
            nr = int(next_line[next_line.index(' ') + 1:len(next_line)])
        elif 'nColumns' in next_line:
            nc = int(next_line[next_line.index('=') + 1:len(next_line)])
        elif 'nRows' in next_line:
            nr = int(next_line[next_line.index('=') + 1:len(next_line)])

        next_line = file_id.readline()
        header.append(next_line)

    # process column labels
    next_line = file_id.readline()
    if next_line.isspace() is True:
        next_line = file_id.readline()

    labels = next_line.split()

    # get data
    data = []
    for i in range(1, nr + 1):
        d = [float(x) for x in file_id.readline().split()]
        data.append(d)

    file_id.close()

    return header, labels, data


def index_containing_substring(list_str, pattern):
    """For a given list of strings finds the index of the element that
    contains the substring.

    Parameters
    ----------
    list_str: list of str

    pattern: str
         pattern


    Returns
    -------
    indices: list of int
         the indices where the pattern matches

    """
    indices = []
    for i, s in enumerate(list_str):
        if pattern in s:
            indices.append(i)

    return indices


###############################################################################
# main

# OpenSim's JRA results
os_header, os_labels, os_data = readMotionFile(os_jra_file)
os_data = np.array(os_data)
joint_index = index_containing_substring(os_labels, joint)
assert(joint_index != -1)

# get all files in the directory
jra_files = os.listdir(jra_results_dir)
# remove files that are not joint reactions
jra_files = [e for e in jra_files if 'ReactionLoads' in e]

if collect:
    # collect simulation data
    print('Processing joint reaction analyses ...')
    # allocate the necessary space
    solutions_to_keep = len(jra_files)
    simulationData = np.empty([solutions_to_keep,
                               os_data.shape[0],
                               os_data.shape[1]],
                              dtype=float)
    # collect data
    for i, f in enumerate(tqdm(jra_files)):
        if i == solutions_to_keep:
            break

        header, labels, data = readMotionFile(jra_results_dir + f)
        simulationData[i, :, :] = np.array(data)


heel_strike_right = [0.65, 1.85]
# heel_strike_right = [0.45, 1.7]
toe_off_right = [0.15, 1.4]
heel_strike_left = [0.0, 1.25]
toe_off_left = [0.8, 2]
# heel_strike_right = [0.65]
# toe_off_right = [1.4]
# heel_strike_left = [1.25]
# toe_off_left = [0.8]
if '_l' in joint:
    heel_strike = heel_strike_left
    toe_off = toe_off_left
else:
    heel_strike = heel_strike_right
    toe_off = toe_off_right

# plot data min/max reactions vs OpenSim JRA
fig, ax = plt.subplots(nrows=1, ncols=joints, figsize=(12, 4))
for i in range(0, joints):
    # plot feasible reaction loads
    min_reaction = np.min(
        simulationData[1:, 1:, joint_index[i]] / body_weight, axis=0)
    max_reaction = np.max(
        simulationData[1:, 1:, joint_index[i]] / body_weight, axis=0)
    ax[i].fill_between(os_data[1:, 0], min_reaction, max_reaction, color='b',
                       alpha=0.2, label='Feasible Reactions')
    # plot OpenSim reaction loads
    ax[i].plot(os_data[1:, 0], os_data[1:, joint_index[i]] / body_weight,
               '-.r', label='OpenSim JRA')
    # annotate the heel strike and toe off regions
    min_min = np.min(min_reaction)
    max_max = np.max(max_reaction)
    ax[i].vlines(x=heel_strike, ymin=min_min, ymax=max_max,
                 color='c', linestyle='--', label='HS')
    ax[i].vlines(x=toe_off, ymin=min_min, ymax=max_max,
                 color='m', linestyle=':', label='TO')
    # figure settings
    ax[i].set_title(os_labels[joint_index[i]])
    # ax[i].set_xlim([heel_strike[0], heel_strike[1]])
    ax[i].set_xlabel('time (s)')
    ax[i].set_ylabel('reaction / body weight')
    if i == joints - 1:
        ax[i].legend()


fig.tight_layout()
fig.savefig(figures_dir + joint + '.pdf',
            format='pdf', dpi=300)
fig.savefig(figures_dir + joint + '.png',
            format='png', dpi=300)
# fig.savefig(figures_dir + joint + '.tif',
#             format='tif', dpi=300)
