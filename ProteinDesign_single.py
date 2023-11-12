import rosetta.core.pack.task
from rosetta.protocols.scoring import Interface
from toolbox import mutate_residue
from rosetta import *
import os
import glob
init()

pdb_file = 'PhaC1437_phenyllactyl-coa.pdb'
pose = pose_from_pdb(pdb_file)

partners = 'A_X'
dock_jobs = 100

# The range to be mutated near the active site of PhaC was designated in advance.
core_mut_range = [(36,58),(72,99),(112,116),(118,271),(273,299),(301,317)]

def make_resfile():
    try:
        os.mkdir('./resfile')
    except:
        pass
    core_res_num = []
    for res_range in core_mut_range:
        start = res_range[0]
        end = res_range[1]+1
        for num in range(start, end):
            core_res_num.append(str(num))
    res_dict = {}
    with open(pdb_file, 'r') as fr:
        lines = fr.read().splitlines()
        for li in lines:
            try:
                keywd = li[0:4].strip()
                if keywd == 'ATOM':
                    res_name = li[17:20].strip()
                    res_num = li[22:26].strip()
                    res_dict[res_num] = res_name
            except:
                pass
    print(core_res_num)
    print(len(core_res_num))
    for i in core_res_num:
        mut_res = res_dict[i]
        with open('./resfile/resfile_'+i, 'w') as fw:
            fw.write('NATRO\n')
            fw.write('USE_INPUT_SC\n')
            fw.write('start\n')
            fw.write('%5s%3s%7s\n' %(i, 'A', 'ALLAA'))

def ProtDesign(pose):
    test_pose = Pose()
    test_pose.assign(pose)

    dock_jump = 1
    movable_jumps = Vector1([dock_jump])
    setup_foldtree(test_pose, partners, movable_jumps)

    scorefxn = get_fa_scorefxn()
    scorefxn(test_pose)

    try:
        os.mkdir('./rosetta_results_wt')
    except:
        pass

    count = 0
    res_files = glob.glob('./resfile/*')
    for res_file in res_files:
        # setup the design PackerTask, use the generated resfile
        test_pose.assign(pose)
        pose_design = standard_packer_task(test_pose)
        rosetta.core.pack.task.parse_resfile(test_pose, pose_design, res_file)

        count += 1
        # perform design
        designmover = PackRotamersMover(scorefxn, pose_design)
        #print '\nPre-design score: ',scorefxn(test_pose)
        #print 'Pre-design sequence: ...' + test_pose.sequence() + '...'

        designmover.apply(test_pose)
        #print '\nPost-design score: ', scorefxn(test_pose)
        #print 'Post-design sequence: ...' + test_pose.sequence() + '...'

        test_pose.dump_pdb('./rosetta_results_wt/designed.pdb')
        pdb_filename = './rosetta_results_wt/designed.pdb'

        res_name = os.path.basename(res_file).split('.')[0]
        dock_output = './rosetta_results_wt/dock_output_wildtype_%s_%s' %(res_name,count)
        ligand_interface(pdb_filename, partners, dock_jobs, dock_output, count)

def ligand_interface(pdb_filename, partners, dock_jobs, dock_output, mut_cnt):

    # 1. creates a pose from the PDB file and introduce ligand parameters
    pose = Pose()
    pose_from_file(pose, pdb_filename)

    # 2. setup the docking FoldTree
    dock_jump = 1
    setup_foldtree(pose, partners, Vector1([dock_jump]))

    # 3. create a copy of the pose for testing
    test_pose = Pose()
    test_pose.assign(pose)

    # 4. create ScoreFunctions for centroid and fullatom docking
    scorefxn = create_score_function('ligand')    ### ligand_soft   scoring functions

    # 5. setup the high resolution (fullatom) docking protocol (DockMCMProtocol)
    docking = DockMCMProtocol()
    docking.set_scorefxn(scorefxn)

    # 6. setup the PyJobDistributor
    jdc = PyJobDistributor(dock_output, dock_jobs, scorefxn)

    # 7. perform protein-ligand docking
    while not jdc.job_complete:
        docking.apply(test_pose)
        jdc.output_decoy(test_pose)


####################
##Perform simulation
####################
#make_resfile()
ProtDesign(pose)

