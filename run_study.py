#!/usr/bin/env python
#SBATCH -J EICSTUDY    # job name
#SBATCH -o logs/eicstudy-%A_%a.out
#SBATCH -e logs/eicstudy-%A_%a.out
#SBATCH -n 1
#SBATCH -c 1
#SBATCH --mem=6G
#SBATCH -p htc             # queue (partition) -- batch, parallel, etc.  parallel-medium
#SBATCH -t 12:00:00        # run time (hh:mm:ss)
#SBATCH -D .               # Directory where executable will be run
#SBATCH --mail-user=ssekula@smu.edu
#SBATCH --mail-type=fail    # email me when the job finishes


# This script was written originally for use on a SLURM-based batch system,
# like the one available at SMU ("ManeFrame II"). It can be run on the command-line
# alone; if the script doesn't detect the requested SLURM environment variables,
# it will ask you to specify them. For instance, to run the first variation in a
# study,
#
# SLURM_ARRAY_TASK_ID=0 ./run_study.py -t ./bfield_study_template.tcl -p "['PARAM_BFIELD']" -v "[[1.0,1.5,2.0,2.5,3.0,3.5]]" -n bfield
#
# The above will manually select the first variation in the array passed to "-v" and run it.
#
# General idea:
#
# A template script is executed by the Delphes simulation executable. Before execution,
# placeholder parameters are filled in by the script and a copy of the template is 
# generated (with the parameter values set) for actual execution. For example:
#
# SLURM_ARRAY_TASK_ID=0 ./run_study.py -t ./bfield_study_template.tcl -p "['PARAM_BFIELD']" -v "[[2.0]]" -n bfield
# 
# will open ./bfield_study_template.tcl, find all instanced of PARAM_BFIELD, and replace them with the number "2.0".
# This allows command-line variation of one or more parameters to perform a study of changes to simulation.
#

import subprocess
import math
import os
import sys
import shutil
import glob
import re
import ast

import argparse
parser = argparse.ArgumentParser()

parser.add_argument("-n", "--name", type=str,
                    help="name for this study set (e.g. bfield)")
parser.add_argument("-t", "--template", type=str,
                    help="template TCL file for study")
parser.add_argument("-c", "--commands", type=str,
                    help="template command file for study")
parser.add_argument("-p", "--params", type=ast.literal_eval, default=[],
                    help="environment variables to set for study")
parser.add_argument("-v", "--values", type=ast.literal_eval, default=[],
                    help="parameter values for study")
parser.add_argument("-f", "--force", default=False, action='store_true',
                    help="force-overwrite existing output")

args = parser.parse_args()



# Create the task superdirectory

if not os.path.exists(args.name):
    try:
        os.makedirs(args.name)
    except os.OSError:
        print("%s already exists... continuing..." % (args.name))


SLURM_ARRAY_TASK_ID="0"

try:
    SLURM_ARRAY_TASK_ID=os.environ["SLURM_ARRAY_TASK_ID"]
except:
    print("Please set the SLURM_ARRAY_TASK_ID environment variable to a number (e.g. 0) before running this script.")
    sys.exit()
    pass

print "Task ID requested: %d" % (int(SLURM_ARRAY_TASK_ID))



value_index = int(SLURM_ARRAY_TASK_ID)

# Handle random number seed issues
random_seed = 0
if "RANDOM_SEED" not in args.params:
    args.params.append("RANDOM_SEED")
    random_seed = abs(hash(args.name)) % (10 ** 8) + value_index
    args.values.append(random_seed)


# Execute the study

taskdir="%s/%d" % (args.name, value_index)
if os.path.exists(taskdir) and not args.force:
    print("Skipping this task directory --- it already exists. Cleanup before overwriting!")
    print(taskdir)
else:
    if not os.path.exists(taskdir):
        os.makedirs(taskdir)

    # Replace parameter placeholders with values for this study
    template_files = [args.template, args.commands]

    for template_file in template_files:
        with open(template_file, "rt") as input_template:
            with open("%s/%s" % (taskdir, template_file), "wt") as output_template:
                for line in input_template:
                    for iparam in range(len(args.params)):
                        param = args.params[iparam]
                        values = args.values[iparam]
                        value = None
                        if type(values) == list:
                            value = values[value_index]
                        else:
                            value = values
                        line = line.replace(param, str(value))
                
                    output_template.write(line)

    # Copy files to the task directory before execution
    copy_files = ["delphes_card_EIC.tcl"]
    for a_file in copy_files:
        subprocess.call("cp %s %s" % (a_file, taskdir), shell=True)
    # Write the random number seed to disk
    rndm_file = open(taskdir+"/random_seed.dat", "w")
    rndm_file.write(str(random_seed))
    rndm_file.close()

    # Execute the study
    subprocess.call("DelphesPythia8 {0[taskdir]}/{0[template]} {0[taskdir]}/{0[commands]} {0[taskdir]}/out.root".format({'taskdir': taskdir, 'template': args.template, 'commands': args.commands}), shell=True)
