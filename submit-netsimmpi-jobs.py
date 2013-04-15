# this is a Python 3.2 program
#
# submit-netsimmpi-jobs.py
#
#

import argparse
import logging
import os.path
import random
import re
import subprocess
import sys
import time

# this configuration only logs to the file
# I wish logging output were going to console and file
logging.basicConfig(
  stream=sys.stdout,
  #filename='log.txt',
  #filemode='w',
  format='%(name)s %(levelname)s: %(message)s',
  level=logging.DEBUG
)

log = logging.getLogger()

def submit(argv):
  cmdLineArgParser = argparse.ArgumentParser()
  cmdLineArgParser.add_argument('jobPrefix')
  cmdLineArgParser.add_argument('outputParentDirPath')
  cmdLineArgParser.add_argument('simulationTimeSeconds')

  torqueOptionGroup = cmdLineArgParser.add_argument_group('Torque options', 'Options for the Torque job scheduler.')
  torqueOptionGroup.add_argument(
    "--wt", "--walltime", dest="walltime", default="16:00:00", help="maximum walltime"
  )
  torqueOptionGroup.add_argument(
    "--ng", "--node-group", dest="nodeGroup", default="all", help="named compute node group"
  )
  torqueOptionGroup.add_argument(
    "--ppn", "--processors-per-node", type=int, dest="ppn", default=8, help="processors per node"
  )
  torqueOptionGroup.add_argument(
    "--ns", "--no-submit", dest="submitJobScripts", action="store_false", default=True, help="do not submit jobs to TORQUE"
  )
  torqueOptionGroup.add_argument(
    "--nd", "--no-delete", dest="deleteJobScripts", action="store_false", default=True, help="do not delete job scripts"
  )

  cmdLineArguments = cmdLineArgParser.parse_args()

  log.info('options: {0}'.format(cmdLineArguments))

  class Object:
    pass

  parameters = Object()
  parameters.jobPrefix = cmdLineArguments.jobPrefix

  # protect me from myself by creating output directories if they do not exist
  parameters.resultsDirPath = os.path.join(cmdLineArguments.outputParentDirPath, cmdLineArguments.jobPrefix, 'results')
  log.info('creating results directory: {0.resultsDirPath}'.format(parameters))
  subprocess.check_call(['mkdir', '-p', parameters.resultsDirPath])

  parameters.stdeoDirPath = os.path.join(cmdLineArguments.outputParentDirPath, cmdLineArguments.jobPrefix, 'stdeo')
  log.info('creating log directory {0.stdeoDirPath}'.format(parameters))
  subprocess.check_call(['mkdir', '-p', parameters.stdeoDirPath])

  jobCount = 0
  submitCount = 0
  submittedJobList = []
  for nodeCount in (1,2):
    for np in (1,2,3,4):
      for m in range(10):
        jobCount = jobCount + 1

        parameters.nodeCount = nodeCount
        parameters.np = np
        parameters.m = m
        parameters.jobName = '{0.jobPrefix}-nc{0.nodeCount}-np{0.np}-job{0.m}'.format(parameters)
        log.info('jobName: {0.jobName}'.format(parameters))

        parameters.stdOutputFilePath = os.path.join(parameters.stdeoDirPath, '{0.jobName}.stdout'.format(parameters))
        parameters.stdErrorFilePath = os.path.join(parameters.stdeoDirPath, '{0.jobName}.stderr'.format(parameters))

        parameters.ppn = np * 2

        parameters.simulationSeed = random.randint(1,1000000)

        ############################################################################
        # write the job script
        #
        jobFilePath = '{0.jobName}.sh'.format(parameters)
        with open(jobFilePath, 'wt') as jobFile:
          jobFile.write(
  """\
  #!/bin/bash
  #PBS -V
  #PBS -l nodes={parameters.nodeCount}:ppn={parameters.ppn}:{args.nodeGroup}
  #PBS -l walltime={args.walltime}
  #PBS -q default
  #PBS -o {parameters.stdOutputFilePath}
  #PBS -e {parameters.stdErrorFilePath}

  mpirun -np {parameters.np} ~/dev-local/netsimmpi/netsimmpi.out {args.simulationTimeSeconds} {parameters.simulationSeed} {parameters.resultsDirPath}
  """.format(args=cmdLineArguments, parameters=parameters))


        if cmdLineArguments.submitJobScripts:
          log.info('submitting {0}'.format(jobFilePath))
          subprocess.check_call(['qsub', jobFilePath])
          submitCount = submitCount + 1
          submittedJobList.append(jobFilePath)

        if cmdLineArguments.deleteJobScripts:
          log.info('deleting {0}'.format(jobFilePath))
          subprocess.check_call(['rm', '-f', jobFilePath])

  log.info('wrote {0} job scripts'.format(jobCount))
  log.info('submitted {0} job script(s):'.format(submitCount))
  log.info('\n'.join(submittedJobList))


if __name__ == "__main__":
  sys.exit(submit(sys.argv)) 
