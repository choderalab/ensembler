# Create compressed archive for models
# Daniel L. Parton <partond@mskcc.org> - Oct 25 2013

import sys
from mpi4py import MPI
from ast import literal_eval

targets = literal_eval( sys.argv[ sys.argv.index('-targets') + 1 ] )

# get communicator
comm = MPI.COMM_WORLD

# get rank and size
rank = comm.rank
size = comm.size

for target_index in range(rank, len(targets), size):
    import subprocess, os
    target = targets[target_index]
    models_path = os.path.join('models', target)
    tar_path = os.path.join('models', target.lower() + '.tgz')

    subprocess.call(['tar', 'zcf', tar_path, project_path])

