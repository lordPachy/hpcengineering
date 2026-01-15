# Login

sudo gpclient connect gp-dmat-saml.vpn.polimi.it
ssh u11172853@10.78.18.100

# General information

With lsblk, we can get information about the partitions available.\
With chmod -R, we can change the permission for r/w of the files inside our directory.\

There is no backup. The maximum amount of memory allowed is 10 GB, but this is a soft limit (the hard limit is 12.5 GB).\

We can enter our private section only if in our folder there is .ssh folder with our private keys (setup this).\

Work directory have 100 GB soft limit, and ~110 GB of hard limit. There is no backup for it, and it is a single drive.\

Home, software and work folders are actually the same on the whole cluster\

With `ip | grep UP`  we can retrieve some information about the path we are working in.

# Information about the filesystem

Files on the scratch folder are deleted after 40 days from the last edit. Information about this can be found on cat /etc/fstab.\

If speed is not a concern, we can work on our own work directory.
If it is, we should use scratch_local, possibly in a new dir.
If a file is too big, we might want to use scratch_global.
LOOK THIS UP\

/software/ is shared as well. There are containers, and in the /containers/mk subfolder there is a little guide on how to use them.

# Transferring files

scp -r u11172853@10.78.18.100:/home/u11172853/<filename> .

# Resource manager

The scheduler is based on queue. The higher the submits, the lower the priority on it.

qsub -I -l host=gpu01
can be used to query a specific device.

# Submitting jobs
We should use the login node only for login. All the other jobs should be submitted.

We can use the files in the PBS folders as examples of job submissions.\
We can add an option select=7 and then place=scatter.

SOLVE $HOME PROBLEM

The first file checks the presence of the required libraries.\
By typing `pbsnodes -aSj` we can retrieve the number of free and total cpu/gpus.\
With `qsub -I` we can retrieve the state of our currently running jobs.\
Some terminology:
 - ncpus : threads (shared memory)
 - mpiprocs : processes (distributed memory)
 - host: choose a specific host

 With ls /software/ -l we can query the installed software. The command in the third line sets up Spack. In Spack, we need to load the necessary modules.\
 Note that the modules are being loaded on the requested cpu, not on the login node.

```
. /etc/profile.d/pbs.sh
qsub -I -q cpu -l select=1:ncpus=6:mpiprocs=2:host=cpu05
source /software/spack/share/spack/setup-env.sh
spack load gcc@15.2.0
spack load intel-oneapi-tbb@2022.3.0
spack load openmpi@5.0.8
which gcc

# Containers

The container is stored in /storage/containers/nvhcc/container.sif. We could use apptainer shell with the flag -nv.\

We can use qsub to access the contianer, in order to send jobs without staying connected.