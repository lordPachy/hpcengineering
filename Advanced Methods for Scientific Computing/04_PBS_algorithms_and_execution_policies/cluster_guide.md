# General information

With lsblk, we can get information about the partitions available.\
With chmod -R, we can change the permission for r/w of the files inside our directory.\

There is no backup. The maximum amount of memory allowed is 10 GB, but this is a soft limit (the hard limit is 12.5 GB).\

We can enter our private section only if in our folder there is .ssh folder with our private keys (setup this).\

Work directory have 100 GB soft limit, and ~110 GB of hard limit. There is no backup for it, and it is a single drive.\

Home, software and work folders are actually the same on the whole cluster

# Information about the filesystem

Files on the scratch folder are deleted after 40 days from the last edit. Information about this can be found on cat /etc/fstab.\

If speed is not a concern, we can work on our own work directory.
If it is, we should use scratch_local, possibly in a new dir.
If a file is too big, we might want to use scratch_global.
LOOK THIS UP\

/software/ is shared as well. There are containers, and in the /containers/mk subfolder there is a little guide on how to use them.

# Transferring files

scp -r u11172853@10.78.18.100:/home/u1172853/<filename> .

# Resource manager

The scheduler is based on queue. The higher the submits, the lower the priority on it.

qsub -I -l host=gpu01
can be used to query a specific device.

# Submitting jobs
We can use the files in the PBS folders as examples of job submissions.\
We can add an option select=7 and then place=scatter.

SOLVE $HOME PROBLEM

# Containers

The container is stored in /storage/containers/nvhcc/container.sif. We could use apptainer shell with the flag -nv.\

We can use qsub to access the contianer, in order to send jobs without staying connected.