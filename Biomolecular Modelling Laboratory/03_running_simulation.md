# VMD coordinates

With `measure minmax $all`, we can measure the minimum and maximum coordinates of the
molecules contained in `all`. Instead, with `measure center $all`, we measure the center.

# Configuration file
In `../common/sample.conf`, we can find a configuration file. It is the template file of a configuration
for a NAMD simulation, with all the input parameters. Refer to the comments inside the file for guidance.

# Force field file
In `par_all27_prot_lipid.inp`, all the force fields are defined. For example, the bonds are defined with the spring
model parameters. In the third column we can indeed find the elastic constant, and in the fourth the reference length.
This is done for all the possible pairs of atom types.\
On top of this, also angles can be found, again modeled with a spring model; atom triplets are considered in this case.\
The file needs to be brought in the same folder as sample.conf.

# Running the simulation
With `namd3 +p2 sample.conf` we can run the simulation. The parameter `+pn` defines the number of cores for the simulation.\
In the `.dcd` file, which is binary, we can find the collection of the frames of the simulation. By entering
`vmd ionized.psf myoutput.dcd` we can see the output. We can play the dynamic simulation with the play button.\
The periodic box can be visualized with display -> ortographic to avoid distortiong, then in the periodic tab we can see 
what the software is seeing.\
We can see that the molecule is vibrating around the equilibrium position.

# Analysis
In the analysis menu, we can find the RMSD trajectory tool. We can select to do it only on the backbone, for example.
In order to plot it, we should select ALIGN and plot before clicking on RMSD.\
Another analysis that can be done is that on hydrogen bonds. There are two options:
1. One selection: hydrogen bonds inside this selection (called "protein" in our example)
2. Two selections: considering two non-overlapping regions, we will see the bonds between the two.
We can calculate details for all bonds with the dedicated options, in order to have more information on each on them, not just the number. Recall that an hydrogen bond is considered when it falls inside certain distance and angle parameters, that can be also changed.\
In the plugins -> analysis menu, we can find the timeline. It can calculate other properties, such as the secundary structure. In this case, VMD will need to re-run a calculation.
