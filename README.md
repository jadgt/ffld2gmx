# ffld2gmx
Python Script to transform Force Field files from Schrödinger's Maestro (.ffld) to GROMACS (.itp, .top)

### Usage

To use this script one just need to have the `ffld2gmx.py` script in the same directory as your `topology.ffld`. The `.ffld` file will be used as argument for the script:
```
python ffld2gmx.py topology.ffld
```

The output of the script are two files:
1. A residue topology file called `UNK.itp` so one can change the name to the residue name desired.
2. A topology file called `topol.top` that serve as template for your simulation topology file. The residue topology `UNK.itp` is automatically linked with the directive: `#include "UNK.itp"`

### Test file

To start checking the usage one can use the `test.ffld` file included in this repository.

### Final disclaimer:

- The automatic force field generated is a OPLS-AA force field. For other force fields one need to use other tools to transform those.
- Transition metals force constants are not reliable and needs to be computed independently.
- At this moment I am working on improvements, right now the masses of the atoms needs to be manually added. Future fixes will be available soon.