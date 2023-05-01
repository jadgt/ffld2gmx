# ffld2gmx
Python Script to transform Force Field files from Schrödinger's Maestro (.ffld) to GROMACS (.itp, .top)

### Usage

To use this script one just need to have the `ffld2gmx.py` script in the same directory as your `topology.ffld`. The `.ffld` file will be used as argument for the script:
```
python ffld2gmx.py topology.ffld
```
**Important:** The script also requires to have the file `ffnotbonded.itp` file in the same directory so it can read the non-bonded parameters from the OPLS-AA force field tailored to GROMACS.

The output of the script are two files:
1. A residue topology file called `UNK.itp` so one can change the name to the residue name desired.
2. A topology file called `topol.top` that serve as template for your simulation topology file. The residue topology `UNK.itp` is automatically linked with the directive: `#include "UNK.itp"`

### Test file

To start checking the usage one can use the `test.ffld` file included in this repository.

### Final disclaimer:

- The automatic force field generated is a OPLS-AA force field. For other force fields one need to use other tools to transform those.
- Transition metals force constants are not reliable and needs to be computed independently.