## FFLD to itp converter
## Created by Juan de Gracia

import sys
import os

itp_file = 'ffnotbonded.itp'

if len(sys.argv) != 2:
    print(f"Usage: {sys.argv[0]} input_file.ffld")
    sys.exit(1)
input_filename = sys.argv[1]

def read_file(filename):
    with open(filename, 'r') as f:
        file_contents = f.read()
    return file_contents

def parse_ffld(ffld_contents):
    sections = {}
    current_section = None
    for line in ffld_contents.split('\n'):
        line = line.strip()
        if not line:
            continue
        elif line.startswith('atom'):
            current_section = 'atoms'
            sections[current_section] = line + '\n'
        elif line.startswith('Stretch'):
            current_section = 'stretch'
            sections[current_section] = line + '\n'
        elif line.startswith('Bending'):
            current_section = 'bending'
            sections[current_section] = line + '\n'
        elif line.startswith('proper Torsion'):
            current_section = 'torsion'
            sections[current_section] = line + '\n'
        elif line.startswith('improper Torsion'):
            current_section = 'improper'
            sections[current_section] = line + '\n'
        else:
            if current_section is not None:
                sections[current_section] += line + '\n'
    
    # Check if the 'improper Torsion' section is present
    if 'improper' not in sections:
        sections['improper'] = ''
    return sections

## Read and parse 

contents = read_file(input_filename) 
sections = parse_ffld(contents)

# Construct atomtype directive

atoms = sections['atoms']

def transform_atoms_section(atoms_section, itp_file):
    # Read the atomtype information from the itp file
    atomtype_info = {}
    with open(itp_file, 'r') as f:
        for line in f:
            if line.startswith(' opls_'):
                parts = line.strip().split()
                atomtype_number = int(''.join(filter(str.isdigit, parts[0])))  # Extract the number from the atomtype name
                atomtype_info[atomtype_number] = {
                    'mass': float(parts[3]),
                    'charge': float(parts[4]),
                    'sigma': float(parts[6]),
                    'epsilon': float(parts[7])
                }

    new_lines = []
    lines = atoms_section.strip().split("\n")[2:-1]
    for i, line in enumerate(lines):
        parts = line.split()
        atom = f"op_unk.{parts[0]}"
        atomtype = int(parts[1])
        if atomtype in atomtype_info:
            atomtype_data = atomtype_info[atomtype]
            mass = atomtype_data['mass']
            charge = atomtype_data['charge']
            sigma = atomtype_data['sigma']
            epsilon = atomtype_data['epsilon']
        else:
            mass = 0  # Define mass manually
            charge = float(parts[4])
            sigma = float(parts[5]) / 10  # Convert sigma from Angstroms to nanometers
            epsilon = float(parts[6]) * 4.18  # Convert epsilon from kcal/mol to kJ/mol
        new_line = f"   {atom:<10s} {atom:<10s} {mass:>11.6f} {charge:>11.3f}   A   {sigma:>14.6e} {epsilon:>14.6e}"
        new_lines.append(new_line)
    return "\n".join(new_lines)

def generate_atomtypes_section(atoms_section, itp_file):
    # Transform the atoms section and iterate over each line
    atomtype_lines= transform_atoms_section(atoms_section, itp_file).split('\n')

    # Add the section header and column headers to the list for [ atomtypes ]
    atomtype_lines.insert(0, "[ atomtypes ]")
    atomtype_lines.insert(1, "; name mass charge ptype sigma(nm) epsilon (kJ/mol)")

    # Join the [ atomtypes ] lines into a single string with newline characters between each line
    atomtypes_section = "\n".join(atomtype_lines)

    # Return the atomtypes section string
    return atomtypes_section

atomtypes_section = generate_atomtypes_section(atoms, itp_file)

## Generate atom section

def transform_atoms_section2(atoms_section, itp_file):
    # Read the atomtype information from the itp file
    atomtype_info = {}
    with open(itp_file, 'r') as f:
        for line in f:
            if line.startswith(' opls_'):
                parts = line.strip().split()
                atomtype_number = int(''.join(filter(str.isdigit, parts[0])))  # Extract the number from the atomtype name
                atomtype_info[atomtype_number] = {
                    'mass': float(parts[3]),
                    'charge': float(parts[4]),
                    'sigma': float(parts[6]),
                    'epsilon': float(parts[7])
                }

    new_lines = []
    lines = atoms_section.strip().split("\n")[2:-1]
    for i, line in enumerate(lines):
        parts = line.split()
        atom = f"op_unk.{parts[0]}"
        atomtype = int(parts[1])
        if atomtype in atomtype_info:
            atomtype_data = atomtype_info[atomtype]
            mass = atomtype_data['mass']
            charge = atomtype_data['charge']
        else:
            mass = 0  # Define mass manually
            charge = float(parts[4])
        new_line = "    {:3d}   {:<10s}{:>5d} {:<5s} {:<6s} {:>3d}  {:>10.6f}  {:>11.6f}".format(i+1, atom, 1, "UNK", parts[0], 1, float(charge), float(mass))
        new_lines.append(new_line)
    return "\n".join(new_lines)

 
def generate_atoms_section(atoms_section, itp_file):
    # Transform the atoms section and iterate over each line
    atom_lines= transform_atoms_section2(atoms_section, itp_file).split('\n')

    # Add the section header and column headers to the list for [ atomtypes ]
    atom_lines.insert(0, "[ atoms ]")
    atom_lines.insert(1, "; nr type  resnr residue  atom  cgnr charge  mass")

    # Join the [ atomtypes ] lines into a single string with newline characters between each line
    atoms_section = "\n".join(atom_lines)

    # Return the atomtypes section string
    return atoms_section

atoms_section = generate_atoms_section(atoms, itp_file)

## Generate bond section
# Tidy ffld
def transform_stretch_section(stretch_section):
    new_lines = []
    lines = stretch_section.strip().split("\n")[1:]
    for i, line in enumerate(lines):
        parts = line.split()
        atom1 = int(''.join(filter(str.isdigit, parts[0])))
        atom2 = int(''.join(filter(str.isdigit, parts[1])))
        k = parts[2]
        r0 = parts[3]
        new_line = f"{atom1:6d} {atom2:6d} {k:12s} {r0:12s}"
        new_lines.append(new_line)
    return "\n".join(new_lines)

stretch = transform_stretch_section(sections['stretch'])

## Generate bonds section

def generate_bonds_section(stretch):
    # Initialize an empty list to store the bond lines
    bond_lines = []

    # Iterate over each line in the stretch string
    for line in stretch.strip().split("\n"):
        # Extract the atom indices and equilibrium distance from the line
        parts = line.split()
        i = int(parts[0]) - 1  # Subtract 1 to convert from 1-indexed to 0-indexed
        j = int(parts[1]) - 1
        r0 = float(parts[3]) / 10  # Convert from angstroms to nanometers

        # Calculate the force constant from the CHARMM form of the potential energy function
        # k = 0.5 * fc * (r - r0)^2
        fc = float(parts[2])  # Convert from kcal/mol A^2 to kJ/mol nm^2
        k = fc * 4.184 * 200

        # Format the bond line and append it to the list
        bond_line = f"{i+1:6d} {j+1:6d} 1 {r0:.6f} {k:.3f}"
        bond_lines.append(bond_line)

    # Join the bond lines into a single string with newline characters between each line
    bonds_section = "[ bonds ]\n; ai    aj    type     r0 (nm)   fc (kJ/(mol nm2))\n" + "\n".join(bond_lines)

    # Return the bonds section string
    return bonds_section

bonds_section = generate_bonds_section(stretch)

##Bending section
# Tidy up function 

angles = sections['bending']

def transform_bending_section(angles_section):
    new_lines = []
    # Split the section into lines and skip the header
    lines = angles_section.strip().split("\n")[1:]
    for line in lines:
        # Extract the atom indices and equilibrium angle from the line
        parts = line.split()
        atom1 = int(''.join(filter(str.isdigit, parts[0])))
        atom2 = int(''.join(filter(str.isdigit, parts[1])))
        atom3 = int(''.join(filter(str.isdigit, parts[2])))
        k = parts[3]
        theta0 = parts[4]

        # Format the new line
        new_line = f"{atom2:6d} {atom1:6d} {atom3:6d} {k:12s} {theta0:12s}"
        # Append the new line to the list of new lines
        new_lines.append(new_line)

    # Join the new lines with newlines and return the result
    return "\n".join(new_lines)

bending_section = transform_bending_section(angles)

#Write the [ angles ] section

def generate_angles_section(angles):
    # Initialize an empty list to store the angle lines
    angle_lines = []

    # Iterate over each line in the angles string
    for line in angles.strip().split("\n"):
        # Extract the atom indices and equilibrium angle from the line
        parts = line.split()
        i = int(parts[0]) - 1  # Subtract 1 to convert from 1-indexed to 0-indexed
        j = int(parts[1]) - 1
        k = int(parts[2]) - 1
        theta0 = float(parts[4])
        fc = float(parts[3]) * 4.184 * 2  # Convert from kcal/mol rad^2 to kJ/mol rad^2

        # Format the angle line and append it to the list
        angle_line = f"{i+1:6d} {j+1:6d} {k+1:6d} 1 {theta0:.3f} {fc:.3f}"
        angle_lines.append(angle_line)

    # Join the angle lines into a single string with newline characters between each line
    angles_section = "[ angles ]\n;  ai    aj    ak type    theta0 (degr)   fc (kJ/(mol rad2))\n" + "\n".join(angle_lines)

    # Return the angles section string
    return angles_section

angles_section = generate_angles_section(bending_section)

## Proper torsion section
#Tidy up function

torsion = sections['torsion']

def transform_torsion_section(torsion_section):
    new_lines = []
    # Split the section into lines and skip the header
    lines = torsion_section.strip().split("\n")[1:]
    for line in lines:
        # Extract the atom indices from the line
        parts = line.split()
        atom1 = int(''.join(filter(str.isdigit, parts[0])))
        atom2 = int(''.join(filter(str.isdigit, parts[1])))
        atom3 = int(''.join(filter(str.isdigit, parts[2])))
        atom4 = int(''.join(filter(str.isdigit, parts[3])))
        # Extract constants
        c1 = float(parts[4])
        c2 = float(parts[5])
        c3 = float(parts[6])
        c4 = float(parts[7])
        # Format the new line
        new_line = f"{atom1:6d} {atom2:6d} {atom3:6d} {atom4:6d} {c1:.3f} {c2:.3f} {c3:.3f} {c4:.3f}"
        # Append the new line to the list of new lines
        new_lines.append(new_line)

    # Join the new lines with newlines and return the result
    return "\n".join(new_lines)

torsion_section = transform_torsion_section(torsion)

#Create the [ dihedrals ] section

def generate_dihedrals_section(torsion_section):
    # Initialize an empty list to store the dihedral lines
    dihedral_lines = []

    # Iterate over each line in the torsion section string
    for line in torsion_section.strip().split("\n"):
        # Extract the atom indices and Fourier coefficients from the line
        parts = line.split()
        i = int(parts[0]) - 1  # Subtract 1 to convert from 1-indexed to 0-indexed
        j = int(parts[1]) - 1
        k = int(parts[2]) - 1
        l = int(parts[3]) - 1
        c1 = float(parts[4])  * 4.184
        c2 = float(parts[5])  * 4.184
        c3 = float(parts[6])  * 4.184
        c4 = float(parts[7])  * 4.184

        # Format the dihedral line and append it to the list
        dihedral_line = f"{i+1:6d} {j+1:6d} {k+1:6d} {l+1:6d} 5 {c1:.6f} {c2:.6f} {c3:.6f} {c4:.6f}"
        dihedral_lines.append(dihedral_line)

    # Join the dihedral lines into a single string with newline characters between each line
    dihedrals_section = "[ dihedrals ]\n; Type 5 Fourier\n;  ai    aj    ak    al  type     coefficients\n" + "\n".join(dihedral_lines)

    # Return the dihedrals section string
    return dihedrals_section

dihedrals_section = generate_dihedrals_section(torsion_section)

##Improper torsion section
# Tidy up

improper = improper = sections['improper']

def transform_improper_torsion_section(improper_torsion_section):
    new_lines = []
    # Split the section into lines and skip the header
    lines = improper_torsion_section.strip().split("\n")[1:]
    for line in lines:
        # Extract the atom indices and constants from the line
        parts = line.split()
        atom1 = int(''.join(filter(str.isdigit, parts[0])))
        atom2 = int(''.join(filter(str.isdigit, parts[1])))
        atom3 = int(''.join(filter(str.isdigit, parts[2])))
        atom4 = int(''.join(filter(str.isdigit, parts[3])))
        constant = float(parts[4])
        # Format the new line
        new_line = f"{atom1:6d} {atom2:6d} {atom3:6d} {atom4:6d} {constant:.3f}"
        # Append the new line to the list of new lines
        new_lines.append(new_line)

    # Join the new lines with newlines and return the result
    return "\n".join(new_lines)

improper_section = transform_improper_torsion_section(improper)

## Generate the improper [ dihedrals ] section:

def generate_improper_dihedrals_section(improper_torsion_section):
    # Check if the improper torsion section is empty
    if not improper_torsion_section:
        return ""

    # Transform the improper torsion section into a list of lines
    lines = improper_torsion_section.strip().split("\n")

    # Initialize an empty list to store the dihedral lines
    dihedral_lines = []

    # Iterate over each line in the improper torsion section
    for line in lines:
        # Extract the atom indices and constant from the line
        parts = line.split()
        i = int(parts[0]) - 1  # Subtract 1 to convert from 1-indexed to 0-indexed
        j = int(parts[1]) - 1
        k = int(parts[2]) - 1
        l = int(parts[3]) - 1
        constant = float(parts[4]) * 4.184/2

        # Format the dihedral line and append it to the list
        dihedral_line = f"{i+1:6d} {j+1:6d} {k+1:6d} {l+1:6d} 4 180.000 {constant:.6f} 2.000"
        dihedral_lines.append(dihedral_line)

    # Join the dihedral lines into a single string with newline characters between each line
    dihedrals_section = "[ dihedrals ]\n; Periodic improper dihedrals (type 4)\n;  ai    aj    ak    al  type     phi0    fc (kJ/mol)   n  \n" + "\n".join(dihedral_lines)

    # Return the dihedrals section string
    return dihedrals_section

#Store

imp_dihedrals_section = generate_improper_dihedrals_section(improper_section)

## Generate the [ pairs ] section from the dihedrals according to the 1-4 rule

# Parse the dihedrals section
dihedrals = []
for line in dihedrals_section.split("\n")[3:]:
    line = line.strip()
    if not line or line.startswith(";"):
        continue
    cols = line.split()
    dihedrals.append((int(cols[0]), int(cols[1]), int(cols[2]), int(cols[3])))

# Generate the set of pairs
pairs = set()
for i in range(len(dihedrals)):
    for j in range(i+1, len(dihedrals)):
        # Check if the first and fourth atoms of the current dihedral are in the current pair
        if dihedrals[i][0] in (dihedrals[j][1], dihedrals[j][2], dihedrals[j][3]) and dihedrals[i][3] in (dihedrals[j][1], dihedrals[j][2], dihedrals[j][3]):
            continue
        if dihedrals[j][0] in (dihedrals[i][1], dihedrals[i][2], dihedrals[i][3]) and dihedrals[j][3] in (dihedrals[i][1], dihedrals[i][2], dihedrals[i][3]):
            continue
        # Check if the first and fourth atoms of the current pair are in another dihedral
        pair = (min(dihedrals[i][0], dihedrals[i][3]), max(dihedrals[i][0], dihedrals[i][3]))
        skip_pair = False
        for k in range(len(dihedrals)):
            if k == i or k == j:
                continue
            if pair[0] in (dihedrals[k][1], dihedrals[k][2], dihedrals[k][3]) and pair[1] in (dihedrals[k][1], dihedrals[k][2], dihedrals[k][3]):
                skip_pair = True
                break
        if skip_pair:
            continue
        pairs.add(pair)

# Generate the pairs section string
pairs_section = "[ pairs ]\n;  ai    aj       f_qq    qi     qj   sigma (nm)  epsilon (kJ/mol) \n"
for pair in sorted(pairs):
    pairs_section += f"  {pair[0]:5d} {pair[1]:5d}  1\n"

## Print out the .itp file 

def generate_itp_file(atoms_section, bonds_section, angles_section, dihedrals_section, imp_dihedrals_section, pairs_section):
    # Construct the ITP file string
    itp_file = f""";Generated by ffld2gmx.py
;Residue topology file 
;Created by Juan de Gracia


[ moleculetype ]
; Name nrexcl
UNK 3

{atoms_section}

{bonds_section}

{angles_section}

{dihedrals_section}

{imp_dihedrals_section}

{pairs_section}
"""

    # Write the ITP file to disk
    with open("UNK.itp", "w") as f:
        f.write(itp_file)

generate_itp_file(atoms_section, bonds_section, angles_section, dihedrals_section, imp_dihedrals_section, pairs_section)

#Print out the topol.top file

def generate_topol_file(atomtypes):
    # Construct the topol.top file string
    topol_file = f""";Generated by ffld2gmx.py
;General topology file using the OPLS-AA force field
;Created by Juan de Gracia    

#include "oplsaa.ff/forcefield.itp"

{atomtypes}

;Include the residue topology
#include "UNK.itp"

[ system ]
; Name
SYS

[ molecules ]
; Compound        #mols
UNK               1
"""

    # Write the topol.top file to disk
    with open("topol.top", "w") as f:
        f.write(topol_file)

generate_topol_file(atomtypes_section)





