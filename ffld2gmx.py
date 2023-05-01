## FFLD to itp converter
## Created by Juan de Gracia

import sys
import os

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

def transform_atoms_section(atoms_section):
    new_lines = []
    lines = atoms_section.strip().split("\n")[2:-1]
    for i, line in enumerate(lines):
        # Print the line to check if it matches the expected format
        #print(f"Processing line {i}: {line}")
        
        parts = line.split()
        atom = parts[0]
        charge = parts[4]
        sigma = parts[5]
        epsilon = parts[6]
        new_line = f"{i+1}  {atom}   {charge}   {sigma}   {epsilon}"
        new_lines.append(new_line)
    return "\n".join(new_lines)

new_atoms = transform_atoms_section(atoms)

#Tidy ffld will be parsed and transformed to itp

def generate_atomtypes_section(atoms_section):
    # Initialize an empty list to store the atomtype lines
    atomtype_lines = []

    # Add the section header and column headers to the list
    atomtype_lines.append("[ atomtypes ]")
    atomtype_lines.append("; name mass charge ptype sigma(nm) epsilon (kJ/mol)")

    # Transform the atoms section and iterate over each line
    for i, line in enumerate(atoms_section.strip().split('\n')[0:]):
        parts = line.split()
        atom = parts[1]
        charge = parts[2]
        sigma = parts[3]
        epsilon = parts[4]
        mass = 0  # Replace this with a function that returns the mass of the atom
        sigma = float(sigma) / 10.0
        epsilon = float(epsilon) * 4.18 * 2  # Convert from kcal/mol to kJ/mol
        new_line = f"{i+1}  op_unk.{atom}   {mass:.6f}   {charge}   A      {sigma:.6e}    {epsilon:.6e}"
        atomtype_lines.append(new_line)

    # Join the atomtype lines into a single string with newline characters between each line
    atomtypes_section = "\n".join(atomtype_lines)

    # Return the atomtypes section string
    return atomtypes_section

atomtypes = generate_atomtypes_section(new_atoms)

## Generate atom section

def generate_atoms_section(new_atoms):
    # Initialize an empty list to store the new lines
    new_lines = []
    
    # Add the section header and column headers to the new lines list
    new_lines.append("[ atoms ]")
    new_lines.append("; nr type  resnr residue  atom  cgnr charge  mass")
    
    # Transform the new_atoms string and iterate over each line
    for i, line in enumerate(new_atoms.strip().split("\n")):
        atom_nr, atom_type, charge, *_ = line.split()
        res_nr = 1
        residue = "UNK"
        atom = atom_type[0]
        cgnr = 1
        mass = 0  # Replace this with a function that returns the mass of the atom
        new_line = f"{int(atom_nr):6d}   op_unk.{atom_type:6s} {res_nr:5d} {residue:8s} {atom:6s} {cgnr:6d} {charge:9s} {mass:10.6f}"
        new_lines.append(new_line)

    # Join the new lines with newlines and return the result
    return "\n".join(new_lines)
        
atomssection = generate_atoms_section(new_atoms)

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

## Print out the .itp file 

def generate_itp_file(atomtypes, atoms_section, bonds_section, angles_section, dihedrals_section, imp_dihedrals_section):
    # Construct the ITP file string
    itp_file = f""";Generated by ffld2gmx.py
;Created by Juan de Gracia


[ moleculetype ]
; Name nrexcl
UNK 3

{atoms_section}

{bonds_section}

{angles_section}

{dihedrals_section}

{imp_dihedrals_section}
"""

    # Write the ITP file to disk
    with open("UNK_new.itp", "w") as f:
        f.write(itp_file)

generate_itp_file(atomtypes, atomssection, bonds_section, angles_section, dihedrals_section, imp_dihedrals_section)

#Print out the topol.top file

def generate_topol_file(atomtypes):
    # Construct the topol.top file string
    topol_file = f"""

;Generated by ffld2gmx.py
;Created by Juan de Gracia    

#include "oplsaa.ff/forcefield.itp"

{atomtypes}

:Include the residue topology
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

generate_topol_file(atomtypes)





