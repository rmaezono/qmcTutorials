#!/usr/bin/python

 #---------#
 # VERSION #
 #  1.0.4  #
 #---------#

 #---------#
 #  NOTES  #
 #---------#
# For PP calculations, any atom with Z > 18, it is assumed that the PP removes only the noble gas core (i.e., krypton has a valence charge of 18, not 8 as if the d electrons were removed.)
# This script does not currently support sp functions, though I left the space for them if someone is inclined to go add them.  
# See the README file for more detailed information.

import struct
import math
import os
import sys
import re

print ('' + "\n" 'Hello, you are converting a MOLDEN output to a CASINO gwfn.data file.' ) 

 #-----------------------------------------------------------------------------------#
 # Open user specified MOLDEN file and determine which code produced the MOLDEN file #
 #-----------------------------------------------------------------------------------#
output_file = raw_input("" + "\n" + "Enter the name of your MOLDEN file: ")
while os.path.exists(str(output_file)) == False:
	print "" + "\n" + "File not found..."
	output_file = raw_input("" + "\n" + "Enter the name of your MOLDEN file: ")
f = open(str(output_file),"r")

code = int(raw_input( '\n' + 'Enter the NUMBER corresponding to the quantum chemistry code used to produce this MOLDEN file: \n\
\n\
0 -- MOLPRO or TURBOMOLE \n\
1 -- PSI4 \n\
2 -- C4 \n\
3 -- ORCA \n\
\n'))
while code > 3 :
	code = int(raw_input('Sorry,  try again.'))
	
print '\n' , "You have entered NUMBER ", code

#### IF a user would like to try this with ORCA, they only need to remove the 'code == 3' portion of the next line.  
#### BUT USE CAUTION! Always run VMC without a Jastrow factor to be sure that the VMC and HF energies agree. 
if code == 2  or code ==3:
	sys.exit('\nSorry, this is not functioning yet, but it is on the to-do list, so I left it in here.')

 #----------------------------------#
 # Get the nuclear repulsion energy #
 #----------------------------------#
f.seek(0,0) # This resets the seek to the begining of the MOLDEN file.
for lines in f:
	nuc = lines
	if re.search('_ENUC=*', nuc): # Currently, MOLPRO is the only code that dumps this variable, and this is the format it uses. 
		break
nuc3 = nuc.split()
if nuc3[0] != '_ENUC=' :
        nuclear_energy = 'Get this from your quantum chemistry package. It is the total n-n repulsion energy divided by the number of atoms.'
else: 
        nuclear_energy = float(nuc3[1])
print ('\n' + "Nuclear energy = %s" % (nuclear_energy,) + "\n")

 #-----------------------------#
 # Get the number of electrons #
 #-----------------------------#
f.seek(0,0)
for lines in f:
        elec = lines
        if re.search('_NELEC=*', elec):
                break
elec2 = elec.split()
if elec2 [0] != '_NELEC=':
	number_elec = 'Also get this from your quantum chemistry calculation.'
else:
	number_elec = int(float(elec2[1]))
print ("Number of electrons = %s" % (number_elec,) + "\n")

 #---------------------------------------------------------------------#
 # Get the number of atoms, atomic numbers, and atom positions in a.u. # 
 #---------------------------------------------------------------------#
raw_atom_list = [] # a list of lines(as strings), and each line contains the atomic symbol, the ith atom number, the valence charge, and the x y z coords.
atomic_number = [] # a list of atomic numbers.
atom_position = [] # a list of atomic coordinates.
ab = 0.52917721    # 1 Angstrom = 0.52917721 Bohr
valence_charge = []
number_to_remove = 0  # To remove things at the end of the atom list that are unnecessary. 
shift = 0      # This is a variable to tell seek how far to move to be in the first column of the first line of the list of atoms.  
if code == 0:  # Unfortunately, it varies based on which qauntum chemistry code produced the MOLDEN file.
	shift = 12
elif code == 1:
	shift = 12
elif code == 2:
	shift = 12
elif code == 3: 
	shift = 11

f.seek(0,0)
lines = f.read()
u = lines.find("[Atoms]")
if u == -1:
	u = lines.find("[ATOMS]")
	f.seek(u+shift,0)
else: f.seek(u+shift,0)     # Puts you on the line of the first atom. 
for lines in f:
	raw_atom_list.append(lines)
        if re.search('[GTO]', lines ):
		break
if raw_atom_list[0] == '\n':
	raw_atom_list.pop(0)


for i in range(len(raw_atom_list)):  # This loop counts anything that may be at the end of the atom list, and removes it.
	if len(raw_atom_list[i]) < 20:  # I chose 20 becasue no item in the list should be longer than ~10. Admittedly, I may have to change this in the future, depending on what codes get implented.
		number_to_remove = number_to_remove + 1
if number_to_remove >= 1:
	raw_atom_list[-number_to_remove:] = []

number_of_atoms = len(raw_atom_list)  # This gives the total number of atoms in the simulation.
if nuclear_energy == 'Get this from your quantum chemistry package. It is the total n-n repulsion energy divided by the number of atoms.':
	nuclear_nuclear_energy_per_atom = nuclear_energy
else:
	nuclear_nuclear_energy_per_atom = nuclear_energy / number_of_atoms

f.seek(0,0)  # This finds out if the atom coordinates are printed in automic units or in Angstroms, and stores them in a.u. 
lines = f.read()
u = lines.find("[Atoms]")
if u == -1:
	u = lines.find("[ATOMS] AU")
f.seek(u,0)
v = f.readline()
w = v.split()
if w[1] == 'Angs':
	for i in range(0,number_of_atoms): 
		y = raw_atom_list[i]
		atom = y.split()
		atomic_number.append(int(float(atom[2])))
		atom_position.append([float(atom[3])/ab, float(atom[4])/ab, float(atom[5])/ab])
else:
        for i in range(0,number_of_atoms):
                y = raw_atom_list[i]
                atom = y.split()
       	        atomic_number.append(int(float(atom[2])))
               	atom_position.append([float(atom[3]), float(atom[4]), float(atom[5])])


 #----------------------------------------------------------------------#
 # Get the correct valence charge for each atom (dependent on AE or PP) #
 #----------------------------------------------------------------------#
valence_charge = []
v = sum(atomic_number)
pseudo_used = 0 # I decided to ask the user if a PP was used, because I feel that it is more common, and changing all of the valence charges in a large calculation would be tedious.
if number_elec == 'Also get this from your quantum chemistry calculation.':
	pseudo_used = int(raw_input('This script did not detect if a pseudopotential was used.  Was a pseudopotential used for the majoritiy of the atoms? \n\
Please eneter the NUMBER that corresponds to your answer: \n\
\n\
1 = YES, a pseudopotential was used for the majority of the atoms in this calculation. \n\
2 = NO, a pseudopotential was not used for the majority of the atoms in this cacluation. \n\
\n'))
	if pseudo_used == 1:
        	for i in range(0,number_of_atoms):
                	if   atomic_number[i] <= 2:
                        	valence_charge.append(atomic_number[i])
	                elif atomic_number[i] <= 10:
        	                valence_charge.append(atomic_number[i] - 2)
                	elif atomic_number[i] <= 18:
	                        valence_charge.append(atomic_number[i] - 10)
        	        elif atomic_number[i] <= 36:
                	        valence_charge.append(atomic_number[i] - 18)
	                elif atomic_number[i] <= 54:
        	                valence_charge.append(atomic_number[i] - 36)
                	elif atomic_number[i] <=86:
                        	print "If you're crazy enough to try QMC on this, you can update this script."
			atomic_number[i] += 200
	else:
	        for i in range(0,number_of_atoms):
	                valence_charge.append(atomic_number[i])
        	print 'You will need to change the valence charge if you used a PP. This script did not detect a PP in the MOLDEN file. As such, it prints the valence chage as the atomic number. \n'

elif v != number_elec:
        for i in range(0,number_of_atoms):
                if   atomic_number[i] <= 2:
                        valence_charge.append(atomic_number[i])
                elif atomic_number[i] <= 10:
                        valence_charge.append(atomic_number[i] - 2)
                elif atomic_number[i] <= 18:
                        valence_charge.append(atomic_number[i] - 10)
                elif atomic_number[i] <= 36:
                        valence_charge.append(atomic_number[i] - 18)
                elif atomic_number[i] <= 54:
                        valence_charge.append(atomic_number[i] - 36)
                elif atomic_number[i] <=86:
                        print "If you're crazy enough to try QMC on this, you can update this script."
		atomic_number[i] += 200
else:
        for i in range(0,number_of_atoms):
                valence_charge.append(atomic_number[i])
        print 'You will need to change the valence charge if you used a PP. This script did not detect a PP in the MOLDEN file. As such, it prints the valence chage as the atomic number. \n'

 #----------------------------------------------#
 # Get and manipulate the basis set information #
 #----------------------------------------------#
highest_ang_mo = 0
raw_basis_list = []
basis_list = [] # This will keep a list of all basis functions, to be used later when reading the AOs.
AO_list = [] # A list for reconstructing the AOs after the d's and f's (and possibly g's in the future) are transformed and normalized.
shell_counter = 0
shell_code = []
primitives = 0
primitives_per_shell_counter = 0
primitives_per_shell = []
two_over_pi_to_point75 = (2/math.pi)**0.75
rt15 = 15.**0.5
basis_functions = 0
sc = [] # A list for s coefficients for an atom.
pc = []
dc = []
fc = []
spc = []
coefficients = [] # Master list to hold all coefficients.
se = [] # A list for s exponents for an atom.
pe = []
de = []
fe = []
spe = []
exponents = [] # Master list to hold all exponents.
shell_position = []
first_shell_sequence_counter = 1
first_shell_sequence = [1]

f.seek(0,0)
lines = f.read()
u = lines.find("[GTO]")
f.seek(u+6,0)     # Puts you on the line of the first basis set.
for lines in f:   # Loop to read all of the basis information. 
        raw_basis_list.append(lines)
        if re.search('[MO]', lines ):
                break
for i in range(0,len(raw_basis_list)):  # Loop to count the total number of shells and code each shell in CASINO format, as well as normalize the contraction coefficients.
	shell = raw_basis_list[i]       # Get the first line of the raw basis data.
	shell1 = shell.split()   # Split the line into a list to be evaulated.
	if len(shell1) == 0: # A blank line, indicating one atom just ended and a new one will be starting.
		first_shell_sequence.append(first_shell_sequence_counter)
		for i in range(len(sc)):
			coefficients.append(sc[i])
                for i in range(len(pc)):
			coefficients.append(pc[i])
                for i in range(len(dc)):
                	coefficients.append(dc[i])
                for i in range(len(fc)):
	                coefficients.append(fc[i])
                for i in range(len(spc)):
	                coefficients.append(spc[i])
                for i in range(len(se)):
	                exponents.append(se[i])
                for i in range(len(pe)):
        	        exponents.append(pe[i])
                for i in range(len(de)):
                	exponents.append(de[i])
                for i in range(len(fe)):
                        exponents.append(fe[i])
                for i in range(len(spe)):
        	        exponents.append(spe[i])
		del sc[:]	# These clear the lists coefficient and exponent lists, so that for each new atom, the previous basis set isn't repeated.	
	        del sc[:]
                del pc[:]
                del dc[:]
                del fc[:]
                del spc[:]
                del se[:]
                del pe[:]
                del de[:]
                del fe[:]
                del spe[:]
	elif len(shell1) == 2:  # Basis set info or atom number
		shell3 = shell1[1]
		if shell3 != "0": # This is a line that contains basis set info, and needs to be stored and manipulated appropriately.  
			if shell2 == "s":
				ss = list(shell1[0]) # Turns the exponent into a list, so that the D can be replaced with an E for python scientific notation.
				if ss[-4]=='D':
					ss[-4] = 'E'         # Replaces said D.
				ss = ''.join(ss)     # Turns the list back into a string, so it can be stored as a float.
				se.append(float(ss)) # Stores the exponent as a float.
				if code == 0: 
					ss=(two_over_pi_to_point75 ) * (float(ss) ** 0.75) # This gets the new normalization constant for the coefficient. 
				elif code == 1:
					ss = 1   # Psi-4 prints the correctly normalized contraction coefficients in the MOLDEN file, so we don't need to modify the exponent  
				#elif code == 2; 
					#whatever C4 wants
				elif code == 3:
					ss = 1
                                sss = list(shell1[1]) # Same deal as the exponent, but now for the coefficient.
				if sss[-4] =='D':
	                                sss[-4] = 'E'
                                sss = ''.join(sss)
				sc.append(float(sss) * ss) # Stores the coefficient (sss) times the new normalization constant (ss, the exponent times the normalization, above).
				primitives_per_shell_counter = primitives_per_shell_counter + 1
			elif shell2 == "p":
				pp = list(shell1[0])
				if pp[-4] =='D':
					pp[-4] = 'E'
				pp = ''.join(pp)
			        pe.append(float(pp))
				if code == 0:
	                                pp = 2 * ((two_over_pi_to_point75 ) * (float(pp) ** 1.25))
				elif code == 1:
					pp = 1
				elif code == 3:
					pp = 1  
				ppp = list(shell1[1])
				if ppp[-4] =='D':
					ppp[-4] = 'E'
				ppp = ''.join(ppp)
                                pc.append(float(ppp) * pp)
                                primitives_per_shell_counter = primitives_per_shell_counter + 1
                        elif shell2 == "d":           
                                dd = list(shell1[0])
				if dd[-4]=='D':
	                                dd[-4] = 'E'
                                dd = ''.join(dd)
                                de.append(float(dd))
				if code == 0:
					 dd = 2*((two_over_pi_to_point75 ) * (float(dd) ** 1.75)) # This part of the normalization is the same for all six cartesian shells
				elif code == 1:							# The rest of the normalization and transformation of cartesian to harmonic is done on the vecotrs (AOs).
                        	        dd = 1   # This part is the same for all five harmonics, but again, Psi4 already normalizes them.
				#elif code == 2:
					#dd = ((two_over_pi_to_point75 ) * (float(dd) ** 1.75))
				elif code == 3:
					dd = 1 /(3.**0.5)
                                ddd = list(shell1[1])                                   
				if ddd[-4] =='D':
        	                        ddd[-4] = 'E'
                                ddd = ''.join(ddd)
                                dc.append(float(ddd) * dd)
                                primitives_per_shell_counter = primitives_per_shell_counter + 1
                        elif shell2 == "f":  
                                ff = list(shell1[0])
				if ff[-4] == 'D':
	                                ff[-4] = 'E'
                                ff = ''.join(ff)          
                                fe.append(float(ff))
				if code == 0:
	                                ff =  float(ff) ** 2.25
				elif code == 1: 
					ff = 1 
				elif code == 2:
                                	ff =  float(ff) ** 2.25 # Not sure what C4 will want.
				elif code == 3:
					ff = 1 / rt15
                                fff = list(shell1[1])
        			if fff[-4] == 'D':
		                        fff[-4] = 'E'
                                fff = ''.join(fff)
                                fc.append(float(fff) * ff)
                                primitives_per_shell_counter = primitives_per_shell_counter + 1
                        elif shell2 == "sp": 
                                spsp= list(shell1[0])
                                spsp[-4] = 'E'
                                spsp = ''.join(spsp)
                                spe.append(float(spsp))
                                spsp = ((2 * float(spsp)) /math.pi ) ** 0.75 # This is probably wrong, so I'm pretty much giving up on sp shells right here.
                                spspsp = list(shell1[1])
                                spspsp[-4] = 'E'
                                spspsp = ''.join(spspsp)
                                spc.append(float(spspsp) * spsp)
                                primitives_per_shell_counter = primitives_per_shell_counter + 1
			primitives = primitives + 1
		else:
			shell4 = int(float(shell1[0])) -1  # This gives the number of the atom that owns the basis set (minus one because MOLDEN counts starting at 1, python counts starting at 0).
	elif len(shell1) == 3: # Shell definition, i.e. s, p, d, f, etc.
		shell2 = shell1[0]
		if   shell2 == "s":
			shell_counter = shell_counter + 1 
			shell_code.append(1) 
			first_shell_sequence_counter = first_shell_sequence_counter +1
                        shell_position.append(atom_position[shell4])
			basis_functions = basis_functions + 1 
			basis_list.append('s')
			AO_list.append('s')
			if highest_ang_mo < 1:
				highest_ang_mo = 1
			primitives_per_shell.append(primitives_per_shell_counter)
			primitives_per_shell_counter = 0
		elif shell2 == "p":
        	        shell_counter = shell_counter + 1
                	shell_code.append(3)
                        first_shell_sequence_counter = first_shell_sequence_counter +1
                        shell_position.append(atom_position[shell4])
                        basis_functions = basis_functions + 3
                        basis_list.append('p')
                        basis_list.append('p')
                        basis_list.append('p')
                        AO_list.append('p')
                        AO_list.append('p')
                        AO_list.append('p')
			if highest_ang_mo < 2:
				highest_ang_mo = 2
                        primitives_per_shell.append(primitives_per_shell_counter)
                        primitives_per_shell_counter = 0
	        elif shell2 == "d":
        	        shell_counter = shell_counter + 1
                	shell_code.append(4)
                        first_shell_sequence_counter = first_shell_sequence_counter +1
                        shell_position.append(atom_position[shell4])
                        basis_functions = basis_functions + 5
			for i in range(5):
				basis_list.append('d')
			if code == 0 or code == 2:  # TURBOMOLE, MOLPRO, and C4 all print Cartesian MOLDEN files.
				basis_list.append('d')
			for i in range(5):
				AO_list.append('d')  # Only five, because the Cartesians will be transformed to spherical later, and this list is used to keep track of the final AOs. 
			if highest_ang_mo < 3:
				highest_ang_mo = 3
                        primitives_per_shell.append(primitives_per_shell_counter)
                        primitives_per_shell_counter = 0
	        elif shell2 == "f":
        	        shell_counter = shell_counter + 1
                	shell_code.append(5)
                        first_shell_sequence_counter = first_shell_sequence_counter + 1
                        shell_position.append(atom_position[shell4])
                        basis_functions = basis_functions + 7
			for i in range(7):
				basis_list.append('f')  # Seven spherical fs
                        if code == 0 or code == 2:
				basis_list.append('f')  # Ten Cartesian fs.
	                        basis_list.append('f')
        	                basis_list.append('f')
			for i in range(7):
				sys.exit('Sorry, Cartesian f functions are not supported')
				#AO_list.append('f')
			if highest_ang_mo < 4:
				highest_ang_mo = 4
                        primitives_per_shell.append(primitives_per_shell_counter)
                        primitives_per_shell_counter = 0
	        elif shell2 == "sp":
        	        shell_counter = shell_counter + 1
                	shell_code.append(2)
                        first_shell_sequence_counter = first_shell_sequence_counter +1
                        shell_position.append(atom_position[shell4])
if first_shell_sequence[-1] == first_shell_sequence[-2]:
	first_shell_sequence[-1:] = [] # This removes the last item on the list, which repeats the last shell number due to the blank line at the end of the raw_basis_set list.
primitives_per_shell.pop(0) # This removes the first item on the list, which is a 0.
primitives_per_shell.append(primitives_per_shell_counter) # This appends the very last shell counter, which was not done above.

 #-----------------------------------------------------------#
 # reading, storing, and manipulating the eigenvectors (AOs) #
 #-----------------------------------------------------------#
raw_vectors = [] # The list of AOs read from the MOLDEN file.
s_vectors = []
p_vectors = []
d6_vectors = [] # Master list of Cartesian vectors read from MOLDEN file.
d6_vectors2 = [] # Temporary list of 6 Cartesian vectors, to be transformed to spherical.
d5_vectors = [] # The final normalized spherical vectors to be resorted for printing.
d5_vectors1 = [] # Master list of spherical vectors read from MOLDEN file.
d5_vectors2 = [] # Temporary list of 5 d vectors, to be renormalized.
f10_vectors = []
f10_vectors2 = []
f7_vectors = []
f7_vectors1 = []
f7_vectors2 = []
vectors = [] # Final list of all vectors, reconstructed in the proper order so that they correspond to the basis sets above.
shell_tracker = 0 # For counting which shell we are on.
s_tracker = 0
p_tracker = 0
d_tracker = 0
f_tracker = 0
AO_tracker = 0
MOs = 0
### SOME CONSTANTS FOR NORMALIZATION ###
one_over_rt_three = 3**(-0.5)
rt15 = 15**0.5
rt6 = 6**0.5
rt10 = 10**0.5
rt2 = 2**0.5
harm_f_1 = (8./15.0)**0.5
harm_f_2 = 2.0 * (1/(5**0.5))/3.0
harm_f_3 = rt2 / 15.0
harm_f_4 = (2.0**1.5) / 30.0
harm_f_5 = one_over_rt_three/15.0

f.seek(0,0)
lines = f.read()
u = lines.find("[MO]")
f.seek(u+5,0)     # Puts you on the line of the first MO.
for lines in f:
        raw_vectors.append(lines) # Reads in all of the vectors.
	if lines == '[Molden Format]\n':  # For some unknown reason, some times Psi4 will print the whole MOLDEN file multiple times. If this happens, the final printing is the one that is needed.  
		break
if raw_vectors[-1] == '[Molden Format]\n':
	 sys.exit("There is a problem with this MOLDEN file. Psi4 will occasionally print the entire MOLDEN file multiple times in one file.  \n\
If this happens, the final printing should be used, and everything before the final occurance of the line [Molden Format] should be deleted.\n")

for i in range(len(raw_vectors)):  # This will get the vectors and sort them by their shell type.
	vec = raw_vectors[i]
	vec2 = vec.split()
	if vec2[0] == "Sym=" : 
		shell_tracker = 0
	elif vec2[0] == "Ene=":
		shell_tracker = 0
	elif vec2[0] == "Spin=":
		shell_tracker = 0
	elif vec2[0] == 'Occup=':
		shell_tracker = 0
        else:
                if basis_list[shell_tracker] == "s": # An s orbital.
                        s_vectors.append(float(vec2[1]))
		if basis_list[shell_tracker] == "p":
			p_vectors.append(float(vec2[1]))
		if basis_list[shell_tracker] == 'd':
			if code == 0 or code == 2:
				d6_vectors.append(float(vec2[1]))
			else: 
				d5_vectors1.append(float(vec2[1]))
		if basis_list[shell_tracker] == 'f':
                        if code == 0 or code == 2:
                                f10_vectors.append(float(vec2[1]))
                        else:
                        	f7_vectors1.append(float(vec2[1]))
		shell_tracker = shell_tracker + 1

if code == 0 or code == 2:
	for i in range(0,len(d6_vectors), 6):  # This loop will take each of the d cartesian orbitals and transform them to spherical harmonic orbitals, and also normalize them.
		for j in range(6):
			d6_vectors2.append(d6_vectors[j+i])
		yy1 =  one_over_rt_three * (d6_vectors2[0] - d6_vectors2[1])
		zz1 =  (1.0/3.0) * one_over_rt_three * ((2*d6_vectors2[2]) - d6_vectors2[0] - d6_vectors2[1])
		d5_vectors.append(zz1)              # z^2
		d5_vectors.append(2*d6_vectors2[4]) # xz
	        d5_vectors.append(2*d6_vectors2[5]) # yz
        	d5_vectors.append(yy1)              # x^2 - y^2
	        d5_vectors.append(2*d6_vectors2[3]) # xy
		del d6_vectors2 [:]

elif code == 1 or code == 3: 
        for i in range(0,len(d5_vectors1), 5):  # This loop will take each of the d spherical orbitals and normalize them correctly.
                for j in range(5):
                        d5_vectors2.append(d5_vectors1[j+i])
                d5_vectors.append(d5_vectors2[0] * one_over_rt_three)
                d5_vectors.append(d5_vectors2[1] * 2.)
                d5_vectors.append(d5_vectors2[2] * 2.) 
                d5_vectors.append(d5_vectors2[3] ) 
                d5_vectors.append(d5_vectors2[4] * 2.)
                del d5_vectors2 [:]

if code == 0 or code == 2:
        for i in range(0,len(f10_vectors), 10):  # This loop will take each of the f cartesian orbitals and transform them to spherical harmonic orbitals and multiply the m dependent normalization.
                for j in range(10):
                        f10_vectors2.append(f10_vectors[j+i])
                zero    = f10_vectors2[2] - ((3/2)*(f10_vectors2[5] + f10_vectors2[8]))
                plus_1  = 4*f10_vectors2[6] - (f10_vectors2[3] + f10_vectors2[0])
                minus_1 = 4*f10_vectors2[7] - (f10_vectors2[1] + f10_vectors2[4])
                plus_2  = f10_vectors2[5] - f10_vectors2[8]
                minus_2 = 2*f10_vectors2[9]
                plus_3  = f10_vectors2[0] - 3*f10_vectors2[3]
                minus_3 = 3*f10_vectors2[4] -  f10_vectors2[1]
		
                f7_vectors.append(zero  *  0.331990278)
		f7_vectors.append(plus_1 * 0.059894236)
		f7_vectors.append(minus_1* 0.064894235)
		f7_vectors.append(plus_2*  0.059227155)
		f7_vectors.append(minus_2* 0.050227159)
		f7_vectors.append(plus_3 * 0.010915709)
		f7_vectors.append(minus_3* 0.010915709)
		del f10_vectors2 [:]


elif code == 1 or code == 3:
        for i in range(0,len(f7_vectors1), 7):  # Tis loop will take each of the f spherical orbitals and normalize them correctly.
                for j in range(7):
                        f7_vectors2.append(f7_vectors1[j+i])
                f7_vectors.append(f7_vectors2[0] * harm_f_1)
                f7_vectors.append(f7_vectors2[1] * harm_f_2)
                f7_vectors.append(f7_vectors2[2] * harm_f_2)
                f7_vectors.append(f7_vectors2[3] * harm_f_3)
                f7_vectors.append(f7_vectors2[4] * harm_f_4)
                f7_vectors.append(f7_vectors2[5] * harm_f_5)
                f7_vectors.append(f7_vectors2[6] * harm_f_5)
                del f7_vectors2 [:]

MOs = len(s_vectors) + len(p_vectors) + len(d5_vectors) + len(f7_vectors)

for i in range(MOs): # This loop will resort the AO's in their original order for printing in gwfn.data.
	if float(i) % float(len(AO_list)) == 0:
		AO_tracker = 0
	if AO_list[AO_tracker] == 's': # An s orbital.
		vectors.append(s_vectors[s_tracker])
		s_tracker = s_tracker + 1
	if AO_list[AO_tracker] == 'p':  
		vectors.append(p_vectors[p_tracker])
		p_tracker = p_tracker + 1
	if AO_list[AO_tracker] == 'd': 
		vectors.append(d5_vectors[d_tracker])
                d_tracker = d_tracker + 1
        if AO_list[AO_tracker] == 'f':
                vectors.append(f7_vectors[f_tracker])
                f_tracker = f_tracker + 1
	AO_tracker = AO_tracker + 1

 #-----------------------------------#
 # Finally, print the gwfn.data file #
 #-----------------------------------#
g = open('gwfn.data' , 'w')
g.write("TITLE\n\
Insert Your Title Here \n\
BASIC_INFO\n\
----------\n\
Generated by:\n\
MOLDEN CONVERSION\n\
Method:\n\nDFT Functional:\n\
\n\
Periodicity:\n\
0\n\
Spin unrestricted:\n\
Put .false. in here if this was a RHF or RKS calculation, put .true. if UHF or UKS. \n\
nuclear-nuclear repulsion energy (au/atom):\n\
   %s\n\
Number of electrons per primitive cell:\n\
          %s\n\n\
GEOMETRY\n\
--------\n\
Number of atoms:\n\
          %s\n\
Atomic positions (au):" % (nuclear_nuclear_energy_per_atom, number_elec, len(atomic_number)))
for i in range(len(atom_position)):
	k = atom_position[i]
	g.write ("\n")
	for j in range(len(k)):
		g.write ("% .13E" %(k[j]))

g.write ('\nAtomic numbers for each atom:')
for i in range(len(atomic_number)):
	if i % 8 != 0:
		g.write( "%10d" % atomic_number[i])
	else:
		g.write("\n %9d" % atomic_number[i])

g.write ('\nValence charges for each atom:')
for i in range(len(valence_charge)):
	if i % 4 !=0:
		g.write ('%1.13E ' % valence_charge[i])
	else:
		g.write ('\n %1.13E ' % valence_charge[i])
g.write ( '\n\nBASIS SET\n\
---------\n\
Number of Gaussian centres\n\
          %s\n\
Number of shells per primitive cell\n\
         %s\n\
Number of basis functions (\'AO\') per primitive cell\n\
         %s\n\
Number of Gaussian primitives per primitive cell\n\
         %s\n\
Highest shell angular momentum (s/p/d/f... 1/2/3/4...)\n\
           %s\n\
Code for shell types (s/sp/p/d/f... 1/2/3/4/5...)' % (len(atom_position) , shell_counter, basis_functions , primitives , highest_ang_mo))
for i in range(len(shell_code)):
        if i % 8 !=0:
                g.write ('%10d' % shell_code[i])
        else:
                g.write ('\n%10d' % shell_code[i])
g.write ('\nNumber of primitive Gaussians in each shell')
for i in range(len(primitives_per_shell)):
        if i % 8 !=0:
                g.write ('%10d' % primitives_per_shell[i])
        else:
                g.write ('\n%10d' % primitives_per_shell[i])
g.write ('\nSequence number of first shell on each centre')
for i in range(len(first_shell_sequence)):
	if i % 8 !=0:
		g.write ('%10d' % first_shell_sequence[i])
	else:
		g.write ('\n%10d' % first_shell_sequence[i])
g.write ('\nExponents of Gaussian primitives')
for i in range(len(exponents)):
        if i % 4 !=0:
                g.write ('%1.13E ' % exponents[i])
        else:
                g.write ('\n %1.13E ' % exponents[i])
g.write ('\nNormalized contraction coefficients')
for i in range(len(coefficients)):
        if i % 4 !=0:
               	if coefficients[i] >= 0:
			g.write (' %1.13E' % coefficients[i])
		else:
			g.write ('%1.13E' % coefficients[i])
        else:
        	if coefficients[i] >= 0:
		        g.write ('\n %1.13E' % coefficients[i])
		else: 
                        g.write ('\n%1.13E' % coefficients[i])
g.write ('\nPosition of each shell (au)')
for i in range(len(shell_position)):
        k = shell_position[i]
        g.write ("\n")
        for j in range(len(k)):
                g.write ("% .13E" %(k[j]))
g.write ('\n\n\
MULTIDETERMINANT INFORMATION\n\
----------------------------\n\
GS\n\
\n\
ORBITAL COEFFICIENTS\n\
------------------------')
for i in range(len(vectors)):
        if i % 4 !=0:
                if vectors[i] == 0.0:
                        vectors[i] = vectors[i]**2
                        g.write (' %1.13E' % vectors[i])
                elif vectors[i] > 0:
                        g.write (' %1.13E' % vectors[i])
                else:
                        g.write ('%1.13E' % vectors[i])
        else:
                if vectors[i] == 0.0:
                        vectors[i] = vectors[i]**2
                        g.write ('\n %1.13E' % vectors[i])
                elif vectors[i] > 0:
                        g.write ('\n %1.13E' % vectors[i])
                else:
                        g.write ('\n%1.13E' % vectors[i])

g.write ('\n\n')
