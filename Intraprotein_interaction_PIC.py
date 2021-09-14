""" Protein Interactions Calculator (PIC)

----------------------------------------------------------------
Nom : Jamay
Prénom : Théo
Cursus : M2BI
----------------------------------------------------------------


"""

import numpy as np
import math
import sys


def get_argument():
	"""Get argument et test le format du fichier.

	Returns
	-------
	filename
		nom du fichier.

	"""
	if len(sys.argv) != 2:
		sys.exit("Wrong number of argument.")
	filename = sys.argv[1]
	if filename[-4:] != ".pdb":
		sys.exit("{} n'est pas un fichier pdb".format(filename))
	if filename[4:6] != "FH" :
		sys.exit("{} doit d'abord être traité par reduce et le fichier doit être de la forme myfileFH.pdb".format(filename))
	return filename


def PDB_parser(filename) :
	"""Parse fichier PDB
	
	 Parameters
	----------
	filename : str
		Nom du fichier au format PDB

	Returns
	-------
	dict
		key : tuple avec comme argument atom_name, res_name, res_num, chain
		items : vecteur de position avec les coordonnées x,y,z
	"""
	dict_vec_pos = {}
	with open(filename, "r") as pdb_file :
		for line in pdb_file:
			if line.startswith("ATOM") :
				words =  line.split()
				atom_name = words[2]
				res_name = words[3]
				res_num = int(words[5])
				chain = words[4]
				x = float(words[6])
				y = float(words[7])
				z = float(words[8])
				dict_vec_pos[(atom_name, res_name, res_num, chain)] = np.array([x, y, z])
	return dict_vec_pos


def dist(x,y) :
	"""Calcul de distance entre 2 coordonnées 3D
	
	 Parameters
	----------
	vecteurs x et y : np.array([x,y,z])
		coordonnées atomes pour calcul distance

	Returns
	-------
	float
		distance en angstrom entre mes coordonnées x et y
	"""
	v = [y[0]-x[0], y[1]-x[1], y[2]-x[2]]
	distance_x_y = math.sqrt(v[0]**2 + v[1]**2 + v[2]**2)
	return distance_x_y


def angle(x,y,z) :
	"""Calcul de l'angle entre 3 coordonnées 3D
	
	 Parameters
	----------
	vecteurs x, y et z : np.array([x,y,z])
		coordonnées atomes pour calcul angle

	Returns
	-------
	float
		angle en degrés entre mes coordonnées x,y et z
	"""
	v_1 = np.array([x[0]-y[0],x[1]-y[1],x[2]-y[2]])
	v_2 = np.array([z[0]-y[0],z[1]-y[1],z[2]-y[2]])
	magni_1 = np.sqrt([v_1[0]**2+v_1[1]**2+v_1[2]**2])
	norma_1 = np.array([v_1[0]/magni_1,v_1[1]/magni_1,v_1[2]/magni_1])
	magni_2 = np.sqrt(v_2[0]**2+v_2[1]**2+ v_2[2]**2)
	norma_2 = np.array([v_2[0]/magni_2,v_2[1]/magni_2,v_2[2]/magni_2])
	res = norma_1[0]*norma_2[0]+norma_1[1]*norma_2[1]+norma_1[2]*norma_2[2]
	angle_rad = np.arccos(res)
	angle_deg = math.degrees(angle_rad)
	return angle_deg


def matrix_intra(dict_vec_pos) :
	"""Calcul de la matrice de distance sous forme de dictionnaire imbriqué
	
	 Parameters
	----------
	dict : dict_vec_pos
		key : tuple avec comme argument atom_name, res_name, res_num, chain
		items : vecteur de position avec les coordonnées x,y,z

	Returns
	-------
	dict in dict :
		dictionnaire imbriqué
		distance_matrix[atom_pdb1][atom_pdb2] = distance entre atom_pdb1 et atom_pdb2
	"""
	distance_matrix = {}
	for pdbid_a, cords_a in dict_vec_pos.items() :
		distance_matrix[pdbid_a] = {}
		for pdbid_b,cords_b in dict_vec_pos.items() :
			distance_matrix[pdbid_a][pdbid_b] = dist(cords_a, cords_b)
	return distance_matrix


def hydrophobic_interaction(distance_matrix,pdb_name) :
	"""Recherche interaction hydrophobe
	
	 Parameters
	----------
	dict in dict : distance_matrix
		distance_matrix[atom_pdb1][atom_pdb2] = distance entre atom_pdb1 et atom_pdb2
	pdb_name : str
		code PDB pour nommer fichier de sortie
	"""
	list_main_ch = ["N", "CA", "C", "O"]#carbone de la chaine principale
	dict_interaction = {} #dictionnaire pour enregistrer lorsque 2 résidus d'aa hydrophobe intéragissent
	list_hydrophobic = ["ALA", "ILE", "LEU", "MET", "PHE", "TRP", "TYR", "VAL","PRO"]
	with open(pdb_name+"_Intraprotein_Hydrophobic_Interactions.txt", "w") as hydrophobic_file :
		hydrophobic_file.write("Position\tResidue\tChain\tPosition\tResidue\tchain\n")
		for atom_a, value in distance_matrix.items() :
			if atom_a[1] in list_hydrophobic and atom_a[0] not in list_main_ch and atom_a[0][0]!="H":
				for atom_b, distance in value.items():
					if atom_a[3] == atom_b[3] :#vérification si interaction intrachaine
						if atom_b[1] in list_hydrophobic and atom_b[0] not in list_main_ch and distance <= 5\
							and atom_b[2]>atom_a[2] and atom_b[0][0]!="H" :
							if (atom_a[2], atom_b[2],atom_a[3]) not in dict_interaction :#ex : l'interaction (1,2,"A") pour ne
								#pas écire plusieur fois la meme intéraction dans le fichier de sortie
								dict_interaction[(atom_a[2], atom_b[2],atom_a[3])] = 1
								hydrophobic_file.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(atom_a[2],atom_a[1],atom_a[3],\
									atom_b[2],atom_b[1],atom_b[3]))


def disulphite_bridges(dict_distance,pdb_name) :
	"""Recherche pont disulfure
	
	 Parameters
	----------
	dict in dict : distance_matrix
		distance_matrix[atom_pdb1][atom_pdb2] = distance entre atom_pdb1 et atom_pdb2
	pdb_name : str
		code PDB pour nommer fichier de sortie
	"""
	dict_interaction = {}
	with open(pdb_name+"_Intraprotein_Disulphide_bridges.txt", "w") as hydrophobic_file :
		hydrophobic_file.write("Disulphide bridges: Between sulphur atoms of cysteines within 2.2 Angstroms\n")
		hydrophobic_file.write("Position\tResidue\tChain\tPosition\tResidue\tchain\tDistance\n")
		for atom_a, value in dict_distance.items() :
			if atom_a[0] == "SG" :#Atome sulfure de la cysteine
				for atom_b, distance in value.items():
					if atom_a[0] == "SG" and atom_b[2] > atom_a[2] and distance <= 2.2 :
						if (atom_a[2], atom_b[2],atom_a[3]) not in dict_interaction :
							dict_interaction[(atom_a[2], atom_b[2],atom_a[3])] = 1
							hydrophobic_file.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(atom_a[2],atom_a[1],atom_a[3],\
									atom_b[2],atom_b[1],atom_b[3],round(distance,2)))


def main_main_hydrogen_bonds(distance_matrix, dict_vec_pos,pdb_name) :
	with open(pdb_name+"_intraprotein_main_main_chain_hydrogen_bonds.txt", "w") as mm_hbond_file :
		mm_hbond_file.write("\t\tDONOR\t\t\t\t\tACCEPTOR\t\t\tPARAMETERS\t\t\n")
		mm_hbond_file.write("POS\tCHAIN\tRES\tATOM\tPOS\tCHAIN\tRES\tATOM\tDd-a\tDh-a\tdon-h-acc\tacc-don-h\n")
		for atom_a, value in distance_matrix.items() :
			if atom_a[0] == "N" :
				for atom_b, distance in value.items():
					if atom_a[3] == atom_b[3] :
						if atom_b[0] == "O" and distance < 3.5 and atom_b[2] != atom_a[2] and atom_b[2] != atom_a[2]+1\
							and  atom_b[2] != atom_a[2]-1 and atom_a[1] != "PRO": #La proline à une main chain différente
							antecedent = dict_vec_pos["C", atom_b[1], atom_b[2], atom_b[3]]#atome antecedent à l'atome accepteur O
							acceptor = dict_vec_pos[atom_a]
							donor = dict_vec_pos[atom_b]
							h = dict_vec_pos["H", atom_a[1], atom_a[2], atom_a[3]]#vecteur x,y,z de l'atome d'hydrogene lie au donneur N
							dist_h_acceptor = distance_matrix["H", atom_a[1], atom_a[2], atom_a[3]][atom_b[0], atom_b[1], atom_b[2], atom_b[3]]
							angle_1 = angle(donor,h,acceptor)
							angle_2 = angle(acceptor,donor,antecedent)
							mm_hbond_file.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(atom_a[2],\
								atom_a[3],atom_a[1],atom_a[0],atom_b[2],atom_b[3],atom_b[1],\
									atom_b[0],round(distance,2),round(dist_h_acceptor,2),round(angle_1,2),round(angle_2,2)))


def aromatic_aromatic_interaction(distance_matrix, dict_vec_pos,pdb_name) :
	list_aromatic = ["TYR","PHE","TRP"]
	dict_interaction = {}
	with open(pdb_name+"_intraprotein_aromatic_aromatic_interaction.txt", "w") as aromatic_file :
		aromatic_file.write("Aromatic-Aromatic Interactions within 4.5 and 7 Angstroms\n")
		aromatic_file.write("Position\tResidue\tChain\tPosition\tResidue\tChain\tD(Centroid-Sulphur)\n")
		for atom_a, value in distance_matrix.items() :
			if atom_a[1] in list_aromatic :
				for atom_b, distance in value.items():
					if atom_a[3] == atom_b[3] :
						if atom_b[1] in list_aromatic and atom_b[2]>atom_a[2] :
							if (atom_a[2], atom_b[2],atom_a[3]) not in dict_interaction :
								dict_interaction[(atom_a[2], atom_b[2],atom_a[3])] = 1
								a1_name = "CG"
								b1_name = "CZ"
								a2_name = "CG"
								b2_name = "CZ"
								if atom_a[1] == "TRP" : #le tryptophane à une configuration differente des aromatiques classiques
									a1_name = "CE3"
									b1_name = "CZ2"
								if atom_b[1] == "TRP" :
									a2_name = "CE3"
									b2_name = "CZ2"
								mid_1_a = dict_vec_pos[a1_name, atom_a[1], atom_a[2], atom_a[3]]#extraction des coordonnées x,y et z
								mid_1_b = dict_vec_pos[b1_name, atom_a[1], atom_a[2], atom_a[3]]
								mid_2_a = dict_vec_pos[a2_name, atom_b[1], atom_b[2], atom_b[3]]
								mid_2_b = dict_vec_pos[b2_name, atom_b[1], atom_b[2], atom_b[3]]
								centroid_1 = np.array([((mid_1_a[0]+mid_1_b[0])/2), ((mid_1_a[1]+mid_1_b[1])/2), ((mid_1_a[2]+mid_1_b[2])/2)])
								#aromatique = hexagone donc on prend le centre d'une des diagonales pour avoir le centroid 						
								centroid_2 = np.array([((mid_2_a[0]+mid_2_b[0])/2), ((mid_2_a[1]+mid_2_b[1])/2), ((mid_2_a[2]+mid_2_b[2])/2)])
								dist_centroid = dist(centroid_1, centroid_2)					
								if (4.5 <= dist_centroid <= 7) :
									aromatic_file.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(atom_a[2],atom_a[1],atom_a[3],atom_b[2],atom_b[1],\
											atom_b[3],round(dist_centroid,2)))


def aromatic_sulfur_interaction(distance_matrix, dict_vec_pos,pdb_name) :
	list_sulfur = ["SG","SD"]
	list_aromatic = ["TYR","PHE","TRP"]
	dict_interaction = {}
	with open(pdb_name+"_intraprotein_aromatic_sulfur_interaction.txt", "w") as sulfur_file :
		sulfur_file.write("Aromatic-Sulphur Interactions within 5.3 Angstroms\n")
		sulfur_file.write("Position\tResidue\tChain\tPosition\tResidue\tChain\tD(Centroid-Sulphur)\n")
		for atom_a, value in distance_matrix.items() :
			if atom_a[1] in list_aromatic :
				for atom_b, distance in value.items() :
					if atom_a[3] == atom_b[3] :
						if atom_b[0] in list_sulfur and atom_b[2]!=atom_a[2] :
							if (atom_a[2], atom_b[2],atom_a[3]) not in dict_interaction :
								dict_interaction[(atom_a[2], atom_b[2],atom_a[3])] = 1
								if atom_a[1] != "TRP" :
									a2_name = "CG"
									b2_name = "CZ"
								if atom_a[1] == "TRP" :
									a2_name = "CE3"
									b2_name = "CZ2"
								mid_2_a = dict_vec_pos[a2_name, atom_a[1], atom_a[2], atom_a[3]]
								mid_2_b = dict_vec_pos[b2_name, atom_a[1], atom_a[2], atom_a[3]]
								sulfur = dict_vec_pos[atom_b]					
								centroid = np.array([((mid_2_a[0]+mid_2_b[0])/2), ((mid_2_a[1]+mid_2_b[1])/2), ((mid_2_a[2]+mid_2_b[2])/2)])
								dist_centroid = dist(sulfur, centroid)
								if dist_centroid < 5.3 :
									sulfur_file.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(atom_a[2],atom_a[1],atom_a[3],atom_b[2],atom_b[1],\
											atom_b[3],round(dist_centroid,2)))


def cation_pi_interaction(distance_matrix, dict_vec_pos,pdb_name) :
	cation = ["ARG","LYS"]
	aromatic = ["TYR","PHE","TRP"]
	dict_cation = {"LYS":"NZ","ARG":"CZ"}
	dict_interaction = {}
	with open(pdb_name+"_intraprotein_cation_pi_interaction.txt", "w") as pi_file :
		pi_file.write("Cation-Pi Interactions within 6 Angstroms\n")
		pi_file.write("Position\tResidue\tChain\tPosition\tResidue\tChain\tD(cation-Pi)\n")
		for atom_a, value in distance_matrix.items() :
			if atom_a[1] in aromatic :
				for atom_b, distance in value.items() :
					if atom_a[3] == atom_b[3] :
						if atom_b[1] in cation and atom_b[2]!=atom_a[2] :
							if (atom_a[2], atom_b[2],atom_a[3]) not in dict_interaction :
								dict_interaction[(atom_a[2], atom_b[2],atom_a[3])] = 1
								if atom_a[1] != "TRP" :
									a2_name = "CG"
									b2_name = "CZ"
								if atom_a[1] == "TRP" :
									a2_name = "CE3"
									b2_name = "CZ2"
								atom_cation = dict_vec_pos[dict_cation[atom_b[1]], atom_b[1], atom_b[2], atom_b[3]]
								mid_2_a = dict_vec_pos[a2_name, atom_a[1], atom_a[2], atom_a[3]]
								mid_2_b = dict_vec_pos[b2_name, atom_a[1], atom_a[2], atom_a[3]]				
								centroid = np.array([((mid_2_a[0]+mid_2_b[0])/2), ((mid_2_a[1]+mid_2_b[1])/2), ((mid_2_a[2]+mid_2_b[2])/2)])
								dist_centroid = dist(centroid, atom_cation)
								if dist_centroid < 6 :
									pi_file.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(atom_a[2],atom_a[1],atom_a[3],atom_b[2],atom_b[1],\
											atom_b[3],round(dist_centroid,2)))


def ionic_interaction(distance_matrix,pdb_name) :
	aa_ion = ["LYS","ARG","HIS","ASP","GLU"]
	cation = ["LYS","ARG","HIS"]
	anion = ["ASP","GLU"]
	atom_ion = ["ND1","NZ","CZ","OD2","OE2"]
	dict_interaction = {}
	with open(pdb_name+"_intraprotein_ionic_interaction.txt", "w") as ion_file :
		ion_file.write("Ionic Interactions within 6 Angstroms\n")
		ion_file.write("Position\tResidue\tChain\tPosition\tResidue\tChain\n")
		for atom_a, value in distance_matrix.items() :
			if atom_a[1] in aa_ion and atom_a[0] in atom_ion :
				for atom_b, distance in value.items() :
					if atom_a[3] == atom_b[3] :
						if atom_b[1] in aa_ion and atom_b[2]>atom_a[2] and distance < 6 and atom_b[0] in atom_ion :
							if (atom_a[2], atom_b[2],atom_a[3]) not in dict_interaction and atom_b[1] != atom_a[1] :
								if (atom_a[1] in cation and atom_b[1] in anion) or (atom_b[1] in cation and atom_a[1] in anion) :
									dict_interaction[(atom_a[2], atom_b[2],atom_a[3])] = 1
									ion_file.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(atom_a[2],atom_a[1],atom_a[3],atom_b[2],atom_b[1],\
										atom_b[3]))


def side_side_hydrogen_bond(distance_matrix, position_vector,pdb_name) :
	list_donor = ["ARG","LYS","TRP"]
	list_acceptor = ["GLU","ASP"]
	list_mixt = ["GLN","ASN","SER","THR","TYR"]
	atom_donor = ["NE","NH1","NH2","ND2","NE2","NZ","OG","OG1","NE1","OH"]
	atom_acceptor = ["OD1","OD2","OE1","OE2","ND1","NE2","OG","OG1","OH"]
	dict_antecedent = {"NE":"CD","NH1":"CZ","NH2":"CZ","ND2":"CG","OD1":"CG","OD2":"CG","NE2":"CD2","OE1":"CD","OE2":"CD","ND1":"CG",\
		"CE1":"CD2","OG":"CB","OG1":"CB","OH":"CZ"}#antecedent de chaque atome accepteur
	dict_hydrogen = {"NE":["HE"],"NH1":["HH11","HH12"],"NH2":["HH21","HH22"],"ND2":["HD21","HD22"],"NE2":["HE21","HE22"],\
		"NZ":["HZ2"],"OG":["HG"],"NE1":["HE1"],"OH":["HH"],"OG1":["HG1"]}#hydrogène lié au donneur
	with open(pdb_name+"_intraprotein_side_side_chain_hydrogen_bonds.txt", "w") as side_side_hbond :
		side_side_hbond.write("\t\tDONOR\t\t\t\t\tACCEPTOR\t\t\tPARAMETERS\t\t\n")
		side_side_hbond.write("POS\tCHAIN\tRES\tATOM\tPOS\tCHAIN\tRES\tATOM\tDd-a\tDh-a\tdon-h-acc\tdon-acc-ant\th-acc-ant\n")
		for atom_a, value in distance_matrix.items() :
			if atom_a[1] in list_donor or atom_a[1] in list_mixt :
				if atom_a[0] in atom_donor :
					for atom_b, distance_a_b in value.items() :
						if atom_a[3] == atom_b[3] :
							if atom_b[1] in list_mixt or atom_b[1] in list_acceptor :
								if atom_b[0] in atom_acceptor and distance_a_b < 3.9 and atom_b[2] != atom_a[2] and  atom_b[2] != atom_a[2]+1 \
									and atom_b[2] != atom_a[2]-1 :
									antecedent = position_vector[dict_antecedent[atom_b[0]], atom_b[1], atom_b[2], atom_b[3]]
									acceptor = position_vector[atom_a]
									donor = position_vector[atom_b]
									for i in range(len(dict_hydrogen[atom_a[0]])) :
										dict_hydrogen_value = dict_hydrogen[atom_a[0]][i]
										h = position_vector[dict_hydrogen_value, atom_a[1], atom_a[2], atom_a[3]]
										dist_h_acceptor = distance_matrix[dict_hydrogen_value, atom_a[1], atom_a[2], atom_a[3]][atom_b[0], atom_b[1], atom_b[2], atom_b[3]]
										angle_1 = angle(donor,h,acceptor)
										angle_2 = angle(donor,acceptor,antecedent)
										angle_3 = angle(h,acceptor,antecedent)
										if dist_h_acceptor < 2.5 and angle_1 > 90 :
											side_side_hbond.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(atom_a[2],\
												atom_a[3],atom_a[1],atom_a[0],atom_b[2],atom_b[3],atom_b[1],\
													atom_b[0],round(distance_a_b,2),round(dist_h_acceptor,2),round(angle_1,2),round(angle_2,2),round(angle_3,2)))


if __name__ == "__main__" :
	filename = get_argument()
	pdb_name = filename[:-6]
	dict_vec_pos = PDB_parser(filename)
	distance_matrix = matrix_intra(dict_vec_pos)
	hydrophobic_interaction(distance_matrix,pdb_name)
	disulphite_bridges(distance_matrix,pdb_name)
	main_main_hydrogen_bonds(distance_matrix, dict_vec_pos,pdb_name)
	aromatic_aromatic_interaction(distance_matrix, dict_vec_pos,pdb_name)
	aromatic_sulfur_interaction(distance_matrix, dict_vec_pos,pdb_name)
	cation_pi_interaction(distance_matrix, dict_vec_pos,pdb_name)
	ionic_interaction(distance_matrix,pdb_name)
	side_side_hydrogen_bond(distance_matrix, dict_vec_pos,pdb_name)