import os
import re
from rdkit import Chem
from rdkit.Chem import Draw, rdDetermineBonds

class MolBuilder():
	'''
	Class to parse the molecule information from .log file and converted to xyz, smile and 2D files.

	Parameters
	----------
	'''
	def __init__(self, output_path, charge=0):
		self.output_xyz = os.path.join(output_path, 'xyz_molecules')
		self.output_img = os.path.join(output_path, 'img_molecules')
		self.charge = charge

		if not os.path.isdir(self.output_xyz):
			os.makedirs(self.output_xyz)

		if not os.path.isdir(self.output_img):
			os.makedirs(self.output_img)

		self.code = {"1" : "H", "2" : "He", "3" : "Li", "4" : "Be", "5" : "B",
		"6"  : "C", "7"  : "N", "8"  : "O",  "9" : "F", "10" : "Ne",
		"11" : "Na" , "12" : "Mg" , "13" : "Al" , "14" : "Si" , "15" : "P",
		"16" : "S"  , "17" : "Cl" , "18" : "Ar" , "19" : "K"  , "20" : "Ca",
		"21" : "Sc" , "22" : "Ti" , "23" : "V"  , "24" : "Cr" , "25" : "Mn",
		"26" : "Fe" , "27" : "Co" , "28" : "Ni" , "29" : "Cu" , "30" : "Zn",
		"31" : "Ga" , "32" : "Ge" , "33" : "As" , "34" : "Se" , "35" : "Br",
		"36" : "Kr" , "37" : "Rb" , "38" : "Sr" , "39" : "Y"  , "40" : "Zr",
		"41" : "Nb" , "42" : "Mo" , "43" : "Tc" , "44" : "Ru" , "45" : "Rh",
		"46" : "Pd" , "47" : "Ag" , "48" : "Cd" , "49" : "In" , "50" : "Sn",
		"51" : "Sb" , "52" : "Te" , "53" : "I"  , "54" : "Xe" , "55" : "Cs",
		"56" : "Ba" , "57" : "La" , "58" : "Ce" , "59" : "Pr" , "60" : "Nd",
		"61" : "Pm" , "62" : "Sm" , "63" : "Eu" , "64" : "Gd" , "65" : "Tb",
		"66" : "Dy" , "67" : "Ho" , "68" : "Er" , "69" : "Tm" , "70" : "Yb",
		"71" : "Lu" , "72" : "Hf" , "73" : "Ta" , "74" : "W"  , "75" : "Re",
		"76" : "Os" , "77" : "Ir" , "78" : "Pt" , "79" : "Au" , "80" : "Hg",
		"81" : "Tl" , "82" : "Pb" , "83" : "Bi" , "84" : "Po" , "85" : "At",
		"86" : "Rn" , "87" : "Fr" , "88" : "Ra" , "89" : "Ac" , "90" : "Th",
		"91" : "Pa" , "92" : "U"  , "93" : "Np" , "94" : "Pu" , "95" : "Am",
		"96" : "Cm" , "97" : "Bk" , "98" : "Cf" , "99" : "Es" ,"100" : "Fm",
		"101": "Md" ,"102" : "No" ,"103" : "Lr" ,"104" : "Rf" ,"105" : "Db",
		"106": "Sg" ,"107" : "Bh" ,"108" : "Hs" ,"109" : "Mt" ,"110" : "Ds",
		"111": "Rg" ,"112" : "Uub","113" : "Uut","114" : "Uuq","115" : "Uup",
		"116": "Uuh","117" : "Uus","118" : "Uuo"}

	def __get_lines(self, file):
		'''
		Extract the file lines from log file
		'''

		with open(file, 'r') as log:
			lines = log.readlines()

		return lines

	def __find_StandardOrientation(self, lines):
		for i, line in enumerate(lines):
			if 'Standard orientation' in line:
				index = i + 5
				break

		for i, line in enumerate(lines):
			if '-----------' in line:
				end = i
				if end > index:
					break

		orientation_matrix = lines[index:end]

		return orientation_matrix
	
	def __correct_matrix(self, matrix):
		for i, line in enumerate(matrix):
			edited_line = (re.sub(r'\s+', ' ',line)).split(' ')

			if edited_line[0] == '':
				edited_line.remove('')

			matrix[i] = edited_line

		return matrix

	def __gen_xyz(self, matrix):
		NAtoms = len(matrix)

		xyz_file = f'{NAtoms}\n\n'
		
		for line in matrix:
			Atom = self.code[line[1]]
			x = line[3]
			y = line[4]
			z = line[5]

			xyz_line = f'{Atom}\t{x}\t{y}\t{z}'
			
			if line != matrix[-1]:
				xyz_line += '\n'

			xyz_file += xyz_line

		return xyz_file
	
	def __write_xyz(self, xyz_file, output_file, force=True):
		if not os.path.isfile(output_file):
			if not force:
				return 'File Already Exist'
		
		with open(output_file, 'w') as xyz:
			xyz.writelines(xyz_file)

	def build_xyz(self, file, save=True, force=True):
		lines = self.__get_lines(file)
		matrix = self.__find_StandardOrientation(lines)
		matrix_corrected = self.__correct_matrix(matrix)

		xyz_file = self.__gen_xyz(matrix_corrected)

		file_name = file.split('/')[-1]
		xyz_name = file_name.replace('.log', '.xyz')

		output_file = os.path.join(self.output_xyz, xyz_name)

		if save:
			self.__write_xyz(xyz_file, output_file, force=force)

	def __get_mol(self, xyz_file):
		mol = Chem.MolFromXYZFile(xyz_file)
		conn_mol = Chem.Mol(mol)
		rdDetermineBonds.DetermineConnectivity(conn_mol, useHueckel=True, charge=self.charge)

		try:
			rdDetermineBonds.DetermineBonds(conn_mol, useHueckel=True, charge=self.charge)
		except:
			pass

		return conn_mol

	def get_smiles(self, xyz_file):
		conn_mol = self.__get_mol(xyz_file)

		if not conn_mol:
			return None

		smiles = Chem.MolToSmiles(conn_mol)
		
		return smiles

	def get_image(self, xyz_file, show=False, save=True):
		smiles = self.get_smiles(xyz_file)
		conn_mol = Chem.MolFromSmiles(smiles)
		
		if not conn_mol:
			return None

		if show:
			img = Draw.MolToImage(conn_mol)
			img.show('img')

		if save:
			file_name = xyz_file.split('/')[-1]
			jpg_name = file_name.replace('.xyz', '.jpg')
			output_file = os.path.join(self.output_img, jpg_name)
			
			Draw.MolToImageFile(conn_mol, output_file)