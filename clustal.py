def matrix(matrix_file):
	blosum62 = {}
	aa = []
	with open (matrix_file,'r') as file :
		first_line = file.readline().strip()
		aa = first_line.split()
		for line in file :
			row_values = list(map(int, line.strip().split()[1:]))
			row_acid = line.split()[0]
			for col_acid, value in zip(aa, row_values):
				blosum62[(row_acid, col_acid)] = value
	return blosum62


