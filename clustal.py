def matrix(matrix_file):
	blosum62 = {}
	aa = []
	with open (matrix_file,'r') as file :
		aa = file.readlines().strip().split()
		for line in file :
			row_values = list(map(int, line.strip().split()[1:]))
        	row_acid = line.strip().split()[0]
        	for col_acid, value in zip(amino_acids, row_values):
            	blosum62[(row_acid, col_acid)] = value


