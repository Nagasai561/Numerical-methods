import copy

#matrices are represented using generic python lists

def add_mat(matrix1, matrix2):  #returns a new matrix
    if((len(matrix1) != len(matrix2))): #checking if number of rows aren't same
        raise Exception("Can't add two incompatible matrices")
    
    elif(all(isinstance(x, list) for x in matrix1)):  #checking if matrix1 is a 2d matrix
        if(not all(isinstance(x, list) for x in matrix2)):    #Then matrix2 must also be a 2d matrix
            raise Exception("Can't add two incompatible matrices")
        
        elif((len(matrix2[0]) != len(matrix1[0]))):     #checking if number of columns are same or not
            raise Exception("Can't add two incompatible matrices")
        
        m = len(matrix1)
        n = len(matrix1[0])
        new_mat = [[0]*n]*m
        for i in range(m):
            for j in range(n):
                new_mat[i][j] = matrix1[i][j] + matrix2[i][j]
        return new_mat  
    else:
        m = len(matrix1)
        new_mat = [0]*m
        for i in range(m):
            new_mat[i] = matrix1[i]+matrix2[i]
        return new_mat
    
def scale_mat(matrix, scaler):  #modifies the existing one
    m = len(matrix)
    if(all(isinstance(x, list) for x in matrix)):
        for i in range(m):
            for j in range(len(matrix[0])):
                matrix[i][j] *= scaler
    else:
        for i in range(m):
            matrix[i] *= scaler

def row_exchange(matrix, row1, row2):   #modifies the existing one
    m = len(matrix)
    if(not all(isinstance(x, list) for x in matrix)):
        matrix[row1], matrix[row2] = matrix[row2], matrix[row1]
        return
    if((not (0 <= row1 < m)) or (not (0 <= row2 < m))):
        raise Exception("row index are out of bounds")
    for j in range(len(matrix[0])):
        matrix[row1][j], matrix[row2][j] = matrix[row2][j], matrix[row1][j]
    
def row_transform(matrix, row_to_change, row_to_use, scaler):   #modifies the existing one
    m = len(matrix)
    if(not all(isinstance(x, list) for x in matrix)):
        matrix[row_to_change] += scaler*matrix[row_to_use]
        return
    if((not (0 <= row_to_change < m)) or (not (0 <= row_to_use < m))):
        raise Exception("row index are out of bounds")
    
    for j in range(len(matrix[0])):
        matrix[row_to_change][j] += scaler*matrix[row_to_use][j]
    
def gauss_jacobi(matrix, column_mat, initial_guess, max_iter, tolerance):   #make sure that a_ii is not zero!
    m = len(matrix)
    n = len(matrix[0])
    for iteration in range(max_iter):
        new_guess = [0]*m
        for i in range(m):
            new_guess[i] = column_mat[i]
            for k in range(m):
                if(k == i):
                    continue
                new_guess[i] -= initial_guess[k]*matrix[i][k]
            new_guess[i] = new_guess[i]/matrix[i][i]

        #check if our improved sol is satisfying the tolerance condition
        old = max([abs(initial_guess[i]-new_guess[i]) for i in range(m)])
        new = max(abs(x) for x in new_guess)
        if((new != 0) and (old/new) <= tolerance):
            return new_guess
        
        initial_guess = new_guess

    return initial_guess

def gauss_seidal(matrix, column_mat, initial_guess, max_iter):
    m = len(matrix)
    n = len(matrix[0])
    for iteration in range(max_iter):
        for i in range(m):
            initial_guess[i] = column_mat[i]
            for k in range(m):
                if(k == i):
                    continue
                initial_guess[i] -= initial_guess[k]*matrix[i][k]
            initial_guess[i] = initial_guess[i]/matrix[i][i]

    return initial_guess

    
def gauss_elimination(matrix, column_mat):
    mat = copy.deepcopy(matrix)
    col_mat = column_mat.copy()
    m = len(mat)
    n = len(mat[0])

    for i in range(m):
        if(mat[i][i] == 0):     #if a_ii is zero, we will exchange this row with an appropriate row, which avoids this problem.
            first = -1           #if we can't find such a row, then it means it doesn't have unique solution
            for j in range(i+1, m):
                if(mat[j][i] != 0):
                    first = j
                    break
            if(first == -1):
                raise Exception("Unique solution doesn't exist")
            row_exchange(mat, i, first)
            row_exchange(col_mat, i, first)
        
        for k in range(i+1, m):
            fac = mat[k][i]/mat[i][i]
            row_transform(mat, k, i, -fac)
            row_transform(col_mat, k, i, -fac)
        
    #we have turned the matrix in echelon form, now it's time to substitute values
    sol = [0]*m
    for i in range(m-1, -1, -1):
        dot_prod = sum(sol[k]*mat[i][k] for k in range(m))
        sol[i] = (col_mat[i]-dot_prod)/mat[i][i]

    return sol


            
