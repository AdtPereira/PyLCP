import numpy as np

class MatrixDistances:    

    def Distance(syst):

        # y in a Cartesian Coordinate System (Cables depth)       
        h = syst.pos_y   

        # x in a Cartesian Coordinate System       
        x = syst.pos_x            
        
        d = np.zeros((syst.nph, syst.nph))
        D = np.zeros((syst.nph, syst.nph))
        
        for m in range(syst.nph):
            for n in range(syst.nph):
                # Self Parameters
                if m == n:
                    d[m][n] = syst.rf
                    D[m][n] = np.sqrt(4*h[m]**2 + syst.rf**2)
                
                # Mutual Parameters                    
                else:
                    r = np.sqrt((x[m] - x[n])**2 + (h[m] - h[n])**2)
                    d[m][n] = np.sqrt((h[m]-h[n])**2 + r**2)
                    D[m][n] = np.sqrt((h[m]+h[n])**2 + r**2)

        return d, D
        
class MatrixOper: 
            
    def LInv(A, B):
        """ 
        Return C = A\\B  
        Solve a linear matrix equation, or system of linear scalar equations,
        where [A] and [B] are matrices.

        Inputs:    
        A:
        [[1 2]
        [3 4]]

        B:
        [[5 6]
        [7 8]]

        Results:
        C = A\\B:
        [[-3. -4.]
        [ 4.  5.]]

        D = AB:
        [[ 3. -2.]
        [ 2. -1.]] 
        """
        C = np.linalg.solve(A, B)
        return C

    def RInv(A, B):
        """ 
        Return C = A/B  
        Solve a linear matrix equation, or system of linear scalar equations,
        where [A] and [B] are matrices.

        Inputs:    
        A:
        [[1 2]
        [3 4]]

        B:
        [[5 6]
        [7 8]]

        D = A/B:
        [[ 3. -2.]
        [ 2. -1.]] 
        """
        C = A @ np.linalg.inv(B)
        return C
    
class Trigonometric:

    def coth(x):
        # Calcular tanh(x)
        tanh_x = np.tanh(x)
        
        # Verificar se tanh(x) é próximo de 0
        if abs(tanh_x) < 1e-15:
            return 1.0 / tanh_x
        else:
            # Usar expressão alternativa para evitar overflow
            return (np.exp(x) + np.exp(-x)) / (np.exp(x) - np.exp(-x))
     
    def csch(x):
        exp_x = np.exp(x)
        exp_minus_x = np.exp(-x)
        return 2 / (exp_x - exp_minus_x)   
