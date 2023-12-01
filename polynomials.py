import math as m

class polynomial:
    #constructor
    def __init__(self, coefficients):   
        self.coefficients = coefficients
        self.degree = len(coefficients)-1
    
    #evaluated the polynomial at x
    def evaluate(self, x: float) -> float: 
        result = 0
        for i in range(len(self.coefficients)-1, -1, -1):
            result = result*x+self.coefficients[i]
        return result
    
    #returns the derivative of the polynomial
    def derivative(self): 
        coeff = []
        for i in range(1,len(self.coefficients)):
            coeff.append(i*self.coefficients[i])
        return polynomial(coeff)
    
    #how to make sense of print(polynomial)
    def __str__(self): 
        return str(self.coefficients)
    
    #adding two polynomials
    def __add__(self, poly): 
        new_coeff = []
        for i in range(max(self.degree, poly.degree)+1):
            temp = 0
            if(i <= self.degree):
                temp += self.coefficients[i]
            if(i <= poly.degree):
                temp += poly.coefficients[i]
            new_coeff.append(temp)
        return polynomial(new_coeff)
    
    #subtracting two polynomials
    def __sub__(self, poly):    
        new_coeff = []
        for i in range(max(self.degree, poly.degree)+1):
            temp = 0
            if(i <= self.degree):
                temp += self.coefficients[i]
            if(i <= poly.degree):
                temp -= poly.coefficients[i]
            new_coeff.append(temp)
        return polynomial(new_coeff)
    
    
    #multiplying two polynomials
    def __mul__(self, poly): 
        result = [0]*(self.degree+poly.degree+2)
        for i in range(poly.degree+1):
            if(poly.coefficients[i] != 0):
                temp = [0]*i + [poly.coefficients[i]*x for x in self.coefficients]
                p = 0
                for x in temp:
                    result[p] += x
                    p += 1

        while(result[-1] == 0):
            result.pop()
        
        return polynomial(result)


#return lagrange interpolating polynomial
def lagrange_polynomial(x_list, y_list): 
    length = len(x_list)
    result = polynomial([0]*(length))
    for i in range(length):
        new_poly = polynomial([1] + [0]*(length-1))
        factor = 1
        for j in range(length):
            if( i == j):
                continue
            else:
                curr = polynomial([-x_list[j], 1])
                new_poly = new_poly*curr
                factor = (x_list[i]-x_list[j])*factor

        new_poly = new_poly*(polynomial([(1/factor)*y_list[i]]))
        result = result + new_poly
    
    return result


# takes in x-coords and y-coords and tha value 'x' at which you want to approximate
def forward_difference(x_list, y_list, x): 
    coeffs = [] # a list [f(x0), delta f(x0), delta^2 f(x0)]
    temp = y_list.copy()
    n = len(y_list)
    for i in range(n):
        coeffs.append(temp[0])
        for i in range(len(temp)-1):
            temp[i] = temp[i+1]-temp[i]
        temp.pop()
    
    result = 0
    h = x_list[1]-x_list[0]
    s = (x-x_list[0])/h
    for i in range(0, n):
        result +=m.comb(s, i)*coeffs[i]
    
    return result

def backward_difference(x_list, y_list, x):
    return forward_difference(x_list[::-1], y_list[::-1], x)


#where n is the nth legendre polynomial (n > 0)
def legendre_polynomial(n): 
    poly = polynomial([-1, 0, 1])
    for i in range(n-1):
        poly = poly*poly
    for i in range(n):
        poly = poly.derivative()
    k = (2**n)*(m.factorial(n))
    coeff = poly.coefficients
    for i in range(poly.degree+1):
        coeff[i] = coeff[i]/k
    return polynomial(coeff)





