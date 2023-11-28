import sympy as sp
import numpy as np

a,b,c=sp.symbols('a b c');
test=np.array([[a**2,a+b],[a*c+b,b/c]]);

w = np.exp((2*np.pi*1j)/3)
A1 = np.array([[1, 0, 0], [0, w, 0], [0, 0, w**2]])
B1 = np.array([[0, 0, 1], [1, 0, 0], [0, 1, 0]])

B = np.array([[0, a], [a, 0]])
bdot = np.dot(B,B)
print(bdot)

# Calculate the product AB
B2 = np.dot(B1, B1)
A2 = np.dot(A1,A1)
A0 = np.dot(A2,A2)
B22 = np.dot(B2,B2)
