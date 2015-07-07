import random
import LDC3

# Generate random alpha values
alphas=[0 for i in range(3)]
alphas=[random.uniform(0.0,1.0) for i in range(3)]
print 'alpha_h = ',alphas[0]
print 'alpha_r = ',alphas[1]
print 'alpha_t = ',alphas[2]

# Calculate the corresponding LDCs
c=LDC3.forward(alphas)
print 'c_2 = ',c[0]
print 'c_3 = ',c[1]
print 'c_4 = ',c[2]

# Test if the LDCs satisfy the seven analytic criteria
passed=LDC3.criteriatest(0,c)
if passed == 1:
   print 'LDCs satisfy the 7 analytic criteria => physically valid'
else:
   print 'LDCs violate one or more of the 7 analytic criteria'

# Invert LDCs back to alphas
alphas_inv=LDC3.inverse(c)
print 'alpha_h = ',alphas_inv[0]
print 'alpha_r = ',alphas_inv[1]
print 'alpha_t = ',alphas_inv[2]
