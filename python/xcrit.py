from scipy.stats import poisson
import sys

lambda_value_txt = open(sys.argv[1], 'r')
lambda_value = lambda_value_txt.read().splitlines()[8]

print 'Lambda Value (06_poisson_calculation/lambda_value.txt)'
print lambda_value

p = float(sys.argv[2]) # User Defined
print 'P Value'
print p 

x_crit = poisson.ppf(1-p,float(lambda_value))
 
print 'Critical Read/Channel Threshold'
print x_crit
