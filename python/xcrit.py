from scipy.stats import poisson

# print 'Calculate  =  Reads/Channel Threshold for Potential Noise Given Lambda'

lambda_value_txt = open('/<your/working/directory>/06_poisson_calculation/lambda_value.txt', 'r')
lambda_value = lambda_value_txt.read().splitlines()[8]

print 'Lambda Value (06_poisson_calculation/lambda_value.txt)'
print lambda_value

p = 0.05 # User Defined
print 'P Value'
print p

x_crit = poisson.ppf(1-p,float(lambda_value))
 
print 'Critical Read/Channel Threshold'
print x_crit

