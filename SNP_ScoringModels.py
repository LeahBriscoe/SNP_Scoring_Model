#Leah Briscoe ID:304-145-856
import sys
import csv
import itertools
import math

ids = []
commons = []
mutants = []
Q_m_list = []
neginf = float('-inf')

#####FUNCTION DEFINITIONS#####

# function that converts a list of strings to a list of floats with the base-10 log probability adjusted to computer's log
def convert_to_ln(list):
	list1 = []
	for i in range(0,len(list)): 
		list1.append(list[i] * math.log(10))
	return list1	

# Actual commands here		
reader=csv.reader(sys.stdin,delimiter='\t')
try:
	for id,common,mutant in reader:
		ids.append(int(id))
		commons.append(float(common))
		mutants.append(float(mutant))
		
	commons = convert_to_ln(commons)
	
	mutants = convert_to_ln(mutants)

except IOError as e:
	print 'Operation failed: %s' % e.strerror
	
def convert_to_base10(x):
	x = x/math.log(10)
	return x	
# function that returns neginf if x = 0	or the natural log for any float x
def safe_log(x):
	if x==0:
		return neginf
	else:
		return math.log(x) # natural base e log

#function that adds log probabilities in proper fashion		
def log_sum(x,y):
	#if p(X) >= p(Y)
	if math.exp(x) >= math.exp(y):
		return x + safe_log(1 + math.exp(y-x))	 
	#if p(Y) >= p(X)
	else:
		return y + safe_log(1 + math.exp(x-y))

#FUNCTION: Calculate the log odds ratio for the basic score	
def calc_basic_score():
	mutation_present = 0
	theta1 = 0.2
	mutation_absent = 0
	theta2 = 0
	
	#iterate through elements of commons and mutants together
	for x,y in itertools.izip(commons,mutants):
		
		#assumes 20% frequency of mutation in population
		#log sum over sigma, sum over the various observations
		mutation_present = mutation_present + log_sum(y+safe_log(theta1),x+safe_log(1-theta1))
		
		#assumes 0% frequency of mutation in population
		#log sum over sigma, sum over the various observations
		mutation_absent = mutation_absent + log_sum(y+safe_log(theta2),x+safe_log(1-theta2))
		
	mutation_present = convert_to_base10(mutation_present)
	mutation_absent = convert_to_base10(mutation_absent)
	score = mutation_present - mutation_absent #Divide by the log-probability given Theta=0
	return score


#FUNCTION: Calculate cumulative Q_k. 
def Q_sub_k():
	Q_k = 0
	L_k = 0
	Q_kminus1 = 0
	L_kminus1 = 0
	for i in range(100): # The Riemann sum over 100 slices of the probability. i/100 for increments. 
		
		if i == 0: # for the Q_0 case, Q_k = neginf
			Q_k = neginf
			L_k = log_sum_p_group(0)
			Q_m_list.append(L_k)
		
		else: # for every other Q_k case. 
			L_kminus1 = L_k
			L_k = log_sum_p_group(float(i)/100) # Calculating the logP(obs|theta) for a given theta
			Q_m_list.append(L_k)
			Q_kminus1 = Q_k
			Q_k = log_sum(Q_kminus1,log_sum(L_kminus1,L_k) + safe_log(0.005)) # Trapezoidal height deltatheta/2
		
	return Q_k

#FUNCTION: Calculate logP(obs|theta) per individual. Multiply all probabilities of independent events (add in log space) per person
def logP_per_kappa_per_individual(KAPPA,list_commons,list_mutants):
	sum_for_kappa = 0.0
	for x,y in itertools.izip(list_commons,list_mutants):
		sum_for_kappa += log_sum(y+safe_log(float(KAPPA)/2),x+safe_log(1-float(KAPPA)/2)) 
	return sum_for_kappa

#FUNCTION: Calculation for kappa given theta	
def binomial_distribution_for_kappa(choose,THETA):
	THETA = float(THETA)
	if choose==0: # for 2 choose 0
		return math.pow(1-THETA,2) 
		
	if choose==1: # for 2 choose 1
		return 2*THETA*(1-THETA)
	else: # for 2 choose 2
		return math.pow(THETA,2)
	
	
#FUNCTION: Sum up all the individuals logP(obs|theta) in log space (multiplying probabilities)	
def log_sum_p_group(THETA):	
	THETA = float(THETA)
	sum_of_individuals = 0
	for i in range(max(ids)+1): # For each individual in the dataset. max(list_id)+1 represents the number of individuals
		individual_commons = []
		individual_mutants = []
		subset_sum = 0
		for j in range(0,len(ids)): # For all data points with that individuals ID
			
			if ids[j] == i:
				# Create lists for the probabilities of P(obs|sigma) for a specific individual based on ID
				individual_commons.append(commons[j])
				individual_mutants.append(mutants[j])
				
		 
		sum_over_kappa_individual = 0 
		for k in range(3):	#for kappa values of 0,1, and 2
			# For each value of kappa, calculate logP for individual
			sum_per_kappa = logP_per_kappa_per_individual(k,individual_commons,individual_mutants)
			
			#Binomial distribution added in log space(multiply) by sum of probabilities of observations per individual
			if k == 0:
				sum_over_kappa_individual = safe_log(binomial_distribution_for_kappa(k,THETA)) + sum_per_kappa
			else:
				#Log sum of probabilities across all values of kappa
			 	sum_over_kappa_individual = log_sum(sum_over_kappa_individual, safe_log(binomial_distribution_for_kappa(k,float(THETA))) + sum_per_kappa)
		# Geometric sum of the probabilities logP per individual, equivalent to multiplying probabilities for independent events.
		sum_of_individuals += sum_over_kappa_individual
	
	#Sum of log-probabilities (the product of independent events)
	#Cast the geometric sum of individual probabilities to float.
	return float(sum_of_individuals)
	
#FUNCTION: Calculate advanced score by taking the log odds ratio of Q_k over the logP(obs|theta)
def calc_adv_score():
	
	Q_m = convert_to_base10(Q_sub_k())
	L_zero = convert_to_base10(log_sum_p_group(0))
	
	return Q_m - L_zero	
	
	
#FUNCTION: Calculate the log of the posterior probability density for theta 
def calc_logp():
	for i in range(0,len(Q_m_list)):
		print "%s\t%s\t%s" % ("LOG_P_PHI", str(float(i)/100), str(convert_to_base10(Q_m_list[i])-convert_to_base10(Q_sub_k())))

#MAIN CODE: For any user-provided sequence data file, prints basic scores, advanced score, and table of posterior probability density for each theta
try:
	basic_score = calc_basic_score()
	print "BASIC_SCORE= " + str(basic_score)
	advanced_score = calc_adv_score()
	print "ADVANCED_SCORE= " + str(advanced_score)
	calc_logp()
except IOError as e:
	print 'Operation failed: %s' % e.strerror


	
