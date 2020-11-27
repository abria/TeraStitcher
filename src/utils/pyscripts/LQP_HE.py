"""
Linear Quadratic Programming + Heuristics (LQP_HE)

Takes a text file with the description of a tile positioning optimization problem and 
generates a text file with the solution computed minimizing the sum of squares of errors
with respect to computed displacements, weighted with the corresponding reliabilities. 
The optimization problem has linear constraints and a quadratic objective function.
Solutions must be integers.

The algorithm uses a solver with linear constraints and quadratic objective function and
then applies heuristics to find an integer approximation of the optimal solution.

The format of the input file is
- #rows of tile matrix
- #columns of tile matrix
- #constraints of the optimization problem
- #variables of the optimization problem
- #constraints x #variables matrix (coefficient matrix describing constraints in the form A * X = 0)
- search region used by the displacement computation along V
- default EAST and south diplacements along V (#variables integers in total)
- computed EAST and SOUTH diplacements along V (#variables integers in total)
- EAST and SOUTH diplacements reliabilities along V (#variables integers in total)
- search region used by the displacement computation along H
- default EAST and south diplacements along H (#variables integers in total)
- computed EAST and SOUTH diplacements along H (#variables integers in total)
- EAST and SOUTH diplacements reliabilities along H (#variables integers in total)
- search region used by the displacement computation along D
- default EAST and south diplacements along D (#variables integers in total)
- computed EAST and SOUTH diplacements along D (#variables integers in total)
- EAST and SOUTH diplacements reliabilities along D (#variables integers in total)


USAGE: python LQP_HE.py [-s="path of problem description"] [-d="path to save solution"]

command line parameters:

-s=inputpath	'inputpath' is the complete path of file containing problem description 
				(default: './opt_probl.txt')

-d=outputpath	'outputpath' is the complete path of file containing problem solution 
				(dafault: './opt_sol.txt')
"""

from numpy import array, dot, sort, argsort, where, reshape, append
from os import path, system, remove
from sys import argv
from math import floor, ceil
from scipy.optimize import minimize


def read_params():
   """
   Read parameters from input string and from a file
   """
   
   # Get the input list
   params = argv
   params.pop(0)
   
   #Check if quoting is necessary
   params = check_double_quote(params)
   
   # Declare input file
   fin = read_item(params, '-s=', './opt_probl.txt')
   # Find standard output directory
   fout = read_item(params, '-d=', './opt_sol.txt')
   
   return ( fin, fout )


def read_item ( input_arr, item, default, message = True ) :
   """
   Read the value related to "item" from the list "input_arr" and if no item are present set it to "default".
   Please note: The function convert the output to the same type of "default" variable
   Input:
      input_arr = List of strings from imput command line
      item = The item to search
      default = The default value if no item are present
   Output:
      value = Output value for the selected item
   """
   tmp = check_flag(input_arr, item,True)
   if tmp == None :
      value = default
      if message:
        print ("The value for ",item," was not declared. It will be set to", value, "by default.")
   else :
      if isinstance(default,bool) :
         value = bool(tmp)
      elif isinstance(default,int) :
         value = int(tmp)
      elif isinstance(default,float) :
         value = float(tmp)
      else :
         value = tmp
   return value


def check_command_option ( params, string, delete ) :
   """
   Check if a parameter (string) was beeen declared in the line of commands (params) and return the associated value.
   If delete is true the related string will be deleted
   If string is not present, return None
   Input:
      params = list of parameters from original command line
      string = string to be searched
      delete = Boolean variable to check if the selected string must be deleted after copied in value variable
   Output:
      value = parameter associated to the selected string
   """
   i = 0
   value = None
   size = len(string)
   for line in params :
      tmp = line.find(string)
      if tmp != -1:
         start = tmp + size
         sel_string = line[start:]
         if delete :
            params.pop(i)
         value = sel_string
      i += 1
   return value


def check_double_quote ( inpstring ) :
   """ 
   Check if some strings needs of a double quote (if some space are inside the string, it will need to be inside two double quote). E.g.: --sfmt="TIFF (unstitched, 3D)"
   Input:
      inpstring: input string or array of strings
   Output:
      newstring = new string (or array of strings) corrected by quoting if necessary
   """
   if type(inpstring) == list :
      newstring = []
      for index in inpstring :
         tmp1 = index.find(' ')
         if tmp1 != -1 :
            tmp2 =  index.find('"')
            if tmp2 == -1 :
               dummy = index.find('=')
               if dummy != -1:
                  newstring.append(index[0:dummy+1] + '"' + index[dummy+1:] + '"')
               else :
                  newstring.append('"' + index + '"')
            else :
               newstring.append(index)
         else :
            newstring.append(index)
   else :
      tmp1 = inpstring.find(' ')
      if tmp1 != -1 :
         tmp2 = inpstring.find('"')
         if tmp2 == -1 :
            dummy = inpstring.find('=')
            if dummy != -1 :
               newstring = inpstring[0:dummy+1] + '"' + inpstring[dummy+1:] + '"'
            else :
               newstring = '"' + inpstring + '"'
         else :
            newstring = inpstring
      else :
         newstring = inpstring
   return newstring


def check_flag ( params, string, delete ) :
   """
   Check if a parameter (string) was beeen declared in the line of commands (params) and return the associated value.
   If delete is true the related string will be deleted
   If string is not present, return None
   Input:
      params = list of parameters from original command line
      string = string to be searched
      delete = Boolean variable to check if the selected string must be deleted after copied in value variable
   Output:
      value = parameter associated to the selected string
   """
   i = 0
   value = None
   size = len(string)
   for line in params :
      tmp = line.find(string)
      if tmp != -1 :
         start = tmp + size
         sel_string = line[start:]
         if delete :
            params.pop(i)
         value = sel_string
      i += 1
   return value


def objective_fun ( x, D, R ) :
	return sum(R * (x - D)**2)


def sol_cost ( S, D, R ) :
	"""
	return the cost of solution S with respect to computed displacements D with reliability R
	"""
	mask = array([0]* len(R))
	mask[R>0] = 1
	maxerr = max(mask * abs(S - D))
	meanerr = float(sum(mask * abs(S - D)))
	if sum(mask) <> 0 :
		meanerr /= sum(mask)
	
	return (sum(R * (S - D)**2),maxerr,meanerr)
	
	
def check_constraints ( A, S, complete ) :
	"""
	check if the constraints are all satisified
	"""
	
	ok = True
	
	for i in range(len(complete)) :
		if complete[i] :
			if not (dot(A[i],S) == 0) :
				ok = False
				print '\n'
				print '*** warning *** constraint %d not verified' % (i)
				vars_inds = (where(abs(A[i]) == 1))[0]
				print 'variables involved:', vars_inds
				print 'displacements:', S[vars_inds]
				print
				#programPause = raw_input("Press the <ENTER> key to continue...")
				
	return ok				 


def draw_algo_state ( rows, cols, fixed, complete ) :
	"""
	print to console the progress of the heuristic algoritm:
	which displacements have been defined (fixed) and which squares have been completed (complete)
	"""

	layout = array([[' ']* (2*cols-1)]* (2*rows-1))
	for i in (2*array(range(rows))) :
		for j in (2*array(range(cols))) :
			layout[i][j] = '*'
			
# 	for i in (2*array(range(rows-1))+1) :
# 		for j in (2*array(range(cols-1))+1) :
	d = 0
	for i in (2*array(range(rows))) :
		for j in (2*array(range(cols-1))+1) :
			if fixed[d] :
				layout[i][j] = '-'
			d += 1
	for i in (2*array(range(rows-1))+1) :
		for j in (2*array(range(cols))) :
			if fixed[d] :
				layout[i][j] = '|'
			d += 1
	d = 0
	for i in (2*array(range(rows-1))+1) :
		for j in (2*array(range(cols-1))+1) :
			if complete[d] :
				layout[i][j] = 'x'
			d += 1

	for i in range(2*rows-1) :
		print ''.join(layout[i])
	
	print '\n' + ('='* (2*cols-1)) + '\n'
	

def sol_to_integer ( S, D, R, A, rows = None, cols = None, verbose = False ) :
	"""
	return the integer solution derived with an heuristic from the optimal (non integer) 
	solution S
	D, R, A are the computed displacements, their reliabilities and the constraints matrix
	"""

	n_vars = len(S)
	(c0,max0,mean0) = sol_cost(S,D,R)
	
	intS = array([0]* (n_vars))
	for j in range(n_vars) :
		if S[j] < D[j] :
			intS[j] = int(ceil(S[j]))
		else :
			intS[j] = int(floor(S[j]))
						
	(c1,max1,mean1) = sol_cost(intS,D,R)
	

	########## FIRST HEURISTIC #######################################
	
	intS2 = intS.copy() # solution with first heuristics

	# no displacements are initially fixed
	fixed = array([False]* n_vars)
	complete = array([False]* n_constraints)
	
	if verbose :
		draw_algo_state(rows,cols,fixed,complete)	
			
	# indices of displacements in decreasing order of reliability
	sorted_inds = argsort(-R)
	
	# complete remaining square in order of less reliable displacements not already set
	for j in range(n_vars) :
		if not fixed[sorted_inds[j]] :
			constr_inds = (where(abs(A[:,sorted_inds[j]]) == 1))[0] # indices of constraints containing the displacement considered 
			ind = constr_inds[0]                                    # adjust the first constraint (all other constraints have at least another displacement not fixed)
			err = dot(A[ind],intS2)
			vars_inds = (where(abs(A[ind]) == 1))[0]                # indices of displacements in constraint ind
			temp_inds = vars_inds[fixed[vars_inds] == False]
			modi = temp_inds[R[temp_inds] == min(R[temp_inds])][0]  # index of displacement not set yet
			if A[ind][modi] > 0 :
				intS2[modi] -= err
			else :
				intS2[modi] += err
			fixed[vars_inds] = True
			complete[ind] = True

			if verbose :
				draw_algo_state(rows,cols,fixed,complete)	
			
			# complete all squares that miss just one displacement
			stop = False
			while not stop :
				stop = True
				for i in range(n_constraints) :
					if not complete[i] :
						vars_inds = (where(abs(A[i]) == 1))[0]            # indices of displacements in constraint i
						if len(where(fixed[vars_inds] == True)[0]) == 3 : # just one displacement is missing
							err = dot(A[i],intS2)
							modi = vars_inds[(where(fixed[vars_inds] == False))[0][0]] # index of displacement not set yet
							if A[i][modi] > 0 :
								intS2[modi] -= err
							else :
								intS2[modi] += err
							fixed[modi] = True
							complete[i] = True
							
							if verbose :
								draw_algo_state(rows,cols,fixed,complete)	
			
							stop = False # one more displacement has been set: check for other squares to complete
			
	if verbose :
		draw_algo_state(rows,cols,fixed,complete)	
	
	complete = array([True]* n_constraints)
	if check_constraints(A,intS2,complete) :
		#raise ValueError('')
		(c2,max2,mean2) = sol_cost(intS2,D,R)
	else :
		(c2,max2,mean2) = (float('inf'),float('inf'),float('inf'))


	########## SECOND HEURISTIC #######################################

	intS3 = intS.copy() # solution with second heuristics

	# no displacements are initially fixed
	fixed = array([False]* n_vars)

	complete_bool = array([False]* n_constraints)

	incomplete = set(range(n_constraints))                     # constraints indices to be completed
	complete   = set([])                                       # constraints indices completed
	R_squares  = array([0.0]* n_constraints)                   # maximum displacement reliability of constraints 
	
	for i in range(n_constraints) :
		vars_inds = (where(abs(A[i]) == 1))[0]                 # indices of displacements in constraint i
		R_squares[i] = sort(R[vars_inds])[len(R[vars_inds])-1] # set maximum displacement reliability of constraint i
		
	# indices of constraints in decreasing order of reliability
	sorted_inds = argsort(-R_squares)
	
	# complete the square with highest reliability
	err = dot(A[sorted_inds[0]],intS3)
	vars_inds = (where(abs(A[sorted_inds[0]]) == 1))[0]        # indices of displacements in constraint with highest reliability
	temp_inds = vars_inds[fixed[vars_inds] == False]
	modi = temp_inds[R[temp_inds] == min(R[temp_inds])][0]     # index of displacement not set yet
	if A[sorted_inds[0]][modi] > 0 :
		intS3[modi] -= err
	else :
		intS3[modi] += err
	fixed[vars_inds] = True
	complete |= set([sorted_inds[0]])
	incomplete -= set([sorted_inds[0]])

	if verbose :
		complete_bool[list(complete)] = True
		draw_algo_state(rows,cols,fixed,complete_bool)	
			
	while incomplete <> set([]) :
		for k in range(1,n_constraints) :
			sq = sorted_inds[k]                                    # candidate constraint to be completed
			if ( sq in incomplete) :
				vars_inds = (where(abs(A[sq]) == 1))[0]            # indices of displacements in constraint sq
				adjacent = False				
				for j in vars_inds :
					constr_inds = (where(abs(A[:,j]) == 1))[0]     # indices of constraints containing the displacement considered (and adjacent to constraint sq) 
					if (set(constr_inds) & complete) <> set([]) :  # al least one of the adjacent contraint is complete
						adjacent = True
						break
				if adjacent :
					break
		# sq is the next squares adjacent to already completed squares with the highest reliability
		
		# complete sq
		err = dot(A[sq],intS3)
		temp_inds = vars_inds[fixed[vars_inds] == False]
		modi = temp_inds[R[temp_inds] == min(R[temp_inds])][0]     # index of displacement not set yet
		if A[sq][modi] > 0 :
			intS3[modi] -= err
		else :
			intS3[modi] += err
		fixed[vars_inds] = True
		complete |= set([sq])
		incomplete -= set([sq])
		
		if verbose :
			complete_bool[list(complete)] = True
			draw_algo_state(rows,cols,fixed,complete_bool)	
			
		# complete all squares that miss just one displacement
		stop = False
		while not stop :
			stop = True
			for i in range(n_constraints) :
				if i in incomplete :
					vars_inds = (where(abs(A[i]) == 1))[0]            # indices of displacements in constraint i
					if len(where(fixed[vars_inds] == True)[0]) == 3 : # just one displacement is missing
						err = dot(A[i],intS3)
						modi = vars_inds[(where(fixed[vars_inds] == False))[0][0]] # index of displacement not set yet
						if A[i][modi] > 0 :
							intS3[modi] -= err
						else :
							intS3[modi] += err
						fixed[modi] = True
						complete |= set([i])
						incomplete -= set([i])
						
						if verbose :
							complete_bool[list(complete)] = True
							draw_algo_state(rows,cols,fixed,complete_bool)	
		
						stop = False # one more displacement has been set: check for other squares to complete
	
	complete_bool[list(complete)] = True
	if check_constraints(A,intS3,complete_bool) :
		#raise ValueError('')
		(c3,max3,mean3) = sol_cost(intS3,D,R)
	else :
		(c3,max3,mean3) = (float('inf'),float('inf'),float('inf'))
		

	########## THIRD HEURISTICS #######################################
	
	intS4 = intS.copy() # solution with second heuristics
	
	
	# find the order in which the constraints have to be processed
	# strategy: starts from the most reliable constraint and then add constraints around
	#           the set of constraints already selected which always is a rectangle of 
	#           constraints

	incomplete = set(range(n_constraints))                     # constraints indices to be completed
	complete   = set([])                                       # constraints indices completed
	R_squares  = array([0.0]* n_constraints)                   # maximum displacement reliability of constraints 
	
	# constraint reliability is the maximum reliability of its displacements
	for i in range(n_constraints) :
		vars_inds = (where(abs(A[i]) == 1))[0]                 # indices of displacements in constraint i
		R_squares[i] = sort(R[vars_inds])[len(R[vars_inds])-1] # set maximum displacement reliability of constraint i
		
	# indices of constraints in decreasing order of reliability
	sorted_inds = argsort(-R_squares)
	
	# indices of the ost reliable constraint
	i0 = sorted_inds[0]/(cols-1)
	i1 = sorted_inds[0]/(cols-1)
	j0 = sorted_inds[0]%(cols-1)
	j1 = sorted_inds[0]%(cols-1)
	complete |= set([sorted_inds[0]])
	incomplete -= set([sorted_inds[0]])
	constr_order = array([sorted_inds[0]]) # initially constraint are sorted by reliability
	while incomplete <> set([]) :
		constr_inds = array([], dtype=int)
		up    = False
		down  = False
		left  = False
		right = False	
		if i0 > 0 : 
			# add a row of constraint above
			f = (i0-1)*(cols-1) + j0
			l = f + (j1-j0+1)
			constr_inds = append(constr_inds,array(range(f,l)))
			up = True
		if i1 < (rows-2) :
			# add a row of constraint below
			f = (i1+1)*(cols-1) + j0
			l = f + (j1-j0+1)
			constr_inds = append(constr_inds,array(range(f,l)))
			down = True
		if j0 > 0 :
			# add a row of constraint on the left
			f = i0*(cols-1) + j0 - 1
			l = (i1+1)*(cols-1) + j0 - 1
			constr_inds = append(constr_inds,array(range(f,l,(cols-1))))
			left = True
		if j1 < (cols-2) :
			# add a row of constraint on the right
			f = i0*(cols-1) + j1 + 1
			l = (i1+1)*(cols-1) + j1 + 1
			constr_inds = append(constr_inds,array(range(f,l,(cols-1))))
			right = True
		# add corner constraints if any
		if up and left :
			constr_inds = append(constr_inds,array([(i0-1)*(cols-1) + j0 - 1]))
		if up and right :
			constr_inds = append(constr_inds,array([(i0-1)*(cols-1) + j1 + 1]))
		if down and left :
			constr_inds = append(constr_inds,array([(i1+1)*(cols-1) + j0 - 1]))
		if down and right :
			constr_inds = append(constr_inds,array([(i1+1)*(cols-1) + j1 + 1]))
		# extends the rectangle of selected constraints
		if up :
			i0 -= 1
		if down :
			i1 += 1
		if left :
			j0 -= 1
		if right :
			j1 += 1
		complete |= set(constr_inds)
		incomplete -= set(constr_inds)
		constr_order = append(constr_order,constr_inds[argsort(-R_squares[constr_inds])])


	# process the contraints

	# no displacements are initially fixed
	fixed = array([False]* n_vars)

	complete_bool = array([False]* n_constraints)

	incomplete = set(range(n_constraints))                     # constraints indices to be completed
	complete   = set([])                                       # constraints indices completed

	# complete the square with highest reliability
	err = dot(A[constr_order[0]],intS4)
	vars_inds = (where(abs(A[sorted_inds[0]]) == 1))[0]        # indices of displacements in constraint with highest reliability
	temp_inds = vars_inds[fixed[vars_inds] == False]
	modi = temp_inds[R[temp_inds] == min(R[temp_inds])][0]     # index of displacement not set yet
	if A[constr_order[0]][modi] > 0 :
		intS4[modi] -= err
	else :
		intS4[modi] += err
	fixed[vars_inds] = True
	complete |= set([constr_order[0]])
	incomplete -= set([constr_order[0]])

	if verbose :
		complete_bool[list(complete)] = True
		draw_algo_state(rows,cols,fixed,complete_bool)	
			
	while incomplete <> set([]) :
		for k in range(1,n_constraints) :
			sq = constr_order[k]                                    # candidate constraint to be completed
			if ( sq in incomplete) :
				vars_inds = (where(abs(A[sq]) == 1))[0]            # indices of displacements in constraint sq
				break
		# sq is the next squares adjacent to already completed squares with the highest reliability
		
		# complete sq
		err = dot(A[sq],intS4)
		temp_inds = vars_inds[fixed[vars_inds] == False]
		modi = temp_inds[R[temp_inds] == min(R[temp_inds])][0]     # index of displacement not set yet
		if A[sq][modi] > 0 :
			intS4[modi] -= err
		else :
			intS4[modi] += err
		fixed[vars_inds] = True
		complete |= set([sq])
		incomplete -= set([sq])
		
		if verbose :
			complete_bool[list(complete)] = True
			draw_algo_state(rows,cols,fixed,complete_bool)	
			
		# complete all squares that miss just one displacement
		stop = False
		while not stop :
			stop = True
			for i in range(n_constraints) :
				if i in incomplete :
					vars_inds = (where(abs(A[i]) == 1))[0]            # indices of displacements in constraint i
					if len(where(fixed[vars_inds] == True)[0]) == 3 : # just one displacement is missing
						err = dot(A[i],intS4)
						modi = vars_inds[(where(fixed[vars_inds] == False))[0][0]] # index of displacement not set yet
						if A[i][modi] > 0 :
							intS4[modi] -= err
						else :
							intS4[modi] += err
						fixed[modi] = True
						complete |= set([i])
						incomplete -= set([i])
						
						if verbose :
							complete_bool[list(complete)] = True
							draw_algo_state(rows,cols,fixed,complete_bool)	
		
						stop = False # one more displacement has been set: check for other squares to complete
	
	complete_bool[list(complete)] = True
	if check_constraints(A,intS4,complete_bool) :
		#raise ValueError('')
		(c4,max4,mean4) = sol_cost(intS4,D,R)
	else :
		(c4,max4,mean4) = (float('inf'),float('inf'),float('inf'))

	########## FINAL RESULT #######################################
	
	if c3 < c2 :
		if c3 < c4 :
			#print 'The second heuristics gives the best result'
			intS2 = intS3
			(c2,max2,mean2) = (c3,max3,mean3)
		else :
			#print 'The third heuristics gives the best result'
			intS2 = intS4
			(c2,max2,mean2) = (c4,max4,mean4)
	elif c4 < c2 :
		#print 'The third heuristics gives the best result'
		intS2 = intS4
		(c2,max2,mean2) = (c4,max4,mean4)
	else :
		#print 'The first heuristics gives the best result'
		pass
				
	return (c0,max0,mean0,c1,max1,mean1,intS2,c2,max2,mean2)
		

if __name__ == '__main__' :

	(fin,fout) = read_params()
		
	din = open(fin)
	dat = din.read().split()
	din.close()
		
	rows = int(dat[0])
	cols = int(dat[1])
		
	n_constraints = int(dat[2])
	n_vars = int(dat[3])
		
	A = array([0]* (n_constraints * n_vars))
	A = reshape(A,(n_constraints,n_vars))
	h = 4
	for i in range(n_constraints) :
		for j in range(n_vars) :
			A[i][j] = int(dat[h])
			h += 1

	dout = open(fout,'w')

	for d in ['V', 'H', 'D'] :

		delay = int(dat[h])
		h += 1
		default = array([0]* (n_vars))
		for j in range(n_vars) :
			default[j] = int(dat[h])
			h += 1								
		D = array([0]* (n_vars))
		for j in range(n_vars) :
			D[j] = int(dat[h])
			h += 1								
		R = array([0.0]* (n_vars))
		for j in range(n_vars) :
			R[j] = float(dat[h])
			h += 1								
		
		bnds = [(default[i]-delay,default[i]+delay) for i in range(n_vars)]
		constr = ({'type': 'eq', 'fun' : lambda x: dot(A,x)})
			
		res = minimize(objective_fun,default,args=(D,R),method='SLSQP',jac=False,bounds=bnds,constraints=constr)
	
		S = res.x
				
		try :
			(c0,max0,mean0,c1,max1,mean1,intS2,c2,max2,mean2) = sol_to_integer(S,D,R,A,rows,cols)
		except :
			(c0,max0,mean0,c1,max1,mean1,intS2,c2,max2,mean2) = (0,0,0,0,0,0,default,0,0,0)

		for j in range(n_vars) :
			dout.write(str(intS2[j]))
			dout.write(' ')
		dout.write('\n')

	dout.close()
