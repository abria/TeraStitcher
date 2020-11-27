#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
This program uses a dispatcher/worker/launcher approach to consume a queue 
of elaborations using teraconverter
Copyright (c) 2016:
Massimiliano Guarrasi (1), Giulio Iannello (2), Alessandro Bria (2)
(1): CINECA
(2): University Campus Bio-Medico of Rome
The program was made in the framework of the HUMAN BRAIN PROJECT.
All rights reserved.


RELEASE 3.2.3


EXAMPLE of usage (X is the major version, Y is the minor version, Z is the patch):

For align step:
mpirun -np XX python ParastitcherX.Y.Z.py -2 --projin=xml_import_file --projout=xml_displcomp_file [--sV=VV] [--sH=HH] [--sD=DD] [--imin_channel=C] [ ... ] 

where:
- XX is the desided level of parallelism plus 1 (for the dispatcher process)
- VV, HH, DD are the half size of the NCC map along V, H, and D directions, respectively
- C is the input channel to be used for align computation


For fusion step:
mpirun -np XX python ParastitcherX.Y.Z.py -6 --projin=xml_import_file --volout=destination_folder --volout_plugin=format_string [--slicewidth=WWW] [--sliceheight=HHH] [--slicedepth=DDD] [--resolutions=RRR] [ ... ] 

where:
- format_string is one the formats: "TIFF (series, 2D)", "TIFF, tiled, 2D", "TIFF, tiled, 3D", "TIFF, tiled, 4D", 
- DDD, HHH, WWW are the values used to partition the image for parallel execution
- RRR are the requested resolutions (according to the convention used by teraconverter)
See teraconverter documentation for more details


*******************************
*        Change Log           *
*******************************

2020-07-20. Giulio. @CHANGED the names of some finctions: introduced the terminology dispatcher/workers/launchers
2019-08-01. Giulio. @FIXED a bug in the distributed termination algorithm in 'dispatcher_step2'
2018-09-20. Giulio  @CHANGED adapted to run also on Python 3 interpreters: changed <> to !=, use parentheses in print commands, converted 'keys' and 'values' of dictionaries to lists
2018-09-05. Giulio. @CHSNGED on non-Windows platforms 'prefix' is automatically switched to './' if executables are not in the system path
2018-08-16. Giulio. @CHANGED command line interface: parameters for step 6 are the same than the sequential implementation
2018-08-16. Giulio. @ADDED debug control
2018-08-07. Giulio. @CREATED from parastitcher2.0.3.py and paraconverter2.3.2.py

*******************************

terastitcher -2 --projin=/Users/iannello/Home/Windows/myTeraStitcher/TestData/Ailey/blending/test_00_01_02_03/xml_import_org.xml --projout=/Users/iannello/Home/Windows/myTeraStitcher/TestData/Ailey/blending/test_00_01_02_03/xml_displcomp_seq.xml

mpirun -np 3 python /Users/iannello/Home/Windows/paratools/parastitcher2.0.3.py -2 --projin=/Users/iannello/Home/Windows/myTeraStitcher/TestData/Ailey/blending/test_00_01_02_03/xml_import_org.xml --projout=/Users/iannello/Home/Windows/myTeraStitcher/TestData/Ailey/blending/test_00_01_02_03/xml_displcomp_par2.xml
mpirun -np 3 python /Users/iannello/Home/Windows/paratools/Parastitcher3.0.0.py -2 --projin=/Users/iannello/Home/Windows/myTeraStitcher/TestData/Ailey/blending/test_00_01_02_03/xml_import_org.xml --projout=/Users/iannello/Home/Windows/myTeraStitcher/TestData/Ailey/blending/test_00_01_02_03/xml_displcomp_par2.xml

teraconverter --sfmt="TIFF (unstitched, 3D)" -s=/Users/iannello/Home/Windows/myTeraStitcher/TestData/Ailey/blending/test_00_01_02_03/xml_merging.xml --dfmt="TIFF (series, 2D)" -d=/Users/iannello/Home/Windows/myTeraStitcher/TestData/temp/result_p1 --resolutions=012 --depth=256 --width=256 --height=256

mpirun -np 3 python /Users/iannello/Home/Windows/paratools/paraconverter2.3.2.py --sfmt="TIFF (unstitched, 3D)" -s=/Users/iannello/Home/Windows/myTeraStitcher/TestData/Ailey/blending/test_00_01_02_03/xml_merging.xml --dfmt="TIFF (tiled, 3D)" -d=/Users/iannello/Home/Windows/myTeraStitcher/TestData/temp/result_p1 --resolutions=012 --depth=256 --width=256 --height=256
mpirun -np 3 python /Users/iannello/Home/Windows/paratools/Parastitcher3.0.0.py -6 --sfmt="TIFF (unstitched, 3D)" -s=/Users/iannello/Home/Windows/myTeraStitcher/TestData/Ailey/blending/test_00_01_02_03/xml_merging.xml --dfmt="TIFF (tiled, 3D)" -d=/Users/iannello/Home/Windows/myTeraStitcher/TestData/temp/result_p1 --resolutions=012 --depth=256 --width=256 --height=256

"""


import os
import sys
import shutil
import time
import datetime
import operator
import platform

from math import *
from glob import glob
from mpi4py import MPI
from collections import deque
from subprocess import *

import os.path
import pickle


"""
The script needs to find the executables of terastitcher (step 2, Align) or teraconverter (step 6, Merge)
Current directory is first searched for the executables, if they are not there the system path is searched
If the excutables are not in the system path an error is raised
"""
prefix = ''
#prefix = './'


"""
Debug level is contrlled by an integer value. Any level includes the previous one
0: no debug additional information
1: the output of instances of executables are saved in files named output_XX.out
   where XX is the ID of the instance that generated it 
"""
debug_level = 0


"""
*************************
* PARAMETERS            *  
*************************
"""
default_tile_size = 256
"""
*************************
"""


"""
*************************
* SUSPEND/RESUME option *  WARNING: this feature has not fully debugged and occcasionally fails
*************************

invert comments if you want to enable suspend/resume behavior
suspend/resume may reduce achievable speedup since it may introduce delays when status is saved
after every instance of TeraConverter terminates
to reduce this effect, status can be saved in a permament RAM disk of limited size (tens of Kbytes)
"""
resume_status_fname = 'para_resume_status.bin'

#suspend_resume_enabled = True  # suspend/resume mechanism enabled
suspend_resume_enabled = False # suspend/resume mechanism disabled

"""
if a special drive (e.g. a RAM disk) should be used to save the status for sispend/resume 
the variable 'save_status_prefix' should be assinged with the path where the status should be 
saved 
if 'save_status_prefix' is the empty string the status is saved in the destination directory
(i.e. where TeraConverter instances store their own suspend/resume status
"""
save_status_prefix = ''
#save_status_prefix = '/media/giannello/ramdisk/'
"""
*************************
"""


###############################################
# functions
###############################################

def partition ( m, n, N ) : 
	"""
	return the number of partitions along V and H, respectively that are optimal to 
	partition a block of size m_V x n_H in at least P sub-blocks
	
	m: block size along V
	n: block size along H
	N:   number of required partitions

	return: 
	p_m, p_n: the  number of partitions along V and H, respectively 
	
	PRE: 
	"""
	
	m = float(m)
	n = float(n)
	N = float(N)
	p_m_min = sqrt(N*m/n)
	p_n_min = sqrt(N*n/m)
	c_min = m * p_n_min + n * p_m_min
	
	# test point p
	p_m = ceil(p_m_min)
	p_n = floor(p_n_min)
	if p_n * p_m < N : # p does no satisfy the conditions
		# test point p'
		p_m = floor(p_m_min)
		p_n = ceil(p_n_min)
		if p_n * p_m < N : # p' does no satisfy the conditions
			# set (p_m, p_n) = p''
			p_m = ceil(p_m_min)
	if 1 <= p_m and p_m <= m/2.0 and 1 <= p_n and p_n <= n/2.0 : # (p_m, p_n) is a candidate solution
		# save p_m0 and set the current cost
		p_m0 = p_m
		c = m * p_n + n * p_m
	else : # the problem does not have an acceptable solution
		return (-1,-1)

	# optimal cost is in interval [c_min, c]
	
	# look for a better solution by increasing p_m
	p_m_cur = p_m + 1 # test a new candidate p_m
	p_n_cur = (c - n * p_m_cur) / m
	while p_m_cur <= m/2.0 and (p_m_cur * p_n_cur) >= N : # there may be a better solution 
		if 1 <= floor(p_n_cur) and (p_m_cur * floor(p_n_cur)) >= N : # this is such a better solution
			# update solution
			p_m = p_m_cur
			p_n = floor(p_n_cur)
			c   = m * p_n + n * p_m
		# look for a better solution further increasing p_m
		p_m_cur = p_m_cur + 1
		p_n_cur = (c - n * p_m_cur) / m
			
	# look for a better solution by decreasing p_m 
	p_m_cur = p_m0 - 1 # start from the initial candidate p_m
	p_n_cur = (c - n * p_m_cur) / m
	while floor(p_n_cur) <= p_m_cur and (p_m_cur * p_n_cur) >= N : # there may be a better solution
		if floor(p_n_cur) <= n/2.0 and (p_m_cur * floor(p_n_cur)) >= N : # this is such a better solution
			# update solution
			p_m = p_m_cur
			p_n = floor(p_n_cur)
			c   = m * p_n + n * p_m
		# look for a better solution further decreasing p_m
		p_m_cur = p_m_cur - 1
		p_n_cur = (c - n * p_m_cur) / m
	
	return (int(p_m), int(p_n))


###############################################################################
#                                                                             #
# ----------------- Common functions ---------------------------------------- #       
#                                                                             #
###############################################################################

def extract_params():
   """
   Extract parameter from line of commands.
   Output: 
      params = list of parameters from original command line
   """
   params = (sys.argv)
   return params


def check_flag(params, string, delete):
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
   for line in params:
      tmp = line.find(string)
      if tmp != -1:
         start = tmp + size
         sel_string = line[start:]
         if delete :
            params.pop(i)
         value = sel_string
      i += 1
   return value

 
def pop_left(dictionary):
   """
   Cuts the first element of dictionary and returns its first element (key:value)
   Input/Output: 
     dictionary = Dictionary of string containing the command lines to use. After reading the dictionary the first element is deleted from the dictionary.
   Output:
     first_el = first element (values) of the dictionary
   """	
   if len(dictionary) > 0:
      first_el ={list(dictionary.keys())[0] : list(dictionary.values())[0]}
      dictionary.pop(list(dictionary.keys())[0])
   else:
      first_el = None
   return first_el

 
def launcher(input_file):
   """
   Perform elaboration for each element of the queue.
   Input/Output
      input_file = command to be executed
   """
   myrank = comm.Get_rank()

   t1 = time.time()
   # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Stringa da Modificare!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   #print ("Terastitcher n. ", list(input_file.keys())[0], " is executed by rank: ", myrank)  
   print ("Scheduled job n. ", list(input_file.keys())[0], " is executed by rank: ", myrank)  
   # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   if debug_level > 0 :
      execution_string = prefix + list(input_file.values())[0] + " > " + "output_" + str(list(input_file.keys())[0]) + ".out"
   else :
      execution_string = prefix + list(input_file.values())[0]
   print (execution_string)
   # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Sleep da cancellare!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   #time.sleep((5.+100./(list(input_file.keys())[0]+1.))/100.)
   # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   os.system(execution_string)
   t2 = time.time()
   print (" ---> Processor ",  myrank, " has calculated for " , t2-t1) 
   return input_file


def worker():
   """
   Worker process.
   """
   myrank = comm.Get_rank()
   WORKTAG = 1   # constants to use as tags in communications
   DIETAG = 2
   while True:
      status = MPI.Status()
      input_name = comm.recv(source=0, tag=MPI.ANY_TAG, status=status)
      # check the tag of the received message
      if status.tag == DIETAG:
         end_signal = ['Exit cpu n. ', myrank]
         print (end_signal[0], end_signal[1])
         comm.send(end_signal, dest=0, tag=1)
         return

      # do the work
      result = launcher(input_name)
      comm.send(result, dest=0, tag=0)

 
###############################################################################
#                                                                             #
# ----------------- Function for step 2 parallelization --------------------- #       
#                                                                             #
###############################################################################

def read_input(inputf,nline=0):
    """
    Reads the file included in inputf at least up to line number nline (if declared).

    """
    i = 0
    data = [] 
    f = open(inputf, 'r')
    for line in f:
    	line = line.strip()
    	l = line.split(' ', 1)
    	data.append(l)
    	if (nline != 0) and (i > nline):
    		break
    	i += 1
    f.close()
    return data


# # suggested by raghu.erapaneedi@uni-wuerzburg.de
# def extract_np(inputf):
#     dom = minidom.parse(inputf)
#     dimensions = dom.getElementsByTagName('dimensions')[0]
#     nrows = int(dimensions.getAttribute('stack_rows'))
#     ncols = int(dimensions.getAttribute('stack_columns'))
#     nslices = int(dimensions.getAttribute('stack_slices'))
            
def extract_np(inputf):
    """
    extract the number of slices along z from the input xml file.

    """
    data = read_input(inputf, 8)
    for tmp_string in data :
    	#tmp_string = data[8][1]
    	if tmp_string[1].find('stack_slices="') != -1 :
    		start = tmp_string[1].find('stack_slices="')+14
    		break
    end = tmp_string[1].find('" />')
    sel_string = tmp_string[1][start:end]
    nslices = int(sel_string)
    
    start = tmp_string[1].find('stack_rows="')+12
    end = start + tmp_string[1][start:].find('"')
    sel_string = tmp_string[1][start:end]
    nrows = int(sel_string)
    
    start = tmp_string[1].find('stack_columns="')+15
    end = start + tmp_string[1][start:].find('"')
    sel_string = tmp_string[1][start:end]
    ncols = int(sel_string)
    
    return (nrows, ncols, nslices)


def find_last_slash(string):
	"""
	Search for / in a string. If one or more / was found, divide the string in a list of two string:
	the first containf all the character at left of the last / (included), 
	and the second contains the remanent part of the text.
	If no / was found, the first element of the list will be set to ''
	"""
	len_string = len(string)
	check = 0
	index = []
	i = 0
	for chara in string:
		if chara == '/' or chara == '\\':
			index.append(i)
			check =1
		i += 1
	if check == 1:
		last_slash = max(index)
		output_string = [string[0:last_slash+1],string[last_slash+1:]]
	else:
		output_string = ['',string]
	
	return output_string


def add_chars(params):
	string =['volin_plugin=','imin_plugin']
	i= 0
	for line in params:
		for local_string in string:
			tmp = line.find(local_string)
			if tmp != -1:
				size = len(local_string)
				sel_string = line.split('=')
				input_string = sel_string[1]
				mod_string = '"'+input_string+'"'
				tot_mod_string = sel_string[0]+'='+mod_string
				params[i] = tot_mod_string
		i += 1
	return params
	

def dispatcher_step2(queue):
    """
    dispatch the work among processors

    queue is a list of job input

    """
    # @ADDED by Giulio Iannello 2019-08-01 to fix distributed termination algorithm
    n_tasks = len(queue)
    done = 0
    
    # constants to use as tags in communications
    WORKTAG = 1   
    DIETAG = 2

    nprocs = comm.Get_size()
    # queue = deque(queue) # deque to enable popping from the left 
    # seed the workers by sending work to each proc
    for rank in range(1, min(len(queue)+1,nprocs)):
        # get the next input to work on 
        input_file = pop_left(queue)
        # send the next input to each rank
        comm.send(input_file, dest=rank, tag=WORKTAG)
        
    print('DISPATCHER: first loop terminated')
    
    # Loop until there's no more work to do. If queue is empty skips the loop.
    while (queue):
        input_file = pop_left(queue)
        # receive result from worker
        status = MPI.Status()
        flag = comm.recv(source=MPI.ANY_SOURCE, tag=MPI.ANY_TAG, status=status)
        # @ADDED by Giulio Iannello 2019-08-01 to fix distributed termination algorithm
        if status.tag != 0 :
            raise ValueError("Wrong tag: a message signalling a finished task expected")
        done += 1 
        # send to the same worker new work
        comm.send(input_file, dest=status.source, tag=WORKTAG) 
    print('DISPATCHER: second loop terminated')

    while done < n_tasks :
        status = MPI.Status()
        flag = comm.recv(source=MPI.ANY_SOURCE, tag=MPI.ANY_TAG, status=status)
        # @ADDED by Giulio Iannello 2019-08-01 to fix distributed termination algorithm
        if status.tag != 0 :
            raise ValueError("Wrong tag: a message signalling a finished task expected")
        done += 1 

    print('DISPATCHER: third loop terminated')

    # there's no more work to do, so receive all the results from the workers
    status = MPI.Status()

    # tell all the workers to exit by sending an empty message with the DIETAG
    for rank in range(1, nprocs):
        comm.send(0, dest=rank, tag=DIETAG)
	
    # Wait all task in order to exit

    for rank in range(1, nprocs):
	    exit_m = comm.recv(source=rank, tag=1, status=status)


def do_additional_partition ( nprocs, nrows, ncols, n_ss ) :
	"""
	All parameters should be float
	"""

	print ('do_additional_partition ( nprocs =', nprocs, ', nrows =', nrows, ', ncols =', ncols, ', nlayers =', n_ss, ')')
	
	if n_ss >= 2*nprocs : # dataset has not to be further partitioned
		return (1,1)
		
	elif floor(nrows/2)*floor(ncols/2)*n_ss < 2*nprocs :
		print ("WARNING: not enough paritions for", nprocs, "processors")
		print ("Dataset partitioned in", int(floor(nrows/2)), "x", int(floor(ncols/2)), "x", n_ss, "=", int(floor(nrows/2))*int(floor(ncols/2))*n_ss, "partitions")
		return (max(1,floor(nrows/2)), max(1,floor(ncols/2)))
		
	else :                 # dataset has to be further partitioned
		m = max(nrows,ncols)
		n = min(nrows,ncols)

		(p_m, p_n) = partition ( m, n, ceil(2*nprocs/n_ss) ) # n_ss partitions are already available along Z
		
		if nrows < ncols : # the greater dimensione is the second: p_m and p_n have to exchanged
			temp = p_m
			p_m  = p_n
			p_n  = temp
			
		return (p_m, p_n)
				
# 		#initialize global variables (tables)
# 		initTables(m,n)
# 		
# 		c = costo(m,n,n_ss,2*nprocs) # each tile is already partitioned in n_ss layers
# 						
# 		return partition(nrows,ncols,n_ss,2*nprocs) # each tile is already partitioned in n_ss layers

	
###############################################################################
#                                                                             #
# ----------------- Function for step 6 parallelization --------------------- #       
#                                                                             #
###############################################################################

def score_function(params):
   """
   Assigns a score value with the formula:
         score = 100*N_of_voxel/max(N_of_voxel)
   Input:
      params =  dictionary containing {input_name : [Nx,Ny,Nz]}
   Output: 
      scores = dictionary containing {input_name : score}
   """
   tmp_scores = {}
   scores = {}
   imp_key = params.keys()
   for i in imp_key:
      tmp = params[i]
      npoints = tmp[0]*tmp[1]*tmp[2]
      tmp_scores[i] = npoints
   den = max(tmp_scores.values())
   for i in tmp_scores:
      scores[i] = 100.*tmp_scores[i]/den
   return scores


def sort_elaborations(scores):
   """
   Create a list of input_name sorted by score
   Input:
     scores = dictionary of the form  {input_name : score}
   Output:
     scored = a list of input_name sorted by score
   """
   scored = sorted(scores, key=scores.__getitem__, reverse=True)
   return scored


def sort_work(params,priority):
   """
   Returns a dictionary as params but ordered by score
   Input:
      params = dictionary of the form  {input_name : value}
      priority = the list of input_name ordered by score calculated by score_function
   Output:
      sorted_dict = the same dictionary as params but ordered by score
   """	
   sorted_dict = {}
   i = 0
   for index in priority:
      sorted_dict.update({i:params[index]})
      i = i + 1 
   return sorted_dict


def dispatcher_step6(queue,rs_fname):  # 2017-02-06. Giulio. @ADDED rs_fname
   """
   Dispatch the work among processors.
   Input:
      queue = list of job inputs
   """

   # support for suspend/resume
   # @ADDED by Giulio Iannello 2017-02-06
   n_tasks = len(queue)
   # @CHANGED by Giulio Iannello 2017-03-12
   if suspend_resume_enabled :
      rs_file = open(rs_fname, 'rb')
      done = pickle.load(rs_file)
      rs_file.close()
   else :
      done = []
      
   # constants to use as tags in communications
   WORKTAG = 1   
   DIETAG = 2
   nprocs = comm.Get_size()
   #queue = deque(queue) # deque to enable popping from the left 
   # seed the workers by sending work to each proc
   for rank in range(1, min(len(queue)+1,nprocs)):
      # get the next input to work on 
      input_file = pop_left(queue)
      while (input_file != None) and (list(input_file.keys())[0] in done) :
         input_file = pop_left(queue)
      if input_file == None :
         # the queue is empty, no more work to do
         break
      # send the next input to each rank
      comm.send(input_file, dest=rank, tag=WORKTAG) 
      # Loop until there's no more work to do. If queue is empty skips the loop.
      
   print('DISPATCHER: first loop terminated')
      
   while (queue) :
      input_file = pop_left(queue)
      while (input_file != None) and (list(input_file.keys())[0] in done) :
         input_file = pop_left(queue)
      if input_file == None :
         # the queue is empty, no more work to do
         break
      # receive result from worker
      status = MPI.Status()
      flag = comm.recv(source=MPI.ANY_SOURCE, tag=MPI.ANY_TAG, status=status)
      if status.tag != 0 :
         raise ValueError("Wrong tag: a message signalling a finished task expected")
      done.append(list(flag.keys())[0])
      if suspend_resume_enabled :
         # save new resume status
         rs_file = open(rs_fname, 'wb')
         pickle.dump(done,rs_file)
         rs_file.close()
      # send to the same worker new work
      comm.send(input_file, dest=status.source, tag=WORKTAG)
      
   print('DISPATCHER: second loop terminated')
      
   while len(done) < n_tasks :
      # receive result from worker
      status = MPI.Status()
      flag = comm.recv(source=MPI.ANY_SOURCE, tag=MPI.ANY_TAG, status=status)
      if status.tag != 0 :
         raise ValueError("Wrong tag: a message signalling a finished task expected")
      done.append(list(flag.keys())[0])
      if suspend_resume_enabled :
         # save new resume status
         rs_file = open(rs_fname, 'wb')
         pickle.dump(done,rs_file)
         rs_file.close()
         
   print('DISPATCHER: third loop terminated')
      
   # there's no more work to do, so receive all the results from the workers
   status = MPI.Status()
   #flag = comm.recv(source=MPI.ANY_SOURCE, tag=0, status=status)
   # tell all the workers to exit by sending an empty message with the DIETAG
   for rank in range(1, nprocs):
      comm.send(0, dest=rank, tag=DIETAG)
    # Whait all task in order to exit
   for rank in range(1, nprocs):
      exit_m = comm.recv(source=rank, tag=1, status=status)
   if suspend_resume_enabled :
      os.remove(rs_fname) # 2017-02-05. Giulio. @ADDED remove the temporary file containing resume state  


###############################################################################
#                                                                             #
# -------------------- Paraconverter Service Functions ---------------------- #       
#                                                                             #
###############################################################################

def read_params():
   """
   Read parameters from input string and from a file
   Input: 
   Output:
      input_name = Input file
      output_name = Standard output directory
      wb1 = Axprossimative depth for the tiles
      wb2 = Axprossimative height for the tiles
      wb3 = Axprossimative width for the tiles
      sfmt = Source format
      dfmt = Destination format
      iresolutions = List of integer values containing all the desidered values for level of resolution
      max_res = Maximum level of resolution available (integer)
      params = Array containing instruction derived from the remanent part of the imput string
      last_string = Remanent part of the input string
      height = Height of the input immage
      width = Width of the input immage
      depth = Depth of the input immage
   """

   # Get the input list
   params = extract_params()
   params.pop(0)
   #Check if quoting is necessary
   params = check_double_quote(params)
   # 2017-03-13. Giulio. @COMMENTED command option --np 
   # 2016-12-10. Giulio. @ADDED command option --np to initialize parameter gi_np
   #gi_np = read_item(params, '--np=', 1)
   #print("np is : ", gi_np) 
   # Declare input file
   input_name = read_item(params, '-projin=', './input.xml')

   #------------- generate a temporary file with image size information
   
   info_string = 'teraconverter ' + ' --sfmt="TIFF (unstitched, 3D)"' + ' --dfmt="TIFF (tiled, 3D)"' + ' -s="' + input_name + '" -d=/'
   #pwd = check_output("pwd",shell=True)
   #execution_string = prefix+info_string + ' --info="' + pwd.strip("\n") + '/__dims__.txt"'
   execution_string = prefix+info_string + ' --info="' + os.getcwd() + '/__dims__.txt"'
   os.system(execution_string)
   print (execution_string)

   #------------- 
   
   #Read from Internal Input file (--origin=file_in) HEIGHT, WIDTH, DEPTH (file_in is Name of the file containing the size of ouf computational domain in Voxel)
   file_in = read_item(params, '--origin=', os.getcwd() + '/__dims__.txt') # changed by Giulio: the --origin option is not needed any longer
   print("Origin file is: ", file_in)
   input_params_in_file = ['HEIGHT=', 'WIDTH=', 'DEPTH=', 'BYTESxCHAN=', 'DIM_C=', 'VXL_V=', 'VXL_H=', 'VXL_D=']
   params_from_file = search_for_entry(input_params_in_file,file_in)
   os.remove(os.getcwd() + '/__dims__.txt') # added by Giulio to remove the temporary file containing image size information
   # Create variables for HEIGHT, WIDTH and DEPTH
   height = int(params_from_file[0])
   width = int(params_from_file[1])
   depth = int(params_from_file[2])
   bytes_x_chan = int(params_from_file[3])
   n_chans = int(params_from_file[4])
   vxl_V = abs(float(params_from_file[5]))
   vxl_H = abs(float(params_from_file[6]))
   vxl_D = abs(float(params_from_file[7]))

   # Find standard output directory
   output_name = read_item(params, '-volout=', './OUT')
   # Declare axprossimative depth for the tiles
   wb1 = read_item(params, '--slicedepth=', 0)
   # Declare axprossimative height for the tiles
   wb2 = read_item(params, '--sliceheight=', 0)
   # Declare axprossimative width for the tiles
   wb3 = read_item(params, '--slicewidth=', 0)
   # Source format is fixed
   sfmt = '"TIFF (unstitched, 3D)"'
   # Declare destination format
   volout_plugin = read_item(params, '--volout_plugin=', '"TiledXY|2Dseries"')
   #dfmt = read_item(params, '--dfmt=', '"TIFF (tiled, 3D)"')
   if volout_plugin == 'TiledXY|3Dseries' :
      dfmt = '"TIFF (tiled, 3D)"'
   else :
      dfmt = '"TIFF (tiled, 2D)"'
   # Declare resolutions
   resolutions = read_item(params, '--resolutions=', '0')
   iresolutions = [int(resolutions[0])]
   len_res = len(resolutions)
   if len_res > 0:
      for i in range(1,len_res):
         iresolutions.append(int(resolutions[i]))
   max_res = max(iresolutions)
   # check isotropic flag
   isotropic = read_item(params, '--isotropic', 'False')
   if isotropic == '' : # isotropic has been set
      isotropic = True
   else :
      isotropic = False
   
   #if isotropic is set computes the halving rule to apply to dimension D
   if isotropic :
      # computes the number of halving to be performed along D
      vxlsz_Vx2 = 2 * vxl_V
      vxlsz_Hx2 = 2 * vxl_H
      vxlsz_D   = vxl_D
      h = 0
      while h < max_res and max(vxlsz_Vx2,vxlsz_Hx2) < vxlsz_D :
         h += 1
         vxlsz_Vx2 *= 2;
         vxlsz_Hx2 *= 2;
   else :
      # halving along D dimension must be always performed
      h = 0
   max_res_D = max_res - h
   print("vxl_V, vxl_H, vxl_D, isotropic, h, max_res, max_res_D :", vxl_V, vxl_H, vxl_D, isotropic, h, max_res, max_res_D)

   # Convert the remanent list of string input parameters into a single string
   last_string = collect_instructions(params)
   # Return values
   # 2016-12-10. Giulio. @ADDED parameter gi_np
   return (input_name, output_name, wb1, wb2, wb3, sfmt, dfmt, iresolutions, max_res, isotropic, max_res_D, params, last_string, height, width, depth, bytes_x_chan, n_chans)



def read_item(input_arr, item, default,message=True):
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
   if tmp == None:
      value = default
      if message:
        print ("The value for ",item," was not declared. It will be set to", value, "by default.")
   else:
      if isinstance(default,int):
         value = int(tmp)
      elif isinstance(default,float):
         value = float(tmp)
      else:
         value = tmp
   return value



def collect_instructions(inst):
    """
    Collect the remanent part of a list of strings in a unique string
    Input:
      inst = Input list of strings
    Output:
      results = String containing all the elements of inst
    """
    len_inst =len(inst)
    if len_inst > 0:
       for i in range(0,len_inst):
          if i == 0:
            results = str(inst[i])
          else:
            results = results + ' ' + str(inst[i])
    else:
       results = '' # if inst is empty a null string must be returned
    return results



def search_for_entry(string_2_serch,file_in, nline=0):
    """
    Extract from the input file (file_in) up to the line number nline (if declared) the value assigned to string_2_serch.
    Input:
      string_2_serch = string (or list of string) containing the variable to search (e.g. 'HEIGHT=')
      file_in = name of the file containing the information we neeed (e.g: prova.txt or /pico/home/prova.txt)
      nline = optional, number of the final row of the file we need to analyze
    Output:
      Output = value or (list of values) assigned to the variable conteined in string_2_serch
    """
    # Read the file
    i = 0
    data = [] 
    f = open(file_in, 'r')
    for line in f:
       line = line.strip()
       l = line.split(' ', 1)
       data = data + l
       if (nline != 0) and (i > nline):
          break
       i += 1
    f.close()
    # Read the value/s from a table of string
    len_string = len(string_2_serch)
    if len_string <= 0:
       print("No possible options! No values will be created!")
    elif len_string == 1:
       tmp = check_flag(data, string_2_serch[0],True)
       if tmp == None:
          output = '0'
          print ("The name of ", string_2_serch," was not declared. It will be set to",  output, "by default.")
       else:
          output = tmp
    elif len_string > 1:
       ii = 0
       output = []
       for i in string_2_serch:
          tmp = check_flag(data, i,True)
          if tmp == None:
             output.append('0')
             print ("The name of ", i," was not declared. It will be set to",  output[ii], "by default.")
          else:
             output.append(tmp)
          ii = ii + 1
    else:
       print("No possible options! No values will be created!")
    return output



def sort_list (len_1, len_2, len_3):
   """
   Create a list sorting the indexes along three directions:
   Input: 
      len_1 = Number of elements of the array for the first index
      len_2 = Number of elements of the array for the second index
      len_3 = Number of elements of the array for the third index
   Output:
      order = An ordered list containig an a sequence of lists of 3 alements (one for each direction) that identify the position on the local index
   """
   order =[]
   for i in range(0,len_1):
      for j in range(0,len_2):
         for k in range(0,len_3):
            order.append([i,j,k])
   return order



def sort_start_end(start_1, start_2, start_3, end_1, end_2, end_3,size_1, size_2, size_3):
   """
   Sort start points and edn point in two lists of elements
   Input:
      start_1 = Array containing all the starting indexes for the tiles on the Depth direction
      start_2 = Array containing all the starting indexes for the tiles on the Height direction
      start_3 = Array containing all the starting indexes for the tiles on the Width direction
      end_1 = Array containing all the ending indexes for the tiles on the Depth direction
      end_2 = Array containing all the ending indexes for the tiles on the Height direction
      end_3 = Array containing all the ending indexes for the tiles on the Width direction
      size_1 = Array containing the size of the tile in the Depth direction
      size_2 = Array containing the size of the tile in the Height direction
      size_3 = Array containing the size of the tile in the Width direction
   Output:
      order = An ordered list containig an a sequence of lists of 3 alements (one for each direction) that identify the position on the local index 
      start_list = Ordered list of lists of starting points. E.g.: [[width_in[0], height_in[0], depth_in[0]], [width_in[1], height_in[1], depth_in[1]], ... ,[width_in[N], height_in[N], depth_in[N]]]
      end_list = Ordered list of lists of starting points. E.g.: [[width_fin[0], height_fin[0], depth_in[0]], [width_fin[1], height_fin[1], depth_fin[1]], ... ,[width_fin[N], height_fin[N], depth_fin[N]]]
      len_arr = Dictionary containing elements like {index:[size_width(i),size_height(i),size_depth(i)],.....}
   """
   len_1 = len(start_1)
   len_2 = len(start_2)
   len_3 = len(start_3)
   #call sort_list
   order = sort_list (len_1, len_2, len_3)
   len_list = len(order)
   start_list = []
   end_list = []
   len_arr = {}
   for i in range(0,len_list):
      tmp = [start_3[order[i][2]], start_2[order[i][1]],start_1[order[i][0]]]
      start_list.append(tmp)
      tmp = [end_3[order[i][2]], end_2[order[i][1]], end_1[order[i][0]]]
      end_list.append(tmp)
      tmp = [size_3[order[i][2]], size_2[order[i][1]], size_1[order[i][0]]]
      len_arr.update({i:tmp})
   return (order, start_list, end_list, len_arr)



def check_double_quote(inpstring):
   """
   Check if some strings needs of a double quote (if some space are inside the string, it will need to be inside two double quote). E.g.: --sfmt="TIFF (unstitched, 3D)"
   Input:
      inpstring: input string or array of strings
   Output:
      newstring = new string (or array of strings) corrected by quoting if necessary
   """
   if type(inpstring) == list:
      newstring = []
      for index in inpstring:
         tmp1 = index.find(' ')
         if tmp1 != -1:
            tmp2 =  index.find('"')
            if tmp2 == -1:
               dummy = index.find('=')
               if dummy != -1:
                  newstring.append(index[0:dummy+1] + '"' + index[dummy+1:] + '"')
               else:
                  newstring.append('"' + index + '"')
            else:
               newstring.append(index)
         else:
            newstring.append(index)
   else:
      tmp1 = inpstring.find(' ')
      if tmp1 != -1:
         tmp2 = inpstring.find('"')
         if tmp2 == -1:
            dummy = inpstring.find('=')
            if dummy != -1:
               newstring = inpstring[0:dummy+1] + '"' + inpstring[dummy+1:] + '"'
            else:
               newstring = '"' + inpstring + '"'
         else:
            newstring = inpstring
      else:
         newstring = inpstring
   return newstring



def eliminate_double_quote(inpstring) :
   """
   Check if the string is already enclosed by quotes
   Input:
      inpstring: input string or array of strings
   Output:
      newstring = new string (or array of strings) corrected by eliminating enclosing quotes if any
   """
   len_str = len(inpstring)
   if (inpstring[0] == '"') and (inpstring[len_str-1] == '"') or (inpstring[0] == "'") and (inpstring[len_str-1] == "'") :
      newstring = inpstring[1:len_str-1]
   return newstring
   


def generate_first_command(input_name, output_name, wb1, wb2, wb3, sfmt, dfmt, iresolutions, max_res,params, last_string):
   """
   Generate first command line
   Input:
      input_name = Input file
      output_name = Standard output directory
      wb1 = Axprossimative depth for the tiles
      wb2 = Axprossimative height for the tiles
      wb3 = Axprossimative width for the tiles
      sfmt = Source format
      dfmt = Destination format
      iresolutions = List of integer values containing all the desidered values for level of resolution
      max_res = Maximum level of resolution available (integer)
      params = Array containing instruction derived from the remanent part of the imput string
      last_string = Remanent part of the input string
   Output:
      first_string = Command line to preprocess the data 
   """
   first_string = 'teraconverter ' + '--height=' + str(wb2)+ ' --width=' + str(wb3) 
   first_string = first_string + ' --depth=' + str(wb1) + ' --sfmt=' + sfmt + ' --dfmt=' + dfmt
   tmp_res = ''
   for i in iresolutions:
      tmp_res = tmp_res + str(i)
   first_string = first_string + ' --resolutions=' +  tmp_res + ' -s="' + input_name + '" -d="' + output_name + '" ' 
   if (last_string != []):
      first_string = first_string + last_string 
   first_string = first_string + ' --makedirs' 
   return first_string



def generate_final_command(input_name, output_name, wb1, wb2, wb3, sfmt, dfmt, iresolutions, max_res,params, last_string):
   """
   Generate last command line to merge metadata
   Input:
      input_name = Input file
      output_name = Standard output directory
      wb1 = Axprossimative depth for the tiles
      wb2 = Axprossimative height for the tiles
      wb3 = Axprossimative width for the tiles
      sfmt = Source format
      dfmt = Destination format
      iresolutions = List of integer values containing all the desidered values for level of resolution
      max_res = Maximum level of resolution available (integer)
      params = Array containing instruction derived from the remanent part of the imput string
      last_string = Remanent part of the input string
   Output:
      final_string = Command line to merge metadata 
   """
   final_string = 'teraconverter ' + '--height=' + str(wb2)+ ' --width=' + str(wb3) 
   final_string = final_string + ' --depth=' + str(wb1) + ' --sfmt=' + sfmt + ' --dfmt=' + dfmt
   tmp_res = ''
   for i in iresolutions:
      tmp_res = tmp_res + str(i)
   final_string = final_string + ' --resolutions=' +  tmp_res + ' -s="' + input_name + '" -d="' + output_name + '" ' 
   if (last_string != []):
      final_string = final_string + last_string 
   final_string = final_string + ' --metadata' 
   return final_string



def generate_parallel_command(start_list, end_list, input_name, output_name, wb1, wb2, wb3, sfmt, dfmt, iresolutions, max_res,params, last_string):
   """
   Generate the list of parallel command lines
   Input:
      start_list = Ordered list of lists of starting points. E.g.: [[width_in[0], height_in[0], depth_in[0]], [width_in[1], height_in[1], depth_in[1]], ... ,[width_in[N], height_in[N], depth_in[N]]]
      end_list = Ordered list of lists of starting points. E.g.: [[width_fin[0], height_fin[0], depth_in[0]], [width_fin[1], height_fin[1], depth_fin[1]], ... ,[width_fin[N], height_fin[N], depth_fin[N]]]
      input_name = Input file
      output_name = Standard output directory
      wb1 = Axprossimative depth for the tiles
      wb2 = Axprossimative height for the tiles
      wb3 = Axprossimative width for the tiles
      sfmt = Source format
      dfmt = Destination format
      iresolutions = List of integer values containing all the desidered values for level of resolution
      max_res = Maximum level of resolution available (integer)
      params = Array containing instruction derived from the remanent part of the imput string
      last_string = Remanent part of the input string
   Output:
      list_string = Dictionary of strings containing the command lines to process the data. E.G.: {i:command[i]} 
   """
   index = len(start_list)
   list_string = {}
   for i in range(0,index):
      dummy = ''
      dummy = 'teraconverter ' + '--height=' + str(wb2)+ ' --width=' + str(wb3) 
      dummy = dummy + ' --depth=' + str(wb1) + ' --sfmt=' + sfmt + ' --dfmt=' + dfmt
      tmp_res = ''
      for j in iresolutions:
         tmp_res = tmp_res + str(j)
      dummy = dummy + ' --resolutions=' +  tmp_res + ' -s="' + input_name + '" -d="' + output_name + '" ' 
      if(last_string != []):
         dummy = dummy + last_string 
      dummy = dummy + ' --parallel'
      dummy = dummy + ' --H0=' + str(start_list[i][0]) + ' --H1=' +str( end_list[i][0])
      dummy = dummy + ' --V0=' + str(start_list[i][1]) + ' --V1=' + str(end_list[i][1])
      dummy = dummy + ' --D0=' + str(start_list[i][2]) + ' --D1=' + str(end_list[i][2])
      list_string.update({i:dummy})
   return list_string



def opt_algo(D, w, n):
    """
    Solves the tiling problem
    patitioning the interval [0, D-1] into k subintervals of size
    2^n b and one final subinterval of size r = D - k 2^n b
    Input:
      D = dimension of the original array
      w = approximate estimation of value for b
      n = desideres level of refinement (e.g. : n = 0 => maximum level of refinement; n =1 => number of point divided by 2^1=2; n = 2 => number of point divided by 2^2=4;)
    Output:
      arr_sizes = [b, r, k, itera]
         b = normalized size of standard blocks (size of standard blocks = b * 2^n)
         r = rest (if not equal to 0, is the size of the last block)
         k = number of standard blocks
         itera = number of itarations to converge
    """
    # step 1
    h = 0
    b_h = w
    k_h = floor(D/(pow(2.,n) * b_h))
    b = b_h
    r = D % (pow(2.,n) * b_h)
    k = k_h

    itera = 0
    verif = bool(1)
    while verif:
       # step 2
       if (D % (pow(2.,n) * b_h)) == 0:
           b = b_h
           r = 0
           k = k_h
           verif = bool(0) # exit form while cicle
       # step 3
       elif (D % (pow(2.,n) * b_h)) > r:
           b = b_h;
           r = (D % (pow(2.,n) * b_h));
           k = k_h;
       # step 4
       if h == floor(w/2):
           verif = bool(0) # exit form while cicle
       h = min( floor(w/2), h + max(1,floor(b_h - D/(pow(2.,n) * (k_h+1)))) );
       b_h = w - h;
       k_h = floor(D/(pow(2.,n) * b_h));
       itera = itera + 1;
    b = int(b)
    r = int(r)
    k = int(k)
    arr_sizes = [b, r, k, itera]
    return arr_sizes


def prep_array(wb, r, k):
    """
    Create a 1D array containing the number of elements per tile.
    Input: 
         wb = size of standard blocks
         r = rest (if not equal to 0, is the size of the last block)
         k = number of standard blocks
    Output:
       array = A list containing the number of element for every tiles.
    """
    for i in range(0,k):
       if i == 0:
          array = [int(wb)]
       elif i > 0:
          array.append(int(wb))
       else:
          print ('Option not permitted!!!!!! i =',i)
          sys.exit(1)
    if r != 0:
       if k != 0:
          array.append(r)
       else:
         array = [r]
    return array


def create_sizes(size, wb, max_res, norest=False):
   """
   Create a 3D array containing the size for each tile on the desidered direction
   Input: 
      start_wb = Start parameter for b
      size = size (in pixel) of the input immage
      wb = Rough depth for the tiles in the desidered direction
      max_res = Maximum level of resolution available (integer)
      norest = Boolean variable to chech if we need of the last array element (if it is different from the preavious one)
   Output:
      arr = Array containing the size for each tile on the desidered direction
   """
   # Make the partition
   values = opt_algo(size, wb, max_res)
   b = values[0]
   r = values[1]
   k = values[2]
   itera = values[3]
   wb = int(pow(2,max_res)*b)
   arr = prep_array(wb, r, k)
   if norest:
      tmp_len = len(arr)
      if arr[tmp_len-1] != arr[tmp_len-2]:
        print('Attention! : ', arr[tmp_len-1],' points was deleted!')
        arr.pop()
   return (arr)


def create_starts_end (array, start_point=0,open_dx=True):
   """
   Create arrays containing all the starting and ending indexes for the tiles on the desidered direction
   Input:
      array = Array containing the size for each tile on the desidered direction
      start_point = Starting index for the input immage (optional)
      open_dx = If true (the default value) ==> ending indexes = subsequent starting indexes ==> Open end
   Output:
      star_arr = Array containing all the starting indexes for the tiles on the desidered direction
      end_arr = Array containing all the ending indexes for the tiles on the desidered direction
   """
   len_arr = len(array)
   ind_arr = range(0, len_arr)
   start_arr = []
   end_arr = []
   if open_dx:
      dx_pad = 0
   else:
      dx_pad = -1
   for i in ind_arr:
      if i != 0:
         start_point = start_point + array[i-1]
      start_arr.append(start_point)
      end_point = start_point + array[i] + dx_pad
      end_arr.append(end_point)
   return (start_arr, end_arr)



def ctrl_parallelism (sfmt,dfmt) :
   partition_depth  = True  # parallelization can exploit depth (D) dimension
   partition_width  = True  # parallelization can exploit width (H) dimension 
   partition_height = True  # parallelization can exploit width (V) dimension 
   
   if sfmt == "TIFF (3D)" or dfmt == "TIFF (series, 2D)" :
      # cannot read subregions along width or height of source, or cannot generate subregions of slices
      partition_width  = False
      partition_height = False
      
   if sfmt == "TIFF (series, 2D)" :
      # cannot read subregions along width of source
      partition_width = False 
      
   # for input tiled formats and unstitched volumes should check if tiling is small enough 

   return (partition_depth,partition_width,partition_height)
   


def create_commands(gi_np,info=False): # 2016-12-10. Giulio. @ADDED parameter gi_np: desired parallelism
   """
   Create commands to run in parallel
   Input:
   Output:
      first_string = String to initialize parallel computation
      list_string = Dictionary of strings containing the command lines to process the data. E.G.: {i:command[i]}
      len_arr = Dictionary containing elements like {index:[size_width(i),size_height(i),size_depth(i)],.....}
      final_string = String to merge all metadadata
   """
   # Call function read_params
   # 2017-03-12. Giulio. @REMOVED parameter gi_np (it is passed from main)
   (input_name, output_name, wb1, wb2, wb3, sfmt, dfmt, iresolutions, max_res, isotropic, max_res_D, params, last_string, height, width, depth, bytes_x_chan, n_chans) = read_params()
   print("#"*80)
   print("Input file = ", input_name)   
   print("Output directory", output_name)
   # set wb1, wb2, wb3 if needed
   if wb1 == 0 :
      wb1 = default_tile_size
   if wb2 == 0 :
      wb2 = default_tile_size
   if wb3 == 0 :
      wb3 = default_tile_size
   print("Rough depth for the tiles in width direction = ", wb3)
   print("Rough depth for the tiles in height direction = ", wb2)
   print("Rough depth for the tiles in depth direction = ", wb1)
   print("Source Format = ", sfmt)
   print("Destination Format = ", dfmt)
   print ("Resolutions = ", iresolutions)
   print ("Max Resolutions", max_res)
   print ('Width (in voxel) of the immage = ', width)
   print ('Height (in voxel) of the immage = ', height)
   print ('Depth (in voxel) of the immage = ', depth)
   print(params)
   if isotropic :
      last_string = last_string + ' --isotropic'
   print("Last input elements of the original string = ", last_string)
   print("#"*80)
   # Call create_sizes function to create 3 arrays containing the sizes for each tile on the Height, Width, Depth directions.
   size_1 = create_sizes(depth, wb1, max_res_D) # max_res_D = max_res if isotropic is not set of if voxels are isotropic
   size_2 = create_sizes(height, wb2, max_res)
   size_3 = create_sizes(width, wb3, max_res)
   
   # added by Giulio 2016-12-10
   #print("initial len of size1, size2, size3: ", len(size_1), " - ", len(size_2), " - ", len(size_3))     
   #print("initial values of size1, size2, size3: ", size_1, " - ", size_2, " - ", size_3) 
   if dfmt == "HDF5 (BigDataViewer)" or dfmt == "HDF5 (Imaris IMS)" :
      raise ValueError("Paraconverter cannot be used with HDF5 output formats")
   (partition_depth,partition_width,partition_height) = ctrl_parallelism(eliminate_double_quote(sfmt),eliminate_double_quote(dfmt))
   print('--------> ',eliminate_double_quote(sfmt),eliminate_double_quote(dfmt),partition_depth,partition_width,partition_height)
   if len(size_1) >= 2*gi_np or (not partition_width and not partition_height) : # parallelism along D is enough
      size_2 = [height]
      size_3 = [width]
   else : # more parallelism is desirable
      if ((len(size_1) * len(size_2)) >= 2*gi_np) or not partition_width : # add parallelism along V only
         size_3 = [width]
      # else use all available parallelism
   print("number of work units (Depth, Height, Width): ", len(size_1), len(size_2), len(size_3))     
   print("size of work units (Depth, Height, Width): ", size_1, size_2, size_3)     
   # end added by Giulio
      
   if info :
      # return null values
      first_string = ''
      list_string  = ''
      final_string = ''
      len_arr = 0
      # print info
      voxel_num = round(1.1 * gi_np *(size_2[0] * size_3[0] * max(64,pow(2,max_res))) * n_chans * bytes_x_chan / (1024 * 1024 * 1024), 3)
      print("#"*80)
      print('Memory needed for ' + str(gi_np) + ' concurrent processes: ' + str(voxel_num) + ' GBytes') 
      print("#"*80)
   else :    
      # Call create_starts_end function to create 6 arrays containing the starting and ending points for each tile on the Height, Width, Depth directions. 
      (start_3,end_3) = create_starts_end(size_3,0)
      (start_2,end_2) = create_starts_end(size_2,0)
      (start_1,end_1) = create_starts_end(size_1,0)
      # Call sort_start_end to sort start points and end point in two lists of elements
      (order, start_list, end_list, len_arr) = sort_start_end(start_1, start_2, start_3, end_1, end_2, end_3,size_1, size_2, size_3)
      # Generate the string to initialize parallel computation
      first_string = generate_first_command(input_name, output_name, wb1, wb2, wb3, sfmt, dfmt, iresolutions, max_res, params, last_string)
      # Generate the list of parallel command lines
      list_string = generate_parallel_command(start_list, end_list, input_name, output_name, wb1, wb2, wb3, sfmt, dfmt, iresolutions, max_res, params, last_string)
      # Generate the string to merge all metadadata
      final_string = generate_final_command(input_name, output_name, wb1, wb2, wb3, sfmt, dfmt, iresolutions, max_res, params, last_string)

   # Return strings
   return (first_string, list_string, output_name, len_arr, final_string) # 2017-02-06. Giulio. @ADDED output_name



###############################################################################
#                                                                             #
# ------------------------ main (parallel) code ----------------------------- #       
#                                                                             #
###############################################################################

if __name__ == '__main__':
   # Initialize env. variables
   comm = MPI.COMM_WORLD
   nprocs = comm.Get_size()
   myrank = comm.Get_rank()
   # Sincronizzation
   comm.Barrier()
   
   tmp = read_item(sys.argv, '--info', 'no_info')
   if tmp == 'no_info' :
      info = False
   else :
      info = True

   step2 = False
   tmp = read_item(sys.argv, '--displcompute', 'no_step',False)
   if tmp != 'no_step' :
      step2 = True
   tmp = read_item(sys.argv, '-2', 'no_step',False)
   if tmp != 'no_step' :
      step2 = True
   step6 = False
   tmp = read_item(sys.argv, '--merge', 'no_step',False)
   if tmp != 'no_step' :
      step6 = True
   tmp = read_item(sys.argv, '-6', 'no_step',False)
   if tmp != 'no_step' :
      step6 = True      
    
   # setting the right prefix
   if platform.system() != 'Windows' :
      # platform is UNIX-like
      if step2 :
         if prefix == '' : # standard behavior 
            # check if terasticher is in the current directory
            if os.system('which ./terastitcher') == 0 :
               # set 'prefix' to current directory
               prefix = './'
            elif os.system('which terastitcher') != 0 :
               # terastitcher it is not in the system path
               raise ValueError("The executable of terasticher is not reachable")
      elif step6 :
         if prefix == '' : # standard behavior 
            # check if teraconverter is in the current directory
            if os.system('which ./teraconverter') == 0 :
               # set 'prefix' to current directory
               prefix = './'
            elif os.system('which teraconverter') != 0 :
               # terastitcher it is not in the system path
               raise ValueError("The executable of teraconverter is not reachable")
   # else do nothing: Windows search first in the current directory

   if myrank == 0:
           
      # timing
      t1 = time.time() 
      print ('*'*80)
      print (str(datetime.datetime.utcnow()), " -- Calculation started on ", nprocs, "- 1 cores.")
      print ('*'*80)
      
   # Sincronizzation
   comm.Barrier()
   
   if myrank == 0:
      if step2 :

         ########################## STEP 2

         #prefix = prefix + "terastitcher "
         execution_flag = '-2'

         # get the input list
         params = extract_params()
         print (params)
         #params = str("terastitcher -2 --projin=/Users/iannello/Home/Windows/myTeraStitcher/TestData/unstitched_RGB_2D/xml_import.xml --projout=/Users/iannello/Home/Windows/myTeraStitcher/TestData/unstitched_RGB_2D/xml_compdispl.xml --subvoldim=200 --imin_channel=G").split()
         params.pop(0)

         print ("Alignment will be performed")

         # select the minimum number of slices
         tmp = check_flag(params, 'subvoldim=',True)
         if tmp == None:
            nsubstring_min = 200
            print ("Number of substring was not declared. It will be set to",  nsubstring_min, "by default.")
         else:
            tmp = int(tmp)
            nsubstring_min = tmp
	
         # Declare input file
         tmp = check_flag(params, 'projin=',False)
         if tmp == None:
            input_name = 'xml_import.xml'
            print ("Name of the input file was not declared. It will be set to",  input_name, "by default.")
         else:
            input_name = tmp

         # Find standard output file
         tmp = check_flag(params, 'projout=',True)
         if tmp == None:
            output_name = 'xml_compdispl.xml'
            print ("name of the output file was not declared. It will be set to",  output_name, "by default.")
         else:
            output_name = tmp
         len_out = len(output_name) 
         output_name = output_name[0:len_out - 4]

         params = add_chars(params)

         # Reads the size of the tile matrix and the number slices
         (nrows, ncols, nslices) = extract_np(input_name)

         # Calculate the size of slices per task
         n_ss = int(ceil(float(nslices) / float(nsubstring_min))) # Number of substacks
         last_size = int(floor(float(nslices) / float(n_ss)))
         first_size = last_size + 1
         n_of_first_stacks = nslices % n_ss
         n_of_last_stacks = n_ss - n_of_first_stacks

         # look for further partitions		
#         nprocs = 6
#         nrows = 11
#         ncols = 12
         (p_r,p_c) = do_additional_partition(float(nprocs-1),float(nrows),float(ncols),float(n_ss)) # dispatcher should not be counted
         # compute start, end of rows partitions
         s_r = int(floor(nrows / p_r))
         r_r = nrows % p_r
         r_start = [0]
         r_end   = []
         for i in range(1,int(p_r - r_r)) :
            r_end.append(i*s_r)
            r_start.append(i*s_r)
         offset = int((p_r - r_r) * s_r)
         for i in range(int(r_r)) :
            r_end.append(offset + i*s_r)
            r_start.append(offset + i*s_r)
         r_end.append(nrows - 1)
         print (r_start, r_end)
         # compute start, end of columns partitions
         s_c = int(floor(ncols / p_c))
         r_c = ncols % p_c
         c_start = [0]
         c_end   = []
         for i in range(1,int(p_c - r_c)) :
            c_end.append(i*s_c)
            c_start.append(i*s_c)
         offset = int((p_c - r_c) * s_c)
         for i in range(int(r_c)) :
            c_end.append(offset + i*s_c)
            c_start.append(offset + i*s_c)
         c_end.append(ncols - 1)
         print (c_start, c_end)
         # generate temporary directory names
         if len(r_start)>1 or len(c_start)>1 :
            gr_dir_names = []
            for i in range(len(r_start)) :
               for j in range(len(c_start)) :
                  gr_dir_names.append('gr_R['+str(r_start[i])+','+str(r_end[i])+']_C['+str(c_start[j])+','+str(c_end[j])+']/')
	
         # Create some dictionaries containing: 
         # - the number of strips per substak (new_params),
         # - the index of the starting strip for each stack (start_dict)
         # - the index of the ending strip for each stack (start_dict)
         # - the entire command line for each stacks
         cmd_string = {}
         params_str = ' '.join(params)
         #if execution_flag == '-2': # it has already been checked

         tmp_out_name = find_last_slash(output_name)

         # Create temporary directory for temporary xml files 
         tmp_xml_dir = tmp_out_name[0]+'tmp/'
         #execution_string = 'mkdir -p ' + tmp_xml_dir
         #print (execution_string)
         #os.system(execution_string)
         try :
            shutil.rmtree(tmp_xml_dir)
         except OSError :
            pass
         print ('removed directory tree ' + tmp_xml_dir + ' if any')
         os.mkdir(tmp_xml_dir)
         print ('created directory ' + tmp_xml_dir)
         #execution_string = 'rm -f '+tmp_xml_dir+'*'
         #print (execution_string)
         #os.system(execution_string)

         if len(r_start) == 1 and len(c_start) == 1 : # no groups

            # generate commands if there are no groups
            start_dict = {}
            end_dict = {}
            start_tmp = 0
            end_tmp = 0
            new_params = {}
            new_output_name = tmp_xml_dir+tmp_out_name[1]		
            for i in range(n_ss):
               start_dict.update({i:end_tmp})
               if i < n_of_first_stacks:
                  new_params.update({i:first_size})
                  end_tmp += first_size
               else:
                  new_params.update({i:last_size})
                  end_tmp += last_size
               end_dict.update({i:end_tmp-1})
               tmp_string = execution_flag+' '+params_str+' --projout='+new_output_name+'-'+str(start_dict[i]).zfill(6)+'-'+str(end_dict[i]).zfill(6)+'.xml'+' --subvoldim='+str(new_params[i])+' --D0='+str(start_dict[i])+' --D1='+str(end_dict[i])+' --noprogressbar'
               cmd_string.update({i:tmp_string})

         else : # there are groups

            # create additional temporary directories
            for r in range(len(r_start)) :
               for c in range(len(c_start)) :
                  # Create temporary directory for temporary xml files 
                  gr_xml_dir = tmp_out_name[0]+gr_dir_names[r*len(c_start)+c]
                  #execution_string = 'mkdir -p ' + gr_xml_dir
                  #print (execution_string)
                  #os.system(execution_string)
                  try :
                     shutil.rmtree(gr_xml_dir)
                  except OSError :
                     pass
                  print ('removed directory tree ' + gr_xml_dir + ' if any')
                  os.mkdir(gr_xml_dir)
                  print ('created directory ' + gr_xml_dir)
                  #execution_string = 'rm -f '+gr_xml_dir+'*'
                  #print (execution_string)
                  #os.system(execution_string)

                  # generate commands for this group
                  start_dict = {}
                  end_dict = {}
                  start_tmp = 0
                  end_tmp = 0
                  new_params = {}
                  new_output_name = gr_xml_dir+tmp_out_name[1]		
                  for i in range(n_ss):
                     start_dict.update({i:end_tmp})
                     if i < n_of_first_stacks:
                        new_params.update({i:first_size})
                        end_tmp += first_size
                     else:
                        new_params.update({i:last_size})
                        end_tmp += last_size
                     end_dict.update({i:end_tmp-1})
                     tmp_string = execution_flag+' '+params_str+' --projout='+new_output_name+'-'+str(start_dict[i]).zfill(6)+'-'+str(end_dict[i]).zfill(6)+'.xml'+' --subvoldim='+str(new_params[i]).zfill(6)+' --D0='+str(start_dict[i]).zfill(6)+' --D1='+str(end_dict[i]).zfill(6)
                     tmp_string = tmp_string+' --R0='+str(r_start[r])+' --R1='+str(r_end[r])+' --C0='+str(c_start[c])+' --C1='+str(c_end[c])+' --noprogressbar --parallel'
                     if c < (len(c_start) - 1) : # it not the last group of columns: disable displacement computation of last column
                        tmp_string = tmp_string+' --disable_last_col'
                     if r < (len(r_start) - 1) : # it not the last group of columns: disable displacement computation of last column
                        tmp_string = tmp_string+' --disable_last_row'
                     cmd_string.update({(r*len(c_start)+c)*n_ss+i:tmp_string})

         ## start scoring
         #scores = score_function(new_params)

         # Sort tasks by score
         #elaborations = sort_elaborations(scores)
         #work_list = sort_work(cmd_string,elaborations)
         work_list = cmd_string

         # Call the routine to distribute the work between cpus
         dispatcher_step2(work_list)

         if len(r_start) == 1 and len(c_start) == 1 : # no groups

            #Collect all the xml files corresponding to subvolumes in a unique xml file (only for step 2)
            slash_pos = len(tmp_xml_dir) - 1
            if debug_level > 0 :
               execution_string = prefix+'mergedisplacements -d='+tmp_xml_dir[0:slash_pos]+' -o='+tmp_out_name[0]+tmp_out_name[1]+'.xml' + ' > ' + 'ouput_mdispls.out'
            else :
               execution_string = prefix+'mergedisplacements -d='+tmp_xml_dir[0:slash_pos]+' -o='+tmp_out_name[0]+tmp_out_name[1]+'.xml'
            print (execution_string)
            os.system(execution_string)
		
         else : # there are groups

            #Collect xml files of each group in a single xml file
            for dname in gr_dir_names :
               slash_pos = dname.find('/')
               suffix = dname[dname.find('_'):slash_pos]
               if debug_level > 0 :
                  execution_string = prefix+'mergedisplacements -d='+tmp_out_name[0]+dname[0:slash_pos]+' -o='+tmp_out_name[0]+'tmp/'+tmp_out_name[1]+suffix+'.xml' + ' > ' + 'ouput_' + suffix + '.out'
               else :
                  execution_string = prefix+'mergedisplacements -d='+tmp_out_name[0]+dname[0:slash_pos]+' -o='+tmp_out_name[0]+'tmp/'+tmp_out_name[1]+suffix+'.xml'
               print (execution_string)
               os.system(execution_string)
		
            #Collect all the xml files corresponding to groups in a unique xml file (only for step 2)
            slash_pos = len(tmp_xml_dir) - 1
            if debug_level > 0 :
               execution_string = prefix+'mergedisplacements --mgroups -d='+tmp_xml_dir[0:slash_pos]+' -o='+tmp_out_name[0]+tmp_out_name[1]+'.xml' + ' > ' + 'ouput_mgroups.out'
            else :
               execution_string = prefix+'mergedisplacements --mgroups -d='+tmp_xml_dir[0:slash_pos]+' -o='+tmp_out_name[0]+tmp_out_name[1]+'.xml'
            print (execution_string)
            os.system(execution_string)

         #clean up all temporary directories and files
         if len(r_start)>1 or len(c_start)>1 :
            for r in range(len(r_start)) :
               for c in range(len(c_start)) :
                  gr_xml_dir = tmp_out_name[0]+gr_dir_names[r*len(c_start)+c]
                  shutil.rmtree(gr_xml_dir)
                  print ('deleted directory ' + gr_xml_dir + ' and all files in it')
         shutil.rmtree(tmp_xml_dir)
         print ('deleted directory ' + tmp_xml_dir + ' and all files in it')

         #################### END STEP 2

      elif step6 :  

         ##################### STEP 6

         if info :
            # get info and print them
            (first_string, list_string, output_name, len_arr, final_string) = create_commands(nprocs-1,True) # 2017-02-06. Giulio. @ADDED output_name
         else :
            # Generate strings of commands
            (first_string, list_string, output_name, len_arr, final_string) = create_commands(nprocs-1) # 2017-02-06. Giulio. @ADDED output_name

            # 2017-02-08. Giulio. @ADDED support for suspend/resume
            if save_status_prefix == '' :
               save_status_prefix = output_name + '/'
            rs_fname = save_status_prefix + resume_status_fname
            if not os.path.exists(rs_fname) :
               # Execute on proc. n. 0 the string to initialize parallel computation
               if debug_level > 0 :
                  execution_string = prefix+first_string + " > " + "output_first.out"
               else :
                  execution_string = prefix+first_string
               os.system(execution_string)
               print (execution_string) # 2016-12-10. Giulio.
               if suspend_resume_enabled :
                  rs_file = open(rs_fname, 'wb')
                  pickle.dump([],rs_file)
                  rs_file.close()
		 
            # Prepare data for parallel computation
            cmd_string = list_string
            npoints = len_arr
            # Start scoring the task list
            scores = score_function(npoints)
            # Sort tasks by score
            elaborations = sort_elaborations(scores)
            work_list = sort_work(cmd_string,elaborations)
            # Call the routine to dinstibute the work between cpus
            dispatcher_step6(work_list,rs_fname) # 2017-02-06. Giulio. @ADDED parameter
            # Execute the command to merge the metadata
            if debug_level > 0 :
               execution_string = prefix+final_string + " > " + "output_final.out"
            else :
               execution_string = prefix+final_string
            os.system(execution_string)
            print (execution_string) # 2016-12-10. Giulio.
			 
         #################### END STEP 6
         
   else:
      if step2 :
         # Start worker sub.
         prefix = prefix + 'terastitcher '
         worker()
      elif step6 :
         if info :
            # do nothing
            dummy = 0
         else :
            # Start worker sub.
            worker()

   # End program
   comm.Barrier()
    
   if myrank == 0:
      t2 = time.time() 
      print ('*'*80)
      print (str(datetime.datetime.utcnow()), "-- Calculation ended after ", (t2-t1), " seconds")
      print ('*'*80) 
