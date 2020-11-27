#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
This program uses a dispatcher/worker approach to consume a queue 
of elaborations using teraconverter
Copyright (c) 2016:
Massimiliano Guarrasi (1), Giulio Iannello (2), Alessandro Bria (2)
(1): CINECA
(2): University Campus Bio-Medico of Rome
The program was made in the framework of the HUMAN BRAIN PROJECT.
All rights reserved.


RELEASE 2.3.6


EXAMPLE of usage (X is the major version, Y is the minor version, Z is the patch):
mpirun -np XX python paraconverterX.Y.Z.py -s=source_volume -d=destination_path --depth=DD --height=HH --width=WW --sfmt=source_format --dfmt=destinatiopn_format --resolutions=RR 

where:
- XX is the desided level of parallelism plus 1 (for the dispatcher process)
- DD, HH, WW are the values used to partition the image for parallel execution
- source and destination format are allowed formats for teraconverter
- RR are the requested resolutions (according to the convention used by teraconverter)
See teraconverter documentation for more details


*******************************
*        Change Log           *
*******************************

v2.3.6 2020-02-08
- changed names of functions: introduced the termnology dispatcher/workers/launchers

v2.3.5 2020-02-08
- fixed bug (for python3) at line 902: integer division is required

v2.3.4 2020-01-15
- adapted to run also on Python 3 interpreters: changed <> to !=, use parentheses in 
  print commands, converted 'keys' and 'values' of dictionaries to lists

v2.3.3 2019-12-01
- added support for fixed_tiling

v2.3.2 2017-10-07
- added management of --isotropic option in the partition algorithm
- corrected a bug in function 'collect_instructions'

v2.2.2 2017-10-07
- revisted platform dependent instructions

v2.2.1 2017-09-19
- added option --info to display the memory needed in GBytes without performing any 
  conversion

v2.2 2017-03-12
- the suspend/resume mechanism can be disabled by changing the value of variable
  'suspend_resume_enabled' (the mechanism is enebled if True, disabled if False
- changed the policy to manage dataset partition and eliminated additional parameter
  to specify the desired degree of parallelism which is now directly passed by the main

v2.1 2017-02-06
- implemented a suspend/resume mechanism
  the mechanism can slow down parallel execution if the dataset chunks are relatively 
  small to avoid this a ram disk can be used to save the status (substitute the name 
  'output_nae'   at line 953 with the path of the ram disk)

v2.0 2016-12-10
- dataset partitioning takes into account the source format in order to avoid that the 
  same image region is read by different TeraConverter instances; requires an additional 
  parameter in the command line (see EXAMPLE of usage above)
"""


import os
import sys
import time
import datetime
import operator
import math

from glob import glob
from mpi4py import MPI
from collections import deque
from subprocess import *

import os.path
import pickle

"""
Set prefix='./' if your teraconverter executable is present in the running director.
If  the teraconverter executable is included in the $PATH dir select prefix=''
"""
#prefix = './'
prefix = ''

"""
*************************
* SUSPEND/RESUME option *
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




###############################################################################
#                                                                             #
# ----------------- Function to parallelize the code ------------------------ #       
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
   print ("Scheduled job n. ", list(input_file.keys())[0], " is executed by rank: ", myrank)  
   # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   execution_string = prefix+ list(input_file.values())[0]
   print (execution_string)
   # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Sleep da cancellare!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   #time.sleep((5.+100./(input_file.keys()[0]+1.))/100.)
   # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   os.system(execution_string)
   t2 = time.time()
   print (" ---> Processor ",  myrank, " has calculated for " , t2-t1) 
   return input_file


def dispatcher(queue,rs_fname):  # 2017-02-06. Giulio. @ADDED rs_fname
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
# -------------------- Paraconverter Service Functions ---------------------- #       
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
   input_name = read_item(params, '-s=', './input.xml')
   # Find standard output directory
   output_name = read_item(params, '-d=', './OUT')
   # Declare axprossimative depth for the tiles
   wb1 = read_item(params, '--depth=', 256)
   # Declare axprossimative height for the tiles
   wb2 = read_item(params, '--height=', 256)
   # Declare axprossimative width for the tiles
   wb3 = read_item(params, '--width=', 256)
   # Declare source format
   sfmt = read_item(params, '--sfmt=', '"TIFF (unstitched, 3D)"')
   # Declare destination format
   dfmt = read_item(params, '--dfmt=', 'RGB')
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
   # check fixed_tiling flag
   fixed_tiling = read_item(params, '--fixed_tiling', 'False')
   if fixed_tiling == '' : # isotropic has been set
      params.append('--fixed_tiling') # the flag must be maintained in command line
      fixed_tiling = True
   else :
      fixed_tiling = False
   
   #------------- Added by Giulio: generate a temporary file with image size information
   
   info_string = 'teraconverter ' + ' --sfmt=' + sfmt + ' --dfmt=' + dfmt + ' -s="' + input_name + '" -d=/'
   #pwd = check_output("pwd",shell=True)
   #execution_string = prefix+info_string + ' --info="' + pwd.strip("\n") + '/__dims__.txt"'
   execution_string = prefix+info_string + ' --info="' + os.getcwd() + '/__dims__.txt"'
   os.system(execution_string)
   print (execution_string)

   #------------- End added by Giulio
   
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
   return (input_name, output_name, wb1, wb2, wb3, sfmt, dfmt, iresolutions, max_res, isotropic, max_res_D, params, last_string, height, width, depth, fixed_tiling, bytes_x_chan, n_chans)



def read_item(input_arr, item, default, message=True):
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
   tmp = check_flag(input_arr, item, True)
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
    k_h = math.floor(D/(math.pow(2.,n) * b_h))
    b = b_h
    r = D % (math.pow(2.,n) * b_h)
    k = k_h

    itera = 0
    verif = bool(1)
    while verif:
       # step 2
       if (D % (math.pow(2.,n) * b_h)) == 0:
           b = b_h
           r = 0
           k = k_h
           verif = bool(0) # exit form while cicle
       # step 3
       elif (D % (math.pow(2.,n) * b_h)) > r:
           b = b_h;
           r = (D % (math.pow(2.,n) * b_h));
           k = k_h;
       # step 4
       if h == math.floor(w/2):
           verif = bool(0) # exit form while cicle
       h = min( math.floor(w/2), h + max(1,math.floor(b_h - D/(math.pow(2.,n) * (k_h+1)))) );
       b_h = w - h;
       k_h = math.floor(D/(math.pow(2.,n) * b_h));
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


def create_sizes(size, wb, max_res, fixed_tiling=False, norest=False):
   """
   Create a 3D array containing the size for each tile on the desidered direction
   Input: 
      size = size (in pixel) of the input immage
      wb = Rough depth for the tiles in the desidered direction
      max_res = Maximum level of resolution available (integer)
      norest = Boolean variable to chech if we need of the last array element (if it is different from the preavious one)
      fixed_tiling = Boolean variable to force that all tiles except possibly the last one have a size that is equal to wb*2^max_res
   Output:
      arr = Array containing the size for each tile on the desidered direction
   """
   if fixed_tiling :
      # make the partition
      b  = wb
      wb = int(math.pow(2,max_res)*b)
      k  = size // wb # it must be integer division
      r  = size - (wb * k)
      if r < math.pow(2,max_res) :
         r = 0
      arr = prep_array(wb, r, k)
   else :
      # Make the partition
      values = opt_algo(size, wb, max_res)
      b = values[0]
      r = values[1]
      k = values[2]
      itera = values[3]
      wb = int(math.pow(2,max_res)*b)
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
   (input_name, output_name, wb1, wb2, wb3, sfmt, dfmt, iresolutions, max_res, isotropic, max_res_D, params, last_string, height, width, depth, fixed_tiling, bytes_x_chan, n_chans) = read_params()
   print("#"*80)
   print("Input file = ", input_name)   
   print("Output directory", output_name)
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
   size_1 = create_sizes(depth, wb1, max_res_D,fixed_tiling) # max_res_D = max_res if isotropic is not set of if voxels are isotropic
   size_2 = create_sizes(height, wb2, max_res,fixed_tiling)
   size_3 = create_sizes(width, wb3, max_res,fixed_tiling)
   
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
    
   if myrank == 0:
      # timing
      t1 = time.time() 
      print ('*'*80)
      print (str(datetime.datetime.utcnow()), " -- Calculation started on ", nprocs, "- 1 cores.")
      print ('*'*80)
   # Sincronizzation
   comm.Barrier()
   if myrank == 0:
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
         dispatcher(work_list,rs_fname) # 2017-02-06. Giulio. @ADDED parameter
         # Execute the command to merge the metadata
         execution_string = prefix+final_string
         os.system(execution_string)
         print (execution_string) # 2016-12-10. Giulio.
   else:
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
