#! /usr/bin/env python2.7
#***************************************************************************
#*   Copyright (C) 2013 by Edson Borin                                     *
#*   edson@ic.unicamp.br                                                   *
#*                                                                         *
#*   This program is free software; you can redistribute it and/or modify  *
#*   it under the terms of the GNU General Public License as published by  *
#*   the Free Software Foundation; either version 2 of the License, or     *
#*   (at your option) any later version.                                   *
#*                                                                         *
#*   This program is distributed in the hope that it will be useful,       *
#*   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
#*   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
#*   GNU General Public License for more details.                          *
#*                                                                         *
#*   You should have received a copy of the GNU General Public License     *
#*   along with this program; if not, write to the                         *
#*   Free Software Foundation, Inc.,                                       *
#*   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
#**************************************************************************/

import sys
import os.path
import shutil

# Try to import rdt and stats modules, if available.
import sys
# Add the perf_base_dir/scripts/misc to the system path to enable importing misc modules
sys.path.append("@PERFTEST_SCRIPTS_SRC_DIR@/misc_modules/")
try:
	import rdt, stats
	summarize_results=True
except ImportError, e: 
	print "WARNING: I wont summarize results because I could not import modules: ", e
	summarize_results=False

# == Adding new tests ==
#
# In order to add new tests one must:
# 1) import the test module (e.g.: import substruct_tst01.test). The module must
#    contain:
#  a) the functions short_description and long_description which return a string
#     containing a short/long description of the test.
#  b) a function called "run_test(ntimes)" that runs the performance test
#     "ntimes" times and returns (status, rdt_files). So that:
#     - status: the test status: 0 == OK, != 0 == ERROR
#     - rdt_files: a dictionary mapping rdt_ids to pairs (rdt_filename,
#       rdt_description).
#       * rdt_id is a unique string identifier for the rdt file.
#       * rdt_filename is the full path to the rdt file
#       * rdt_description is a short description of the what is being measured
#         by the rdt file.
#
# 2) add the test module to one (or more) of the test lists.

# == Import tests modules ==
import substruct_tst01.test
import substruct_tst02.test
import substruct_tst03.test
import substruct_tst04.test
import substruct_tst05.test
import substruct_tst06.test
import substruct_tst07.test
import substruct_tst08.test
import substruct_tst09.test
import substruct_tst10.test
import substruct_tst11.test
import substruct_tst12.test
import substruct_tst13.test

import skyline_tst01.test  # s (p2) s Chlsk+sfwd+sbck
import skyline_tst02.test  # m (p3) m Chlsk+sfwd+sbck???
import skyline_tst03.test  # m (p4) m Chlsk+sfwd+sbck
import skyline_tst04.test  # s (p2) s LDLt
import skyline_tst05.test  # s (p3) s LDLt
import skyline_tst06.test  # m (p4) m LDLt
import skyline_tst07.test  # s (p2) s multadd
import skyline_tst08.test  # m (p3) m multadd
import skyline_tst09.test  #   (p7) l sor????
import skyline_tst10.test  #   (p3) m sor????
import skyline_tst11.test  #   (p6) l sor????
import skyline_tst12.test  #   (p5) m Chlsk+sfwd+sbck???
# =========================

# == Available test lists ==
# See Usage for the description of the test lists

# Tests with execution time shorter than 1 minute
short_tests = [("substruct_tst01",substruct_tst01.test), 
	       ("substruct_tst02",substruct_tst02.test),
	       ("skyline_tst01",skyline_tst01.test),  # ~1s
	       ("skyline_tst04",skyline_tst04.test),  # ~1s
	       ("skyline_tst05",skyline_tst05.test),  # ~28s
	       ("skyline_tst07",skyline_tst07.test)]  # ~10s

other_tests = [("skyline_tst09",skyline_tst09.test),
               ("skyline_tst10",skyline_tst10.test),
               ("skyline_tst11",skyline_tst11.test),
               ("skyline_tst12",skyline_tst12.test)]

# TODO: How about:
# substruct_tst05
# substruct_tst14
# substruct_tst15
# substruct_tst16

# Tests with execution time between 1 and 10 minutes
medium_tests= [("substruct_tst03",substruct_tst03.test),
	       ("substruct_tst04",substruct_tst04.test),
	       ("substruct_tst07",substruct_tst07.test),
	       ("substruct_tst08",substruct_tst08.test),
	       ("substruct_tst09",substruct_tst09.test),
	       ("substruct_tst11",substruct_tst11.test),
	       ("substruct_tst13",substruct_tst13.test),
	       ("skyline_tst02",skyline_tst02.test),  # ???
	       ("skyline_tst03",skyline_tst03.test),  # ~190s
	       ("skyline_tst06",skyline_tst06.test),  # ~312s
	       ("skyline_tst08",skyline_tst08.test)]  # ~286s

# Tests with execution time longer than 10 minutes
long_tests = [
#	      ("skyline_tst06",skyline_tst06.test),
	      ("substruct_tst06",substruct_tst06.test),
	      ("substruct_tst10",substruct_tst10.test),
              ("skyline_tst09",skyline_tst09.test),
              ("skyline_tst11",skyline_tst11.test)]

short_regression_tests = short_tests
regression_tests       = short_tests + medium_tests
full_regression_tests  = short_tests + medium_tests + long_tests

all_tests = short_tests + medium_tests + long_tests + other_tests
# =========================

# Default value for ntimes
ntimes_dft = 3

def error(message, status):
	sys.stderr.write('ERROR: '+message+'\n')
	if status != 0 : sys.exit(status)

def get_test(v) : 
	for t in all_tests : 
		if t[0] == v : 
			return t
	error("Could not find test " + str(v))

# Functions for stand alone tests
def usage():
	print "\nUsage: test.py -t test_name [-n #times] [-a|-s|-m|-l|-r] [-t test_name] [-h]"
	print "\nARGUMENTS"
	print "\t-n #times    : Run each test #times times. (default = ",ntimes_dft,")"
	print "\t-a           : Run all tests."
	print "\t-s           : Run short tests."
	print "\t-m           : Run medium tests."
	print "\t-l           : Run long tests."
	print "\t-r res_dir   : Run daily regression tests and move results to \"res_dir\" directory."
	print "\t-F res_dir   : Run full regression tests and move results to \"res_dir\" directory."
	print "\t-S res_dir   : Run short regression tests and move results to \"res_dir\" directory."
	print "\t-t test_name : Run test test_name."
	print "\nDESCRIPTION"
	print "\tExecutes a set of performance tests. The following tests are available:"
	print "\tShort tests:"
	for t in short_tests :
		print "\t* ", t[0], ":", t[1].short_description()
	print "\tMedium tests:"
	for t in medium_tests :
		print "\t* ", t[0], ":", t[1].short_description()
	print "\tLong tests:"
	for t in long_tests :
		print "\t* ", t[0], ":", t[1].short_description()
	print "\tRegression tests:"
	for t in regression_tests :
		print "\t* ", t[0], ":", t[1].short_description()
	print "\tFull regression tests:"
	for t in full_regression_tests :
		print "\t* ", t[0], ":", t[1].short_description()
	print "\tShort regression tests:"
	for t in short_regression_tests :
		print "\t* ", t[0], ":", t[1].short_description()
	sys.exit(1)

# Main - for stand alone tests only
if __name__ == "__main__":
	import getopt
	results_dir=0
	ntimes=ntimes_dft
	tests_to_run = {}
	# Process arguments
	try :
		opts, extra_args = getopt.getopt(sys.argv[1:], 't:n:hasmlr:S:F:')
	except getopt.GetoptError, e:
		error(str(e), 1)
	for f, v in opts:
		if   f == '-a': 
			for t in short_tests      : tests_to_run[t[0]] = t
			for t in medium_tests     : tests_to_run[t[0]] = t
			for t in long_tests       : tests_to_run[t[0]] = t
			for t in regression_tests : tests_to_run[t[0]] = t
			for t in full_regression_tests : tests_to_run[t[0]] = t
			for t in short_regression_tests : tests_to_run[t[0]] = t
		elif f == '-n': ntimes=int(v)
		elif f == '-t': tests_to_run[v] = get_test(v)
		elif f == '-s':
			for t in short_tests      : tests_to_run[t[0]] = t
		elif f == '-m':
			for t in medium_tests     : tests_to_run[t[0]] = t
		elif f == '-l':
			for t in long_tests       : tests_to_run[t[0]] = t
		elif f == '-r':
			for t in regression_tests : tests_to_run[t[0]] = t
			results_dir = v
		elif f == '-F':
			for t in full_regression_tests : tests_to_run[t[0]] = t
			results_dir = v
		elif f == '-S':
			for t in short_regression_tests : tests_to_run[t[0]] = t
			results_dir = v
		elif f == '-h': usage()

	all_results={}

	# Run tests
	for f, t in tests_to_run.iteritems() :
			obj=t[1]
			test_name=t[0]
			try:
				status,rdt_files = obj.run_test(ntimes)
				all_results[test_name]=(status,obj.short_description(),obj.long_description(),rdt_files)
			except:
				error('Could not run test '+test_name,0)
				all_results[test_name]=(-1,obj.short_description(),obj.long_description(),{})

	# Move/Print results
	for k, v in all_results.iteritems() :
		status = v[0]
		s_desc = v[1]
		l_desc = v[2]
		rdt_files = v[3]
		test_name = k
		print '** ' + test_name + ' **'
		print 'desc: ', s_desc
		# Print results
		if status != 0: 
			print "Status [FAILED] ("+str(status)+")"
		else :
			print "Status [OK]"
			print "Results summary ----------------------------"
			
			for rdt_id,v in rdt_files.iteritems() : 
				if summarize_results :
					try: 
						fn=v[0]
						rdt_d=rdt.read(fn)
						elapsed_list=rdt.get_column_values(rdt_d,"ELAPSED")
						try:
							av=stats.average(elapsed_list)
						except stats.StatsError, e:
							print "WARNING: Could not compute average for results at", fn, "(", e, ")"
							av=0.0
						try:
							ci=stats.conf_int(elapsed_list, 95.0)
						except stats.StatsError, e:
							print "WARNING: Could not compute confidence interval for results at", fn, "(", e, ")"
							ci=0.0
					except rdt.RdtError, e:
						print "WARNING: error when summarizing results for", fn, "(", e, ")"
						av=0.0
						ci=0.0
					print '{0:15s} : {1:>16f} +- {2:<16f} : {3:s}'.format(rdt_id, av, ci, v[1])
				else:
					print '{0:15s} : {1:s} : {2:s}'.format(rdt_id,v[0],v[1])
			print "--------------------------------------------"
		# Move results
		if results_dir != 0 :
			# Record results to results_dir. For each test, create a directory,
			# record the test status and copy the rdt result files to the directory.
			result_dir = os.path.join(results_dir,test_name)
			if not os.path.isdir(result_dir) :
				try:    
					os.makedirs(result_dir)
				except os.error, e: 
					warning(str(e))
					continue
			result_info = os.path.join(results_dir,test_name+".info")
			try:
				f = open(result_info, 'w+')
				f.write("test_name : "+test_name+"\n")
				f.write("test_desc : "+s_desc+"\n")
				f.write("test_status : "+str(status)+"\n")
				# Copy rdt files
				for rdt_id,v in rdt_files.iteritems() : 
					f.write(rdt_id+" : "+v[1]+"\n")
				f.close()
			except IOError, e:
				warning(str(e))
				continue
			# Copy rdt files
			for rdt_id,v in rdt_files.iteritems() : 
				rdt_fn = v[0]
				rdt_dsc = v[1]
				# Copy rdt file (rdt_fn) to result_dir
				try: 
					shutil.copy2(rdt_fn,result_dir)
				except (IOError, os.error) as why :
					warning(str(why))
					continue
