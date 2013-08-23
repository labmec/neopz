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
import substruct_tst1.test
import substruct_tst2.test

ntimes_dft = 3

# Available tests
short_tests = [("substruct_tst1",substruct_tst1.test), 
	       ("substruct_tst2",substruct_tst2.test)]
medium_tests = []
long_tests = []

tests = short_tests + medium_tests + long_tests

def error(message, status):
	sys.stderr.write('ERROR: '+message+'\n')
	if status != 0 : sys.exit(status)

# Functions for stand alone tests
def usage():
	print "\nUsage: test.py -t test_name [-a] [-h] [-n #times]\n"
	print "\nARGUMENTS"
	print "\t-a : Run all tests."
	print "\t-t test_name : Run test test_name."
	print "\t-n #times    : Run each test #times times. (default = ",ntimes_dft,")"
	print "\t-s           : Run all short tests."
	print "\t-m           : Run all medium tests."
	print "\t-l           : Run all long tests."
	print "\nDESCRIPTION"
	print "\tExecute a set of performance tests. The following tests are available:"
	print "\tShort tests:"
	for t in short_tests :
		print "\t* ", t[0], ":", t[1].description
	print "\tMedium tests:"
	for t in medium_tests :
		print "\t* ", t[0], ":", t[1].description
	print "\tLong tests:"
	for t in long_tests :
		print "\t* ", t[0], ":", t[1].description
	sys.exit(1)

# Main - for stand alone tests only
if __name__ == "__main__":
	import getopt
	allt=0
	ntimes=ntimes_dft
	test = {}
	# Process arguments
	try :
		opts, extra_args = getopt.getopt(sys.argv[1:], 't:an:hsml')
	except getopt.GetoptError, e:
		error(str(e), 1)
	for f, v in opts:
		if   f == '-a': allt=1
		elif f == '-n': ntimes=int(v)
		elif f == '-t': test[v] = True
		elif f == '-s':
			for t in short_tests : test[t[0]] = True
		elif f == '-m':
			for t in medium_tests : test[t[0]] = True
		elif f == '-l':
			for t in long_tests : test[t[0]] = True
		elif f == '-h': usage()

	all_results={}

	# Run tests
	for t in tests :
		if allt or (t[0] in test) :
			obj=t[1]
			try:
				status,results = obj.run_test(ntimes)
				all_results[t[0]]=(status,obj.description,results)
			except:
				error('Could not run test '+t[0],0)
				all_results[t[0]]=(-1,obj.description,{})

	# Print all results
	for k, v in all_results.iteritems() :
		status = v[0]
		desc = v[1]
		res = v[2]
		print '** ' + k + ' **'
		print 'desc: ', desc
		if status != 0: 
			print "Status [FAILED]"
		else :
			print "Status [OK]"
			print "Results summary ----------------------------"
			for rk,rv in res.iteritems() : 
				print '{0:10s} : {1:>16f} +- {2:<16f}'.format(rk, rv[0], rv[1])
			print "--------------------------------------------"
