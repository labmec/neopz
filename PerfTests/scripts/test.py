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

# Functions for stand alone tests
def usage():
	print "\nUsage: test.py -t test_name [-a] [-h] [-n #times]\n"
	print "\nARGUMENTS"
	print "\t-a : Run all tests."
	print "\t-t test_name : Run test test_name."
	print "\t-n #times    : Run each test #times times."
	print "\t-s           : Run all short tests."
	print "\t-m           : Run all medium tests."
	print "\t-l           : Run all long tests."
	print "\nDESCRIPTION"
	print "\tExecute the following performance tests:"
	print "\t TODO -- describe them "
	sys.exit(1)

# Main - for stand alone tests only
if __name__ == "__main__":
	import getopt
	allt=0
	ntimes=1
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
			# Add short tests here
			test["substruct_tst1"] = True
		#elif f == '-m':
			# Add medium tests here
		#elif f == '-l':
			# Add long tests here
		elif f == '-h': usage()

	# Run tests
	if allt or ("substruct_tst1" in test) :
		try:
			status,results = substruct_tst1.test.run_test(ntimes)
			if status == 0: print "Execution [OK]"
			else          : print "Execution [FAILED]"
			print "Results summary ----------------------------"
			for k,v in results.iteritems() : print '{0:10s} : {1:>16f} +- {2:<16f}'.format(k, v[0], v[1])
			print "--------------------------------------------"
		except:
			error("Could not run test substruct-tst1")
