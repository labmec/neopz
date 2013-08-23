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

class RdtError(Exception):
	def __init__(self, msg):
		self.msg = msg
	def __str__(self):
		return self.msg

# Handles RDT (raw data table) files
# TODO: Improve the description...
#
# The RDT file is a CSV (comma separated value) file.
# The first row contains the headers of the columns.
# The remaining rows contains the data for each row. Example:
# RUN, ELAPSED
# 1, 12.2
# 2, 43.2
# 3, 87.3

# Read a RDT file into a dictionary.
# The dictionary is indexed by the column index.
# The value is a 2-entry tuple describing the column header and the list of
#  values at the column
def read(filename):
	# Check configuration file
	if not os.path.isfile(filename): 
		raise RdtError(filename+' is not a valid rdt file.')
	d={}
	with open(filename) as f:
		# Process the first row to setup the headers
		hlist=f.readline().split(',')
		hi=1
		for h in hlist:
			header = h.strip()
			d[hi] = (header,[])
			hi=hi+1
		# Process the remaining rows
		for line in f:
			vlist = line.split(',')
			first = vlist[1]
			if first.strip() != "#":
				hi=1
				for v in vlist:
					d[hi][1].append(float(v.strip())) 
					hi=hi+1
	return d

def get_column_values(rdt_d,field):
	for k, v in rdt_d.iteritems():
		if v[0] == field : 
			return v[1]
	return []

# Functions for stand alone tests
def usage():
	print "\nUsage: rdt.py -i input.rdt [-a] [-f field] [-h]\n"
	print "\nARGUMENTS"
	print "\t-i input.rdt: input rdt file."
	print "\t-f field: print selected column."
	print "\t-a : print all columns."
	print "\t-h : this help :-)."
	print "\nDESCRIPTION"
	print "\tReads the contents of a RDT file and prints its contents."
	sys.exit(1)

def error(message, status):
	sys.stderr.write('ERROR: '+message+'\n')
        sys.exit(status)

# Main - for stand alone tests only
if __name__ == "__main__":
	import getopt
	inputfn=0
	field=0
	printa=False
	opts, extra_args = getopt.getopt(sys.argv[1:], 'i:f:ah')
	for f, v in opts:
		if   f == '-i': inputfn=v
		elif f == '-f': field=v
		elif f == '-a': printa=True
		elif f == '-h': usage()

	if inputfn == 0:
		error("An input file name must be provided.", 1)

	try:
		rdt_d = read(inputfn)
	except RdtError, e:
		error(str(e), 1)

	if printa :
		print "-- All values ---------"
		print rdt_d
		print "-----------------------"
		
	if field != 0 :
		values = get_column_values(rdt_d,field)
		print field, " : ", values
