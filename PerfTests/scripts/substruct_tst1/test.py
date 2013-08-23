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
import shlex, subprocess
import resource

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

# Variables to be defined by cmake
builddir="@PERFTEST_BASE_DIR@"
datadir="@PERFTEST_SMALL_DATA_DIR@"

assfn=("ass", "ass.rdt")
crefn=("cre", "cre.rdt")
prefn=("pre", "pre.rdt")
solfn=("sol", "sol.rdt")
totfn=("tot", "tot.rdt")
tpzdohrassfn=("dohrass", "tpzdohrass.rdt")
tpzdohrdecfn=("dohrdec", "tpzdohrdec.rdt")

# List of rdt files produced by the test
rdtfiles_l=[assfn, crefn, prefn, solfn, totfn, tpzdohrassfn, tpzdohrdecfn]

def error(message, status):
	sys.stderr.write('ERROR (test.py): '+message+'\n')
        sys.exit(status)

# Setup the command line
def setup_cmd():
	# Check build directory
	if not os.path.isdir(builddir) :
		error(builddir+' is an invalid build directory.', 1)
	# Check run directory
	rundir = os.path.join(builddir,'scripts','substruct_tst1')
	if not os.path.isdir(rundir) :
		error(rundir+' is an invalid run directory.', 1)
	if not os.path.isdir(builddir) :
		error(builddir+' is an invalid build directory.', 1)
	# Check executable
	executable=os.path.join(builddir,"progs","substruct", "substruct-perf")
	if not os.path.isfile(executable) :
		error(executable+' is an invalid executable file name.', 1)
	# Check input file
	inputfn = os.path.join(datadir,"substruct","inputs","Cubo986.msh")
	if not os.path.isfile(inputfn) :
		error(inputfn+' is an invalid input file name.', 1)	
	# Put the arguments together
        arguments = ' -mc '+inputfn
	arguments = arguments + ' -nsub 64'
	arguments = arguments + ' -nt_a 1' 
	arguments = arguments + ' -nt_d 1' 
	arguments = arguments + ' -nt_m 1' 
	arguments = arguments + ' -nt_sm 1' 
	arguments = arguments + ' -p 2' 
	arguments = arguments + ' -ass_rdt ' + assfn[1]
	arguments = arguments + ' -cre_rdt ' + crefn[1]
	arguments = arguments + ' -pre_rdt ' + prefn[1]
	arguments = arguments + ' -sol_rdt ' + solfn[1]
	arguments = arguments + ' -tot_rdt ' + totfn[1]
	arguments = arguments + ' -tpz_dohr_ass ' + tpzdohrassfn[1]
	arguments = arguments + ' -tpz_dohr_dec ' + tpzdohrdecfn[1]
	# TODO: Add arguments to enforce output checking!
	return rundir, executable+arguments

# Limits for this test
limits = { "cpu"   : (resource.RLIMIT_CPU,     60, "Max CPU time in seconds"), 
	   "nofile": (resource.RLIMIT_NOFILE,   7, "The maximum number of open file descriptors for the current process."),
#	   "rss"   : (resource.RLIMIT_RSS,   1024, "The maximum resident set size that should be made available to the process"),
#	   "fsize" : (resource.RLIMIT_FSIZE,    1, "Max size of a file which the process may create"),
#	   "data"  : (resource.RLIMIT_DATA,  1024, "The maximum size (in bytes) of the process's heap"),
#	   "nproc" : (resource.RLIMIT_NPROC,    0, "The maximum number of processes the current process may create")
	 }

# Set the rlimits of the chidren process (see limits above)
# TODO: Improve the handling of sandboxing limits
def setlimits():
	print "Setting resource limit in child"
	for k, v in limits.iteritems() : 
		resource.setrlimit(v[0], (v[1],v[1])) 
		#print k, " : ", v[0], " => ", v[1]

# Sumarizes the RDT (Raw data table) results (generated by the run)
def sumarize_rdt_results(rundir) :
	results = {}
	# Compute average and confidence interval for each rdt file.
	for f in rdtfiles_l : 
		k =f[0] # Step name
		fn=f[1] # RDT file name
		if summarize_results :
			try: 
				rdtfn=os.path.join(rundir,fn)
				rdt_d=rdt.read(rdtfn)
				elapsed_list=rdt.get_column_values(rdt_d,"ELAPSED")
				av=stats.average(elapsed_list)
				ci=stats.conf_int(elapsed_list, 95.0)
			except stats.StatsError, e:
				print "WARNING: error when summarizing results for", fn, "(", e, ")"
				av=0.0
				ci=0.0
		else :
			av=0.0
			ci=0.0
		results[k]=(av,ci)
	return results

# Execute the test.
def run_test(ntimes):
	rundir,cmd=setup_cmd()
	args = shlex.split(cmd)
	sout = None
	serr = None
	for i in range(ntimes) : 
		p = subprocess.Popen(args, preexec_fn=setlimits, stdout=sout, stderr=serr, cwd=rundir)
		p.wait()
		if (p.returncode != 0) : 
			return p.returncode, {}
	results = sumarize_rdt_results(rundir)
	return 0, results

# Functions for stand alone tests
def usage():
	print "\nUsage: test.py -r [-h]\n"
	print "\nARGUMENTS"
	print "\t-r : Run the experiment."
	print "\nDESCRIPTION"
	print "\tExecute the substruct tool collecting statistics for the following steps:"
	print "\t ", assfn[0], ": assembling the system (serial) -- results at", assfn[1]
	print "\t ", tpzdohrassfn[0], ": assembling (ass part) the system (serial) -- results at", tpzdohrassfn[1]
	print "\t ", tpzdohrdecfn[0], ": assembling (dec part) the system (serial) -- results at", tpzdohrdecfn[1]
	print "\t ", crefn[0], ": creating the sytem (serial) -- results at", crefn[1]
	print "\t ", prefn[0], ": pre-processing (serial) -- results at", prefn[1]
	print "\t ", solfn[0], ": solver (serial) -- results at", solfn[1]
	print "\t ", totfn[0], ": total -- results at", totfn[1]
	sys.exit(1)

# Main - for stand alone tests only
if __name__ == "__main__":
	import getopt
	run=0
	ntimes=1
	# Process arguments
	try :
		opts, extra_args = getopt.getopt(sys.argv[1:], 'rn:h')
	except getopt.GetoptError, e:
		error(str(e), 1)
	for f, v in opts:
		if   f == '-r': run=1
		elif f == '-n': ntimes=int(v)
		elif f == '-h': usage()

	# Run test
	if run == 1: 
		status,results = run_test(ntimes)
		if status == 0: print "Execution [OK]"
		else          : print "Execution [FAILED]"
		print "Results summary ----------------------------"
		for k,v in results.iteritems() : print '{0:10s} : {1:>16f} +- {2:<16f}'.format(k, v[0], v[1])
		print "--------------------------------------------"
	else:
		print "WARNING: No options provided. (use -h for help)"
