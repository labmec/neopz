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

def error(message, status):
	sys.stderr.write('ERROR (test.py): '+message+'\n')
        sys.exit(status)

#  (rdt_id, rdt_opt, rdt_filename, rdt_description)
assfn=("ass", "-ass_rdt", "ass.rdt", "Assemble: dohrstruct->Assemble(...). Assemble element matrices and decompose matrices.")
crefn=("cre", "-cre_rdt", "cre.rdt", "Create: dohrstruct->Create()")
solfn=("sol", "-sol_rdt", "sol.rdt", "Solver: cg.Solve(...)")
totfn=("tot", "-tot_rdt", "tot.rdt", "Total: all the steps")
tpzdohrassfn=("dohrass", "-tpz_dohr_ass", "tpzdohrass.rdt", "Assemble element matrices")
tpzdohrdecfn=("dohrdec", "-tpz_dohr_dec", "tpzdohrdec.rdt", "Decompose matrices")
# List of rdt files produced by the test
rdtfiles_l=[assfn, crefn, solfn, totfn, tpzdohrassfn, tpzdohrdecfn]

# Setup the command line
def setup_cmd():
	# Check build directory
	if not os.path.isdir(builddir) :
		error(builddir+' is an invalid build directory.', 5)
	# Check run directory
	rundir = os.path.join(builddir,'scripts','substruct_tst14')
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
	arguments = arguments + ' -p 2' 
	#NUMA aware Dohrman Assembly List thread work objects re-allocation.
	arguments = arguments + ' -naDALora'
	#NUMA aware Dohrman Assembly List thread work objects re-allocation threshold.
	#arguments = arguments + ' -naDALorat 1835008' # 2/2MB(l2) + 6/8MB(l3)
	#NUMA aware (node round-robin) Dohrman Assembly List thread work scheduling.
	arguments = arguments + ' -naDALtws' 
	for rdtarg in rdtfiles_l :
		arguments = arguments + ' ' + rdtarg[1] + ' ' + rdtarg[2]
	# TODO: Add arguments to enforce output checking!
	return rundir, executable+arguments

# Limits for this test
limits = { "cpu"   : (resource.RLIMIT_CPU, 3600, "Max CPU time in seconds"), 
#	   "nofile": (resource.RLIMIT_NOFILE,     7, "The maximum number of open file descriptors for the current process."),
#	   "rss"   : (resource.RLIMIT_RSS,     1024, "The maximum resident set size that should be made available to the process"),
#	   "fsize" : (resource.RLIMIT_FSIZE,      1, "Max size of a file which the process may create"),
#	   "data"  : (resource.RLIMIT_DATA,    1024, "The maximum size (in bytes) of the process's heap"),
#	   "nproc" : (resource.RLIMIT_NPROC,      0, "The maximum number of processes the current process may create")
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
		fn=f[2] # RDT file name
		if summarize_results :
			try: 
				rdtfn=os.path.join(rundir,fn)
				rdt_d=rdt.read(rdtfn)
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
		else :
			av=0.0
			ci=0.0
		results[k]=(av,ci)
	# erase time files
	for arq in rdtfiles_l:
		os.remove(arq[1])
	return results

# Sumarizes the RDT (Raw data table) files information
def sumarize_rdt_files(rundir) :
	results = {}
	for f in rdtfiles_l : 
		rdt_id  = f[0]   # Step name
		rdt_fn  = os.path.join(rundir,f[2]) # RDT file name
		rdt_dsc = f[3]   # Description
		results[rdt_id] = (rdt_fn, rdt_dsc)
	return results

def short_description() : return "substructure -- cubo986.msh -- parallel --  different number of substructures and number of threads with Numa naDALora and naDALtws"

# Execute the test.
def run_test(ntimes, nsub, thread):
	rundir,cmd=setup_cmd()
        print short_description() 
	print "CMD:",cmd
	# substructures
	cmd = cmd + " -nsub " + nsub
	# threads
	cmd = cmd + ' -nt_a ' + thread 
	cmd = cmd + ' -nt_d ' + thread
	cmd = cmd + ' -nt_m ' + thread
	cmd = cmd + ' -nt_sm ' + thread
	#print cmd
	args = shlex.split(cmd)
	sout = subprocess.PIPE # redirect output (before: None)
	#sout = None
	serr = None
	for i in range(ntimes) : 
		p = subprocess.Popen(args, preexec_fn=setlimits, stdout=sout, stderr=serr, cwd=rundir)
		p.wait()
		if (p.returncode != 0) : 
			return p.returncode, {}
	#results = sumarize_rdt_files(rundir)
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

	subs = [ "4", "8", "12", "16", "20", "24", "28", "32", "36", "48", "64", "128"]
	numthreads = [ "1" , "6", "12", "18", "24", "32", "48" , "64"]
	# Run test
	if run == 1: 
		for thr in numthreads:
			print "# Results Summary"
			print '{0:s}; {1:s}; '.format("nsub", "threads"),
			for arq in rdtfiles_l:
				print'{0:s}; {1:s};'.format(arq[0],"err_"+arq[0]),
			for item in subs:
				status,results = run_test(ntimes, item, thr)
				#if status == 0: print "Execution [OK]"
				if status != 0 : print "Execution [FAILED] (status = ", status, ")"
				print '\n{0:s}; {1:s};'.format(item, thr ),
				for k,v in results.iteritems() : print '{0:.2f}; {1:.2f};'.format(v[0], v[1]),
			print "\n--------------------------------------------"
	else:
		print "WARNING: No options provided. (use -h for help)"
