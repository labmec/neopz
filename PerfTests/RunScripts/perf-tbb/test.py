#! /usr/bin/env python2.7
################################################################################
#   Copyright (C) 2014 by:                                                     #
#   Gilvan Vieira (gilvandsv@gmail.com)                                        #
#                                                                              #
#   This program is free software; you can redistribute it and/or modify       #
#   it under the terms of the GNU General Public License as published by       #
#   the Free Software Foundation; either version 2 of the License, or          #
#   (at your option) any later version.                                        #
#                                                                              #
#   This program is distributed in the hope that it will be useful,            #
#   but WITHOUT ANY WARRANTY; without even the implied warranty of             #
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the              #
#   GNU General Public License for more details.                               #
#                                                                              #
#   You should have received a copy of the GNU General Public License          #
#   along with this program; if not, write to the                              #
#   Free Software Foundation, Inc.,                                            #
#   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.                  #
################################################################################

# libraries
import sys
import os.path
import shlex, subprocess
import resource
import getopt
import math # sqrt

# program information
executable='./Perf-TBB'
perfdata='dec.rdt'
inputdata='cube1.txt'
rundir='.'

# limits for this test
limits = { "cpu"   : (resource.RLIMIT_CPU, 3600, "Max CPU time in seconds") }

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
        # process the remaining rows
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
    
def average(vlist) :
    sz = len(vlist)
    if sz == 0 :
        return 0.0
    else :
        s=reduce(lambda x,y : float(x)+float(y), vlist, 0.0)
    return s/sz

def variance(vlist) :
    size=len(vlist)
    if size == 1 or size == 0:
        raise StatsError("Cannot calculate variance for a sample of size "+str(size))
    sum1=reduce(lambda x,y : float(x)+float(y), vlist, 0.0)
    sum2=reduce(lambda x,y : float(x)+(float(y)*float(y)), vlist, 0.0)
    v = (sum2-(pow(sum1,2)/float(size)))/float(size-1)
    if v < 0.0 :
        if v > -1.e-15 :  # tolerate floating point rounding errors
            return 0.0
        else :
            raise StatsError("Something went wrong while calculating the variance (result < 0.0)")
    return v

def stdev(vlist) :
    var = variance(vlist)
    return math.sqrt(var)
    
# set the rlimits of the chidren process (see limits above)
def setlimits():
    print "Setting resource limit in child"
    for k, v in limits.iteritems() :
        resource.setrlimit(v[0], (v[1],v[1]))
def show_results(res_data, ths):
    # read and show stats
    stat_res = read(res_data)
    elapsed_list = get_column_values(stat_res,  "ELAPSED")
    self_list = get_column_values(stat_res,     "SELF_RU_UTIME")
    lst = 0
    print '{0:7s} \t {1:7s} \t {2:7s} \t {3:7s} \t {4:7s}'.format('threads', 'elapsed', 'user', 'elapsed_error', 'user_error')
    for t in ths:
        print '{0:7s} \t {1:5.2f} \t {2:5.2f} \t '.format(t, average(elapsed_list[lst:lst+ntimes]), average(self_list[lst:lst+ntimes])),
        if ntimes > 1:
            print '{0:5.2f} \t {1:5.2f} '.format(stdev(elapsed_list[lst:lst+ntimes]), stdev(self_list[lst:lst+ntimes]))
        else:
            print '{0:5.2f} \t {1:5.2f} '.format(0, 0)
        lst = lst + ntimes
def run_exp(nthreads, ntimes, res_data, t):
    # setup the command line
    arguments = executable + ' -p 3'
    arguments = arguments  + ' -perf ' + res_data
    arguments = arguments  + ' -mc '   + inputdata
    arguments = arguments  + ' -nt '   + nthreads
    arguments = arguments  + ' -type ' + t
    # execute
    print arguments
    args = shlex.split(arguments)
    sout = None
    serr = None
    for i in range(ntimes) :
        p = subprocess.Popen(args, preexec_fn=setlimits, stdout=sout, stderr=serr, cwd=rundir)
        p.wait()
        if (p.returncode != 0) :
            print "Execution [FAILED]"
            
# main - principal execution
if __name__ == "__main__":
    # init arguments
    ntimes = 1
    threads = [ '2', '4', '6', '8', '12' ]
    # process arguments
    try :
        opts, extra_args = getopt.getopt(sys.argv[1:], 't:n:')
    except getopt.GetoptError, e:
        error(str(e), 1)
    for f, v in opts:
        if f == '-n':
            # number of times to repeat the experiments
            ntimes=int(v)
    # running serial
    res = 'res_data_1.rdt'
    run_exp('1', ntimes, res, '1')

    # running experiments parallel
    for t in range(2,6):
        res = 'res_data_' + str(t) + '.rdt'
        for i in threads:
             run_exp(i, ntimes, res, str(t))
    print '# Execution Type : 1 #'
    show_results(res, ['0'])
    for t in range(2,6):
        res = 'res_data_'+str(t)+'.rdt'
        print '# Execution Type : ' + str(t) + ' #'
        show_results(res, threads)
    