#! /usr/bin/env python2.7

# import libraries
import sys
import os.path
import resource
import getopt
import math #sqrt

# read a RDT file into a dictionary
# the dictionary is indexed by the column header
# and the list of values at the column
def read(filename):
    # check if the file exists
    if not os.path.isfile(filename):
        raise rdt_error(filename + 'is not a valid file.')
    d={}
    with open(filename) as f:
        # process the first row to setup headers
        hlist = f.readline().split(',')
        hi=1
        for h in hlist:
            header = h.strip()
            d[hi] = (header,[])
            hi = hi + 1
        # process the remaining rows
        for line in f:
            vlist = line.split(',')
            first = vlist[1]
            if first.strip() != "#":
                hi = 1
                for v in vlist:
                    d[hi][1].append(float(v.strip()))
                    hi = hi + 1
        return d
# find the column by field and return
# the list of values
def get_column_values(rdt_d, field):
    for k, v in rdt_d.iteritems():
        if v[0] == field:
            return v[1]
    return []
# calculates the average of the values
# in the vlist
def average(vlist):
    sz = len(vlist)
    if sz == 0:
        return 0.0
    else:
        s=reduce(lambda x, y : float(x) + float(y), vlist, 0.0)
    return s/sz
# calculates the variance of the values
# in the vlist
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
# calculate the standard deviation
# using the function variance of the values
# int the vlist
def stdev(vlist) :
    var = variance(vlist)
    return math.sqrt(var)
# main
if __name__ == "__main__":
    # process arguments
    try:
        opts, extra_args = getopt.getopt(sys.argv[1:], 'f:')
    except:
        error(str(e), 1)
    for f, v in opts:
        if f == '-f':
            filename = v
    stats_res = read(filename)
    elapsed_lst = get_column_values(stats_res, "ELAPSED")
    self_lst = get_column_values(stats_res, "SELF_RU_UTIME")
    print '{0:7s} \t {1:7s} \t {2:7s} \t {3:7s}'.format('User Time', 'User Error', 'Wall Time', 'Wall Error')
    elapsed_avg = average(elapsed_lst)
    self_avg = average(self_lst)
    elapsed_err = stdev(elapsed_lst)
    self_err = stdev(self_lst)
    print '{0:5.2f} \t {1:5.2f} \t {2:5.2f} \t {3:5.2f}'.format(self_avg, self_err, elapsed_avg, elapsed_err)
