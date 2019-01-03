from scipy.stats import ttest_ind as ttest
import sys
from sys import argv

import numpy as np
import scipy as sp
import scipy.stats

def welch_ttest(x1, x2):
#custom t test code https://rosettacode.org/wiki/Welch%27s_t-test#Python
    n1 = x1.size
    n2 = x2.size
    m1 = np.mean(x1)
    m2 = np.mean(x2)
    v1 = np.var(x1, ddof=1)
    v2 = np.var(x2, ddof=1)
    t = (m1 - m2) / np.sqrt(v1 / n1 + v2 / n2)
    df = (v1 / n1 + v2 / n2)**2 / (v1**2 / (n1**2 * (n1 - 1)) + v2**2 / (n2**2 * (n2 - 1)))
    p = 2 * sp.stats.t.cdf(-abs(t), df)
    return t, df, p




#readin and split and transform to float

print >> sys.stderr, "start to readin points"

id = {1:"AA",2:"AB",3:"BA",4:"BB",5:"U"}
data = {}

fh = open(argv[1],"r")
cnt = 0
for line in fh:
    line = line.strip()
    box = line.split(",")
    box.pop()
    box_float = map (float, box)
    cnt += 1
    data[id[cnt]] = box_float


print >> sys.stderr, "done"

#do t.test

print >> sys.stderr, "do t.test"
if data.has_key("AA") and data.has_key("AB"):
    t_stat, p_val = ttest(data["AA"],data["AB"],equal_var = False)
    t, df, p = welch_ttest(np.array(data["AA"]), np.array(data["AB"])) #customed
    print "AA vs AB is : ",t_stat,p_val,t, df, p
else: 
    print >> sys.stderr, "error data has no AA or AB"

if data.has_key("BA") and data.has_key("BB"):
    t_stat, p_val = ttest(data["BA"],data["BB"],equal_var = False)
    t, df, p = welch_ttest(np.array(data["BA"]), np.array(data["BB"]))
    print "BA vs BB is : ",t_stat,p_val,t, df, p
else:
    print >> sys.stderr, "error data has no BA or BB"


if data.has_key("AB") and data.has_key("BA"):
    t_stat, p_val = ttest(data["AB"],data["BA"],equal_var = False)
    t, df, p = welch_ttest(np.array(data["AB"]), np.array(data["BA"]))
    print "AB vs BA is : ",t_stat,p_val,t, df, p
else:
    print >> sys.stderr, "error data has no AB or BA"


print >> sys.stderr, "done"



