#!/usr/bin/env python

import sys




# input comes from STDIN (standard input)
for line in sys.stdin:
	# remove leading and trailing whitespace
	line = line.strip()
	# split the line into words
	words = line.split(',')

	year = words[0]
	origin = words[16]
	cancelled = words[21]

	print  '%s,%s\t%s'% (year, origin, cancelled)