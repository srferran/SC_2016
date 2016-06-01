#!/usr/bin/env python
import sys

d = {}
# input comes from STDIN (standard input)
for line in sys.stdin:
	# remove leading and trailing whitespace
	# split the line into words
	words = line.strip().split(',')
	
	# read year, airport, delay and cancelled status
	year = words[0]
	airport = words[16]
	delay = words[15]
	cancelled = words[21] == "1"
	
	
	# Only for not cancelled flights
	# cast delay from string to float
	if not cancelled:
		try:
			delay = float(delay)
		except ValueError:
			continue
		# check if key exists at dict
		# If not exist make new key & init to [current delay value, number of processed items = 1 ]
		if (year, airport) not in d:
			d[(year, airport)] = [delay, 1]
		#if key exists add delay and increase processed items by 1
		else:
			d[(year, airport)][0] += delay
			d[(year, airport)][1] += 1

# Output from mapper  
# flush dict
for y, a in d:
	print "%s,%s\t%f\t%d" % (y, a, d[(y, a)][0], d[(y, a)][1])