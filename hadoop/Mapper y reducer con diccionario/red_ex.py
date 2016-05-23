#!/usr/bin/env python

from operator import itemgetter
import sys

current_airport = None
current_counts = 0
current_delay = 0
# input comes from STDIN
for line in sys.stdin:
	# remove leading and trailing whitespace
	#line = line.strip()
	# parse the input we got from mapper.py
	key, delay, counts = line.split('\t')
	year, airport = key.split(',')
	
	try:
		delay = float(delay)
		counts = int(counts)
	except ValueError:
		continue
	
	if current_airport == airport:
		current_delay += delay
		current_counts += counts

	else:
		if current_airport:
			# write result to STDOUT
			mean = current_delay/current_counts
			print '%s\t%s\t%f' % (year,current_airport, mean)  

		current_delay = delay
		current_counts = counts
		current_airport = airport

if current_airport == airport:
	mean = current_delay/current_counts
	print '%s\t%s\t%f' % (year,current_airport, mean)