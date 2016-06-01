#!/usr/bin/env python

from operator import itemgetter
import sys

current_airport = None
current_counts = 0
current_delay = 0
current_year = 0
airport = ""
# input comes from STDIN
for line in sys.stdin:
	# remove leading and trailing whitespace
	#line = line.strip()
	# parse the input we got from mapper.py
	# split into key =(year, airport) delay = delay, counts = processed items for this (year,airport) pair

	key, delay, counts = line.split('\t')
	# split key into year, airport
	year, airport = key.split(',')
	
	try:
		delay = float(delay)
		counts = int(counts)
	except ValueError:
		continue

	# while input has same airport add values
	if current_airport == airport:
		current_delay += delay
		current_counts += counts
		
		# if new airport it's loaded just output current stored data
	else:
		# avoid initial case or no items processed. ( current_airport = None)
		if current_airport:
			
			# write result to STDOUT
			mean = current_delay/current_counts
			print '%s\t%s\t%f' % (current_year,current_airport, mean)  

		current_delay = delay
		current_counts = counts
		current_airport = airport
		current_year = year
		
# last iteration
if current_airport == airport:
	mean = current_delay/current_counts
	print '%s\t%s\t%f' % (current_year,current_airport, mean)