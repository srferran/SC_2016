#!/usr/bin/env python

from operator import itemgetter
import sys

current_year = None
current_origin = None
current_count = 0
year = None
origin = None
total = 0
# input comes from STDIN
for line in sys.stdin:
	# remove leading and trailing whitespace
	line = line.strip()

	# parse the input we got from mapper.py
	key, cancelled = line.split('\t')
	year, origin = key.split(',')

	# convert count (currently a string) to int
	try:
		cancelled = int(cancelled)
	except ValueError:
		# count was not a number, so silently
		# ignore/discard this line
		continue

	# this IF-switch only works because Hadoop sorts map output
	# by key (here: word) before it is passed to the reducer
	if current_year == year and current_origin == origin:
		current_count += cancelled
	else:
		if current_year and current_origin:
			# write result to STDOUT
			print '%s\t%s\t%s'% (current_year, current_origin, current_count)

		current_count = cancelled
		current_year = year
		current_origin = origin

# do not forget to output the last word if needed!
if current_year == year and current_origin == origin:
	print '%s\t%s\t%s'% (current_year, current_origin, current_count)

	