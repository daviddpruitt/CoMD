

from sys import stdin
within = 0
for line in stdin.readlines():
    line = line[:-1]            # strip lf
    if "Timings for Rank" in line:
        within = 1
    if line == "Timing Statistics Across 1 Ranks:":
        within = 1
    if line == "Network Results:":
        within = 0
    if within:
        print line

