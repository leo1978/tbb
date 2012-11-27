#!/usr/bin/python
import subprocess

#---------------------------------------------------------
# Test everything related to coarse grained multithreadin
#---------------------------------------------------------
def ResultLine(blob):
  lines = blob.strip().split('\n')
  return lines[-1]

paralleltests = [
  ("./pcontour --numThreads 4", "2856 points, 5504 cells 5504 cell data values"),
  ("./pcontour --numThreads 4", "2856 points, 5504 cells 5504 cell data values"),
  ("./pcontour --filter cutter --numThreads 4","331 points, 600 cells 600 cell data values"),
  ("./contouramr --numThreads 4", "2876 polygons in 8 blocks"),
  ("./threaded_streamtrace --numThreads 4","257437 points, 200 cells")
  ]

for (cmd,expected) in paralleltests:
  print cmd
  out = subprocess.check_output(cmd.split())
  assert(ResultLine(out)==expected)

#---------------------------------------------------------
# Test everything related to speeding up data insertion
#---------------------------------------------------------

arraytests = ["./arraytest --op insert",
              "./arraytest --op copy --helper --check",
              "./arraytest --op copy --lazy --check"]

for cmd in arraytests:
  print cmd
  assert(subprocess.call(cmd.split())==0)

print "All tests passed"
