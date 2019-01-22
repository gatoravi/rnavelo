import sys
import loompy

loompy.combine(sys.argv[1:-1], sys.argv[-1])
# usage: script.py file1.loom file2.loom file3.loom merged.loom
