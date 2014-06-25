![alt tag](https://raw.github.com/sqaxomonophonen/ptks/master/ptks-screenshot.png)

Attempt at real-time path tracing on the CPU. Yet another one of those problems
that could use a 1000x faster CPU. The path tracer is simple; no bidirectional
path tracing nor metropolis light transport. Initially inspired by smallpt. It
traces against a Build engine-like 2.5D structure (so guess where "ks" came
from).

The above 192x108 screenshot renders at 7 frames per second using all 4 cores on my `Intel(R) Core(TM) i5-3230M CPU @ 2.60GHz` laptop CPU (according to `/proc/cpuinfo`). That's about 7 million ray traces per second, each with up to 4 bounces.
