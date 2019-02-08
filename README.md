Ray tracer assignment for intro to graphics course taught at the university of waterloo. Definitely a mess, just archiving. Cool pictures can be found in `Assets/` though. Maybe one day I'll clean it up.

Previous readme
```
Compilation ---------------------
NOTE: always run make clean before switching configurations, ie if you just compiled with config=debug, then run make clean before recompiling with config=nosupersample

Run `make` to compile with supersampling (and bounding volumes).

The following other configurations exist:
  make config=drawboundingvol
    No supersampling. Draws the bounding vol of mesh objects
  make config=nosupersample
    No supersampling.
  make config=nobb
    No supersampling. No bounding volume. Used to show speedup gained from the bounding vol

To run, cd to the Assets dir and run ../A4 <name of lua file>
for example, to test the sample.lua run "../A4 sample.lua" from the Assets dir.

Manual --------------------------

File extensions:
-bb.png --- drawn bounding volumes
-noss.png --- no supersampling

Notes:
- the program will note progress
- sample.lua contains the custom scene. It is supposed to be a helix w/ some other primitives floating around
- the background is a closeup of the mandelbrot set
- THE EXTRA FEATURE IS SUPERSAMPLING: done with 16 samples per pixel
- ignored bounding volumes if mesh name is "plane.obj" (so the plane is still visible in this mode)


proof of speedup from bounding volumes: (running /usr/bin/time -f '%Uu %Ss %Er')

running: ../A4 macho-cows.lua

config=nobb (no bounding volumes):
118.99u 0.02s 1:59.70r

config=nosupersample (bounding volumes):
9.68u 0.00s 0:09.69r

a 12x speedup!


```
