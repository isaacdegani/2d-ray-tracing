# RayTracing
Code for ray tracing through simple optical systems
Requires Scipy.

AngularIntensity.json decribes the source beam.

An optical system is built up of classes that inherent from the abstract "Element" Class.
The system is simulated as lines depicting rays in a matplotlib plot.
The beam originates as a point source from the origin and ends at an 
imaging element. Intensity graphs of the source beam and image are generated.

User definable code is found at the top of the main function in RayTracer.py
