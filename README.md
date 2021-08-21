# RayTracing
Code for ray tracing through simple (2D) optical systems.
Requires Matplotlib and Scipy.

AngularIntensity.json decribes the angular intensity of a source beam.

An optical system is built up of classes that inherent from the abstract "Element" Class.
The system is simulated as lines depicting rays in a matplotlib plot.
The beam originates as a point source from the origin and ends at an 
imaging element. Intensity graphs of the source beams and image are generated.

User definable code is all in the main function of RayTracer.py
