import raytracer as rt
import math


def main():
    """
    2D ray-tracing in simple optical system
    Create system built of optical elements
    Add graphics manually to 'rt.ax_sim' matplotlib axes element if needed
    """
    print("Constructing optical system...")

    """Place two lambertian point sources at the edges of an object"""
    image_distance = 200  # Distance from center of bed to sensor

    object_length = 100

    beam_width = math.radians(80)  # Arbitrary; just limits the number of rays to propagate
    # First source
    obj_x1, obj_y1 = 0, object_length/2
    point_source1 = rt.LambertianPointSource(obj_x1, obj_y1, 0, -beam_width/2, beam_width/2, math.radians(0.1), color='blue')
    # Second source
    obj_x2, obj_y2 = 0, -object_length/2
    point_source2 = rt.LambertianPointSource(obj_x2, obj_y2, 0, -beam_width/2, beam_width/2, math.radians(0.1), color='blue')

    # Draw object
    rt.ax_sim.plot([obj_x1, obj_x2], [obj_y1, obj_y2], linewidth=5, color='gray')

    # Combine the sources into one element
    object_sources = [point_source1, point_source2]

    """Model of optics"""
    # Aperture
    aperture_radius = 3
    aperture_start = image_distance - 30.0
    aperture = rt.Aperture(aperture_start, -aperture_radius, aperture_radius)  # Start of sensor column

    # Define lens geometry
    lens_thickness = 2.5
    lens_diameter = 12
    lens_curvature = 20
    n_bk7 = 1.51
    lens_start_distance = image_distance - 25.0

    lens = rt.SphericalLens(lens_start_distance, lens_diameter,lens_curvature,lens_thickness, 0, n_bk7)

    image = rt.Image(image_distance, -10.0, image_distance, 10.0)  # Sensor die perpendicular to lens/apertures

    """Simulate!"""
    # System elements need to be in order (source --> image)
    system = rt.System(object_sources, [aperture, lens], image)
    system.run()


if __name__ == "__main__":
    main()
