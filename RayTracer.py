import numpy as np
import matplotlib.pyplot as plt
from matplotlib import patches
import json
import math
from scipy import optimize
from abc import ABC, abstractmethod
# TODO: add typing, error handling


class Path:
    """ Collection of Rays"""
    def __init__(self, x, y, theta, color='blue'):
        self.rays = [[x, y, theta]]
        self.curr_ray = 0
        self.color = color
        self.blocked = False
        self.image = False  # Set to true if the path makes it to the image plane

    def propagate_to(self, element):
        if not self.blocked:
            # Move to next element; calculate new ray (x,y)
            x, y = element.get_intersect(self.rays[self.curr_ray])
            # Update direction of new ray add the new ray to path
            self.rays.append([x, y, element.update_direction(x, y, self.rays[self.curr_ray][2])])
            element.update_blocking(x, y, self)
            self.curr_ray += 1

    def plot(self):
        # Convert to numpy array
        self.rays = np.array(self.rays)
        # Slice and plot
        ax_sim.plot(self.rays[:, 0], self.rays[:, 1], color=self.color, zorder=0)


class Source(ABC):
    @abstractmethod
    def __init__(self):
        pass

    @abstractmethod
    def get_paths(self):
        """returns all the paths created by a source"""


class PointSourceFromFile(Source):
    """A set of rays originating from a single point. Intensity from JSON"""
    def __init__(self, x, y, file_path):
        super().__init__()
        self.file_path = file_path
        self.paths = []

        with open(self.file_path) as f:
            angularIntensity = json.load(f)

        def get_intensity(theta, angularIntensity):
            min_angle = angularIntensity[0][0]
            max_angle = -min_angle
            increment = max_angle / ((len(angularIntensity) - 1) / 2)
            return angularIntensity[int(theta / increment + min_angle / increment)][1]

        """
        # Plot Intensity vs. theta
        intensity_numpy = np.array(angularIntensity)
        ax_intensity.plot(intensity_numpy[:, 0], intensity_numpy[:, 1])
        """

        # Define initial ray-thetas in density derived from angular intensity
        theta_bins = range(-15, 15, 1)
        initial_thetas = []
        for theta in theta_bins:
            rayDensity = get_intensity(theta + .5, angularIntensity)
            num_rays = int(rayDensity * 30)
            if num_rays > 0:
                spacing = 1.0 / num_rays
                for ray in range(num_rays):
                    initial_thetas.append(math.radians(theta + ray * spacing))

        # Create ray objects
        for theta in initial_thetas:
            self.paths.append(Path(x, y, theta))

    def get_paths(self):
        return self.paths


class LambertianPointSource(Source):
    """A set of rays originating from a single point with Lambertian intensity"""
    def __init__(self, x, y, normal_theta, min_theta, max_theta, resolution, color='blue'):
        super().__init__()
        self.paths = []

        # Define initial ray-thetas in density derived from angular intensity of lambertian emmiter
        theta_bins = np.arange(min_theta, max_theta, step=resolution)
        theta_ray_density = np.ndarray(len(theta_bins))
        initial_thetas = []
        for i, theta in enumerate(theta_bins):
            rayDensity = math.cos(theta + normal_theta)
            num_rays = int(rayDensity * 20)  # second number is intensity resolution
            if num_rays > 0:
                spacing = resolution / num_rays
                for ray in range(num_rays):
                    initial_thetas.append(theta + ray * spacing)
            theta_ray_density[i] = num_rays


        # Plot Intensity vs. theta
        ax_intensity.plot(theta_ray_density, theta_bins)


        # Create ray objects
        for theta in initial_thetas:
            self.paths.append(Path(x, y, theta, color))

    def get_paths(self):
        return self.paths


class CompoundSource(Source):
    def __init__(self, sources):
        self.paths = list()
        self.sources = list()
        self.sources.extend(sources)
        for source in self.sources:
            self.paths.extend(source.get_paths())

    def append(self, source):
        """Add one source element to Sources"""
        self.paths.append(source.get_paths())

    def get_paths(self):
        return self.paths


class Element(ABC):
    @abstractmethod
    def __init__(self):
        pass

    @abstractmethod
    def get_intersect(self, ray):
        """returns the point that the ray intersects the element"""

    @abstractmethod
    def update_blocking(self, x, y, path):
        """set 'blocked' status of path"""

    @abstractmethod
    def update_direction(self, x, y, initial_theta):
        """returns the updated theta when a ray crosses an element"""


class Aperture(Element):
    def __init__(self, x, min_y, max_y):
        super().__init__()
        self.center = x
        self.min_y = min_y
        self.max_y = max_y

        ax_sim.plot([x, x], [min_y, max_y], linewidth=5, color='red')  # plots aperture as line

    def get_intersect(self, ray):
        x = self.center
        y = (self.center-ray[0])*math.tan(ray[2]) + ray[1]
        return x, y

    def update_blocking(self, x, y, path):
        if y <= self.min_y or y >= self.max_y:
            path.blocked = True

    def update_direction(self, x, y, initial_theta):
        new_theta = initial_theta
        return new_theta


class LensTransition(Element):
    """Spherical or flat lens surface"""
    def __init__(self, x_intersect, diameter, radius_of_curvature, index_ratio):
        """x is defined as the point where the transition surface intersects the x axis
        radius_of_curvature > 0 for convex (center of curvature after the interface)
        index_ratio is n1/n2
        This element assumes all non-blocked rays pass through it with a non-zero normal angle"""
        super().__init__()
        self.x_intersect = x_intersect
        self.r = radius_of_curvature
        self.n = index_ratio
        self.diameter = diameter
        self.block_current_path = False

        self.is_flat = False
        if self.r == 0:
            self.is_flat = True
            # Plot a line to represent flat lens transition
            ax_sim.plot([x_intersect, x_intersect], [-diameter/2, diameter/2], linewidth=5, color='red', zorder=2)  # plots diameter as line
            ax_sim.plot([x_intersect, x_intersect], [-diameter/2, diameter/2], 'ro', zorder=2)  # plot opening points for clarity
        else:
            # Plot an arc to represent curved lens transition
            start_angle = math.degrees(math.asin(diameter/np.abs(radius_of_curvature)))/2
            arc = patches.Arc((x_intersect + radius_of_curvature, 0), radius_of_curvature*2, radius_of_curvature*2,
                              theta1=180-start_angle, theta2=180+start_angle, linewidth=5, fill=False, zorder=2, color='red')
            ax_sim.add_patch(arc)

    def get_intersect(self, ray):
        tan_theta = math.tan(ray[2])
        if np.abs(ray[1]) < self.diameter/2:  # Check that ray actually instersects lens

            if not self.is_flat:  # Flat transistions only require simple snell's law application
                if self.r > 0:  # Check positive/negative curvature
                    brentq_bracket = [self.x_intersect, self.x_intersect + self.r/4]
                else:
                    brentq_bracket = [self.x_intersect + self.r/4, self.x_intersect]

                x_equation = lambda x: (x - self.x_intersect - self.r)**2 + (
                        (x - ray[0])*tan_theta + ray[1])**2 - self.r**2
                sol = optimize.root_scalar(
                    x_equation, bracket=brentq_bracket, method='brentq')
                x = sol.root
                y = (self.x_intersect - ray[0]) * tan_theta + ray[1]
            else:
                x = self.x_intersect
                y = (self.x_intersect - ray[0]) * tan_theta + ray[1]
        else:
            x = ray[0]
            y = ray[1]

        return x, y

    def update_blocking(self, x, y, path):
        """Deal with total internal reflection"""
        if self.block_current_path:
            path.blocked = True
            self.block_current_path = False

    def update_direction(self, x, y, initial_theta):
        if np.abs(y) < self.diameter/2:  # Check that ray actually instersects lens

            if self.is_flat:
                angle_of_incidence = initial_theta
                sin_theta2 =self.n * math.sin(angle_of_incidence)  # Snell's law
                if sin_theta2 > 1 or sin_theta2 < - 1:  # total internal reflection
                    self.block_current_path = True
                    new_theta = 0  # Path is blocked so this is arbitrary
                else:
                    new_theta = math.asin(sin_theta2)
            else:
                normal_theta = math.asin(y/-self.r)
                angle_of_incidence = normal_theta - initial_theta

                sin_theta2 = self.n * math.sin(angle_of_incidence)  # Snell's law
                if sin_theta2 > 1 or sin_theta2 < - 1:  # total internal reflection
                    self.block_current_path = True
                    new_theta = 0  # arbitrary, path is blocked
                else:
                    angle_of_exit = math.asin(sin_theta2)
                    delta_theta = angle_of_exit - angle_of_incidence
                    new_theta = initial_theta - delta_theta

        else:
            new_theta = initial_theta

        return new_theta


class Image(Element):
    """Image plane, defined by 2 points in 2D"""
    def __init__(self, x1, y1, x2, y2):
        super().__init__()
        self.x1 = x1
        self.x2 = x2
        self.y1 = y1
        self.y2 = y2
        ax_sim.plot([x1, x2], [y1, y2], linewidth=5, color='red')

    def get_intersect(self, ray):
        tan_theta = math.tan(ray[2])
        if self.x1 == self.x2:
            x = self.x1
            y = (x-ray[0])*tan_theta +ray[1]
        else:
            m = (self.y2 - self.y1)/(self.x2-self.x1)
            x = (m*self.x1 - self.y1 -ray[0]*tan_theta
                 + ray[1])/(m - tan_theta)
            y = (x-ray[0])*tan_theta +ray[1]
        return x, y

    def update_blocking(self, x, y, path):
        path.blocked = True
        """Only for image plane, update if the path impinges the image plane"""
        # Vertical image plane
        if (self.x2 == self.x1):
            if self.y2 >= self.y1:
                if(y >= self.y1 and y <= self.y2):
                    path.image = True
            if self.y2 <= self.y1:
                if(y >= self.y2 and y <= self.y1):
                    path.image = True
        # Any other angle, use x to check intersection
        if self.x2 >= self.x1:
            if x >= self.x1 and x <= self.x2:
                path.image = True
        if self.x1 >= self.x2:
            if x >= self.x2 and x <= self.x1:
                path.image = True

    def update_direction(self, x, y, initial_theta):
        new_theta = initial_theta
        return new_theta

    def plot_image_intensity(self, paths):
        ray_intersects = []
        for path in paths:
            if path.image:
                x = path.rays[-1][0]
                y = path.rays[-1][1]
                if self.y2 > self.y1:
                    image_x = math.sqrt((x-self.x1)**2 + (y-self.y1)**2)
                else: image_x = math.sqrt((x-self.x2)**2 + (y-self.y2)**2)
                ray_intersects.append(image_x)
        """Plot histogram"""
        hist_points = []
        bins = 50
        image_length = image_x = math.sqrt((self.x2-self.x1)**2 + (self.y2-self.y1)**2)
        bin_length = image_length/bins
        curr_bin = 0
        curr_location = 0
        while curr_location < image_length:
            counts = 0
            for intersect in ray_intersects:
                if intersect >= curr_location and intersect < curr_location + bin_length:
                    counts += 1
            hist_points.append([curr_location, counts])
            # interate
            curr_bin+=1
            curr_location = curr_bin * bin_length
        hist_points = np.array(hist_points)
        average_bins = 5 # Should be odd so average is symmetric around center bin
        clip = int((average_bins - 1)/2)
        hist_points[clip:-clip,1] = np.convolve(
                            hist_points[:,1], np.ones((average_bins))/average_bins, mode='valid')
        hist_points[0:clip-1,1] = hist_points[clip-1,1]
        print(hist_points)
        hist_points[-clip,1] = hist_points[-clip -1,1]
        hist_points[:,1] = hist_points[:,1]/np.sum(a = hist_points, axis = 1) # Normalize
        ax_image.plot(hist_points[:, 0], hist_points[:, 1])


class System:
    """The collection of source, image plane, and all optical elements to simulate"""
    def __init__(self, sources, elements, image):
        self.sources = sources
        self.image = image
        self.elements = elements
        self.elements.append(image)  # Images are just elements anyway

    def run(self):
        paths = self.sources.get_paths()  # all paths start with the source rays
        # Calculate each ray propagating through the system
        for path in paths:
            for element in self.elements:
                path.propagate_to(element)
            path.plot()

        # Plot the intensity resulting in the image plane
        self.image.plot_image_intensity(paths)

        # Show the generated simulation
        plt.show()


# Set up plotting
fig_intensity, ax_intensity = plt.subplots()
ax_intensity.set_title('Source Intensity vs. Theta')
fig_sim, ax_sim = plt.subplots()
ax_sim.set_title('Optical System')
fig_image, ax_image = plt.subplots()
ax_image.set_title('Image Intensity vs. Y location')


def main():
    """2D ray-tracing in simple optical system"""

    """**Start USER CODE**"""
    # Create system built of optical elements
    point_source1 = LambertianPointSource(0, 0, 0, math.radians(-89), math.radians(89), math.radians(1), color='blue')
    point_source2 = LambertianPointSource(0, 5, 0, math.radians(-89), math.radians(89), math.radians(1), color='green')

    #point_source3 = PointSourceFromFile(0, 10, 'AngularIntensity.json')
    compound_source = CompoundSource([point_source1, point_source2])

    lens_thickness = 2.5
    lens_diameter = 25
    lens_start_distance = 29.1
    lens_start = LensTransition(lens_start_distance, lens_diameter, 0, 1/4.01)
    lens_end = LensTransition(lens_start_distance + lens_thickness, lens_diameter, -45.1, 4.01/1)

    aperture_start = Aperture(33, -1.4, 1.4)
    aperture_end = Aperture(34, -1.4, 1.4)

    """Image extents"""
    bed_center_distance = 343.9
    bed_len = 165
    bed_theta = math.radians(90 - 19.83)
    bed_x1, bed_y1 = bed_center_distance - math.cos(bed_theta)*bed_len/2, -bed_len/2*math.sin(bed_theta)
    bed_x2, bed_y2 = bed_center_distance + math.cos(bed_theta)*bed_len/2, bed_len/2*math.sin(bed_theta)

    bed = Image(bed_x1, bed_y1, bed_x2, bed_y2)

    system = System(compound_source, [lens_start, lens_end, aperture_start, aperture_end], bed)
    system.run()
    """**END USER CODE**"""


if __name__ == "__main__":
    main()
