import numpy as np
import matplotlib.pyplot as plt
import json
import math
from scipy import optimize
from abc import ABC, abstractmethod
# TODO: add typing, error handling


class Path:
    """ Collection of Rays"""
    def __init__(self, x, y, theta):
        self.rays = [[x, y, theta]]
        self.curr_ray = 0
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
        ax_sim.plot(self.rays[:, 0], self.rays[:, 1])

        # show elements
        ax_sim.plot(self.rays[:, 0], self.rays[:, 1], 'ro')


class Source(ABC):
    @abstractmethod
    def __init__(self):
        pass

    @abstractmethod
    def get_paths(self):
        """returns all the paths created by a source"""


class PointSource(Source):
    """A set of rays originating from a single point"""
    def __init__(self, x, y, angular_i='lambertian', num_rays=15, path_to_angular_i=None):
        super().__init__()
        self.angular_i = angular_i
        self.path_to_angular_i = path_to_angular_i
        self.num_rays = num_rays
        self.paths = []

        # TODO: make this a case switch with new type
        if angular_i == 'lambertian':
            angularIntensity = [[]]
        elif angular_i == 'from_path':  # Load intensity vs theta data
            with open(self.path_to_angular_i) as f:
                angularIntensity = json.load(f)
        elif angular_i == 'iso':
            angularIntensity = [[]]
        else:
            print('What gives?')

        def get_intensity(theta, angularIntensity):
            min_angle = angularIntensity[0][0]
            max_angle = -min_angle
            increment = max_angle / ((len(angularIntensity) - 1) / 2)
            return angularIntensity[int(theta / increment + min_angle / increment)][1]

        # Plot Intensity vs. theta
        intensity_numpy = np.array(angularIntensity)
        ax_intensity.plot(intensity_numpy[:, 0], intensity_numpy[:, 1])

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
        ax_sim.plot([x, x], [min_y, max_y], 'bo')  # plot endpoints

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
    def __init__(self, x_intersect, radius_of_curvature, index_ratio):
        """x is defined as the point where the transition surface intersects the x axis
        radius_of_curvature > 0 for convex (center of curvature after the interface)
        index_ratio is n1/n2
        This element assumes all non-blocked rays pass through it with a non-zero normal angle"""
        super().__init__()
        self.x_intersect = x_intersect
        self.r = radius_of_curvature
        self.n = index_ratio
        self.block_current_path = False

        self.is_flat = False
        if self.r == 0:
            self.is_flat = True

    def get_intersect(self, ray):
        tan_theta = math.tan(ray[2])
        if not self.is_flat:
            if self.r > 0:
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

        return x, y

    def update_blocking(self, x, y, path):
        """Deal with total internal reflection"""
        if self.block_current_path:
            path.blocked = True
            self.block_current_path = False

    def update_direction(self, x, y, initial_theta):
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

        return new_theta


class Image(Element):
    """Image plane, defined by 2 points in 2D"""
    def __init__(self, x1, y1, x2, y2):
        super().__init__()
        self.x1 = x1
        self.x2 = x2
        self.y1 = y1
        self.y2 = y2
        ax_sim.plot([x1, x2], [y1, y2])
        ax_sim.plot([x1, x2], [y1, y2], 'bo')  # plot endpoints

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
        average_bins = 5 # Should be odd
        clip = int((average_bins - 1)/2)
        hist_points[clip:-clip,1] = np.convolve(
                            hist_points[:,1], np.ones((average_bins))/average_bins, mode='valid')
        hist_points[0:clip-1,1] = hist_points[clip,1]
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
fig_sim, ax_sim = plt.subplots()
fig_image, ax_image = plt.subplots()


def main():
    """2D ray-tracing in simple optical system"""

    """**Start USER CODE**"""
    # Create system built of optical elements
    point_source1 = PointSource(0, 0, angular_i='from_path', path_to_angular_i='AngularIntensity.json')
    point_source2 = PointSource(0, 10, angular_i='from_path', path_to_angular_i='AngularIntensity.json')
    compound_source = CompoundSource([point_source1, point_source2])

    lens_thickness = 2.5
    lens_start_distance = 29.1
    lens_start = LensTransition(lens_start_distance, 0, 1/4.01)
    lens_end = LensTransition(lens_start_distance + lens_thickness, -45.1, 4.01/1)

    aperture_start = Aperture(33, -1.4, 1.4)
    aperture_end = Aperture(34, -1.4, 1.4)

    """Bed extents"""
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
