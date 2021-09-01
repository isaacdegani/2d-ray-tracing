import numpy as np
import matplotlib.pyplot as plt
from matplotlib import patches
import json
import math
from scipy import optimize
from abc import ABC, abstractmethod
from collections.abc import Iterable
from inspect import isclass
# TODO: add typing, error handling


class Path:
    """ Collection of Rays"""
    def __init__(self, x, y, theta, visible=True, color='blue'):
        self.rays = [[x, y, theta]]
        self.curr_ray = 0
        self.color = color
        self.visible = visible
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
        if self.visible:
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
            rayDensity = math.cos(theta - normal_theta)
            num_rays = int(rayDensity * 500)  # second number is intensity resolution
            if num_rays > 0:
                spacing = resolution / num_rays
                for ray in range(num_rays):
                    initial_thetas.append(theta + ray * spacing)
            theta_ray_density[i] = num_rays

        # Plot Intensity vs. theta
        ax_intensity.plot(theta_ray_density, theta_bins)

        # Create ray objects
        for i, theta in enumerate(initial_thetas):
            # Show fewer rays for visibility
            if i % 1000 == 0:
                self.paths.append(Path(x, y, theta, color=color))
            else:
                self.paths.append(Path(x, y, theta, visible=False, color=color))

    def get_paths(self):
        return self.paths


class CompoundSource(Source):
    def __init__(self, sources):
        super().__init__()
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

    def get_sources(self):
        return self.sources


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

    @abstractmethod
    def plot_element(self):
        """draw the optical element"""


class CompoundElement():
    def __init__(self, elements=None):  # Initialize with optional list of elements
        super().__init__()
        self.elements = list()
        if elements is not None:
            self.elements.extend(elements)

    def append(self, element):
        """Add one source element to Sources"""
        self.elements.append(element)

    def get_elements(self):
        return self.elements

    def plot_element(self):
        for element in self.elements:
            element.plot_element()


class Aperture(Element):
    def __init__(self, x, min_y, max_y):
        super().__init__()
        self.center = x
        self.min_y = min_y
        self.max_y = max_y

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

    def plot_element(self):
        ax_sim.plot([self.center, self.center], [self.min_y, self.min_y*4], linestyle='dashdot', linewidth=5, color='black')
        ax_sim.plot([self.center, self.center], [self.max_y, self.max_y*4], linestyle='dashdot', linewidth=5, color='black')


class Wall(Element):
    """Opaque line, defined by 2 points in 2D"""
    def __init__(self, x1, y1, x2, y2):
        super().__init__()
        self.x1 = x1
        self.x2 = x2
        self.y1 = y1
        self.y2 = y2
        ax_sim.plot([x1, x2], [y1, y2], linewidth=5, color='black')

    def intersect_in_range(self, x, y):
        in_range = False
        # Vertical wall
        if (self.x2 == self.x1):
            if self.y2 >= self.y1:
                if(y >= self.y1 and y <= self.y2):
                    in_range = True
            if self.y2 <= self.y1:
                if(y >= self.y2 and y <= self.y1):
                    in_range = True
        # Any other angle, use x to check intersection
        if self.x2 >= self.x1:
            if x >= self.x1 and x <= self.x2:
                in_range = True
        if self.x1 >= self.x2:
            if x >= self.x2 and x <= self.x1:
                in_range = True
        return in_range

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
        if self.intersect_in_range(x, y):
            return x, y
        else:
            return ray[0], ray[1]

    def update_blocking(self, x, y, path):
        """Any paths that intersect are blocked"""
        path.blocked = self.intersect_in_range(x, y)

    def update_direction(self, x, y, initial_theta):
        new_theta = initial_theta
        return new_theta

    def plot_element(self):
        ax_sim.plot([self.x1, self.x2], [self.y1, self.y2], linewidth=5, color='black')


class LensTransition(Element):
    """Spherical or flat lens surface"""
    def __init__(self, x_intersect, diameter, radius_of_curvature, index_ratio):
        """x is defined as the point where the transition surface intersects the x axis
        radius_of_curvature > 0 for convex (center of curvature after the interface)
        index_ratio is n1/n2"""
        super().__init__()
        self.x_intersect = x_intersect
        self.r = radius_of_curvature
        self.n = index_ratio
        self.diameter = diameter
        self.block_current_path = False
        self.is_flat = False

        if self.r == 0:
            self.is_flat = True

        if not self.is_flat:
            self.delta_x = self.r - math.sqrt(self.r**2 - (self.diameter/2)**2)  # The x location of the top/bottom end of the lens transition
        else:
            self.delta_x = 0

    def get_delta_x(self):
        return self.delta_x

    def get_intersect(self, ray):
        # Calculate intersect
        tan_theta = math.tan(ray[2])
        y_at_delta_x = (self.x_intersect - ray[0] + self.delta_x) * tan_theta + ray[1]
        if not self.is_flat and np.abs(y_at_delta_x) < self.diameter/2:  # Check if a curved interface is encountered
                try:
                    if self.r > 0:  # Check positive/negative curvature
                        brentq_bracket = [self.x_intersect - 0.01, self.x_intersect + self.r + 0.01]
                    else:
                        brentq_bracket = [self.x_intersect + self.r - 0.01, self.x_intersect + 0.01]

                    x_equation = lambda x: (x - self.x_intersect - self.r)**2 + (
                            (x - ray[0])*tan_theta + ray[1])**2 - self.r**2
                    sol = optimize.root_scalar(
                        x_equation, bracket=brentq_bracket, method='brentq')
                    x = sol.root
                    y = (self.x_intersect - ray[0]) * tan_theta + ray[1]
                except ValueError:
                    print("Could not propagate ray:")
                    print(ray)
                    print("Y intercept:")
                    print(y_at_delta_x)
                    print("bracket")
                    print(brentq_bracket)
                    x = ray[0]
                    y = ray[1]
                    self.block_current_path = True
        else:
            x = self.x_intersect
            y = (self.x_intersect - ray[0] + self.delta_x) * tan_theta + ray[1]

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

    def plot_element(self):
        if self.is_flat:
            # Plot a line to represent flat lens transition
            ax_sim.plot([self.x_intersect, self.x_intersect], [-self.diameter/2, self.diameter/2], linewidth=5, color='red', zorder=2)  # plots optical surface as line
        else:
            # Plot an arc to represent curved lens transition
            start_angle = math.degrees(math.asin((self.diameter/2)/np.abs(self.r)))
            arc = patches.Arc((self.x_intersect + self.r, 0), self.r*2, self.r*2,
                              theta1=180-start_angle, theta2=180+start_angle, linewidth=5, fill=False, zorder=2, color='red')
            ax_sim.add_patch(arc)


class SphericalLens(CompoundElement):
    """Dual Spherical or Plano Spherical Lens Compound Element"""
    def __init__(self, x_intersect, diameter, r_of_curvature1, thickness, r_of_curvature2, index):
        super().__init__()
        lens_start = LensTransition(x_intersect, diameter, r_of_curvature1, 1/index)
        lens_upper_wall = Wall(x_intersect + lens_start.get_delta_x() + 0.01, diameter/2,
                            x_intersect + thickness - 0.01, diameter/2)
        lens_lower_wall = Wall(x_intersect + lens_start.get_delta_x() + 0.01, -diameter/2,
                            x_intersect + thickness - 0.01, -diameter/2)
        lens_end = LensTransition(x_intersect + thickness, diameter, r_of_curvature2, index/1)
        self.elements = [lens_start, lens_upper_wall, lens_lower_wall, lens_end]


class Image(Element):
    """Image plane, defined by 2 points in 2D"""
    def __init__(self, x1, y1, x2, y2):
        super().__init__()
        self.x1 = x1
        self.x2 = x2
        self.y1 = y1
        self.y2 = y2

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
        hist_points[-clip,1] = hist_points[-clip -1,1]
        hist_points[:,1] = hist_points[:,1]/np.sum(a = hist_points, axis = 1) # Normalize
        ax_image.plot(hist_points[:, 0], hist_points[:, 1])

    def plot_element(self):
        ax_sim.plot([self.x1, self.x2], [self.y1, self.y2], linewidth=5, color='blue')


class System:
    """The collection of source, image plane, and all optical elements to simulate"""
    def __init__(self, sources, elements, image):
        self.sources = [sources]
        self.image = image
        self.elements = [elements]
        self.elements.append(image)  # Images are just elements anyway

    def run(self):

        # Helper function to decompose nested lists and compound elements
        def flatten(el_list):
            for el in el_list:
                if isinstance(el, Iterable):
                    yield from flatten(el)
                elif isinstance(el, CompoundElement) or (isclass(el) and issubclass(el, CompoundElement)):
                    yield from flatten(el.get_elements())
                elif isinstance(el, CompoundSource) or (isclass(el) and issubclass(el, CompoundSource)):
                    yield from flatten(el.get_sources())
                else:
                    yield el

        # All paths start with the source rays
        print("Generating source rays...")
        self.simple_sources = list(flatten(self.sources))  # Unpack sources
        self.paths = [path
                      for source in self.simple_sources
                      for path in source.get_paths()]  # Add up paths

        # Calculate each ray propagating through the system
        print("Propagating rays through system...")
        self.elements = list(flatten(self.elements))  # Unpack elements
        for element in self.elements:  # Draw elements in system
            element.plot_element()

        for path in self.paths:
            for element in self.elements:
                path.propagate_to(element)
            path.plot()

        # Plot the intensity resulting in the image plane
        self.image.plot_image_intensity(self.paths)

        # Show the generated simulation
        x_scale = np.abs(ax_sim.get_xlim()[1] - ax_sim.get_xlim()[0])
        ax_sim.set_ylim(-x_scale/2, x_scale/2)  # Equal aspect ratio, limit y range
        print("Generating plots...")
        plt.show()


# Set up plotting
#TODO: Move this to run() function and get rid of global variables
fig_sim, ax_sim = plt.subplots()  # This is the main system plot
ax_sim.set_title('Optical System')

fig_intensity, ax_intensity = plt.subplots()  # Source characterization
ax_intensity.set_title('Source Intensity vs. Theta')

fig_image, ax_image = plt.subplots()  # Image
ax_image.set_title('Image Intensity vs. Y location')
