# Angular-Radial-Position-Distribution

Description:

Groups a spots population of germ cells into binned populations based on their positions relative to the ovary's dorsal and ventral sides. The statistics output includes a user defined binned histogram of the radial and angular distribution of the spots population.

Description of Algorithm:

Transforming coordinate systems to quantify angular and radialdistributions of cells within the mouse ovary.   The 2D cartesian coordinates of the Imaris spots objects, representing the cells, were recorded. Additionally, in cartesian coordinates, the center of geometry and the end points (top and bottom) of the Imaris surface representing the ovarian outer surface were recorded. Using the the center of geometry and two end points, the radius of curvature was calculated, and all spots object coordinates were offset so as to shift the circle described by the radius of curvature to the cartesian origin.&nbsp;These coordinates were then transformed from cartesian to polar coordinates.     Left Ovary:  If analyzing the left ovary, defined as having the dorsal on the left side, all polar coordinates were rotated by −180◦ so that the center of the ovary is placed at an angle of 0◦. To match the orientation of the right ovary, such that top of the left ovary is in the positive angular direction, the y-component of the cartesian coordinates are multiplied by −1 before transforming in the polar coordinate system.   Right Ovary:  If analyzing the right ovary, defined as having ventral on the right side, no rotations or inverses are necessary.
