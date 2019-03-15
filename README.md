# Angular-Radial-Position-Distribution

Description:

Groups a spots population of germ cells into binned populations based on their positions relative to the ovary's dorsal and ventral sides. The statistics output includes a user defined binned histogram of the radial and angular distribution of the spots population.

Description of Algorithm:

Transforming coordinate systems to quantify angular and radialdistributions of cells within the mouse ovary. The 2D cartesian coordinates of the Imaris spots objects, representing the cells, were recorded. Additionally, in cartesian coordinates, the center of geometry and the end points (top and bottom) of the Imaris surface representing the ovarian outer surface were recorded. Using the the center of geometry and two end points, the radius of curvature was calculated, and all spots object coordinates were offset so as to shift the circle described by the radius of curvature to the cartesian origin.&nbsp;These coordinates were then transformed from cartesian to polar coordinates.     




Left Ovary: If analyzing the left ovary, defined as having the dorsal on the left side, all polar coordinates were rotated by −180◦ so that the center of the ovary is placed at an angle of 0◦. To match the orientation of the right ovary, such that top of the left ovary is in the positive angular direction, the y-component of the cartesian coordinates are multiplied by −1 before transforming in the polar coordinate system.   





Right Ovary: If analyzing the right ovary, defined as having ventral on the right side, no rotations or inverses are necessary.


Dependencies:

-Imaris 8  
-Matlab 2015a  
  
Installation:
  
To make these function available to Imaris, copy the 'spotsRadialDist_right.m' and the 'spotsRadialDist_left.m ' file into the Imaris XT Matlab directory. After copying this file to the XTensions folder and restarting Imaris, you can find this function in the Spots Function menu, as well as in the Image Processing menu under the Spots Functions group.   
  
Executing the XTension:

1) Before execution the user will need to create a Spots population of the germ cells and a Surface object of the ovary.  
2) Upon execution the first dialog box will ask the user if the mouse is older than 13.5 days. If the mouse is older than 13.5 days, the algorithm that computes the radius of curvature of the ovary will fail and a follow up dialog box will appear asking the user to manually input the extremum coordinates of the top and bottom of the ovary.  
3) The output of the measurements will be stored as histograms and the next dialog box will ask the user the number of desired bins. 
4) The number of bins determined by the user will create a sub-populations of the original Spots population.    

Example Output:

1) The output produced includes three histograms which measure the angular position of the germ cells relative to the center of the ovary, the arc distance of the germ cells relative to the center of the ovary, and the user defined binned histogram of the dorsal versus ventral distribution. Additionally, there is a polar plot produced which shows how the ovary is binned.  
2) Excel files of the histograms are also produced.
