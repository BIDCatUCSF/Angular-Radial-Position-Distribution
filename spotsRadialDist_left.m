%  Installation:
%
%  - Copy this file into the XTensions folder in the Imaris installation directory
%  - You will find this function in the Image Processing menu
%
%    <CustomTools>
%      <SurpassTab>
%        <SurpassComponent name="bpSpots">
%          <Item name="Spots Radial Distribution Left" icon="Matlab" tooltip="Find spots close to surface.">
%            <Command>MatlabXT::spotsRadialDist_left(%i)</Command>
%          </Item>
%        </SurpassComponent>
%        <SurpassComponent name="bpSurfaces">
%          <Item name="Spots Radial Distribution Left" icon="Matlab" tooltip="Find spots close to surface.">
%            <Command>MatlabXT::spotRadialDist_left(%i)</Command>
%          </Item>
%        </SurpassComponent>
%      </SurpassTab>
%    </CustomTools>
% 
%
%  Description:
%   
%   Bins a spots population into subpopulations based on the radius of curvature of
%   a surface encasing the spots. Left corresponds to the left facing ovary.
%
%   Authors: Adam Fries and Bikem Soygur 2021
% 	
%	Modification History:
% 	February 2021 -  distribution mask can now be applied to a second Spots object 

function spotsRadialDist_left(aImarisApplicationID)


 warning( 'off', 'MATLAB:xlswrite:AddSheet' ) ;


% connect to Imaris interface
if ~isa(aImarisApplicationID, 'Imaris.IApplicationPrxHelper')
    javaaddpath ImarisLib.jar
    vImarisLib = ImarisLib;
    if ischar(aImarisApplicationID)
        aImarisApplicationID = round(str2double(aImarisApplicationID));
    end
    vImarisApplication = vImarisLib.GetApplication(aImarisApplicationID);
else
    vImarisApplication = aImarisApplicationID;
end

% the user has to create a scene with some spots and surface
vSurpassScene = vImarisApplication.GetSurpassScene;
if isequal(vSurpassScene, [])
    msgbox('Please create some Spots and Surface in the Surpass scene!')
    return
end
vNumChans = vImarisApplication.GetDataSet.GetSizeC;
% get the spots and the surface object
vSpots = vImarisApplication.GetFactory.ToSpots(vImarisApplication.GetSurpassSelection);
vSurfaces = vImarisApplication.GetFactory.ToSurfaces(vImarisApplication.GetSurpassSelection);

vSpotsSelected = ~isequal(vSpots, []);
vSurfaceSelected = ~isequal(vSurfaces, []);
if vSpotsSelected
    vParent = vSpots.GetParent;
elseif vSurfaceSelected
    vParent = vSurfaces.GetParent;
else
    vParent = vSurpassScene;
end

% get the spots and surfaces
vSpotsSelection = 1;
vSurfaceSelection = 1;
vNumberOfSpots = 0;
vNumberOfSurfaces = 0;
vSpotsList = [];
vSurfacesList = [];
vSpotsName = {};
vSurfacesName = {};
for vIndex = 1:vParent.GetNumberOfChildren
    vItem = vParent.GetChild(vIndex-1);
    if vImarisApplication.GetFactory.IsSpots(vItem)
        vNumberOfSpots = vNumberOfSpots + 1;
        vSpotsList(vNumberOfSpots) = vIndex;
        vSpotsName{vNumberOfSpots} = char(vItem.GetName);
        
        if vSpotsSelected && isequal(vItem.GetName, vSpots.GetName)
            vSpotsSelection = vNumberOfSpots; 
        end
    elseif vImarisApplication.GetFactory.IsSurfaces(vItem)
        vNumberOfSurfaces = vNumberOfSurfaces + 1;
        vSurfacesList(vNumberOfSurfaces) = vIndex;
        vSurfacesName{vNumberOfSurfaces} = char(vItem.GetName);
        
        if vSurfaceSelected && isequal(vItem.GetName, vSurfaces.GetName)
            vSurfaceSelection = vNumberOfSurfaces;
        end
    end
end

% check that both spots and surfaces have been created
if min(vNumberOfSpots,vNumberOfSurfaces) == 0
    msgbox('Please create some spots AND a surface object!')
    return
end

% drop down menu for choosing surface objects
if vNumberOfSpots>1
    [vSpotsSelection,vOk] = listdlg('ListString',vSpotsName, ...
        'InitialValue', vSpotsSelection, 'SelectionMode','multiple', ...
        'ListSize',[300 300], 'Name','Find Spots Close To Surface', ...
        'PromptString',{'Please select the spots:'});
    if vOk<1, return, end
end

% get the number of surfaces and spots objects selected by the user
if vNumberOfSurfaces>1
    [vSurfaceSelection,vOk] = listdlg('ListString',vSurfacesName, ...
        'InitialValue', vSurfaceSelection, 'SelectionMode','multiple', ...
        'ListSize',[300 300], 'Name','Find Spots Close To Surface', ...
        'PromptString',{'Please select the surface:'});
    if vOk<1, return, end
end

% get the number of surfaces and spots objects selected by the user
vNumberOfSurfacesSelected = numel(vSurfaceSelection);
vNumberOfSpotsSelected = numel(vSpotsSelection);

% list of colors to rotate through for display
colors = {'ffff00', 'ff00ff', '00ffff', 'ff0000', '00ff00', '00ffff', 'ffff00', ...
    'ff00ff', 'ff0000', '00ff00', '00ffff', 'ff0000', 'ff00ff', 'ffff00'};

% loop through the surfaces, for this application there is just the one surface
%	but we may expand the imaging for multiple surfaces in the future
for vSurfaceIndex = 1:vNumberOfSurfacesSelected
    vItem = vParent.GetChild(vSurfacesList( ...
        vSurfaceSelection(vSurfaceIndex)) - 1);
    vSurface = vImarisApplication.GetFactory.ToSurfaces(vItem);

	% surface vertices stats via GetVertices were discontinued for Imaris versions later
	%	than 8.3.1
    vSurfaceVertices = [];
    for vIndex = 0:vSurface.GetNumberOfSurfaces - 1
      vSurfaceVertices = [vSurfaceVertices; vSurface.GetVertices(vIndex)];
    end
    vNumberOfVertices = size(vSurfaceVertices, 1);
    
    % limit the memory usage to 3*10000*vNumberOfSpots for each block
    vBlockLimit = 10000;
    vBlockIndices = 1:vBlockLimit:vNumberOfVertices;
    vNumberOfBlock = size(vBlockIndices,2);
 
    % obtain the center of geometry for the surface object, 
	%	when the user wants an automated calculation for the 
	% 	radius of curvature, this value is used as one of the 
	%	3 points
	 surf_com = vSurface.GetCenterOfMass(0);
    
	% loop through the spots objects selected and grab the relavent stats		
    for vSpotsIndex = 1:vNumberOfSpotsSelected
        vItem = vParent.GetChild(vSpotsList( ...
            vSpotsSelection(vSpotsIndex)) - 1);
        vSpots = vImarisApplication.GetFactory.ToSpots(vItem);
        vSpotsStats = vSpots.GetStatistics;
        vSpotsValues = vSpotsStats.mValues;
        vSpotsIDs = vSpots.GetSelectedIds;
		vSpotsPosition = vSpots.GetPositionsXYZ;
		vSpotsTime = vSpots.GetIndicesT;
		vSpotsRadius = vSpots.GetRadiiXYZ;
        vNumberOfSpots = size(vSpotsPosition, 1);
        chars2comp = length('Intensity Mean');
        spotNames = cell(vSpotsStats.mNames);
        
        %% dorsal ventril histogram
        zzzpos = vSpotsPosition(:,3);
        zzzbins = 3;
        [uu, vv] = hist(zzzpos, zzzbins);
        figure(10)
        
        vDistancesMin = [];
       
		% grab the file name 
        fname = strcat(char(vSurface.GetName), char(vSpots.GetName));
        
        % from the vertices grab the coords for automatic radius of curvature calculation
        topind = (vSurfaceVertices(:,2)==max(vSurfaceVertices(:,2)));
        ptopy = vSurfaceVertices(topind,2);
        ptopx = vSurfaceVertices(topind,1);
        ptopy = ptopy(1);
        ptopx = ptopx(1);
        
        botind = (vSurfaceVertices(:,2)==min(vSurfaceVertices(:,2)));
        pboty = vSurfaceVertices(botind,2);
        pbotx = vSurfaceVertices(botind,1);
        pboty = pboty(1);
        pbotx = pbotx(1);
        
        midind = (vSurfaceVertices(:,1)==min(vSurfaceVertices(:,1)));
        
		% user input for automatic or manual radius of curvature calculation
		%	mice older than 13.5 days work better with manual measurements of 
		% 	extents of the surface 
        yn = questdlg('Is the mouse older than 13.5 days?', 'Custom Ends?','Yes', 'No', 'No');
         
        if strcmp(yn, 'Yes')
            prompt = {'Input Anterior (Top) x-coordinate', 'Input Top y-cooridnate', ...
                'Input Posterior (Bottom) x-coordinate', 'Input bottom y-coordinate'};
            dlg_title = 'Anterior, Posterior Input';
            numlines = 1;
            defaultans = {'0', '0', '0', '0'};
            ansinput = inputdlg(prompt, dlg_title, numlines, defaultans);
            indata = str2double(ansinput);
            ptopx = indata(1);
            ptopy = indata(2);
            pbotx = indata(3);
            pboty = indata(4);
        end 
        
        % begin center of mass check
        pmidx = surf_com(1);
        pmidy = surf_com(2);
        
        eps =0.1;
        
        myzs = vSurfaceVertices(:,3);
        surf_ind = (myzs >= (median(myzs) - eps)) & (myzs <= (median(myzs) + eps));
        
      	% radius of cruvature 
        k = kfrom3pts([ptopx pmidx pbotx],[ptopy pmidy pboty]);
        radcurv = 1/k;
        
        surf_outlinex = vSurfaceVertices(surf_ind, 1) - radcurv - pmidx;
        surf_outliney = vSurfaceVertices(surf_ind, 2) - pmidy;
        fsize = 24;
        
		% bins for creating the histogram
        numbinsans = inputdlg('Input Number of Bins', 'Histogram Bins', 1, {'7'});
        numbinsans = str2double(numbinsans);
       
		% go around the circle formed by the surface's radius of curvature
		%	grab the polar coordinates of the spots
        ntheta = 360;
        thetamax = 2*pi;
        dtheta = 2*pi/ntheta;
        totalchunks = zeros(1,ntheta);
        % 4-quadrant arctan
        xcen = pmidx + radcurv;
        ycen = pmidy;
        xxs = vSpotsPosition(:,1) - xcen;
        yys = vSpotsPosition(:,2) - ycen;
        rrs = sqrt(xxs.^2+yys.^2);
        
       	% figure 1: histogram of the radial distribution 
        figure(1)
        tbins = calcnbins(rrs);
        tbins = numbinsans;
        [py,px] = hist(rrs, tbins);
        plot(px, py, '-o');
        axis([min(px),(max(px)),0,max(py)+10]);
        set(gca,'FontSize',fsize, 'FontName', 'Times')
        xlabel('Distance from Polar Origin [\mum]')
        ylabel('Number/Bin')
        title(char(vSpots.GetName))
        rads = 0:0.05:2.1*pi;
        
        grey = [0.5 0.5 0.5];
        thetas = atan2(yys, xxs) + pi;
        
        % figure 2: begin polar plot
        figure(2)
        pp = polar(thetas, rrs, 'o');
        set(pp, 'markersize', 3);
        hold on
        greyrho = zeros(1, length(rads));
        binwidth = px(2)-px(1);
        for l=1:tbins
            gg(l) = polar(rads,greyrho + px(l)-binwidth/2, '--');
            set(gg(l), 'color', grey);
            hold on
        end
        ggl = polar(rads, greyrho + px(tbins)+binwidth/2, '--');
        set(ggl, 'color', grey);
        hold on
        
		theta_surf = atan2(surf_outliney, surf_outlinex) + pi;
        rr_surf = sqrt(surf_outlinex.^2 + surf_outliney.^2);
        hh = polar(theta_surf, rr_surf, '.r');
        set(hh,'markersize',2)
        hold on
        argmid = atan2(pmidy, pmidx);
        rmid = sqrt(pmidx.^2 + pmidy.^2);
        polar(argmid, rmid, '-r');
        hold on
        
		curv_ptsx = [ptopx pmidx pbotx] - xcen;
        curv_ptsy = [ptopy pmidy pboty] - ycen;
        argcen = atan2(curv_ptsy, curv_ptsx) + pi;
        argr = sqrt(curv_ptsx.^2 + curv_ptsy.^2);
        txtA = 'Anterior';
        txtP = 'Posterior';
        xP = double(argr(1)*cos(argcen(1)));
        yP = double(argr(1)*sin(argcen(1)));
        yP = yP + 0.3*yP;
        xA = double(argr(3)*cos(argcen(3)));
        yA = double(argr(3)*sin(argcen(3)));
        yA = yA + 0.3*yA;
        text(xA, yA, txtA);
        text(xP, yP, txtP);
        
        set(gca,'FontSize',fsize, 'FontName', 'Times');
        thetas = thetas + (thetas < 0)*2*pi;
       
		% loop through angles and count spots within the arclength of the section 
        adamangs =[];
        for i = 0:ntheta
            ang = i*2*pi/ntheta;
           
            thetachunk = thetas>ang&thetas<=(ang+dtheta);
            totalchunks(i+1) = sum(thetachunk);
            angs = ones(sum(thetachunk),1)*ang*180/pi;
            adamangs = [adamangs, angs']; 
            
        end
        
        adamgt180 = adamangs > 180;
        adamangs(adamgt180) = adamangs(adamgt180) - 360;
    
		% create histogram
        [nn, xouty] = hist(adamangs, 7);
        
        % figure 3: angular distrubution from the center of the surface
        figure(3)

        xchunk = 0:dtheta:thetamax;
        nonzers = totalchunks>0;  
        xchunknz = xchunk(nonzers)*180/pi;
        totalchunksnz = totalchunks(nonzers);
        indsgt180 = xchunknz > 180;
        xchunknz(indsgt180) = xchunknz(indsgt180) - 360;
        [xchunknz, xinds] = sort(xchunknz);
        totalchunksnz = totalchunksnz(xinds);
        chunklen = length(xchunknz);
        plot(xchunknz, totalchunksnz, '-o');
        axis([min(xchunknz), max(xchunknz), 0, max(totalchunksnz)+10]);
        set(gca,'FontSize',fsize, 'FontName', 'Times')
        xlabel('Angle from center [\circ]')
        ylabel('Number/Bin')
        title(char(vSpots.GetName)) 
      
		% figure 4: arcdistance distribution from center
        figure(4)
        newthetas = thetas;
        indsbt180 = newthetas > pi;
        newthetas(indsbt180) = newthetas(indsbt180) - 2*pi;
        arcss = newthetas.*rrs;
        nbins4 = chunklen;
        [yarc, xarc] = hist(arcss, nbins4);
        plot(xarc, yarc, '-o');
        axis([min(xarc),max(xarc),0,max(yarc)+10]);
        set(gca,'FontSize',fsize, 'FontName', 'Times')
        xlabel('Arc Distance from Center [\mum]')
        ylabel('Number/Bin')
        title(char(vSpots.GetName))
        
        % write the statistics files using xlswrite
        xlswrite(strcat(fname, '_polardistance.xls'), [px;py]', 'polar_distance');
        xlswrite(strcat(fname, '_anteriorangles.xls'), [xchunknz;totalchunksnz]', 'anterior_angles');
        xlswrite(strcat(fname, '_anteriorangles_binned.xls'), nn;xouty]', 'anterior_angles_binned');
        xlswrite(strcat(fname, '_DV.xls'), [vv;uu]', '3binnedbyz');
              
        % loop through the bin centers and filter primary spots population, create new spots 
		%	groups	 
        for jj = 1:length(px)
            if jj == length(px)
                vSpotsrbin = rrs >= (px(jj)-binwidth/2) & rrs <= (px(jj)+binwidth/2);
            else
                vSpotsrbin = rrs >= (px(jj)-binwidth/2) & rrs < (px(jj)+binwidth/2);
            end
            if any(vSpotsrbin)
                vNewSpotsrBin = vImarisApplication.GetFactory.CreateSpots;
                vNewSpotsrBin.Set(vSpotsPosition(vSpotsrbin, :), ...
                    vSpotsTime(vSpotsrbin), zeros(sum(vSpotsrbin),1));
                vNewSpotsrBin.SetRadiiXYZ(vSpotsRadius(vSpotsrbin,:));
                vNewSpotsrBin.SetColorRGBA(hex2dec(colors(jj)));
                vNewSpotsrBin.SetName(sprintf('%s from [%.2f um] to [%.2f um]', ...
                    char(vSpots.GetName), px(jj)-binwidth/2, px(jj)+binwidth/2));
                vParent.AddChild(vNewSpotsrBin, -1);
            end
        end
    end
end

% get the 2nd spots population for which to apply the borders defined by 
%  the first spots population
if vNumberOfSpots>1
    [vSpotsSelection,vOk] = listdlg('ListString',vSpotsName, ...
        'InitialValue', vSpotsSelection, 'SelectionMode','multiple', ...
        'ListSize',[300 300], 'Name','Find Spots Close To Surface', ...
        'PromptString',{'Please select the Secondary dependent spots:'});
    if vOk<1, return, end
end

vItem = vParent.GetChild(vSpotsList(vSpotsSelection(vSpotsIndex)) - 1);
vSpots = vImarisApplication.GetFactory.ToSpots(vItem);
vSpotsPosition = vSpots.GetPositionsXYZ;

xxs = vSpotsPosition(:,1) - xcen;
yys = vSpotsPosition(:,2) - ycen;
rrs = sqrt(xxs.^2+yys.^2);

% loop through the secondary spots population and filter based on 
%  primary spots binning
sec_spots = zeros(length(px), 1);

for jj = 1:length(px)
            if jj == length(px)
                vSpotsrbin = rrs >= (px(jj)-binwidth/2) & rrs <= (px(jj)+binwidth/2);
                sec_num = sum(vSpotsrbin);
            else
                vSpotsrbin = rrs >= (px(jj)-binwidth/2) & rrs < (px(jj)+binwidth/2);
                sec_num = sum(vSpotsrbin);
            end
            if any(vSpotsrbin)
                vNewSpotsrBin = vImarisApplication.GetFactory.CreateSpots;
                vNewSpotsrBin.Set(vSpotsPosition(vSpotsrbin, :), ...
                    vSpotsTime(vSpotsrbin), zeros(sum(vSpotsrbin),1));
                vNewSpotsrBin.SetRadiiXYZ(vSpotsRadius(vSpotsrbin,:));
                vNewSpotsrBin.SetColorRGBA(hex2dec(colors(jj)));
                vNewSpotsrBin.SetName(sprintf('%s from [%.2f um] to [%.2f um]', ...
                    char(vSpots.GetName), px(jj)-binwidth/2, px(jj)+binwidth/2));
                vParent.AddChild(vNewSpotsrBin, -1);
                sec_spots(jj) = sec_num;
            end
end

% write the stats from the secondary spots objects            
xlswrite(strcat(fname, '_2ndary_polardistance.xls'), sec_spots', '2ndary polar_distance');     
end

function nbins = calcnbins(x, method, minimum, maximum)
% Calculate the "ideal" number of bins to use in a histogram, using a
% choice of methods.
% 
% NBINS = CALCNBINS(X, METHOD) calculates the "ideal" number of bins to use
% in a histogram, using a choice of methods.  The type of return value
% depends upon the method chosen.  Possible values for METHOD are:
% 'fd': A single integer is returned, and CALCNBINS uses the
% Freedman-Diaconis method,
% based upon the inter-quartile range and number of data.
% See Freedman, David; Diaconis, P. (1981). "On the histogram as a density
% estimator: L2 theory". Zeitschrift fr Wahrscheinlichkeitstheorie und
% verwandte Gebiete 57 (4): 453-476.

% 'scott': A single integer is returned, and CALCNBINS uses Scott's method,
% based upon the sample standard deviation and number of data.
% See Scott, David W. (1979). "On optimal and data-based histograms".
% Biometrika 66 (3): 605-610.
% 
% 'sturges': A single integer is returned, and CALCNBINS uses Sturges'
% method, based upon the number of data.
% See Sturges, H. A. (1926). "The choice of a class interval". J. American
% Statistical Association: 65-66.
% 
% 'middle': A single integer is returned.  CALCNBINS uses all three
% methods, then picks the middle (median) value.
% 
% 'all': A structure is returned with fields 'fd', 'scott' and 'sturges',
% each containing the calculation from the respective method.
% 
% NBINS = CALCNBINS(X) works as NBINS = CALCNBINS(X, 'MIDDLE').
% 
% NBINS = CALCNBINS(X, [], MINIMUM), where MINIMUM is a numeric scalar,
% defines the smallest acceptable number of bins.
% 
% NBINS = CALCNBINS(X, [], MAXIMUM), where MAXIMUM is a numeric scalar,
% defines the largest acceptable number of bins.
% 
% Notes: 
% 1. If X is complex, any imaginary components will be ignored, with a
% warning.
% 
% 2. If X is an matrix or multidimensional array, it will be coerced to a
% vector, with a warning.
% 
% 3. Partial name matching is used on the method name, so 'st' matches
% sturges, etc.
% 
% 4. This function is inspired by code from the free software package R
% (http://www.r-project.org).  See 'Modern Applied Statistics with S' by
% Venables & Ripley (Springer, 2002, p112) for more information.
% 
% 5. The "ideal" number of depends on what you want to show, and none of
% the methods included are as good as the human eye.  It is recommended
% that you use this function as a starting point rather than a definitive
% guide.
% 
% 6. The wikipedia page on histograms currently gives a reasonable
% description of the algorithms used.
% See http://en.wikipedia.org/w/index.php?title=Histogram&oldid=232222820
% 
% Examples:     
% y = randn(10000,1);
% nb = calcnbins(y, 'all')
%    nb = 
%             fd: 66
%          scott: 51
%        sturges: 15
% calcnbins(y)
%    ans =
%        51
% subplot(3, 1, 1); hist(y, nb.fd);
% subplot(3, 1, 2); hist(y, nb.scott);
% subplot(3, 1, 3); hist(y, nb.sturges);
% y2 = rand(100,1);
% nb2 = calcnbins(y2, 'all')
%    nb2 = 
%             fd: 5
%          scott: 5
%        sturges: 8
% hist(y2, calcnbins(y2))
% 
% See also: HIST, HISTX
% 
% $ Author: Richard Cotton $		$ Date: 2008/10/24 $    $ Version 1.5 $

% Input checking
narginchk(1, 4); %#ok<*NCHKN>

if ~isnumeric(x) && ~islogical(x)
    error('calcnbins:invalidX', 'The X argument must be numeric or logical.')
end

if ~isreal(x)
   x = real(x);
   warning('calcnbins:complexX', 'Imaginary parts of X will be ignored.');
end

% Ignore dimensions of x.
if ~isvector(x)
   x = x(:);
   warning('calcnbins:nonvectorX', 'X will be coerced to a vector.');
end

nanx = isnan(x);
if any(nanx)
   x = x(~nanx);
   warning('calcnbins:nanX', 'Values of X equal to NaN will be ignored.');
end

if nargin < 2 || isempty(method)
   method = 'middle';
end

if ~ischar(method)
   error('calcnbins:invalidMethod', 'The method argument must be a char array.');
end

validmethods = {'fd'; 'scott'; 'sturges'; 'all'; 'middle'};
methodmatches = strmatch(lower(method), validmethods);
nmatches = length(methodmatches);
if nmatches~=1
   error('calnbins:unknownMethod', 'The method specified is unknown or ambiguous.');
end
method = validmethods{methodmatches};

if nargin < 3 || isempty(minimum)
   minimum = 1;
end

if nargin < 4 || isempty(maximum)
   maximum = Inf;
end
   
% Perform the calculation
switch(method)
   case 'fd'
      nbins = calcfd(x);
   case 'scott'
      nbins = calcscott(x);
    case 'sturges'
      nbins = calcsturges(x);
   case 'all'
      nbins.fd = calcfd(x);    
      nbins.scott = calcscott(x);
      nbins.sturges = calcsturges(x);
   case 'middle'
      nbins = median([calcfd(x) calcscott(x) calcsturges(x)]);
end

% Calculation details
   function nbins = calcfd(x)
      h = diff(prctile0(x, [25; 75])); %inter-quartile range
      if h == 0
         h = 2*median(abs(x-median(x))); %twice median absolute deviation
      end
      if h > 0
         nbins = ceil((max(x)-min(x))/(2*h*length(x)^(-1/3)));
      else
         nbins = 1;
      end
      nbins = confine2range(nbins, minimum, maximum);
   end

   function nbins = calcscott(x)
      h = 3.5*std(x)*length(x)^(-1/3);
      if h > 0 
         nbins = ceil((max(x)-min(x))/h);
      else 
         nbins = 1;
      end
      nbins = confine2range(nbins, minimum, maximum);
   end

   function nbins = calcsturges(x)
      nbins = ceil(log2(length(x)) + 1);
      nbins = confine2range(nbins, minimum, maximum);
   end

   function y = confine2range(x, lower, upper)
      y = ceil(max(x, lower));
      y = floor(min(y, upper));
   end

   function y = prctile0(x, prc)
      % Simple version of prctile that only operates on vectors, and skips
      % the input checking (In particluar, NaN values are now assumed to
      % have been removed.)
      lenx = length(x);
      if lenx == 0
         y = [];
         return
      end
      if lenx == 1
         y = x;
         return
      end
      
      function foo = makecolumnvector(foo)
         if size(foo, 2) > 1 
            foo = foo';
         end
      end
         
      sortx = makecolumnvector(sort(x));
      posn = prc.*lenx/100 + 0.5;
      posn = makecolumnvector(posn);
      posn = confine2range(posn, 1, lenx);
      y = interp1q((1:lenx)', sortx, posn);
   end
end


function k = kfrom3pts(xs,ys)

%KFROM3POINTS Calculate curvature of a circle given 3 points
%
% K = KFROM3POINTS(XS,YS)
% Where
% XS holds 3 x values
% YS holds 3 y values
% Then
% K is the curvature of the circle through the points

xs=xs(:); % columnize xs
ys=ys(:); % columnize ys

os = ones(3,1);
ss = xs.^2+ys.^2; % sum of squares
a = det([xs ys os]); % Eq. 31
d = -det([ss ys os]); % Eq. 32
e = det([ss xs os]); % Eq. 33
f = -det([ss xs ys]); % Eq. 34
r = sqrt((d^2+e^2)/(4*a^2)-(f/a)); % Eq. 30
k=1/r; % curvature 

end



