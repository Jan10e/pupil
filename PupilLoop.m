%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% README

% author: Jantine Broek
% date: May 2017
% for: McCormick lab, Yale University, New Haven, USA

% Pre-req:
% - fit_ellipse from https://www.mathworks.com/matlabcentral/fileexchange/3215-fit-ellipse
% Conic Ellipse representation = a*x^2+b*x*y+c*y^2+d*x+e*y+f=0
%  (Tilt/orientation for the ellipse occurs when the term x*y exists (i.e. b ~= 0)) 
%   EDIT: made changes in this function to create plots
% - ellipse from https://www.mathworks.com/matlabcentral/fileexchange/289-ellipse-m

% Goal
% (1) Using edge detection, the pupil is found and the coordinates of the
% pupil location will be used to estimate the function for (2) ellipse
% fitting. With the formula for ellipse fitting, the (3) regionprops will
% be used to compare the dynamics of the pupil fluctuation.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load data
cd /home/jantine/newnas/Garrett/'Pupil Videos'/010316/Mouse882/1/010316_Mouse882_1_CompressedROIs/010316_Mouse882_1_Eyes/010316_EyeFrames/set/noEllipse
yourFolder = '/home/jantine/newnas/Garrett/Pupil Videos/010316/Mouse882/1/010316_Mouse882_1_CompressedROIs/010316_Mouse882_1_Eyes/010316_EyeFrames/set/noEllipse';
addpath(yourFolder);

filePattern = fullfile(yourFolder, '*.jpeg');
srcFiles = dir(filePattern)
numFiles = length(srcFiles)
if numFiles == 0
	message = sprintf('There are no jpeg files are in folder:\n%s', yourFolder);
	uiwait(warndlg(message));
else
	fprintf('There are %d files in %s:\n', numFiles, yourFolder);
	for k = 1 : numFiles %change 10 to numFiles to see all files
		fprintf('    %s\n', srcFiles(k).name);
	end
end

  
%% Edge detection + Ellipse fit
contents = dir('*.jpeg') 
n=numel(contents);
for k = 1:n
  filename(k) = {contents(k).name};
  I = imread(filename{k});
  
  % Creating title for images
    [pathstr, name, ext] = fileparts(filename{k});
    ix=strfind(name,'_');
    t = name(ix(4)+1:end);
    title(sprintf('file %s', t))
  
  % Canny edge detection
  I = rgb2gray(I);
  BW1 = edge(I, 'Canny'); 
  figure;
  imshow(BW1)

  % get rid of small blobs
  BW1 = bwareaopen(BW1, 210);
  
  % get the boundary and outline over original image
  boundaries = bwboundaries(BW1);
  x = boundaries{1}(:, 2);
  y = boundaries{1}(:, 1);
  % display original image
  imshow(I); title(sprintf('Outline over image %s', t)); % Outline over original image
  % plot boundaries over image
  hold on;
  plot(x, y, 'g-', 'LineWidth', 2);
  hold off
  
  % get centroid, area, axis, etc
  s = regionprops(BW1, 'Centroid', 'Area', 'Eccentricity', ...
    'MajorAxisLength', 'MinorAxisLength', 'Orientation');
  centroids = cat(1, s.Centroid);
%   majAxis = cat(1, s.MajorAxisLength);
%   minAxis = cat(1, s.MinorAxisLength);
%   area = cat(1, s.Area);
%   orientation = cat(1,s.Orientation);
  
  % if there are more than one centroids (centroids(1,:) is one centroid), this means that it didn't find the
  % pupil. So we need to redefine
  if centroids(:,1) > 1 
      centroids = [];
%       majAxis = [];
%       minAxis = [];
%       area = [];
%       orientation = [];
  
      % get the boundary and outline over original image
      boundaries = bwboundaries(BW1);
      x = boundaries{2}(:, 2);
      y = boundaries{2}(:, 1);
      % display original image
      imshow(I); title(sprintf('Outline 2 over image %s', t)); % Outline over original image
      % plot boundaries over image
      hold on;
      plot(x, y, 'g-', 'LineWidth', 2);
      hold off
      
      % get centroid, area, axis, etc
      s = regionprops(BW1, 'Centroid', 'Area', 'Eccentricity', ...
            'MajorAxisLength', 'MinorAxisLength', 'Orientation');
      centroids = cat(1, s.Centroid);
      centroids = centroids(2,:);
%       
%       majAxis = cat(1, s.MajorAxisLength);
%       majAxis = majAxis(2,:);
%       
%       minAxis = cat(1, s.MinorAxisLength);
%       minAxis = minAxis(2,:);
%       
%       area = cat(1, s.Area);
%       area = area(2,:);
%       
%       orientation = cat(1,s.Orientation);
%       orientation = orientation(2,:);
      
  end  
  
  % plot centroid and boundaries in original image
  imshow(I); title(sprintf('Boundaries + centroid of image %s', t));
  hold on
  plot(x, y, 'g-', 'LineWidth', 1);
  plot(centroids(:,1), centroids(:,2), 'r*'); 
  hold off
  
  
  % fit ellipse  
  pupilEllipse = fit_ellipse(x,y);
  %if ellipse found, use short axis as pupil diameter
     if ~isempty(pupilEllipse.long_axis) && pupilEllipse.long_axis~=0
         %record pupil diameter
          pupilD(k) = pupilEllipse.long_axis;
         %record pupil position
          pupilXY(1,k) = pupilEllipse.Y0_in;
          pupilXY(2,k) = pupilEllipse.X0_in;
     end
     
   % plot ellipse  
   imshow(I);title(sprintf('Pupil Ellipse, image %s', t));
   hold on
   %pupil outline
    if ~isempty(pupilEllipse)
%          ellipse(pupilEllipse.b,pupilEllipse.a,pupilEllipse.phi,pupilEllipse.X0_in,pupilEllipse.Y0_in, 'r'); %the X0_in and Y0_in are reversed in Dov's script
       
        %pupil center
         scatter(pupilEllipse.X0_in,pupilEllipse.Y0_in,'r');  
    
    
    cos_phi = pupilEllipse.cos_phi;
    sin_phi = pupilEllipse.sin_phi;
         
    % rotation matrix to rotate the axes with respect to an angle phi
    R = [ cos_phi sin_phi; -sin_phi cos_phi ];
    
    % the axes
    ver_line        = [ [pupilEllipse.X0 pupilEllipse.X0]; pupilEllipse.Y0+pupilEllipse.b*[-1 1] ];
    horz_line       = [ pupilEllipse.X0+pupilEllipse.a*[-1 1]; [pupilEllipse.Y0 pupilEllipse.Y0] ];
    new_ver_line    = R*ver_line;
    new_horz_line   = R*horz_line;
    
    % the ellipse - is a better fit than above
    theta_r         = linspace(0,2*pi);
    ellipse_x_r     = pupilEllipse.X0 + pupilEllipse.a*cos( theta_r );
    ellipse_y_r     = pupilEllipse.Y0 + pupilEllipse.b*sin( theta_r );
    rotated_ellipse = R * [ellipse_x_r;ellipse_y_r];
    
    % draw
    plot( new_ver_line(1,:),new_ver_line(2,:),'g' ); %major axis
    plot( new_horz_line(1,:),new_horz_line(2,:),'g' ); %minor axis
    plot( rotated_ellipse(1,:),rotated_ellipse(2,:),'r' );

    
    hold off
    end     
  
end

 
%% Get diameter measurement
% We can use the MajorAxisLength and MinorAxisLength to compare across the
% pupils to look for changes. As getting the diameter/radii circles
% transforms your data much. Maybe we need to get more data points of the
% pupil ellipse to compare. 

% and compare area