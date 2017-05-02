%% Pre-req:
% - fit ellipse from https://www.mathworks.com/matlabcentral/fileexchange/3215-fit-ellipse
% - ellipse from https://www.mathworks.com/matlabcentral/fileexchange/289-ellipse-m

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

contents = dir('*.jpeg') 
n=numel(contents);
for k = 1:n
  filename(k) = {contents(k).name};
  I = imread(filename{k});
  
  % Canny edge detection
  I = rgb2gray(I);
  BW1 = edge(I, 'Canny'); 
  figure;
  imshow(BW1)
  
  % Creating title for images
  [pathstr, name, ext] = fileparts(filename{k});
  ix=strfind(name,'_');
  t = name(ix(4)+1:end);
  title(sprintf('file %s', t))

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
  majAxis = cat(1, s.MajorAxisLength);
  minAxis = cat(1, s.MinorAxisLength);
  area = cat(1, s.Area);
  orientation = cat(1,s.Orientation);
  
  % if there are more than one centroids (centroids(1,:) is one centroid), this means that it didn't find the
  % pupil. So we need to redefine
  if centroids(:,1) > 1 
      centroids = [];
      majAxis = [];
      minAxis = [];
      area = [];
      orientation = [];
  
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
      
      majAxis = cat(1, s.MajorAxisLength);
      majAxis = majAxis(2,:);
      
      minAxis = cat(1, s.MinorAxisLength);
      minAxis = minAxis(2,:);
      
      area = cat(1, s.Area);
      area = area(2,:);
      
      orientation = cat(1,s.Orientation);
      orientation = orientation(2,:);
      
  end
  
  
  % plot centroid and boundaries in original image
  imshow(I); title(sprintf('Boundaries + centroid of image %s', t));
  hold on
  plot(x, y, 'g-', 'LineWidth', 1);
  plot(centroids(:,1), centroids(:,2), 'r*'); 
  
  % plot major and minor axis
  % The axes length is R, compute R*cos(theta), R*sin(theta), and center those displacements on the centroid.

%   hlen = majAxis/2;
%   xCenter = centroids(1);
%   yCenter = centroids(2);
%   cosOrient = cosd(orientation);
%   sinOrient = sind(orientation);
%   xcoords = xCenter + hlen * [cosOrient -cosOrient];
%   ycoords = yCenter + hlen * [-sinOrient sinOrient];
%   line(xcoords, ycoords);
  
  hlen = minAxis/2;
  xCenter = centroids(1);
  yCenter = centroids(2);
  cosOrient = cosd(orientation);
  sinOrient = sind(orientation);
  xcoords = xCenter + hlen .* [cosOrient -cosOrient];
  ycoords = yCenter + hlen .* [-sinOrient sinOrient];
  line(xcoords, ycoords);

  hold off

  
end

 
%% Get diameter measurement
% We can use the MajorAxisLength and MinorAxisLength to compare across the
% pupils to look for changes. As getting the diameter/radii circles
% transforms your data much. Maybe we need to get more data points of the
% pupil ellipse to compare. 

% and compare area
 
%% 

% dus = regionprops(BW1, 'Image');
% BW2 = [dus(1).Image];
% BW3 = [dus(2).Image];
% 
%       % get the boundary and outline over original image
%       boundaries = bwboundaries(BW3);
%       x = boundaries{1}(:, 2);
%       y = boundaries{1}(:, 1);
%       % display original image
%       imshow(I); title(sprintf('Outline 2 over image %s', t)); % Outline over original image
%       % plot boundaries over image
%       hold on;
%       plot(x, y, 'g-', 'LineWidth', 2);
%       hold off