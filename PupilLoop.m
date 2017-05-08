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
% - natsortfiles from http://www.mathworks.com/matlabcentral/fileexchange/47434-natural-order-filename-sort

% Goal
% (1) Using edge detection, the pupil is found and the coordinates of the
% pupil location will be used to estimate the function for (2) ellipse
% fitting. With the formula for ellipse fitting, the (3) regionprops will
% be used to compare the dynamics of the pupil fluctuation.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load data
cd /home/jantine/newnas/Garrett/'Pupil Videos'/010316/Mouse882/1/010316_Mouse882_1_CompressedROIs/010316_Mouse882_1_Eyes/010316_EyeFrames/set
yourFolder = '/home/jantine/newnas/Garrett/Pupil Videos/010316/Mouse882/1/010316_Mouse882_1_CompressedROIs/010316_Mouse882_1_Eyes/010316_EyeFrames/set';
% cd /Volumes/'Seagate Backup Plus Drive'/Pupil_Garrett/set
% yourFolder = '/Volumes/Seagate Backup Plus Drive/Pupil_Garrett/set';
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
tic
contents = dir('*.jpeg') 
n = natsortfiles({contents.name}); % put files in natural order
[~,ndx] = natsortfiles({contents.name}); % put files in natural order
contents = contents(ndx);
for k = 1:numel(n)
  filename(k) = {contents(k).name};
  I = imread(filename{k});
  
  % Creating title for images
  [pathstr, name, ext] = fileparts(filename{k});
  ix=strfind(name,'_');
  t = name(ix(4)+1:end);
  %title(sprintf('file %s', t));
  
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
  
  % get centroid and area
  s = regionprops(BW1, 'Centroid', 'Area');
  centroids = cat(1, s.Centroid);
  nCentroids = size(centroids,1);
  edgeArea(k) = s.Area;
  
  % if there are more than one centroids (centroids(1,:) is one centroid), this means that it didn't find the
  % pupil. So we need to redefine
  if nCentroids > 1 
      centroids = [];
  
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
      
      % get area
      s = regionprops(BW1, 'Area');
      edgeArea(k) = s.Area;
      
  end  
  
%   % plot centroid and boundaries in original image + NOT NECESSARY AFTER
%   % DEBUGING
%   imshow(I); title(sprintf('Boundaries + centroid of image %s', t));
%   hold on
%   plot(x, y, 'g-', 'LineWidth', 1);
%   plot(centroids(:,1), centroids(:,2), 'r*'); 
%   hold off
  
  
  % fit ellipse  
  pupilEllipse = fit_ellipse(x,y);
  %if ellipse found, collect long and short axis and pupil position LEAVE
  %IF STATEMENT IN FOR FUTURE VIDEO DETECTION
     if ~isempty(pupilEllipse.long_axis) && pupilEllipse.long_axis~=0
         %record pupil axis
          shortAxis(k) = pupilEllipse.short_axis;
          longAxis(k) = pupilEllipse.long_axis;   
         %record pupil position
          pupilXY(1,k) = pupilEllipse.Y0_in;
          pupilXY(2,k) = pupilEllipse.X0_in;
     end
     
         
   % plot ellipse  
   imshow(I);title(sprintf('Pupil Ellipse, image %s', t));
   hold on
   % pupil outline
    if ~isempty(pupilEllipse)
       
     % pupil center
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
      plot( new_ver_line(1,:),new_ver_line(2,:),'g' ); % major axis
      plot( new_horz_line(1,:),new_horz_line(2,:),'g' ); % minor axis
      plot( rotated_ellipse(1,:),rotated_ellipse(2,:),'r' );

    
    hold off
    end       

end
toc

% extract filename for x-axis (this was done for debugging)
filename = filename';
token = strtok(filename,'.');
D = regexp(token, '_', 'split');
D = vertcat(D{:});
filenr_temp = D(:,5);
S = sprintf('%s*', filenr_temp{:});
filenr = sscanf(S, '%f*');
filenr = filenr';


% plot diameter fluctuation
figure;
hold on
subplot(3,1,1);
plot(edgeArea, '-.og')
xticks(1:1:numel(filenr));
set(gca,'XTickLabel',filenr)
xtickangle(45)
xlabel('filenr');
ylabel('area of edge detection (px)');
legend('edgeArea');

subplot(3,1,2);
plot(shortAxis, '-.ob')
xticks(1:1:numel(filenr));
set(gca,'XTickLabel',filenr)
xtickangle(45)
xlabel('filenr');
ylabel('length axis (px)');
legend('shortAxis');

subplot(3,1,3);
plot(longAxis, '-.or')
xticks(1:1:numel(filenr));
set(gca,'XTickLabel',filenr)
xtickangle(45)
xlabel('filenr');
ylabel('length axis (px)');
legend('longAxis');
hold off

%% Collect data

% collect data
data = [filenr; edgeArea; shortAxis; longAxis]; % as the eye and face of the animal is curved, it is better to take the longAxis as
% a diameter and compare this across time
data = data';

% collect data with DataJoint
tuple.filenr = filenr;
tuple.edgeArea = edgeArea;
tuple.shortAxis = shortAxis;
tuple.longAxis = longAxis;

insert(preprocess.PupilJPEG, tuple)
