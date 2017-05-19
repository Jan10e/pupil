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


% filePattern = fullfile(yourFolder, '*.jpeg');
% srcFiles = dir(filePattern)
% numFiles = length(srcFiles)
% if numFiles == 0
% 	message = sprintf('There are no jpeg files are in folder:\n%s', yourFolder);
% 	uiwait(warndlg(message));
% else
% 	fprintf('There are %d files in %s:\n', numFiles, yourFolder);
% 	for k = 1 : numFiles %change 10 to numFiles to see all files
% 		fprintf('    %s\n', srcFiles(k).name);
% 	end
% end


contents = dir('*.jpeg'); 
n = natsortfiles({contents.name}); % put files in natural order
[~,ndx] = natsortfiles({contents.name}); % put files in natural order
contents = contents(ndx);

filename = cell(1,numel(n));
for j = 1:numel(n)
  filename(j) = {contents(j).name};
  I = imread(filename{j});
end

%% Select eye ROI

picEyeAvg = mean(I,3);

%choose region of interest around eye
figure
uiwait(msgbox('Select region of interest -> Right click -> Create Mask '));
maskEye = roipoly(picEyeAvg./255);
close all

%% Choose threshold value 

threshDark = 25;

figure
uiwait(msgbox('Choose threshDark value so that most pixels of pupil are showing. Press ESC to progress.'));
for iTest = 1:15
    randFrame = ceil(rand*size(filename,2));
    disp(['random frame #' num2str(randFrame)])
    
    randImg = imread(filename{randFrame});
    randImg = randImg(:,:,1);
    
    randImg_threshed = randImg<threshDark;
    imshow(randImg_threshed.*maskEye)
    pause
end
close all

  
%% Edge detection + Ellipse fit

%eye may be closed if very few dark pixels found
threshNdarkPix = 200; %default 200
%fudge factor for estimated radius
radiusFudge = 1.4;  %default 1.4
%edge detection threshold fudge factor
edgeFudge = 1;   %default 1 in Canny Edge detection (lower means more edges will be found)
%edge pixels must be at least this close to dark pixels
threshPupilDist = 3;  %default 3
%edge pixels must be continguous with at least this many other pixels
minEdgeSize = 10;     %default 15



tic
for k = 1:numel(n)
  filename(k) = {contents(k).name};
  I = imread(filename{k});
  
  % Creating title for images
  [pathstr, name, ext] = fileparts(filename{k});
  ix=strfind(name,'_');
  t = name(ix(4)+1:end);
  %title(sprintf('file %s', t));
  
  % Display loop progression
  if mod(k,1000)==0
      disp([num2str(k) ' / ' num2str(size(n)) ' image ']);
  end

  % Define image frame and threshold
    im_threshed = (I<threshDark) .* maskEye;
  % Estimate pupil center and size based of number of pixels
    [row1, col1] = find(im_threshed);
  % Check to see if eye is open
    if length(row1) > threshNdarkPix
        pupilCenter_estimate = [median(row1), median(col1)];  % pupil center
        pupilRadius_estimate = sqrt(length(row1)/pi);   % pupil radius (from A=pi*r^2)
        % Delete pixels which are too far away from pupil
        for l = 1:length(row1)
            % Distance from estimated center
            dXY = pupilCenter_estimate - [row1(l), col1(l)];
            d = sqrt(sum(dXY .* dXY));
            % Delete pixels if too far
            if d > pupilRadius_estimate * radiusFudge
                im_threshed(row1(l), col1(l)) = false;
            end
        end
        
   % Canny Edge Detection to find edges
     im = rgb2gray(I);
     [~, threshold] = edge(im, 'canny');
     im_edges = edge(im, 'canny', threshold * edgeFudge).*maskEye;
   % Delete edge pixels which are too far from center
     [row2, col2] = find(im_edges);
     for m = 1:length(row2)
        % distance from estimated center
           dXY = pupilCenter_estimate - [row2(m), col2(m)];
           d = sqrt(sum(dXY .* dXY));
%           % delete if too far - REMOVES ALL PIXELS
%             if d > pupilRadius_estimate * radiusFudge
%                im_edges(row2(m), col2(m)) = false;
%             else
              % delete edge pixels which are too far from dark pixels
                [~,d] = dsearchn([row1, col1], [row2(m), col2(m)]);
                if d > threshPupilDist
                    im_edges(row2(m),col2(m)) = false;
                end
              %end
     end       
        
    
    % Delete edge pixels which are not part of large group
      im_edges2 = bwlabel(im_edges); %label ROIs by unique number
    % Number of ROIs in edges matrix
      nRois = max(im_edges2(:));
    % Clear edge matrix and add back in if ROI is big
      im_edges = false(size(im_edges));
      for iROI = 1:nRois
            roiSize = sum(sum(im_edges2 == iROI));
            if roiSize >= minEdgeSize
                im_edges(im_edges2==iROI) = true;
            end
      end
      
    % Ellipse Fitting if there are sufficient edges pixels
        if sum(im_edges(:)) > minEdgeSize
            % fit ellipse
            [row3, col3] = find(im_edges);
            pupilEllipse = fit_ellipse(row3, col3);
            % if ellipse found, use long axis as pupil diameter
            if ~isempty(pupilEllipse.long_axis) && pupilEllipse.long_axis~=0
                % record pupil axis
                shortAxis(k) = pupilEllipse.short_axis;
                longAxis(k) = pupilEllipse.long_axis; 
                % record pupil position
                pupilXY(1,k) = pupilEllipse.Y0_in;
                pupilXY(2,k) = pupilEllipse.X0_in;
            else
                error_fit(k) = filename(k); % no ellipse found but parabola or hyperbola
            end
        end
    end
    
        
%   % Canny edge detection
%   I = rgb2gray(I);
%   BW1 = edge(I, 'Canny');
% %   figure;
% %   imshow(BW1)
% 
%   % get rid of small blobs
%   BW1 = bwareaopen(BW1, 160); % change threshold number when error occurs for boundaries ... it might be too strict
%   
%   % get the boundary and outline over original image
%   try
%     boundaries = bwboundaries(BW1);
%     x = boundaries{1}(:, 2);
%     y = boundaries{1}(:, 1);
%   catch ME
%       error_edge(k) = filename(k);      
%       fprintf('Data not collected and skip to next iteration %s\n', ME.message);
%       continue; % skips whole procedure for this file and starts new iteration
%   end
% %   % display original image
% %   imshow(I); title(sprintf('Outline over image %s', t)); % Outline over original image
% %   % plot boundaries over image
% %   hold on;
% %   plot(x, y, 'g-', 'LineWidth', 2);
% %   hold off
%    
%   % get centroid and area
%   s = regionprops(BW1, 'Centroid', 'Area');
%   centroids = cat(1, s.Centroid);
%   nCentroids = size(centroids,1);
%   %edgeArea(k) = s.Area;
%   
%   % if there is more than one centroid (centroids(1,:) is one centroid), this means that it didn't find the
%   % pupil. So we need to redefine
%   if nCentroids > 1 
%       centroids = [];
%   
%       % get the boundary and outline over original image
%       boundaries = bwboundaries(BW1);
%       x = boundaries{2}(:, 2);
%       y = boundaries{2}(:, 1);
% %       % display original image
% %       imshow(I); title(sprintf('Outline 2 over image %s', t)); % Outline over original image
% %       % plot boundaries over image
% %       hold on;
% %       plot(x, y, 'g-', 'LineWidth', 2);
% %       hold off
%       
%       % get area
%       s = regionprops(BW1, 'Area');     
%   end  
%   
%   edgeArea(k) = s.Area;
  
%        % fit ellipse  
%   pupilEllipse = fit_ellipse(x,y);
%   %if ellipse found, collect long and short axis and pupil position LEAVE
%   %IF STATEMENT IN FOR FUTURE VIDEO DETECTION
%      if ~isempty(pupilEllipse.long_axis) && pupilEllipse.long_axis~=0
%          %record pupil axis
%           shortAxis(k) = pupilEllipse.short_axis;
%           longAxis(k) = pupilEllipse.long_axis;   
%          %record pupil position
%           pupilXY(1,k) = pupilEllipse.Y0_in;
%           pupilXY(2,k) = pupilEllipse.X0_in;
%      else
%          error_fit(k) = filename(k); % no ellipse found but parabola or hyperbola
%      end    
%          


%    % plot ellipse
% %     figure;
% %     imshow(im);title(sprintf('Pupil Ellipse 1, image %s', t));
% %     hold on
%     % pupil outline
%         if ~isempty(pupilEllipse)
%        
%          % pupil center
%           scatter( pupilEllipse.Y0_in,pupilEllipse.X0_in,'r' );  
% 
%           cos_phi = pupilEllipse.cos_phi;
%           sin_phi = pupilEllipse.sin_phi;
% 
%           % rotation matrix to rotate the axes with respect to an angle phi
%           R = [ cos_phi sin_phi; -sin_phi cos_phi ];
% 
%           % the axes
%           ver_line        = [ [pupilEllipse.X0 pupilEllipse.X0]; pupilEllipse.Y0+pupilEllipse.b*[-1 1] ];
%           horz_line       = [ pupilEllipse.X0+pupilEllipse.a*[-1 1]; [pupilEllipse.Y0 pupilEllipse.Y0] ];
%           new_ver_line    = R*ver_line;
%           new_horz_line   = R*horz_line;
% 
%           % the ellipse - is a better fit than above
%           theta_r         = linspace( 0,2 * pi );
%           ellipse_x_r     = pupilEllipse.X0 + pupilEllipse.a*cos( theta_r );
%           ellipse_y_r     = pupilEllipse.Y0 + pupilEllipse.b*sin( theta_r );
%           rotated_ellipse = R * [ ellipse_x_r; ellipse_y_r ];
% 
%           % draw
%           plot( new_ver_line(2,:),new_ver_line(1,:),'b' ); % major axis
%           plot( new_horz_line(2,:),new_horz_line(1,:),'g' ); % minor axis
%           plot( rotated_ellipse(2,:),rotated_ellipse(1,:),'r' );
%           
%           hold off
%         end
        
        
        %plot pic with fitted ellipse
        figure;
        imshow(im); title(sprintf('Pupil Ellipse 2, image %s', t));
        set(gcf, 'Position', [50, 50, 500, 500]); % Maximize figure.
        %hold on
        if ~isempty(pupilEllipse)
            
              ellipse(pupilEllipse.b,pupilEllipse.a,pupilEllipse.phi,pupilEllipse.Y0_in,pupilEllipse.X0_in,'r');
            
            % pupil center
              hold on
              scatter(pupilEllipse.Y0_in,pupilEllipse.X0_in,'r'); 
            
            % pupil axis
              cos_phi = pupilEllipse.cos_phi;
              sin_phi = pupilEllipse.sin_phi;

              % rotation matrix to rotate the axes with respect to an angle phi
              R = [ cos_phi sin_phi; -sin_phi cos_phi ];

              % the axes
              ver_line        = [ [pupilEllipse.X0 pupilEllipse.X0]; pupilEllipse.Y0+pupilEllipse.b*[-1 1] ];
              horz_line       = [ pupilEllipse.X0+pupilEllipse.a*[-1 1]; [pupilEllipse.Y0 pupilEllipse.Y0] ];
              new_ver_line    = R*ver_line;
              new_horz_line   = R*horz_line;

              % draw
              plot( new_ver_line(2,:),new_ver_line(1,:),'b' ); % major axis
              plot( new_horz_line(2,:),new_horz_line(1,:),'g' ); % minor axis

              hold off
        end



% get filename and give NaN for when file was skipped due to no ellipse
filename_ell(k) = filename(k);

% get edge
%edgeArea_ell(k) = edgeArea(k);

end
toc




% %% Extract data tuples
% % extract session date, session number, filename for x-axis, movie number (first number of filenr), frame
% % number (rest of numbers in filenr), animal id
%filename_ell = filename_ell';
token = strtok(filename,'.');
D = regexp(token, '_', 'split');
D = vertcat(D{:});

%filenr
% ix = cellfun('isempty', filename_ell);
% filename_ell(ix) = {'nan'};
filenr_temp = D(:,5);
G = sprintf('%s*', filenr_temp{:}); % change cell to double
filenr = sscanf(G, '%f*');
filenr = filenr';


% session date and number
session_date_temp = D(:,1);
session_date = datetime(session_date_temp, 'InputFormat', 'ddMMyy', 'Format', 'yyyy-MM-dd');
session_date = string(session_date); % as string
session_date = char(session_date); % as char temporarily for DataJoint
session_date = session_date';

% E = sprintf('%s*', session_date_temp{:}); 
% session_date = sscanf(E, '%f*'); %as double
% 

% Just creating a column with strings to see whether it works
%session_date = repmat('2016-05-15', numel(session_date_temp), 1);


session_number_temp = D(:,3);
F = sprintf('%s*', session_number_temp{:}); % change cell to double
session_number = sscanf(F, '%f*');


% movie number and frame number
for m = 1:numel(filenr)
    G = num2str(filenr(m));
    movie_number(m) = str2double(G(1));
    frame_number(m) = str2double(G(2:end));
end

% extract animal_id column
animal_id_temp = D(:,2);
animal_id_temp = animal_id_temp(:,1);
animal_id_temp = sprintf('%s*', animal_id_temp{:});
animal_id = regexp(animal_id_temp, '\d*', 'Match');
animal_id = str2double(animal_id');


%% Plot variables
% plot diameter fluctuation
figure;
hold on
% subplot(3,1,1);
% plot(edgeArea, '-.og')
% xticks(1:1:numel(filenr));
% set(gca,'XTickLabel',filenr)
% xtickangle(45)
% xlabel('filenr');
% ylabel('area of edge detection (px)');
% legend('edgeArea');

subplot(2,1,1); % (3,1,2)
plot(shortAxis, '-.ob')
xticks(1:1:numel(filenr));
set(gca,'XTickLabel',filenr)
xtickangle(45)
xlabel('filenr');
ylabel('length axis (px)');
legend('shortAxis');

subplot(2,1,2); %(3,1,3)
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
data(data == 0) = NaN;
data = data';
fname = sprintf('pupilAxis_Mouse%d.mat', animal_id(1:1));
save(fname, 'data');

[row, col] = find(isnan(data));
nrNaN = unique(row);

% find files where edge detection didn't work
[~, error_edge_col] = find(~cellfun('isempty', error_edge));
error_edge_filenr = error_edge(error_edge_col)';




% % Error upon analyzing whole dataset
% Error using vertcat
% Dimensions of matrices being concatenated are not consistent.
% 
% Error in PupilLoop (line 262)
% data = [filenr; edgeArea; shortAxis; longAxis]; % as the eye and face of the animal is curved, it is better to take the longAxis
% as

% axis = [shortAxis; longAxis]';

addpath /media/jantine/Data/04_DataJoint/2PE/schemas

% collect data for schema preprocess.EyeROI
tuple.animal_id = animal_id;
tuple.session_date = session_date; % should be in string format "2017-05-12"
tuple.session_number = session_number;
tuple.movie_number = movie_number';
tuple.frame_number = frame_number';
tuple.edge_area = edgeArea';
tuple.short_axis = shortAxis';
tuple.long_axis = longAxis';

% order structure as DataJoint keys
a = preprocess.EyeROI;
fields = a.header.names;
tuple = orderfields(tuple, fields);

insert(preprocess.EyeROI, dj.struct.fromFields(tuple))
