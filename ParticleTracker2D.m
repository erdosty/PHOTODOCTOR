%===============================================================================
% 2D Particle Tracker
%
% This code converts a video file (.avi or .mp4) into mat. file and track particles in 2D plane.
%
% If you use this tool for your work, we kindly ask you to cite the following article for which it was included:
%
% Sridhar, V., Yildiz, E., Rodríguez-Camargo, A., Lyu, X., Yao, L., Wrede, P., Aghakhani, A., Akolpoglu, M.B., Podjaski, F., Lotsch, B.V. and Sitti, M., 2023. 
% Designing Covalent Organic Framework-based Light-driven Microswimmers towards Intraocular Theranostic Applications. arXiv preprint arXiv:2301.13787.
%
% This work is licensed under a Creative Commons Attribution-NonCommercial-NoDerivatives 4.0 International License.
% Copyright (c) 2023, Erdost YILDIZ
% All rights reserved.
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:
% 
%     * Redistributions of source code must retain the above copyright
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright
%       notice, this list of conditions and the following disclaimer in
%       the documentation and/or other materials provided with the distribution
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.
%===============================================================================

clear all
close all
clc

name = 'example2';


fileFolder = '';
dataName = [name,'.avi'];
outFileName = [name,'.mat'];

carobj = VideoReader([fileFolder, dataName]);
nFrames = carobj.NumberOfFrames;
M = carobj.Height;
N = carobj.Width;
data = zeros(M,N,nFrames,'uint8');
frameRate = carobj.FrameRate;

for k= 1 : nFrames
    im= read(carobj,k);
    im=im(:,:,1);
    data(:,:,k)=im;
end

save(outFileName, 'data', 'frameRate', '-v7.3');
warning('off','MATLAB:MKDIR:DirectoryExists')
dataFolder = '';
resultsFolder = '';
dataFiles{1} = [name];
resultsNames{1} = [dataFiles{1},'_calData'];

%%
for fileNum=1:length(dataFiles)
    % Location to Store Results
    dataFile=dataFiles{fileNum};
    resultsName=resultsNames{fileNum};
    disp(sprintf('Processing file %d ...', fileNum));
    tic
    
    %-------------------------------------------------------------------------------
    % FLAGS
    % All flags must be set to either 0 for false or 1 for true
    dropCollisions  = 0; % Whether to detect and remove possibly colliding paths
    is3D            = 0; % Whether to try to calculate 3D trajectories
    tumbleAnalysis  = 0; % Whether to try to calculate tumbling related parameters
    manualTumble    = 0; % Whether to include manual checking of tumble events
    trajJoining     = 0; % Whether to attempt connecting different parts of same path
    manualJoining   = 0; % Whether to attempt to manually join trajectories
    writeVideo      = 1; % Whether to create video showing calculated trajectories
    singleparticle  = 0; % Whether attempting to track only a single particle
    
    flags = struct('dropCollisions', dropCollisions, 'is3D', is3D,...
        'tumbleAnalysis', tumbleAnalysis, 'manualTumble', manualTumble,...
        'trajJoining', trajJoining, 'writeVideo', writeVideo, 'manualJoining', ...
        manualJoining, 'singleparticle', singleparticle);
    clear dropCollisions is3D tumbleAnalysis manualTumble trajJoining writeVideo
    clear manualJoining singleparticle
    
    %-------------------------------------------------------------------------------
    % TRACKING PARAMETERS
    
    scale     = 2.27; % pix/um, for the fluorescent video.
    avgSpace  = 5;  % pix - spatial averaging (median) window width
    avgTime   = 0;  % frames - temporal averaging (median) window width -
    % # of frames in each direction
    threshold = 40; % threshold to pick out likable spots
    minArea   = 10;  % Pixels
    maxArea   = 500; % Pixels
    minLength = 50; % Frames - minimum time length trajectory to consider valid
    maxDist   = 10;  % Pixels
    profDist3 = 100; % Pixels, radial distance to find rings for 3D tracking
    searchWidth = 10; % pixels, width rage for searching particle in the next frame (10)
    minSpeed = 8; % min speed to count as moving swimmer
    
    trackParams = struct('scale', scale, 'avgSpace', avgSpace, 'avgTime', avgTime,...
        'threshold', threshold, 'minArea', minArea, 'maxArea', maxArea,...
        'minLength', minLength, 'maxDist', maxDist, 'profDist3', profDist3, 'searchWidth', searchWidth, 'minSpeed', minSpeed);
    clear scale avgSpace avgTime threshold minArea maxArea minLength maxDist profDist3
    
    %-------------------------------------------------------------------------------
    % ANALYSIS PARAMETERS
    
    velocitySmoothWin = 5; % frames - time over which to smooth paths for v calc
    headingSmoothWin = 5; % frames - time over which to smooth paths for heading calc.
    tumbleType = 'headingThreshold'; % What kind of tumble identification algorithm
    tumbleVelThresh = 1; % Minimum average trajectory speed for which to attempt tumble calc
    tumbleJoinWin = 5; % frames - time over which to consider tumbling to be single event
    angleThresh = 0.05; % Radians, angular difference which indicates a tumble
    
    anlysParams = struct('velocitySmoothWin', velocitySmoothWin, 'tumbleType',...
        tumbleType, 'tumbleVelThresh', tumbleVelThresh, 'tumbleJoinWin',...
        tumbleJoinWin);
    
    %-------------------------------------------------------------------------------
    % LOAD DATA
    
    lData = load([dataFolder, dataFile]);
    data = lData.data;
    clear dataf
    
    back=median(double(data(:,:,1:20:end)),3);
    
    % Get timing information
    times = [];
    vidLength = size(data, 3);
    frameRate = lData.frameRate;
    stdTimes = std(diff(times));
    
    % Calculate background image
    
    

    %-------------------------------------------------------------------------------
    % TRACKING
    paths = trackpart2D(data, vidLength, back, trackParams);   
    tracks = [paths{:}];
    tracks = tracks([tracks.len] > trackParams.minLength);
    
    %-------------------------------------------------------------------------------
    %% joining trajectories
    if flags.trajJoining == 1
        % Determine whether trajectory has ended because edge of screen was reached,
        % because two particle overlapped, or because particle moved out of focus.
        for i = 1:length(tracks)
            
            % Determine trajectory termination--------------------------------------
            xf = round(tracks(i).pathx(tracks(i).len));
            yf = round(tracks(i).pathy(tracks(i).len));
            
            if tracks(i).active == 1
                tracks(i).endType = 'video';
                tracks(i).active = 0;
                
            elseif ((xf - trackParams.maxDist < 0) || ...
                    (xf + trackParams.maxDist > size(data, 2))) && ...
                    ((yf - trackParams.maxDist < 0) || ...
                    (yf + trackParams.maxDist > size(data, 1)))
                tracks(i).endType = 'wall'; % a little bit not sure(Jiang)
                
            else
                %trackingFigure(data, tracks, i, xf, yf);
                tracks(i).endType = 'unknown';
            end
            
            % Determine trajectory beginning----------------------------------------
            xi = round(tracks(i).pathx(1));
            yi = round(tracks(i).pathy(1));
            
            if tracks(i).start == 1
                tracks(i).startType = 'video';
                
            elseif ((xi - 3*trackParams.maxDist < 0) || ...
                    (xi + 3*trackParams.maxDist > size(data, 2))) && ...
                    ((yi - 3*trackParams.maxDist < 0) || ...
                    (yi + 3*trackParams.maxDist > size(data, 1)))
                tracks(i).startType = 'wall';
                
            else
                tracks(i).startType = 'unknown';
            end
        end
        
        % Begin automatically joining trajectories - if speed, angle, and position
        % are approximately the same, join two separte trajectories together,
        % filling in the spaces with linear interpolation.
        for i = 1:length(tracks)
            if strcmp(tracks(i).endType, 'unknown')
                for j = i+1:length(tracks)
                    tracks = joinTrajectories(tracks, i, j, 60, 20, ...
                        trackParams, anlysParams);
                end
            end
        end
        tracks = tracks([tracks.active] == 0);
    end
       
    %%
    %-------------------------------------------------------------------------------
    % SCALING
    
    % Convert x and y components of path to microns from pixels (z component already
    % correct
    for i = 1:length(tracks)
        tracks(i).pathx_um = tracks(i).pathx/trackParams.scale;
        tracks(i).pathy_um = tracks(i).pathy/trackParams.scale;
        tracks(i).x_smooth=smooth(tracks(i).pathx, velocitySmoothWin)';
        tracks(i).y_smooth=smooth(tracks(i).pathy, velocitySmoothWin)';
    end
    
    %%
    %-------------------------------------------------------------------------------
    % VELOCITY ANALYSIS
    v_mean2D = zeros(length(tracks), 1);
    
    for i = 1:length(tracks)
        x = smooth(tracks(i).pathx_um, velocitySmoothWin);
        y = smooth(tracks(i).pathy_um, velocitySmoothWin);
        
        velocity2D = sqrt(diff(x).^2 + diff(y).^2)*frameRate;
        
        tracks(i).velocity2D = velocity2D;
        tracks(i).corrVelocity2D = ...
            velocity2D(velocitySmoothWin:(length(velocity2D) - velocitySmoothWin));
        tracks(i).corrLen = length(tracks(i).corrVelocity2D);
        
        v_mean2D(i) = mean(velocity2D);
        tracks(i).v_mean2D=v_mean2D(i);
    end
    
    %% Heading
    for i = 1:length(tracks)
        x = smooth(tracks(i).pathx_um, headingSmoothWin);
        y = smooth(tracks(i).pathy_um, headingSmoothWin);
        dx = x((headingSmoothWin+1):length(x)) - x(1:(length(x)-headingSmoothWin));
        dy = y((headingSmoothWin+1):length(y)) - y(1:(length(y)-headingSmoothWin));
        
        heading2D = atan2(dy, dx);
        
        % Correct for turns around 360
        headOffset = 0;
        for j = 1:(length(heading2D) - 1)
            heading2D(j) = heading2D(j) + headOffset ;
            if (heading2D(j+1) + headOffset - heading2D(j)) > 1.5*pi
                headOffset = headOffset - 2*pi;
            end
            if (heading2D(j+1) + headOffset - heading2D(j)) < -1.5*pi
                headOffset = headOffset + 2*pi;
            end
        end
        
        tracks(i).heading2D = heading2D;
    end
    
    %%
    %-------------------------------------------------------------------------------
    
    tracksTumb = tracks(v_mean2D>0); % Discard tracks with very low average speeds
    for i=1:length(tracksTumb)
       tracksTumb(i).x_smooth=smooth(tracksTumb(i).pathx, velocitySmoothWin)';
       tracksTumb(i).y_smooth=smooth(tracksTumb(i).pathy, velocitySmoothWin)';
    end
    
    if (isempty(tracksTumb) == 0)
        v_mean=mean([tracksTumb.v_mean2D]); % mean 2D velocity of all trajectories
    end
    
    tracks = tracks([tracks.v_mean2D] > trackParams.minSpeed);
    smoothers = ones(1, length(tracks));
    for i = 1 : length(tracks)
        if (std(tracks(i).velocity2D) > 8)
                smoothers(i) = 0;
        end
    end
        tracks = tracks(smoothers == 1);
    
    %-------------------------------------------------------------------------------
    % FIGURES and RESULTs
    resultsDirectory=[resultsFolder, resultsName];
    mkdir(resultsDirectory);
    
    results = struct('numTraj', length(tracks), 'numFiltTraj', length(tracksTumb),...
        'v_mean2D', mean(v_mean2D(v_mean2D > 1)),...
        'v_std2D', std(v_mean2D(v_mean2D > 1)));
    
    %%  
    save([resultsDirectory, '/resultData'], 'tracks', 'trackParams','anlysParams');
    
    %%
    %-------------------------------------------------------------------------------
    % VIDEO RESULTS
    
    % Write video file showing trajectories in tracksTumb
    if flags.writeVideo == 1
        
        % Create video object and format figure
        writerObj = VideoWriter([resultsDirectory, '/trackingVideo'], 'Motion JPEG AVI');
        writerObj.FrameRate = frameRate;
        open(writerObj);
        frame = figure(1);
        %set(frame, 'OuterPosition', [100, 100, size(data, 2), size(data, 1)])
        
        set(gcf, 'Color', [1,1,1]);
        set(gca,'position',[0 0 1 1],'units','normalized')
        axis equal
        for i = 1:vidLength
            figure(1)
            imagesc(double(data(:,:,i)),[0,255]);
            colormap gray
            hold on
            daspect([1, 1, 1]);
            axis off
            
            % Plot calculated positions and trajectories
            for j = 1:length(tracks)
                if (i >= tracks(j).start) &&...
                        (i < tracks(j).start + tracks(j).len)
                    
                    x = tracks(j).x_smooth(1:(i - tracks(j).start + 1));
                    y = tracks(j).y_smooth(1:(i - tracks(j).start + 1));
                    
                    plot(x, y, 'r-')
                    
                    plot(tracks(j).x_smooth(i - tracks(j).start + 1),...
                        tracks(j).y_smooth(i - tracks(j).start + 1), 'go')
                    text(tracks(j).x_smooth(i - tracks(j).start + 1) + 2,...
                        tracks(j).y_smooth(i - tracks(j).start + 1), ...
                        int2str(j),'color','w')
                end
            end
            text(920, 20, int2str(i),'color','w');
            % Copy frame to video writer
            f1 = getframe(frame);
            writeVideo(writerObj, f1);
            hold off
        end
        close(writerObj)
        close(figure(1))
    end
    toc
    disp(sprintf('End file %d.', fileNum));
end

% The end of the code