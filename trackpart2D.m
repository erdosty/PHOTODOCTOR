 function paths = trackpart2D(data, vidLength, back, trackParams)
% trackparticle2D, tracks mobile particles in 2D

pCounter = 0; % Highest numbered path in paths
paths = {};

for i = 1:vidLength
    
    frame = double(data(:,:,i)) - back;
    % Take the median in space
    frame = abs(frame);
    frame(frame < trackParams.threshold) = 0;
    framebw = im2bw(frame);
%     framebw=bwmorph(framebw,'erode');
    framebw = bwareaopen(framebw,trackParams.minArea,4);
    framebw = bwmorph(framebw,'dilate',1);
    % get particle clusters info
    particle = regionprops(framebw);
    particle = particle([particle.Area] < trackParams.maxArea);
    
    if pCounter == 0
        j = 1;
    else
        for j = 1:pCounter
            if paths{j}.active == 1
                x = paths{j}.pathx(length(paths{j}.pathx));
                y = paths{j}.pathy(length(paths{j}.pathy));
                foundNext = 0;
                
                mDist = trackParams.maxDist;
                for k = 1:length(particle)
                    dist = sqrt((particle(k).Centroid(1) - x)^2 + ...
                        (particle(k).Centroid(2) - y)^2);
                    if dist < mDist
                        mDist = dist;
                        foundNext = 1;
                        index = k;
                    end
                end
                
                if foundNext == 0
                    paths{j}.active = 0;   
                else
                    paths{j}.pathx = [paths{j}.pathx, particle(index).Centroid(1)];
                    paths{j}.pathy = [paths{j}.pathy, particle(index).Centroid(2)];
                    paths{j}.len = paths{j}.len + 1;
                    paths{j}.area = [paths{j}.area, particle(index).Area];
                    particle = [particle(1:(index-1)); particle((index+1):length(particle))];
                end
            end
        end
    end
    for m = 1:length(particle)
        pCounter = pCounter + 1;
        paths{pCounter}  = struct('active', 1, 'start', i,...
            'len', 1,...
            'pathx', [particle(m).Centroid(1)], 'pathy', [particle(m).Centroid(2)],...
            'area', [particle(m).Area], 'collision', 0, 'v_mean2D', 0);
    end  
end
end

