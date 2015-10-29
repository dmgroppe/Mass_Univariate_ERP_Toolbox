% spatial_neighbors() - Given a set of electrode coordinates, this function 
%                       returns a binary 2D matrix that indicates which 
%                       electrodes are close enough to be considered 
%                       neighbors. For use with cluster-based permutation 
%                       tests.
%
% Usage:
%  >>chan_hood=spatial_neighbors(chanlocs,max_dist);
%
% Required Inputs:
%   chanlocs - An EEGLAB chanlocs structure (e.g., EEG.chanlocs from an EEG
%              variable)
%   max_dist - All electrodes within max_dist of another electrode are
%              considered spatial neighbors.  Max_dist is in whatever units
%              your EEGLAB chanlocs coordinates are in.  If your chanlocs 
%              coordingates are on an idealized sphere with unit radius 
%              then you can convert max_dist into centimeters by measuring 
%              the circumference of a participant's head in centimeters and
%              using the following formulas:
%                 radius=circumference/(2*pi);
%                 radius*max_dist=max_dist in units of cm
%
% Optional Inputs:
%   head_radius - The radius of the head in whatever units the Cartesian
%                 coordinates in chanlocs are in. This is used to
%                 convert scalar values of chan_hood into centimeters.
%                 {default: estimated from chanlocs by assuming center of
%                 head is at 0,0,0}
%
% Outputs:
%   chan_hood - A symmetric binary matrix indicating which channels are
%               neighbors. If chan_hood(a,b)=1, then Channel A and Channel
%               B are nieghbors.
%
% Notes:
% -This function outputs an estimate of max_dist (in cm) on the command line
% by assuming a circumference of 56 cm and EEGLAB chanloc coordinates based 
% on a spherical head with unit radius.  It also summarizes the number of 
% neighbors per channel using basic measures of central tendency and 
% dispersion.
% -You will have more statistical power to detect effects at electrodes
% that have more neighbors.  Thus you may have significantly less power to
% detect effects on the edge of the montage.  Thanks to Manish Saggar for
% bringing this to my attention.
%
% Author:
% David Groppe
% Kutaslab, 5/2011

%%%%% Future Work %%%%%%
% -Perhaps allow max_dist to be specified in centimeters

function chan_hood=spatial_neighbors(chanlocs,max_dist,head_radius)

n_chan=length(chanlocs);

if nargin<3,
    head_radius=[];
end

if isempty(head_radius)
    fprintf('Estimating the radius of the head by assuming that the center of the head is [0,0,0] in Cartesian coordinates.\n');
    dst_from_origin=zeros(1,n_chan);
    for a=1:n_chan,
        dst_from_origin(a)=sqrt(sum([chanlocs(a).X chanlocs(a).Y chanlocs(a).Z].^2));
    end
    mn=min(dst_from_origin);
    mx=max(dst_from_origin);
    md=median(dst_from_origin);
    fprintf('Min/Max electrode distance from origin (in chanlocs units): %f/%f\n',mn,mx);
    uni=unique(dst_from_origin);
    if ((mx-mn)/md)>.001,
        fprintf('WARNING: It appears that all electrodes are NOT the same distance from the origin!!!\n');
        fprintf('Your electrodes'' Cartesian coordinates either are not spherical or are not centered on [0 0 0].\n');
        fprintf('This function will still work properly but its estimate of the radius of the head of will not be correct.\n');
    end
    head_radius=median(dst_from_origin);
    fprintf('Radius of head (in chanlocs units) is estimated to be %f\n',head_radius);
else
    fprintf('Using provided head radius of %f (in chanlocs units)\n',head_radius);
end

circumference=56; %very rough estimate based on the average circumference of 10 Kutaslab participants
max_dist_cm=max_dist*circumference/(2*pi*head_radius); % Radius=Circumference/(2*pi)

fprintf('max_dist value of %g corresponds to an approximate distance of %.2f cm (assuming\n',max_dist,max_dist_cm);
fprintf('  a 56 cm great circle circumference head and that your electrode coordinates are based on an idealized\n');
fprintf('  spherical head with radius of %f).\n',head_radius);
    
chan_hood=zeros(n_chan,n_chan);
n_neighbors=zeros(1,n_chan);
chan_dist=zeros(1,n_chan*(n_chan-1)/2);
ct=0;
for c=1:n_chan,
    coordA=[chanlocs(c).X chanlocs(c).Y chanlocs(c).Z];
    for d=c:n_chan,
        coordB=[chanlocs(d).X chanlocs(d).Y chanlocs(d).Z];
        dstnce=sqrt(sum((coordA-coordB).^2));
        if dstnce<=max_dist,
            chan_hood(c,d)=1;
            chan_hood(d,c)=1;
        end
        
        if c~=d
            %don't count channels with themselves
            ct=ct+1;
            chan_dist(ct)=dstnce;
        end
    end
    n_neighbors(c)=sum(chan_hood(c,:))-1;
end

fprintf('Min/Max distances between all pairs of channels (in chanlocs units): %f/%f\n', ...
    min(chan_dist),max(chan_dist));
fprintf('Median (semi-IQR) distance between all pairs of channels (in chanlocs units): %f (%f)\n', ...
    median(chan_dist),iqr(chan_dist)/2);
fprintf('Mean (SD) # of neighbors per channel: %.1f (%.1f)\n',mean(n_neighbors), ...
    std(n_neighbors));
fprintf('Median (semi-IQR) # of neighbors per channel: %.1f (%.1f)\n',median(n_neighbors), ...
    iqr(n_neighbors)/2);
fprintf('Min/max # of neighbors per channel: %d to %d\n',min(n_neighbors), ...
    max(n_neighbors));