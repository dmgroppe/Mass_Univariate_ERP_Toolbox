% spatial_nearest_neighbors() - Given a set of electrode coordinates, this 
%                       function returns a binary 2D matrix that indicates  
%                       which electrodes are close enough to be considered 
%                       neighbors. Specifically, the k nearest electrodes 
%                       to each electrode are defined as neighbors. For use 
%                       with cluster-based permutation tests. This method 
%                       was suggested by Manish Saggar is useful when your
%                       electrode coordinates are different for each
%                       subject or if you want to have relatively equal
%                       statistical power across all electrodes. Defining
%                       spatial neighbors by a simple distance threshold 
%                       (using spatial_neighbors.m) means that electrodes 
%                       on the edge of the montage will have fewer 
%                       neighbors than others and will thus have less 
%                       opportunities to join clusters.
%
% Usage:
%  >>chan_hood=spatial_nearest_neighbors(chanlocs,k);
%
% Required Inputs:
%   chanlocs - An EEGLAB chanlocs structure (e.g., EEG.chanlocs from an EEG
%              variable)
%   k        - The number of electrodes that will be defined as a neighbor
%              of each electrode. The nearest k electrodes are chosen as
%              neighbors. In case of ties, the lowest indexed channels are 
%              chosen over higher indexed channels. If this happens, you
%              should probably use a different value of k or manually
%              adjust chan_hood to avoid the issue.
%
% Outputs:
%   chan_hood - A symmetric binary matrix indicating which channels are
%               neighbors. If chan_hood(a,b)=1, then Channel A and Channel
%               B are nieghbors.
%
% Notes:
% -See also the function spatial_neighbors.m
% -This function won't work if it produces an asymmetric neighborhood
% index, which appears to be hard to avoid.
%
% Authors:
% Idea from Manish Saggar
% Written by David Groppe 
% Feinstein Institute for Medical Research, 3/2013

%%%%% Future Work %%%%%%
% figure out how to deal with asymmetric 

function chan_hood=spatial_nearest_neighbors(chanlocs,k)

n_chan=length(chanlocs);

chan_hood=zeros(n_chan,n_chan);
n_neighbors=zeros(1,n_chan);
chan_dist=zeros(1,n_chan*(n_chan-1)/2);
ct=0;
for c=1:n_chan,
    coordA=[chanlocs(c).X chanlocs(c).Y chanlocs(c).Z];
    dists=zeros(1,n_chan);
    for d=1:n_chan,
        coordB=[chanlocs(d).X chanlocs(d).Y chanlocs(d).Z];
        dists(d)=sum((coordA-coordB).^2); % Ok to keep distance in units squared
    end
    [vals ids]=sort(dists);
    chan_hood(c,ids(1:k+1))=1; %k+1 is used because each channel should be a neighbor with itself
end

% Check for symmetry
is_sym=1;
for a=1:n_chan,
    for b=a+1:n_chan,
        if chan_hood(a,b)~=chan_hood(b,a)
           is_sym=0;
           break;
        end
    end
end

if ~is_sym,
    error('Spatial neighborhood matrix is not symmetric. You need to change the value of k, define chan_hood manually, or use spatial_neighbors.m');
end
