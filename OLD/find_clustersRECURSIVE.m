% find_clusters() - Given a 2D matrix of t-scores, this function bundles
%                   neighboring above threshold t-scores into clusters 
%                   and returns a 2D matrix of cluster assignments. For
%                   use with cluster-based permutation tests.
%             
% Usage:
%  >>clust_membership=find_clusters(tscores,thresh,chan_hood,thresh_sign,recursion_limit);
%
% Required Inputs:
%   tscores     - A 2D matrix of t-scores (channel x time point). Note,
%                 t-scores SHOULD be signed so that negative and positive
%                 deviations from the null hypothesis are NOT clustered
%                 together.
%   thresh      - The thershold for cluster inclusion (i.e., only t-scores
%                 more extreme than the threshold will be included in clusters). 
%                 If thresh is positive, clusters of positive t-scores will 
%                 be returned.  Otherwise clusters of negative t-scores will 
%                 be returned.
%   chan_hood   - A symmetric 2d matrix indicating which channels are 
%                 neighbors with other channels.  If chan_hood(a,b)=1,
%                 then Channel A and B are neighbors. This is produced by the
%                 function spatial_neighbors.m.
%   thresh_sign - If greater than zero, t-scores greater than thresh will be
%                 included in clusters.  Otherwise, t-scores less than thresh
%                 will be included in clusters.
%
% Optional Inputs:
%   recursion_limit - The number of recursive embeddings you want MATLAB to
%                     allow. MATLAB limits this and, since clusters are
%                     found using recursion, it is possible that you will
%                     exceed MATLAB's limit with large clusters. The
%                     default behavior is to estimate the largest number of
%                     embeddings, set this as the limit, and then re-set
%                     the original limit after the clustering is done.
%                     However, if this doesn't work, use a higher number
%                     (MATLAB's default is 500 in our lab).
%
% Outputs:
%   clust_membership - A 2D matrix (channel x time point) indicating which
%                      cluster each channel/time point pair belongs to. For 
%                      example, if clust_membership(2,10)=3, then the t-score 
%                      at the 2nd channel and 10th time point belongs to the 
%                      3rd cluster.  0 indicates that the t-score was not 
%                      included in any cluster.
%   n_clust          - The number of clusters found.
%   used_rec_lim     - The recursion limit used by the function. This may
%                      differ from the specified recursion limit since this 
%                      function will automatically double the recursion
%                      limit whenever it is too small.  This output is
%                      useful when you need to re-run this function on a
%                      similarly sized data set (e.g., when one is
%                      performing a cluster-based permutation test) since
%                      having to automatically increase the recursion limit
%                      wastes time.
%
% Notes:
% -A global variable, clust_ids, is created to keep track of which cluster
% each channel/time point pair belongs to.  This is done to minimize memory
% demands during the recusion. It is erased after all the clusters are
% established.
%
% Author:
% David Groppe
% Kutaslab, 5/2011

%%%%%%%%%%%%%%%% REVISION LOG %%%%%%%%%%%%%%%%%
% ?/?/??-

function [clust_membership, n_clust, used_rec_lim]=find_clusters(tscores,thresh,chan_hood,thresh_sign,recursion_limit)

[n_chan, n_tpt]=size(tscores);

global clust_ids; % channel x time point matrix indicated which cluster each "voxel" belongs to

clust_ids=zeros(n_chan,n_tpt);

time_matrix=repmat(1:n_tpt,n_chan,1);
chan_matrix=repmat([1:n_chan]',1,n_tpt);

if thresh_sign>0
    %looking for positive clusters
    above_thresh_ids=find(tscores>=thresh);
else
    %looking for negative clusters
    above_thresh_ids=find(tscores<=thresh);
end
above_thresh_times=time_matrix(above_thresh_ids);
above_thresh_chans=chan_matrix(above_thresh_ids); 

%Clear chan & time matrix to save memory (just in case)
clear chan_matrix time_matrix

n_above=length(above_thresh_ids);

orig_limit=get(0,'RecursionLimit');
if nargin<5 || isempty(recursion_limit),
    if orig_limit<n_above,
        set(0,'RecursionLimit',n_above);
    end
else
    set(0,'RecursionLimit',recursion_limit);
end

n_clust=0;
for a=1:n_above,
    voxel_id=above_thresh_ids(a);
    if ~clust_ids(voxel_id)
        %this "voxel" isn't in a cluster yet
        n_clust=n_clust+1;
        clust_ids(voxel_id)=n_clust;
        
        %go through all the remaining voxels and find all the above
        %threshold voxels this voxel is neighbors with and give them this
        %cluster #
        did_it=0;
        clust_ids_copy=clust_ids;
        while ~did_it,
            try
                follow_clust(a+1,n_above,a,n_clust,above_thresh_ids,above_thresh_times,above_thresh_chans,chan_hood,n_chan);
                did_it=1;
            catch me
                did_it=0;
                clust_ids=clust_ids_copy; %reset clust_ids to what it was before the attempt
                if strcmpi(me.identifier,'MATLAB:recursionLimit')
                    current_limit=get(0,'RecursionLimit');
                    %fprintf('Doubling RecursionLimit to %d\n',current_limit*2);
                    set(0,'RecursionLimit',current_limit*2);
                end
            end
        end
    end
end

%reset RecursionLimit to its original value
used_rec_lim=get(0,'RecursionLimit');
set(0,'RecursionLimit',orig_limit);
clust_membership=clust_ids;
clear global clust_ids

%% End of Main Function

function follow_clust(start_voxel,n_above,current_voxel_id,current_clust_num,above_thresh_ids,above_thresh_times,above_thresh_chans,chan_hood,n_chan)
%Function for recursively finding all the members of cluster based on
%single "voxel" seed
%
% start_voxel        - All voxels above this might belong to the cluster
%                      (i.e., all voxels before this have already been 
%                      assigned to clusters)
% n_above            - The total number of above threshold voxels
% current_voxel_id   - The index of the seed voxel into
%                      above_thresh_ids, above_thresh_chans, above_thresh_times
% current_clust_num  - The # of the cluster the seed voxel belongs to
% above_thresh_ids   - A vector of indices of all the above threshold
%                      voxels into the original 2d data matrix (channel x
%                      time point)
% above_thresh_times - The time points corresponding to above_thresh_ids
% above_thresh_chans - The time points corresponding to above_thresh_chans
% chan_hood          - A symmetric 2d matrix indicating which channels are 
%                      neighbors with other channels.  If chan_hood(a,b)=1,
%                      then Channel A and B are neighbors.
% n_chan             - The total number of channels
%

global clust_ids; % channel x time point matrix indicated which cluster each "voxel" belongs to

new_members=zeros(1,n_chan*3); %pre-allocate memory
new_members_ct=0;
for b=start_voxel:n_above,
    if ~clust_ids(above_thresh_ids(b))
        temp_dist=abs(above_thresh_times(b)-above_thresh_times(current_voxel_id));
        if above_thresh_chans(current_voxel_id)==above_thresh_chans(b)
            %voxels are at same channel
            chan_dist=0;
        elseif chan_hood(above_thresh_chans(current_voxel_id),above_thresh_chans(b))
            %channels are neighbors
            chan_dist=1;
        else
            %voxels aren't spatially compatible
            chan_dist=2;
        end
        if (temp_dist+chan_dist)<=1,
            %if voxels are at same time point and neighboring channels OR
            %if voxels are at same channel and neighboring time points,
            %merge them into the same cluster
            
            clust_ids(above_thresh_ids(b))=current_clust_num;
            %keep track of which other voxels are joined to this
            %cluster
            new_members_ct=new_members_ct+1;
            new_members(new_members_ct)=b;
        end
    end
end

%recursively search for above threshold neighbors of the other
%voxels that were just joined to this cluster
for c=1:new_members_ct,
    follow_clust(start_voxel,n_above,new_members(c),current_clust_num,above_thresh_ids,above_thresh_times,above_thresh_chans,chan_hood,n_chan);
end

