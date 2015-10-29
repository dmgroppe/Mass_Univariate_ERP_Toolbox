% find_clusters() - Given a 2D matrix of t-scores, this function bundles
%                   neighboring above threshold t-scores into clusters
%                   and returns a 2D matrix of cluster assignments. For
%                   use with cluster-based permutation tests.
%
% Usage:
%  >>clust_membership=find_clusters(tscores,thresh,chan_hood,thresh_sign);
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
%
% Outputs:
%   clust_membership - A 2D matrix (channel x time point) indicating which
%                      cluster each channel/time point pair belongs to. For
%                      example, if clust_membership(2,10)=3, then the t-score
%                      at the 2nd channel and 10th time point belongs to the
%                      3rd cluster.  0 indicates that the t-score was not
%                      included in any cluster.
%   n_clust          - The number of clusters found.
%
%
% Notes:
% -A global variable, clust_ids, is created to keep track of which cluster
% each channel/time point pair belongs to.  This speeds things up because
% of the subfunction this variable has to be passed to.
%
% Authors:
%  Amy Guthormsen
%  Applied Modern Physics Group
%  Los Alamos National Laboratory
%
%  David M. Groppe
%  Kutaslab
%  University of California, San Diego

%%%%%%%%%%%%%%%% REVISION LOG %%%%%%%%%%%%%%%%%
% 12/11/2011-Made reursiveless by Amy Guthormsen

function [clust_membership, n_clust]=find_clusters(tscores,thresh,chan_hood,thresh_sign)

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
        voxels_not_checked = ones(size(above_thresh_ids));
        check_me = zeros(size(above_thresh_ids));
        
        check_me(a) = 1;
        
        while sum(check_me>0)
            first = find(check_me,1);
            new = follow_clust(n_above,first,n_clust,above_thresh_ids,above_thresh_times,above_thresh_chans,chan_hood,n_chan);
            check_me(new)=1;
            voxels_not_checked(first)=0;
            check_me = check_me&voxels_not_checked;
        end
    end
end

clust_membership=clust_ids;
clear global clust_ids

%% End of Main Function
    function [new_members] = follow_clust(n_above,current_voxel_id,current_clust_num,above_thresh_ids,above_thresh_times,above_thresh_chans,chan_hood,n_chan)
        %Function for finding all the members of cluster based on
        %single "voxel" seed. If it finds new members of a cluster, it indicates
        %which cluster they are a member of in the global variable clust_ids.
        %
        % Inputs:
        %   clust_ids          - Matrix indicating which cluster each
        %                        voxel belongs to (if any)
        %   n_above            - The total number of above threshold voxels
        %   current_voxel_id   - The index of the seed voxel into
        %                      above_thresh_ids, above_thresh_chans, above_thresh_times
        %   current_clust_num  - The # of the cluster the seed voxel belongs to
        %   above_thresh_ids   - A vector of indices of all the above threshold
        %                      voxels into the original 2d data matrix (channel x
        %                      time point)
        %   above_thresh_times - The time points corresponding to above_thresh_ids
        %   above_thresh_chans - The time points corresponding to above_thresh_chans
        %   chan_hood          - A symmetric 2d matrix indicating which channels are
        %                      neighbors with other channels.  If chan_hood(a,b)=1,
        %                      then Channel A and B are neighbors.
        %   n_chan             - The total number of channels
        %
        %
        % Output:
        %  new_members - above_thresh_ids indices to new members of the cluster
        
        %   global clust_ids; % channel x time point matrix indicated which cluster each "voxel" belongs to
        
        new_members=zeros(1,n_chan*3); %pre-allocate memory
        new_members_ct=0;
        
        for b=current_clust_num:n_above,
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
        
        new_members=nonzeros(new_members);
    end
end