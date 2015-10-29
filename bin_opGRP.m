% bin_opGRP() - Computes the difference waves from two bins or the mean of 
%               two bins of a GRP variable and stores the resultant ERPs as 
%               a new bin in the GRP variable.
%             
% Usage:
%  >> GRP=bin_opGRP(GRP,op,groupA,binA,groupB,binB,neo_bindesc,neo_ccode,neo_ccdescr)
%
% Required Inputs:
%   GRP    - A GRP structure variable.  You can create a GRP variable 
%            from files of GND variables (i.e., *.GND files)   
%            using GNDs2GRP.m.  See Mass Univariate ERP Toolbox documentation
%            for detailed information about the format of a GRP variable. 
%   op     - [string] One of the following three "bin operators":
%              'A-B'    = new bin will be the difference wave produced by
%                         subtracting the bin from Group B from the bin from 
%                         Group A. Grand average t-scores are based on  
%                         independent samples t-statistics.
%              'B-A'    = new bin will be the difference wave produced by
%                         subtracting the bin from Group A from the bin from 
%                         Group B. Grand average t-scores are based on two 
%                         independent samples t-statistics.
%              '(A+B)/n'= new bin will be the mean of the bins from Group
%                         A and B.  The means will be weighted by the
%                         number of participants in each bin.  Thus the
%                         results are equivalent to what you would derive
%                         if all participants were in the same group.
%                         Grand average t-scores are based on one sample
%                         t-statistics.
%   groupA - [string] The group descriptor for Group A (e.g., 'patients') 
%            Use the function headinfo.m to see what groups are available.
%   binA   - [integer] A bin index from Group A's GND variable. Note this is
%            NOT a bin in the GRP variable. Use the function headinfo.m on 
%            Group A's GND variable to see what bins are available
%   groupB - [string] Same as groupA but for Group B.
%   binB   - [integer] Same as binA but for Group B.
%
% Optional Inputs:
%   neo_bindesc - [string] The bin descriptor for the new bin being
%                 created. {default: '$, Bin #-$$ Bin ##', where # is the 
%                 value of binA, ## is the value of binB, $ is Group A's
%                 descriptor, and $$ is Group B's descriptor}
%   neo_ccode   - [integer] The condition code of the new bin being
%                 created. {default: biggest condition code}
%   neo_ccdescr - [string] The condition code descriptor for the new bin being
%                 created. {default: existing descriptor for that condition
%                 code or 'Not Cals' if no condition codes exist for this 
%                 GRP variable}
%
% Notes:
% -binA and binB should index bins starting at 1 (i.e., Bin 0 is assumed to
% contain cal pulses in Kutaslab data and is stored apart from bins in GND 
% variables)
%
% Author:
% David Groppe
% Kutaslab, 3/2010

%%%%%%%%%%%%%%%% REVISION LOG %%%%%%%%%%%%%%%%%
% 5/6/2010-Fixed bug with default bin descriptor when 'op'='(A+B)/n' and
% made combining bins a weighted average (weighted by the number of
% participants in each bin).  Function was simple taking the mean of both
% bins (treating each group equally)
%
% 11/1/2010-Fixed bug when GRP variable initially has no bins.
%
% 5/24/2011-Function wouldn't work if the Group A wasn't the first group in
% GRP variable.  That problem is fixed.

%%%%%%%%%%%%%%%% FUTURE WORK %%%%%%%%%%%%%%%%%
% - add verblevel? 


function GRP=bin_opGRP(GRP,op,groupA,binA,groupB,binB,neo_bindesc,neo_ccode,neo_ccdescr)

if nargin<7,
   neo_bindesc=[]; 
end

if nargin<8,
   neo_ccode=[]; 
end

if nargin<9,
   neo_ccdescr=[]; 
end

if strcmpi(groupA,groupB),
   error('groupA and groupB are the same (%s).  For within-subject difference waves use bin_dif.m.',groupA); 
end

if strcmpi(op,'A-B') || strcmpi(op,'B-A'),
    op_sign='-';
elseif strcmpi(op,'(A+B)/n')
    op_sign='+';
else
    error('Argument "op" must be ''A-B'',''B-A'',or ''(A+B)/n''.');
end

n_chans=length(GRP.chanlocs);
n_pts=length(GRP.time_pts);
n_binsGRP=length(GRP.bin_info);
neo_bin=n_binsGRP+1;
n_groups=length(GRP.group_desc);

%Get the GRP variable indices of Group A and Group B 
grpA_id=0;
grpB_id=0;
for a=1:n_groups,
    if strcmpi(GRP.group_desc{a},groupA),
        if grpA_id,
            error('groupA (%s) matches multiple group descriptors in this GRP variable.',groupA);
        else
            grpA_id=a;
        end
    elseif strcmpi(GRP.group_desc{a},groupB),
        if grpB_id,
            error('groupB (%s) matches multiple group descriptors in this GRP variable.',groupB);
        else
            grpB_id=a;
        end
    end
end
if ~grpA_id,
    error('groupA (%s) not found in this GRP variable.  Enter "headinfo(GRP)" to see possible names.',groupA);
end
if ~grpB_id,
    error('groupB (%s) not found in this GRP variable.  Enter "headinfo(GRP)" to see possible names.',groupB);
end

load(GRP.GND_fnames{grpA_id},'-MAT');
n_subA=GND.sub_ct;
%n_subA=sum(GRP.indiv_bin_ct{grpA_id}(:,binA)>0);
%n_subB=sum(GRP.indiv_bin_ct{grpB_id}(:,binB)>0);
load(GRP.GND_fnames{grpB_id},'-MAT');
n_subB=GND.sub_ct;
if isempty(neo_bindesc),
    if op_sign=='-',
        GRP.bin_info(neo_bin).bindesc=sprintf('%s, Bin %d%c%s, Bin %d',groupA,binA, ...
            op_sign,groupB,binB);
    else
        GRP.bin_info(neo_bin).bindesc=sprintf('((%d*%s, Bin %d)%c(%d*%s, Bin %d))/%d', ...
            n_subA,groupA,binA,op_sign,n_subB,groupB,binB,n_subA+n_subB);
    end
else
    GRP.bin_info(neo_bin).bindesc=neo_bindesc;
end
if isempty(neo_ccode)
    if isempty(GRP.condesc),
        fprintf('No condition code specified for new bin and no condition codes exist in GRP variable.\n');
        fprintf('Creating new condition code 1: "Not Cals" for the new bin being created.\n');
        GRP.bin_info(neo_bin).condcode=1;
        GRP.condesc{1}='Not Cals';
    else
        neo_ccode=length(GRP.condesc);
        fprintf('No condition code specified for new bin.\n');
        fprintf('Using greatest condition code (%d: %s) for the new bin being created by default.\n', ...
            neo_ccode,GRP.condesc{neo_ccode});
        GRP.bin_info(neo_bin).condcode=neo_ccode;
    end
else
    %Check to see if the new condition code greater than the current possible
    cur_ccodes=length(GRP.condesc);
    if neo_ccode>cur_ccodes,
        if isempty(neo_ccdescr),
            error('You specified a condition code of %d, but this GRP variable only has %d existing conditions and you did not specify a new condition code descriptor.', ...
                neo_ccode,cur_ccodes);
        end
        GRP.bin_info(neo_bin).condcode=neo_ccode;
        GRP.condesc{neo_ccode}=neo_ccdescr;
    else
        GRP.bin_info(neo_bin).condcode=neo_ccode;
    end
end
GRP.bin_info(neo_bin).groupA=grpA_id;
GRP.bin_info(neo_bin).groupB=grpB_id;
GRP.bin_info(neo_bin).source_binA=binA;
GRP.bin_info(neo_bin).source_binB=binB;


%preallocate memory
if strcmpi(op,'(A+B)/n'),
    all_subs=zeros(n_chans,n_pts,n_subA+n_subB); %third dimension is participant
else
    sums=zeros(n_chans,n_pts,2); %third dimension is group (1=A, 2=B)
    means=sums;
    ss=sums;
end
bins=[binA binB];
n_subs=[0 0];

%get data from Group A & B
for grp_ct=1:2,
    % Group count of 1=Group A, 2=Group B
    if grp_ct==1,
        grp_id=grpA_id;
    else
        grp_id=grpB_id;
    end
    
    load(GRP.GND_fnames{grp_id},'-MAT');
    if bins(grp_ct)>length(GND.bin_info)
        error('The GND variable from %s only has %d bins but you''re trying to import Bin %d.', ...
            GRP.GND_fnames{grp_id},length(GND.bin_info),bins(grp_ct));
    end
    
    %Get the set of channels to use
    use_chans=zeros(1,length(GND.chanlocs));
    for c=1:length(GND.chanlocs),
        for d=1:n_chans, %n_chans is the number of channels in the GRP variable
           if strcmpi(GND.chanlocs(c).labels,GRP.chanlocs(d).labels),
              use_chans(c)=1;
              break;
           end
        end
    end
    %Make sure channel information is compatible
    use_chans=find(use_chans);
    if ~isequal(GND.chanlocs(use_chans),GRP.chanlocs),
        error('The channel location information in the GND variable from file %s differs from that in the GRP variable.', ...
            GRP.GND_fnames{grp_id});
    end
    
    %Get the set of subjects to use
    use_subs=find(GND.indiv_bin_ct(:,bins(grp_ct)))';
    n_subs(grp_ct)=length(use_subs);
    if grp_ct==1,
        % # of participants who contributed to this bin from Group A
        GRP.bin_info(neo_bin).n_subsA=n_subs(grp_ct); 
    else
        % # of participants who contributed to this bin from Group B
        GRP.bin_info(neo_bin).n_subsB=n_subs(grp_ct);
    end
    
    %get individual participant trial counts
    indiv_bin_ct{grp_ct}=GND.indiv_bin_ct(:,bins(grp_ct));
    indiv_bin_raw_ct{grp_ct}=GND.indiv_bin_raw_ct(:,bins(grp_ct));
    
    data=squeeze(GND.indiv_erps(use_chans,:,bins(grp_ct),use_subs));

    %make sure squeezed data has proper dimensions
    if length(size(data))~=3,
       error('Contributing ERPs from file %s have only one participant, time point, and/or channel.  I can''t deal with that.', ...
           GRP.GND_fnames{grp_id});
    end

    if strcmpi(op,'(A+B)/n'),
         if grp_ct==1,
             all_subs(:,:,1:n_subs(grp_ct))=data;
         else
             all_subs(:,:,n_subs(1)+1:end)=data;
         end
    else
        sums(:,:,grp_ct)=sum(data,3);
        means(:,:,grp_ct)=sums(:,:,grp_ct)/n_subs(grp_ct);
        ss(:,:,grp_ct)=sum(data.^2,3)-(squeeze(sums(:,:,grp_ct)).^2)/n_subs(grp_ct);
    end
end
GRP.bin_info(neo_bin).op=op;

for a=1:n_groups,
    n_subs_grp=length(GRP.indiv_subnames{a});
    if a==grpA_id,
        GRP.indiv_bin_ct{a}(1:n_subs_grp,neo_bin)=indiv_bin_ct{1};
        GRP.indiv_bin_raw_ct{a}(1:n_subs_grp,neo_bin)=indiv_bin_raw_ct{1};
    elseif a==grpB_id,
        GRP.indiv_bin_ct{a}(1:n_subs_grp,neo_bin)=indiv_bin_ct{2};
        GRP.indiv_bin_raw_ct{a}(1:n_subs_grp,neo_bin)=indiv_bin_raw_ct{2};
    else
        %just zeros for this groups since they don't contribute to the new
        %bin
        GRP.indiv_bin_ct{a}(1:n_subs_grp,neo_bin)=zeros(1:n_subs,1);
        GRP.indiv_bin_raw_ct{a}(1:n_subs_grp,neo_bin)=zeros(1:n_subs,1);
    end
end

%new ERPs
if strcmpi(op,'A-B'),
    GRP.grands(:,:,neo_bin)=means(:,:,1)-means(:,:,2);
elseif strcmpi(op,'B-A'),
    GRP.grands(:,:,neo_bin)=means(:,:,2)-means(:,:,1);
else
    % weighted mean of both bins (the right thing to do for a one-sample
    % t-test)
    [p_values, t_scores, mn, stder]=fast_t1(all_subs,0,0);
    GRP.grands(:,:,neo_bin)=mn;
    GRP.grands_t(:,:,neo_bin)=t_scores;
    GRP.grands_stder(:,:,neo_bin)=stder;
end
    
if ~strcmpi(op,'(A+B)/n'),
    %pooled sum of squares
    ssP=(ss(:,:,1)+ss(:,:,2))/(sum(n_subs)-2);
    GRP.grands_stder(:,:,neo_bin)=sqrt( ssP*sum(n_subs)/(n_subs(1)*n_subs(2)) );
    %t-score
    GRP.grands_t=GRP.grands./GRP.grands_stder;
end

fprintf('<<New bin successfully created>>\n');
cc=GRP.bin_info(neo_bin).condcode;
fprintf('Condition Code %d: %s\n',cc,GRP.condesc{cc});
fprintf('Bin %d: %s\n',neo_bin,GRP.bin_info(neo_bin).bindesc);
hist_cmd=sprintf('GRP=bin_opGRP(GRP,''%s'',''%s'',%d,''%s'',%d,''%s'',%d,''%s'');', ...
    op,groupA,binA,groupB,binB,GRP.bin_info(neo_bin).bindesc,cc,GRP.condesc{cc});
GRP.history{length(GRP.history)+1}=hist_cmd;
GRP.saved='no';


