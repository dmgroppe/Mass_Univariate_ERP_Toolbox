% bin_dif() - Computes the difference waves from two bins and stores the 
%             difference as a new bin in a GND variable. The difference
%             wave will be binA-binB.
%             
% Usage:
%  >> GND=bin_dif(GND,binA,binB,dif_bindesc,dif_ccode,dif_ccdescr);
%
% Required Inputs:
%   GND    - A GND structure variable.  To create a GND variable 
%            from Kutaslab ERP files (e.g., *.mas files) use avgs2GND.m.  
%            To do the same from EEGLAB *.set files use sets2GND.m.
%            See Mass Univariate ERP Toolbox documentation for detailed  
%            information about the format of a GND variable. 
%   binA   - [integer] A bin index. Use the function headinfo.m to see what 
%            bins are available.
%   binB   - [integer] A bin index. The difference wave will be binA-binB.
%
% Optional Inputs:
%   dif_bindesc - [string] The bin descriptor for the new bin being
%                 created. {default: 'Bin #-Bin ##', where # is the value
%                 binA and ## is the value of binB}
%   dif_ccode   - [integer] The condition code of the new bin being
%                 created. Condition codes are specific to Kutaslab data
%                 and can be ignored if your lab doesn't support them.
%                 {default: condition code of binA}
%   dif_ccdescr - [string] The condition code descriptor for the new bin being
%                 created. Condition codes are specific to Kutaslab data 
%                 and can be ignored if your lab doesn't support them.
%                 {default: existing descriptor for that condition code or 
%                 'Not Specified'}
%
% Notes:
% -binA and binB should index bins starting at 1 (i.e., Bin 0 is assumed to
% contain cal pulses and is stored apart from bins in GND variables)
%
% Example:
% >> GND=bin_dif(GND,2,1,'Targets-Standards'); 
%
% Author:
% David Groppe
% Kutaslab, 3/2010

%%%%%%%%%%%%%%%% REVISION LOG %%%%%%%%%%%%%%%%%
% 3/15/2010-bin_dif command is added to GND variable history
%
% 12/9/2010-Original GND variable returned in case of failure (as opposed
% to NaN)

%%%%%%%%%%%%%%%% FUTURE WORK %%%%%%%%%%%%%%%%%
% - add verblevel? 

function GND=bin_dif(GND,binA,binB,dif_bindesc,dif_ccode,dif_ccdescr)

if nargin<4,
   dif_bindesc=sprintf('Bin %d-Bin %d',binA,binB); 
end
if nargin<5,
   dif_ccode=[]; 
end
if nargin<6,
   dif_ccdescr=[]; 
end

GND_copy=GND; %copy original variable in case something fails

[n_chans, n_pts, n_bins, n_subs]=size(GND.indiv_erps);
neo_bin=n_bins+1;
GND.indiv_erps(:,:,neo_bin,:)=zeros(n_chans,n_pts,1,n_subs);
use_subs=[];
for sub=1:n_subs,
    if sum(sum(isnan(GND.indiv_erps(:,:,binA,sub))))
        fprintf('Sub %d does not have ERPs for Bin %. and will be ignored',sub,binA);
        GND.indiv_erps(:,:,neo_bin,sub)=NaN;
        GND.indiv_bin_ct(sub,neo_bin)=0;
        GND.indiv_bin_raw_ct(sub,neo_bin)=0;
    else
        if sum(sum(isnan(GND.indiv_erps(:,:,binB,sub))))
            fprintf('Sub %d does not have ERPs for Bin %. and will be ignored',sub,binB);
            GND.indiv_erps(:,:,neo_bin,sub)=NaN;
            GND.indiv_bin_ct(sub,neo_bin)=0;
            GND.indiv_bin_raw_ct(sub,neo_bin)=0;
        else
            GND.indiv_erps(:,:,neo_bin,sub)=GND.indiv_erps(:,:,binA,sub)- ...
                GND.indiv_erps(:,:,binB,sub);
            GND.indiv_bin_ct(sub,neo_bin)=-1;
            GND.indiv_bin_raw_ct(sub,neo_bin)=-1;
            use_subs=[use_subs sub];
        end
    end
end
if isempty(use_subs),
    GND=GND_copy;
    error('No subjects have ERPs in both Bin %d and Bin %d.  No changes were made to GND variable.\n',binA,binB);
else
    GND.sub_ct(neo_bin)=length(use_subs);
    GND.grands(:,:,neo_bin)=mean(GND.indiv_erps(:,:,neo_bin,use_subs),4);
    GND.grands_stder(:,:,neo_bin)=std(GND.indiv_erps(:,:,neo_bin,use_subs),0,4)/sqrt(GND.sub_ct(neo_bin));
    GND.grands_t(:,:,neo_bin)=GND.grands(:,:,neo_bin)./GND.grands_stder(:,:,neo_bin);
    GND.bin_info(neo_bin).bindesc=dif_bindesc;
    ccode1=GND.bin_info(binA).condcode;
    ccode2=GND.bin_info(binB).condcode;
    if isempty(dif_ccode)
        if ccode1==ccode2,
            GND.bin_info(neo_bin).condcode=ccode1;
        else
            watchit(sprintf(['Bin %d and %d belong to different condition codes. ', ...
                'Condition code for new bin not specified.  Using Condition Code of %d by default.\n'], ...
                binA,binB,ccode1));
            GND.bin_info(neo_bin).condcode=ccode1;
        end
    else
        GND.bin_info(neo_bin).condcode=dif_ccode;
        if length(GND.condesc)>=dif_ccode,
            exist_desc=1;
        else
            exist_desc=0;
        end
        if isempty(dif_ccdescr),
            if ~exist_desc,
               watchit(sprintf('There is no condition code descriptor for Condition Code %d.  ', ... 
                   dif_ccode));
               GND.condesc{dif_ccode}='Not Specified';
            end
        else
            if exist_desc,
                if ~strcmpi(dif_ccdescr,GND.condesc{dif_ccode}),
                    fprintf('Existing descriptor for Condition Code %d is "%s".\n', ...
                        dif_ccode,GND.condesc{dif_ccode});
                    resp=input(sprintf('Do you want me to overwrite existing descriptor with new descriptor of "%s"? (y or n)', ...
                        dif_ccdescr),'s');
                    if strcmpi(resp,'y') || strcmpi(resp,'yes'),
                        GND.condesc{dif_ccode}=dif_ccdescr;
                    end
                end
            else
                GND.condesc{dif_ccode}=dif_ccdescr;
            end
        end
    end
end

fprintf('<<New bin successfully created>>\n');
cc=GND.bin_info(neo_bin).condcode;
fprintf('Condition Code %d: %s\n',cc,GND.condesc{cc});
fprintf('Bin %d: %s\n',neo_bin,GND.bin_info(neo_bin).bindesc);
hist_cmd=sprintf('GND=bin_dif(GND,%d,%d,''%s'',%d,''%s'');',binA,binB, ...
    GND.bin_info(neo_bin).bindesc,cc,GND.condesc{cc});
GND.history{length(GND.history)+1}=hist_cmd;
GND.saved='no';

