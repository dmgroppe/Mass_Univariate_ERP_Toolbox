function [justpath, justname]=pathNname(full_file)
%function [justpath, justname]=pathNname(full_file)
%
% Splits a string specifiying a file name and its path into two
% strings (one for just the filename and one for just the
% path). For some reason EEGLAB wants the two separate when
% importing data.
%
% Input:
%  full_file = a string specifying a file name and its path
%
% Outputs:
%  justpath = a string containing just the path part of full_file
%  justname = a string containing just the file name part of full_file
%
% Example: 
% >> setfile='/homes/dgroppe/MK2MAT/tryRT.set';
% >> [jpath, jname]=pathNname(setfile);
%


justname='';
justpath='';
n=length(full_file);

%Figure out which path convention is being used
if ismac || isunix
    %Mac or Unix OS
    slash='/';
else
    %Windows OS
    slash='\';
end

flg=0;
for a=n:-1:1,
  if (flg==0) && (full_file(a)~=slash)
    justname=[full_file(a) justname];
  else
    if full_file(a)==slash
      flg=1;
    end
    justpath=[full_file(a) justpath];
  end
end

