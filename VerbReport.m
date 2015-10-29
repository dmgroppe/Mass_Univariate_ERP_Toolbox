function VerbReport(report,verbtag,VERBLEVEL)
% VerbReport() - Outputs messages if they exceed a desired level of importance.  
%
% Usage:
%  >> VerbReport(report,verbtag,VERBLEVEL)
%
% Inputs:
%   report    = a string that is some message to the user
%   verbtag   = an intger specifiying the importance of report
%   VERBLEVEL = an integer specifiying a threshold of importance
%               for displaying reports. If verbtag is less than VERBLEVEL, the
%               report will be displayed..
%
% Author: Tom Urbach
% Kutaslab

if nargin<3
  tmpVERBLEVEL = 3;
elseif isempty(VERBLEVEL)
  tmpVERBLEVEL = 3;
else
  tmpVERBLEVEL = VERBLEVEL;
end;

if verbtag <= tmpVERBLEVEL
	if ischar(report)
		fprintf('%s\n',report);
	else
		fprintf('%d\n',report);
	end;
end;
