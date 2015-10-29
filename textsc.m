% textsc() - places text in screen coordinates and places
%            a title at the top of the figure.
%
% Usage:
%   H = textsc(X,Y,TXT) places the text string, TXT
%   at the normalized coordinates X and Y.  H is the
%   handle to the text object.
%
%   H = textsc(TXT,'title') places a title at the top
%   of the figure window.  This is useful when you
%   want to place a single title above multiple
%   subplots.
%
% Notes: textsc creates an invisible AXES which occupies
% the entire FIGURE.  The units of the AXES are
% normalized (range from 0 to 1).  textsc checks
% all the children of the current FIGURE to determine
% if an AXES exist which meets these criteria already
% exist.  If one does, then it places the text relative
% to that AXES.
%
% Author: John L. Galenski, January 21, 1994

% Written by John L. Galenski III
% All Rights Reserved  January 21, 1994
% LDM031695jlg

% $Log: textsc.m,v $
% Revision 1.2  2007/01/26 18:07:16  arno
% *** empty log message ***
%
% Revision 1.1  2002/04/05 17:36:45  jorn
% Initial revision
%
% June 2015, Made compatible with Matlab 2014

function H = textsc(x,y,txt);

% Basic error checking
if nargin < 2
    error('TEXTSC requires at least 2 inputs')
end

% Check to see if AXES already exist
ch = get(gcf,'Children');
ax = findobj(gcf,'Type','axes','Tag','TEXTSC');
if isempty(ax)
    ax = axes('Units','Normal','Position',[0 0 1 1], ...
        'Visible','Off','Tag','TEXTSC');
else
    axes(ax)
end

% Place the text
if nargin == 2 & isstr(x) & strcmp(lower(y),'title')  % Subplot title
    txt = x;
    x = .5;
    tmp = text('Units','normal','String','tmp','Position',[0 0 0]);
    ext = get(tmp,'Extent');
    delete(tmp)
    H = ext(4);
    y = 1 - .60*H;
end
h = text(x,y,txt,'VerticalAlignment','Middle', ...
    'HorizontalAlignment','Center','interpreter', 'none' );

% Make the original AXES current
if ~isempty(ch)
    if verLessThan('matlab', '8.4')
        set(gcf,'CurrentAxes',min(ch))
    else
        a=1;
        axSwitch=0;
        while ~axSwitch
            try
                set(gcf,'CurrentAxes',ch(a)); %needed to be compatible with M2014
                axSwitch=1;
                % colorbar axes cannot be made current axes
            catch
                a=a+1;
                if a>length(ch),
                   axSwitch=1; %no axes will work 
                end
            end
        end
    end
end

% Check for output
if nargout == 1
    H = h;
end
