function autolabelfigures()

% AUTOLABELFIGURES auto generate figure titles in titlebar
%
% Useful function if you make a lot of figures and want 
% to automatically Name them (so the title shows up in the Taskbar 
% for easy identification) based on Titles or Y Axis labels. 
% Could use some enhancement: it doesn't work so well if there are subplots

% Peter Adamczyk

kids = get(0,'Children');
for jj = 1:length(kids)
    figure(kids(jj));
    if isempty(get(gcf,'Name'))
        if ~isempty(get(get(gca,'Title'),'String'))
            set(gcf,'Name',get(get(gca,'Title'),'String'))
        else
            set(gcf,'Name',get(get(gca,'YLabel'),'String'))
        end
    end
end