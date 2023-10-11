function h = rplot(x,y,varargin)

ncolores = size(y,2);
if ncolores>1
    %rainbow_colors(ncolores);
    colores = [[0 100 255]; % blue
    [3 201 183];
    [33 255 0];
    [235 255 2];
    [255 113 0];
    [239 1 42]]/255; % red
else
    colores = [0 0 1];
end

indToRemove = [];
for i=1:length(varargin)
    if ischar(varargin{i}) && strcmp(varargin{i},'colorType')
        colores = rainbow_colors(ncolores,'colorType',varargin{i+1});
        indToRemove = [indToRemove i,i+1];
    elseif ischar(varargin{i}) && strcmp(varargin{i},'colors')
        colores = varargin{i+1};
        %remove from varargin
        indToRemove = [indToRemove i,i+1];
    end
end

if not(isempty(indToRemove))
    varargin(indToRemove) = [];
end


hold on
if not(isempty(varargin))
    if not(isempty(x))
        for i=1:ncolores
            h(i) = plot(x,y(:,i),'color',colores(i,:),'markerfacecolor',colores(i,:),varargin{:});
        end
    else
        for i=1:ncolores
            h(i) = plot(y(:,i),'color',colores(i,:),'markerfacecolor',colores(i,:),varargin{:});
        end
    end
else
    if not(isempty(x))
        for i=1:ncolores
            h(i) = plot(x,y(:,i),'color',colores(i,:));
        end
    else
        for i=1:ncolores
            h(i) = plot(y(:,i),'color',colores(i,:));
        end
    end
end

format_figure(gcf);
set(gcf,'Color','w')
set(gca,'FontSize',16)
set(gca,'TickDir','out')

set(gca,'TickLength',[0.015 0.015])

