function setTScm(varargin)

temp_cmap = [1.0000         0    1.0000;
             0.6000         0    0.6000;
             0.3000         0    0.5000;
             0              0    0.4000;
             0              0    0.7000;
             0              0    1.0000;
             0              0.3333    1.0000;
             0              0.6667    1.0000;
             0              1.0000    1.0000;
             0              0.8000         0;
             1.0000         1.0000         0;
             1.0000         0.7500         0;
             1.0000    0.5000         0;
             1.0000         0         0;
             0.8000         0         0;
             0.6000         0         0;
             0.4000         0         0];
         
temp_levels = [ -2.0000    -1.0000   -0.5000   -0.2500  0 ...
                0.2500    0.5000    0.7500    1.0000    1.5000 ...
                2.0000    3.0000    4.0000    5.0000    6.0000 ...
                8.0000   10.0000   12.0000];

            
sal_cmap = [0.3000         0    0.5000;
            0         0    0.4000;
            0         0    0.6000;
            0         0    0.8000;
            0         0    1.0000;
            0    0.3333    1.0000;
            0    0.6667    1.0000;
            0    1.0000    1.0000;
            0    0.8000         0;
            1.0000    1.0000         0;
            1.0000    0.5000         0;
            1.0000         0         0;
            0.7667         0         0;
            0.5333         0         0;
            0.3000         0         0;
            0.5000    0.2500         0;
            0.7000    0.5000         0];
            
sal_levels = [31.0000   33.0000    34.0000   34.2000   34.4000...
              34.6000   34.8000   34.8500   34.8800   34.9000...
              34.9100   34.9200   34.9400   34.9700   35.0000...
              35.1000   35.2000   35.3000];
    
          
h = varargin{1};
if strcmp(varargin{2},'temp')
    cmap = temp_cmap;
    levels = temp_levels;
elseif strcmp(varargin{2},'sal')
    cmap = sal_cmap;
    levels = sal_levels;
else
    return
end
% load ~/Documents/projects/kogur/code/cmaps_toBen.mat


if(isa(h,'matlab.graphics.chart.primitive.Contour'))
    levs = get(h,'levelList');
    a = cat(2,h.FacePrims.ColorData);
    cm2 = uint8([cmap';ones(1,size(cmap,1))]*255);
    for i = 1:size(a,2)
        ci = find(levels==levs(i));
        h.FacePrims(i).ColorData = cm2(:,ci);
    end
    
else
    levs = str2num(h.YTickLabel);
    a = cat(1,h.Children.FaceColor);
    if(size(a,1)~=size(levs,1))
        levs = levs(1:end-1);
    end
    cm2 = uint8(cmap*255);
    for i = 1:size(a,1)
        j = size(a,1)-i+1;
        ci = find(levels==levs(i));
        h.Children(j).FaceColor = cm2(ci,:); 
    end
end
% 
% if nargin > 1
%     minC = min(Cdata(:));
%     if(minC<-2)
%         minC = -2;
%     end
%     minCPT = temp_levels-minC;
%     ci = find(minCPT<=0,1,'last');
% end

% for i = 1:size(a,2)
%     vec = get(blocks(i),'Cdata');
%     if(~isempty(vec))
%         if(i==length(blocks))
%             minval = temp_levels-vec;
%             val = find(minval<=0,1,'last');
%         else
%             val = findnear(temp_levels,vec);
%         end
%         val(val>size(temp_cmap,1)) = size(temp_cmap,1);
%         set(blocks(i),'facecolor',temp_cmap(val,:));
%     else 
%         set(blocks(length(blocks)-i+1),'facecolor',temp_cmap(ci+i-1,:));
%     end
%         
% end







