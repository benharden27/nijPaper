function [h,cm2] = cmapmaker(varargin)

h = varargin{1};
if nargin>1
    switch(nargin)
        case 2
            if(ishandle(varargin{2}))
                hc = varargin{2};
            else
                cmch = varargin{2};
            end
        case 3
            if(ishandle(varargin{2}))
                hc = varargin{2};
                cmch = varargin{3};
            else
                hc = varargin{3};
                cmch = varargin{2};
            end
    end
end

if(~exist('cmch','var'))
    cmch = 'midw';
end

levs = get(h,'levelList');
a = cat(2,h.FacePrims.ColorData);
nlev = length(levs);

switch(cmch)
    
    case 'ds_temp'
    cmraw = [160,0,160;
        144,0,160;
        88,0,160;
        0,0,120;
        0,140,255;
        0,187,255;
        0,220,255;
        119,255,135;
        183,255,71;
        247,255,7;
        255,220,0;
        255,175,0;
        255,135,0;
        255,71,0;
        255,7,0;
        193,0,0;
        139,0,0;
        100,0,0]/255;
    
    case 'midw'
        cmraw = [103,0,31;...
            178,24,43;...
            214,96,77;...
            244,165,130;...
            253,219,199;...
            209,229,240;...
            146,197,222;...
            67,147,195;...
            33,102,172;...
            5,48,97]/255;
        cmraw = cmraw(end:-1:1,:);

        
    case 'zis'
        cmraw = [59,154,178;...
            120 183 197;...
            235 204 42;...
            225 175 0;...
            242 26 0]/255;

    case 'midy'
        cmraw = [165,0,38;...
            215,48,39;...
            244,109,67;...
            253,174,97;...
            254,224,144;...
            255,255,191;...
            171,217,233;...
            116,173,209;...
            69,117,180;...
            49,54,149]/255;
        cmraw = cmraw(end:-1:1,:);
        
    case 'brix'
        cmraw = [166,206,227;...
            31,120,180;...
            178,223,138;...
            51,160,44;...
            251,154,153;...
            227,26,28]/255;
    case 'bri'
        cmraw = [5,48,97;...
            31,120,180;...
            166,206,227;...
%             116,196,118;...
            43,129,86;...
%             253,191,111;...
%             254,217,118;...
            254,243,118;...
            255,127,0;...
            227,26,28]/255;
        
    case 'gray'
        cmraw = gray*1.1;
        cmraw(cmraw>1) = 1;
        
    case 'nice'
        cmraw = [158,1,66
            213,62,79
            244,109,67
            253,174,97
            254,224,139
            255,255,191
            230,245,152
            171,221,164
            102,194,165
            50,136,189
            94,79,162]/255;
        cmraw = cmraw(end:-1:1,:);

end


cmvec1 = linspace(0,1,length(cmraw(:,1)));
cmvec2 = linspace(0,1,nlev);
for n = 1:3
    cm(:,n) = interp1(cmvec1,cmraw(:,n),cmvec2);
end    
cm2 = uint8([cm';ones(1,size(cm,1))]*255);

for n = 1:size(a,2)
%    val = find(levs==get(a(n),'Cdata'),1);
%    if(isempty(val))
%        val = find(levs>get(a(n),'Cdata'),1)-1;
%    end
% %    val
% %    if(val>nlev)
% %        val = val-1;
% %    end
%    if(~isempty(val))
%        set(a(n),'facecolor',cm(val,:)); 
%    end
    h.FacePrims(n).ColorData = cm2(:,n);
end

if(exist('hc','var'))
a = get(hc,'children');
val1 = min(get(a(end),'Ydata'));
val1 = find(levs==val1);
% 
% for n = 1:length(a)
%    vals = get(a(n),'Ydata');
%    val = find(levs==min(vals));
%    if(val>nlev)
%        val = val-1;
%    end
%    if(~isempty(val))
%    set(a(n),'facecolor',cm(val,:)); 
%    end
% end
% val1
% length(a)
for n = 1:length(a)
%     n
   if(length(get(a(length(a)-n+1),'Ydata'))>3)
      set(a(length(a)-n+1),'facecolor',cm(val1+n-1,:));
   else
      set(a(length(a)-n+1),'facecolor',cm(end,:));
   end
end


end
    
   