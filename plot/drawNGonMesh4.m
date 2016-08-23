function h = drawNGonMesh4( VK, EToV, lcolor, OP1, VAL1, elcolor, OP2, VAL2, ndcolor )
% Plot DG mesh (wireframe).
% Note:
%    VK- Coordinate
%    EToV     - mesh connectivity column one is number of polygon sides
%    xndoes   - DG Nodes


[NN ~] = size(VK) ;
[Nel ~] = size(EToV)  ;


elems = [1:Nel];
nodes = [1:NN];





h = findobj ;
if ( findobj == 0 )
    figure ;
    hold ;
end
OPD{1,1} = 'ElNum' ;
OPD{1,2} = 'off' ;
OPD{2,1} = 'NodeNum' ;
OPD{2,2} = 'off' ;
% OPD{3,1} = 'NodeNum' ;
% OPD{3,2} = 'off' ;

switch nargin
    case 5
        Opt{1,1} = OP1 ;
        Val{1,1} = VAL1 ;
    case 9
        Opt{1,1} = OP1 ;
        Opt{2,1} = OP2 ;
        Val{1,1} = VAL1 ;
        Val{2,1} = VAL2 ;
    otherwise
        Opt = [] ;
end
if ( length(Opt) > 0 )
    for i = 1: length(Opt)
        for l = 1: 2
            if ( strcmp((lower(Opt{i,1})),lower(OPD{l,1})) )
                OPD{l,2} = Val{i,1} ;
                break ;
            end
        end
    end
end

% Nel = length(VK1(:,1)) ;

%lcolor = 'k' ;
for i = 1: Nel
    
%     if ismember(i,elems)
    
    % Plot polygonal %
    ngon = EToV(i,1) ;
    if ngon == 9
       ngon = 8; 
    end
    vx  = VK(EToV(i,2:ngon+1),:) ;
    for ig = 1: ngon
        ib = ig ;
        ie = mod(ig,ngon) + 1 ;  
        plot( [vx(ib,1) vx(ie,1)], ...
              [vx(ib,2) vx(ie,2)], 'Color',lcolor ) ;
    end
%     end
end


if ( strcmpi(OPD{1,2},'on') )
    for i = 1:Nel
        ngon = EToV(i,1) ;
        
        xp = VK(EToV(i,2:ngon+1),1) ;
        yp = VK(EToV(i,2:ngon+1),2) ;
        
        x = sum(xp)/ngon ;
        y = sum(yp)/ngon ;
        
        if ismember(i,elems) 
          txt = sprintf('%d',i) ;
          text(x,y,txt,'Color',elcolor,'FontSize',6) ;
        end
    end
end

if ( strcmpi(OPD{2,2},'on') )
    for i = 1: length(VK)
        x = VK(i,1) ;
        y = VK(i,2) ;
        
        if ismember(i,nodes)
          txt = sprintf('%d',i) ;
          text(x,y,txt,'Color',ndcolor,'FontSize',6) ;
        end
    end
end
    
h = gcf ; 