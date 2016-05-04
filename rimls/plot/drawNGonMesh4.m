function h = drawNGonMesh4( VK, EToV, OP1, VAL1, OP2, VAL2, lcolor,lwidth )
% Plot DG mesh (wireframe).
% Note:
%    VK- Coordinate
%    EToV     - mesh connectivity column one is number of polygon sides
%    xndoes   - DG Nodes
%
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
    case 7
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
% ncol = length(VK1(1,:)) ;
Nel = length(EToV)  ;
%lcolor = 'k' ;
for i = 1: Nel
    
    % Plot polygonal %
    ngon = EToV(i,1) ;
    vx  = VK(EToV(i,2:ngon+1),:) ;
    for ig = 1: ngon
        ib = ig ;
        ie = mod(ig,ngon) + 1 ;  
        plot( [vx(ib,1) vx(ie,1)], ...
              [vx(ib,2) vx(ie,2)], 'Color',lcolor,'LineWidth',lwidth ) ;
    end
    
    xp(1:ngon,i) = vx(:,1) ;
    yp(1:ngon,i) = vx(:,2) ;
   
    x = sum(xp(1:ngon,i))/ngon ;
    y = sum(yp(1:ngon,i))/ngon ;
    
    if ( strcmpi(OPD{1,2},'on') )
        txt = sprintf('%d',i) ;
        text(x,y,txt,'Color',lcolor,'FontSize',6) ;
    end
end

if ( strcmpi(OPD{2,2},'on') )
    for i = 1: length(VK)
        x = VK(i,1) ;
        y = VK(i,2) ;
        txt = sprintf('%d',i) ;
        text(x,y,txt,'Color',lcolor,'FontSize',6) ;
    end
end
    
h = gcf ; 