function h = drawNGonMesh4_bou( VK, EToV,boudat,lcolor, OP1, VAL1, elcolor, OP2, VAL2, ndcolor )
% Plot DG mesh (wireframe).
% Note:
%    VK- Coordinate
%    EToV     - mesh connectivity column one is number of polygon sides
%    xndoes   - DG Nodes


[NN ~] = size(VK) ;
[Nel ~] = size(EToV)  ;


% elems = [1:Nel];
elems = [113489 113490];
% elems = [3934
%          8721
%         11519
%         12164
%         12170
%         12397
%         13345
%         13762
%         14113
%         16611
%         18260
%         20303
%         22508
%         23298
%         26491
%         26521
%         26562
%         28735
%         37625
%         63934
%         67396
%         71049
%        111640
%        129836
%        156119
%        186856
%        216766
%        220760
%        230667
%        238085
%        261041
%        279486
%        281644
%        318004
%        378101
%        406123
%        439630
%        440140
%        460746];

% nodes = [1:NN];
%nodes = [994 73376 73384];


plot_flag = zeros(Nel,1);


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
    case 6
        Opt{1,1} = OP1 ;
        Val{1,1} = VAL1 ;
    case 10
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

    for j = 1:ngon
       if ismember(EToV(i,j+1),boudat)
          plot_flag(i) = 1;
          break
       else
          plot_flag(i) = 0;       
       end
    end
    
    if plot_flag(i) == 1
        vx  = VK(EToV(i,2:ngon+1),:) ;
        for ig = 1: ngon
            ib = ig ;
            ie = mod(ig,ngon) + 1 ;
            plot( [vx(ib,1) vx(ie,1)], ...
                [vx(ib,2) vx(ie,2)], 'Color',lcolor ) ;
        end
    end
%     end
end


if ( strcmpi(OPD{1,2},'on') )
    for i = 1:Nel
        if ismember(i,elems)
%         if plot_flag(i) == 1
            ngon = EToV(i,1) ;
            
            xp = VK(EToV(i,2:ngon+1),1) ;
            yp = VK(EToV(i,2:ngon+1),2) ;
            
            x = sum(xp)/ngon ;
            y = sum(yp)/ngon ;
            
            
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