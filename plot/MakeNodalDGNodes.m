function [G,ElIdDG,xDG,vKN,vKI] = MakeNodalDGNodes( p, EToV, VX )
% Read Adcirc grid
%
%
% 
[r,s] = Nodes2D( p ) ;
[r,s] = xytors( r, s ) ;

% No. nodes in each element
Np = (p + 1)*(p + 2)/2 ;

% No. of element
Nel = length(EToV(:,1)) ;
G.nEl = Nel  ;

pmax = 0 ;
pmin = 100 ;

idbeg = 0 ;
for iel = 1: Nel
    pmax = max(pmax,p) ; 
    pmin = min(pmin,p) ;    
    
    G.GE(iel).p = p ; % order of polynomail
    G.GE(iel).Np = Np ; % no. of nodes in each element
    
    % vertices
    G.GE(iel).xVi = VX(EToV(iel,:)',1:2) ;
    G.GE(iel).idVi = EToV(iel,:) ; % Global index
    
    % DG ID
    G.GE(iel).idDG = zeros(Np,1) ;
    for idx = 1: Np
        G.GE(iel).idDG(idx) = idbeg + idx ;
    end
    idbeg = idbeg + Np ;
    
    % Local index of the vertices
    G.GE(iel).idViL = [1, p + 1, Np] ;
    
    % DG: Global index of the vertices 
    G.GE(iel).idViG = G.GE(iel).idDG(G.GE(iel).idViL) ;
end
G.nDGNodes = idbeg ;
G.NodeCoord = zeros(G.nDGNodes,2) ;

npm = (pmax - pmin + 1) ;

%------------
for ip = 1: npm
    RObj.p(ip) = pmin + (ip - 1) ;
   
    [rr,ss] = Nodes2D( pmin + (ip - 1) ) ;
    xx = [rr ss] ;
   
    TRI = delaunay(xx(:,1),xx(:,2)) ; 
    
    % Check the quality of TRI
    % and single a flat element 
    J = zeros(length(TRI(:,1)),1) ;
    for l = 1: length(TRI(:,1))
        xv = xx(TRI(l,:),:) ; 
        
        xr = (xv(2,:) - xv(1,:))/2 ;
        xs = (xv(3,:) - xv(1,:))/2 ;
        
        J(l) = xr(1)*xs(2) - xs(1)*xr(2) ;
    end
    Jmax = max(abs(J)) ; 
    idx = find( J/Jmax > 1.0e-10 ) ; 
      
    TRI = TRI(idx,:)' ;
    
    RObj.SubEL(ip).Tri = TRI ;  
end
%------------

for iel = 1: Nel
    [x,y] = rstoxy( r, s, G.GE(iel).xVi ) ;
    
    % Coordinates of nodes
    G.NodeCoord(G.GE(iel).idDG,1:2) = [x y] ;
end

% connectivity
sk = 0 ;
for i = 1: Nel
    ii = find( RObj.p == G.GE(i).p ) ;
    
    TRI = RObj.SubEL(ii).Tri ;
    
    G.GE(i).TriDat = TRI ;
    
    nn = size(TRI,2) ;
    
    ibeg = sk + 1 ;
    iend = ibeg + nn - 1 ;
    G.GE(i).SubElmID = ibeg:iend ;
    
    sk = sk + nn ;
end

% For visualization
ElIdDG = zeros(Nel,3) ;
for iel = 1: Nel
    ElIdDG(iel,1) = iel ;
    ElIdDG(iel,2) = G.GE(iel).idDG(1) ;
    ElIdDG(iel,3) = G.GE(iel).idDG(Np) ;
end

% X
xDG = zeros(G.nDGNodes,4) ;
m = 0 ;
for iel = 1: Nel
    idDG = G.GE(iel).idDG ;
    for n = 1: Np
        m = m + 1 ;
        xDG(m,1:2) = G.NodeCoord(idDG(n),1:2) ;
        xDG(m,3) = idDG(n) ;
        xDG(m,4) = n ;
    end
end

% VKN
m = 0 ;
for iel = 1: Nel
    for n = 1: length(G.GE(iel).idViG)
        m = m + 1 ;
        vKN(m,1) = G.GE(iel).idViG(n) ;
    end
end

% VKI
vKI = zeros(iel,3) ;
m = 0 ;
for iel = 1: Nel
    vKI(iel,1) = iel ;
    me = m + length(G.GE(iel).idViG) ;
    
    vKI(iel,2) = m + 1 ;
    vKI(iel,3) = me ;
    
    m = me ;
end
