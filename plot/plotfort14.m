% plotmesh

% direc = '/home/sbrus/data-drive/galveston/SL18+TX33_cutout/';
% grid = 'galveston_0cuttoff_boundary_jetty.grd';


direc = '/home/sbrus/data-drive/galveston/dgswe/tri/';
grid = 'galveston_tri.grd';

figure
set(gcf,'PaperPositionMode','auto') ;

[EToV,VX,B,nelnds,opedat,boudat,title] = readfort14([direc,grid]) ;
p = 1 ;

pdeplot( VX', [], EToV','xydata',B,'colormap', 'jet', 'mesh','on' ) ;
axis image


% % draw mesh
% [G,ElIdDG,xDG,vKN,vKI] = MakeNodalDGNodes( p, EToV, VX ) ;
% %drawNGonMesh3( vKI, vKN, ElIdDG, xDG, 'ElNum', 'off', 'Nodes', 'on') ;
% 
% DEToV = zeros(length(EToV(:,1)),4) ;
% DEToV(:,1) = 3 ;
% DEToV(:,2:4) = EToV ;
% drawNGonMesh4( VX, DEToV, 'ElNum', 'off', 'NodeNum', 'off' )
% 
% % BT = zeros(3*length(EToV),1) ; 
% % for i = 1: length(EToV)
% %    ibeg = 3*(i - 1) + 1 ;
% %    iend = 3*i ; 
% %    BT(ibeg:iend) = -B(EToV(i,:)) ;
% %end
% % % draw bathymetry
% %VisG = CreateDGVisOBJ( vKI, vKN, ElIdDG, xDG ) ;
% %[xPatch,yPatch,TriDAT] = MakeDGPatchDATA( VisG ) ;
% %cdat = MakeDGZDATA( VisG, TriDAT, BT ) ;
% %figure ; 
% %draw2DDGsolLV0( xPatch, yPatch, cdat, ...
% %   'EdgeColor', 'black' ) ;
% 
% 
% %figure ;
% %trisurf( EToV, VX(:,1), VX(:,2), -B ) ; 