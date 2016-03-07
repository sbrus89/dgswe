restoredefaultpath
clear all
close all
clc


grd_direc = {'~/data-drive/galveston/dgswe/tri/';
             '~/data-drive/galveston/dgswe/quad/';
             '~/data-drive/galveston/dgswe/quad/'};
sol_direc = {'~/data-drive/galveston/dgswe/tri/';
             '~/data-drive/galveston/dgswe/quad/';
             '~/data-drive/galveston/dgswe/quad_curve/'};
grd_name = {'galveston_tri.grd';
            'galveston_quad.grd';
            'galveston_quad.grd'};

nsnap = 60;

nsol = length(sol_direc);
ne = zeros(nsol,1);
nn = zeros(nsol,1);
mndof = zeros(nsol,1);

for sol = 1:nsol
    fid = fopen([grd_direc{sol},grd_name{sol}]) ;    
    agrid = fgetl(fid) ;   
    disp(agrid)
    N = fscanf(fid,'%g %g',2) ;
    ne(sol) = N(1);
    nn(sol) = N(2);
    fclose(fid);
end

mne = max(ne);
mnn = max(nn);

EToV = zeros(mne,4,nsol);
NX = zeros(mnn,2,nsol);
HB = zeros(mnn,nsol);
nelnds = zeros(mne,nsol);

for sol = 1:nsol
    
    [EToV(1:ne(sol),1:4,sol),VX(1:nn(sol),1:2,sol),HB(1:nn(sol),sol),nelnds(1:ne(sol),sol),~,~,~] = readfort14([grd_direc{sol},grd_name{sol}]);
    
    
    fid = fopen([sol_direc{sol},'modal2nodal.d']);
    
    nel_type = 4;
    
    nvert = zeros(1,nel_type);
    ndof = zeros(1,nel_type);
    m2n = zeros(4,36,nel_type);
    for et = 1:nel_type
        data  = fscanf(fid, '%g %g',2);
        ndof(et) = data(2);
        nvert(et) = data(1);
        m2n(1:nvert(et),1:ndof(et),et) = fscanf(fid,' %g ',[ndof(et) nvert(et)])';
    end
    fclose(fid);
    
    mndof(sol) = max(ndof);
    
end

mne = max(ne);

Hv = zeros(4,mne,nsnap,nsol);
Qxv = zeros(4,mne,nsnap,nsol);
Qyv = zeros(4,mne,nsnap,nsol);
zv = zeros(4,mne,nsnap,nsol);

for sol = 1:nsol
    
    fid_H = fopen([sol_direc{sol},'solution_H.d']);
    fid_Qx = fopen([sol_direc{sol},'solution_Qx.d']);
    fid_Qy = fopen([sol_direc{sol},'solution_Qy.d']);
    
    line1 = fgetl(fid_H);
    line2 = fgetl(fid_Qx);
    line3 = fgetl(fid_Qy);
    
    snap = 0;
    while ~feof(fid_H) && snap < nsnap
        
        snap = snap + 1;
        
        th = fscanf(fid_H,' %g ', 1); % read in time
        H = fscanf(fid_H,' %g ', [ne(sol) mndof(sol)])'; % read in H solution at time t
        
        tqx = fscanf(fid_Qx,' %g ', 1); % read in time
        Qx = fscanf(fid_Qx,' %g ', [ne(sol) mndof(sol)])'; % read in Qx solution at time t
        
        t(snap) = fscanf(fid_Qy,' %g ', 1); % read in time
        Qy = fscanf(fid_Qy,' %g ', [ne(sol) mndof(sol)])'; % read in Qy solution at time t
        
        for el = 1:ne(sol)
            if nelnds(el,sol) == 3
                Hv(1:3,el,snap,sol) = m2n(1:3,1:ndof(1),1)*H(1:ndof(1),el);
                Qxv(1:3,el,snap,sol) = m2n(1:3,1:ndof(1),1)*Qx(1:ndof(1),el);
                Qyv(1:3,el,snap,sol) = m2n(1:3,1:ndof(1),1)*Qy(1:ndof(1),el);
                zv(1:3,el,snap,sol) = Hv(1:3,el,snap)-HB(EToV(el,1:3,sol),sol);
            elseif nelnds(el,sol) == 4
                Hv(1:4,el,snap,sol) = m2n(1:4,1:ndof(2),2)*H(1:ndof(2),el);
                Qxv(1:4,el,snap,sol) = m2n(1:4,1:ndof(2),2)*Qx(1:ndof(2),el);
                Qyv(1:4,el,snap,sol) = m2n(1:4,1:ndof(2),2)*Qy(1:ndof(2),el);
                zv(1:4,el,snap,sol) = Hv(1:4,el,snap)-HB(EToV(el,1:4,sol),sol);
            end
        end
    end    
    
    
end
nsnap = snap;

vel = sqrt((Qxv./Hv).^2 + (Qyv./Hv).^2);
velscale = zeros(nsnap,1);

for snap = 1:nsnap
    velscale(snap) = max(max(max(vel(:,:,snap,:))));
end

for sol = 1:nsol
    save([sol_direc{sol},'velscale.mat'],'velscale')
end