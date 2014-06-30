restoredefaultpath
clear all
close all
clc

grd_direc = '~/dgswe/grids/';
sol_direc = '~/dgswe/output/';
grd_name = 'inlet1_quad.grd';
nsnap = 100;




[EToV,VX,HB,nelnds,~,~,~] = readfort14([grd_direc,grd_name]);
[ne,~] = size(EToV);
[nn,~] = size(VX);





fid = fopen([sol_direc,'modal2nodal.d']);

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

mndof = max(ndof);

fid_H = fopen([sol_direc,'solution_H.d']);
fid_Qx = fopen([sol_direc,'solution_Qx.d']);
fid_Qy = fopen([sol_direc,'solution_Qy.d']);

line = fgetl(fid_H);
line = fgetl(fid_Qx);
line = fgetl(fid_Qy);

Hv = zeros(4,ne,nsnap);
Qxv = zeros(4,ne,nsnap);
Qyv = zeros(4,ne,nsnap);

snap = 0;
while ~feof(fid_H) && snap < nsnap
    
    snap = snap + 1;
    
    th = fscanf(fid_H,' %g ', 1); % read in time
    H = fscanf(fid_H,' %g ', [ne mndof])'; % read in H solution at time t
    
    tqx = fscanf(fid_Qx,' %g ', 1); % read in time
    Qx = fscanf(fid_Qx,' %g ', [ne mndof])'; % read in Qx solution at time t
    
    t(snap) = fscanf(fid_Qy,' %g ', 1); % read in time
    Qy = fscanf(fid_Qy,' %g ', [ne mndof])'; % read in Qy solution at time t
    
    for el = 1:ne
        if nelnds(el) == 3
            Hv(1:3,el,snap) = m2n(1:3,1:ndof(1),1)*H(1:ndof(1),el);
            Qxv(1:3,el,snap) = m2n(1:3,1:ndof(1),1)*Qx(1:ndof(1),el);
            Qyv(1:3,el,snap) = m2n(1:3,1:ndof(1),1)*Qy(1:ndof(1),el);
        elseif nelnds(el) == 4
            Hv(1:4,el,snap) = m2n(1:4,1:ndof(2),2)*H(1:ndof(2),el);
            Qxv(1:4,el,snap) = m2n(1:4,1:ndof(2),2)*Qx(1:ndof(2),el);
            Qyv(1:4,el,snap) = m2n(1:4,1:ndof(2),2)*Qy(1:ndof(2),el);
        end
    end
end


vel = sqrt((Qxv./Hv).^2 + (Qyv./Hv).^2);

Hmax = max(max(max(Hv)));
Qxmax = max(max(max(Qxv)));
Qymax = max(max(max(Qyv)));
velmax = max(max(max(vel)));

Hmin = min(min(min(Hv)));
Qxmin = min(min(min(Qxv)));
Qymin = min(min(min(Qyv)));
velmin = min(min(min(vel)));


for tsnap = 1:snap
    
    figure
    axis equal

    disp(['Time snap: ',num2str(tsnap),'/',num2str(snap)])
    [day,hr,minute,sec] = s2dhms(t(tsnap));
    
    hold on
    velmin = min(min(vel(:,:,tsnap)));
    velmax = max(max(vel(:,:,tsnap)));
    caxis([velmin velmax])
    for el = 1:ne
        if nelnds(el) == 3
            fill(VX(EToV(el,1:3),1),VX(EToV(el,1:3),2),vel(1:3,el,tsnap))
        elseif nelnds(el) == 4
            fill(VX(EToV(el,1:4),1),VX(EToV(el,1:4),2),vel(1:4,el,tsnap))
        end
    end
    hold off
    
    ttext = ['Velocity solution: t = ',num2str(t(tsnap)),' (Day:  ',num2str(day),', Hour:  ',num2str(hr),', Minute:  ',num2str(minute),', Second:  ',num2str(sec),')'] ;
    title(ttext)
    xlabel('x')
    ylabel('y')
    colorbar
    axis image
    
    pause(.01)
end