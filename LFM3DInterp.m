% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %
% %          Interpolation code
% %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NI=53; % number of i-shells
% NJ=50; % number of j-shells
% NK=65; % number of k-shells
%
% % load mhd positions, change the directory when you are working on other
% % systems
% %
% posname='/thayerfs/research/lfm/bzhang/Poynting/para_clock/EB/mhdpos.txt';
% %posname='Q:\lfm\bzhang\Poynting\para_clock\EB\mhdpos.txt'
% d=load(posname);
%
% filename=['/thayerfs/research/lfm/bzhang/Poynting/m10c/BE-',int2str(200),'.txt']
% %filename=['Q:\lfm\bzhang\Poynting\m10c\BE-',int2str(200),'.txt']
% input=load(filename);
%
% Bx=input(:,1);
% By=input(:,2);
% Bz=input(:,3);
%
% x=d(:,1);
% y=d(:,2);
% z=d(:,3);
%
% x=reshape(x,NI,NJ,NK);    % the positions and B-fields need to be converted to (i,j,k)
% y=reshape(y,NI,NJ,NK);    % this nice "reshape" matlab function can do this job easily!
% z=reshape(z,NI,NJ,NK);
% bx=reshape(Bx,NI,NJ,NK);
% by=reshape(By,NI,NJ,NK);
% bz=reshape(Bz,NI,NJ,NK);
%
% data=bz;

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% First determine which two k-shells is the point in between
%
% x0=-5; % this is just a random point for testing
% y0=3;
% z0=0;

function d=LFM3DInterp(x,y,z,data,x0,y0,z0)

[NI NJ NK]=size(x);

for k=1:NK-1
    
    theta1=pi-sign(z(5,5,k))*acos(-y(5,5,k)/sqrt(y(5,5,k)^2+z(5,5,k)^2));         % theta1 and theta2 are angle between k-shells and x-y plane
    theta2=pi-sign(z(5,5,k+1))*acos(-y(5,5,k+1)/sqrt(y(5,5,k+1)^2+z(5,5,k+1)^2)); % theta is the angle of the test point... for a k value (say k=KK)
    % when theta1<theta<theta2, we say that this point is between KK and KK+1 k-shells
    if (z0==0)                                                                    % **should be really careful when KK=64, which is the last k-shell
        theta=pi/2-pi/2*sign(y0);                                                 % for example, for the k=64 shell theta=6.22
    else                                                                          % then when we go up one more shell, k=65 shell which is the same with k=1shell,
        theta=pi-sign(z0)*acos(-y0/sqrt(y0^2+z0^2));                              % the theta value for k=65 shell is 0.05. Thus We have to add 2pi to this theta
    end                                                                           % ****THINK ABOUT THIS: WHY CANNOT JUST USE theta=acos(y/sqrt(y^2+z^2))?*****
    
    KK=64;
    if ((theta-theta1)*(theta-theta2)<0) % This condition is to decide if theta is between (theta1,theta2)
        KK=k;                            % KK is the index, which means the point is between KK and KK+1 k-shells
        break;
    end
end

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%figure
%plot_2kshells(x,y,z,zeros(size(x))+1,KK,KK+1,'surf');hold on;h1=colorbar;   % these functions are for plotting.
%plot_jshells(x,y,z,zeros(size(x)),15);
%plot_ishells(x,y,z,zeros(size(x))-1,5);
%caxis([-2 2]);shading faceted;view(90,0)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

rho0=sqrt(y0^2+z0^2);  % this part calculates the two points mapped onto the two k-shells: k=KK and KK+1
y1=rho0*cos(theta1);   % now we have (x0,y1,z1) on k=KK shell and (x0,y2,z2) on k=KK+1 shell for interpolating
z1=rho0*sin(theta1);   % Since we found that k-shells are planes, we can reduce the 3D interpolation to 2D interpolation
% This is done in the kshell_tri_interp.m function
y2=rho0*cos(theta2);   % when these two interpolations are done, we can use theta1,theta2 and theta to construct a 1D interpolation
z2=rho0*sin(theta2);   % which will give us the value on that point

% plot3(x0,y0,z0,'*')
% plot3(x0,y1,z1,'r*')
% plot3(x0,y2,z2,'r*')
% xlabel('x'),ylabel('y'),zlabel('z')
%%

if (KK==65)
    d1=kshell_tri_interp(x,y,z,data,x0,y1,z1,KK);         % do the interpolation on k=KK shell, the point is (x0,y1,z1)
    d2=kshell_tri_interp(x,y,z,data,x0,y2,z2,1);       % do the interpolation on k=KK+1 shell, the point is (x0,y2,z2)
else
    d1=kshell_tri_interp(x,y,z,data,x0,y1,z1,KK);         % do the interpolation on k=KK shell, the point is (x0,y1,z1)
    d2=kshell_tri_interp(x,y,z,data,x0,y2,z2,KK);
end

if (KK<64)
    d=(d2-d1)*(theta-theta1)/(theta2-theta1)+d1;       % this is the simple 1D interpolation
else                                                  % note that when KK=64, theta1=6.22 and theta2=0.05, because k=65 shell is
    theta1=theta1-2*pi;                               % the same with k=1 shell. Thus we need to add or subtract 2pi to one of them
    d=(d2-d1)*(theta-theta1)/(theta2-theta1)+d1;       % to make the relation linear
end

%griddata3(x,y,z,data,x0,y0,z0)