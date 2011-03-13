%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This code calculates the cusp center position based on the maximum dipole
% depression. The LFM3DInterp function only works in the SM LFM grid rather
% than GSM grid. If the grid is GSM, then LFM3DInterp might be in trouble and
% OBrambles sucks.
%                             ####Algorithm####:
% STEP 1: load grid position (x,y,z) and reshape to xyz(NI,NJ,NK)
% STEP 2: load b field from LFM and reshape, too. Only taking the first 10
%         ishells data to save a little bit computational time.
% STEP 3: find the grid point with greatest dipole depression on i=8 shell.
% STEP 4: trace the field line from that grid point to R=2.5 RE, then
%         dipole map it down to the ionosphere to calculate the cusp latitu
% STEP 5: go to STEP 2 for the next cusp position
% STEP 6: Oliver Brambles smells like cheese.
%
% INPUT: grid position x,y,z(in SM coordinate)
%        b field bx, by, bz(in SM coordinate)
% OUTPUT: cusp center latitude lat2
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;clc;
positions=load('positions.txt'); % load positions

x=positions(:,1);
y=positions(:,2);
z=positions(:,3);

x1=reshape(x,53,50,65);
y1=reshape(y,53,50,65);
z1=reshape(z,53,50,65);

x1(5,5,1)
y1(5,5,1)
z1(5,5,1)

NI=1:53; % first 10 i-shells
%NJ=20:50; % ~ night side
NJ=1:50;
NK=1:65; % k=1:32, northern hemisphere

x=x1(NI,NJ,NK);
y=y1(NI,NJ,NK);
z=z1(NI,NJ,NK);
count=1;
dl=0.3;
start=300;
%finish=600;
finish=301;



for jj=start:finish
    jj
    filename=['',int2str(jj)]
    
    %bfield=load(filename);
    % data=load('/thayerfs/research/lfm/entropy/testdata3.txt');
    data=load(filename);
    pressure=data(:,1);
    bx=data(:,2);
    by=data(:,3);
    bz=data(:,4);
    
    pressure1=reshape(pressure,53,50,65);
    bx1=reshape(bx,53,50,65);
    by1=reshape(by,53,50,65);
    bz1=reshape(bz,53,50,65);
    
    pressure=pressure1(NI,NJ,NK);
    bx=bx1(NI,NJ,NK);
    by=by1(NI,NJ,NK);
    bz=bz1(NI,NJ,NK);

    [X Y]= meshgrid(-5:-2:-60,-25:2:25);
    
    [Nx Ny]=size(X);
    
    volumetotal = zeros(Nx,Ny);
    pressuretotal = zeros(Nx,Ny);
    entrop = zeros(Nx,Ny);
    bzfailure = zeros(Nx,Ny)
    
    for j=1:Nx
        for k=1:Ny
            [j k]
            xstart=X(j,k);
            ystart=Y(j,k);
            zstart=0.1;
            
            xtrack=zeros(1,10);
            ytrack=zeros(1,10);
            ztrack=zeros(1,10);
            
            xtrack(1)=xstart;
            ytrack(1)=ystart;
            ztrack(1)=zstart;
            xtrack2(1)=xstart;
            ytrack2(1)=ystart;
            ztrack2(1)=zstart;
            
            bztest=LFM3DInterp(x,y,z,bz,xstart,ystart,zstart);
            if (bztest > 0);
                
                volume=0;
                
                press=LFM3DInterp(x,y,z,pressure,xstart,ystart,zstart);
                for i=1:1000
                    
                    bxtrack=LFM3DInterp(x,y,z,bx,xtrack(i),ytrack(i),ztrack(i));
                    bytrack=LFM3DInterp(x,y,z,by,xtrack(i),ytrack(i),ztrack(i));
                    bztrack=LFM3DInterp(x,y,z,bz,xtrack(i),ytrack(i),ztrack(i));
                    
                    btrack=sqrt(bxtrack^2+bytrack^2+bztrack^2);
                    
                    xtrack(i+1)=bxtrack/btrack*dl+xtrack(i);
                    ytrack(i+1)=bytrack/btrack*dl+ytrack(i);
                    ztrack(i+1)=bztrack/btrack*dl+ztrack(i);
                    
                    volume=dl/btrack+volume;
                    
                    radi=sqrt((xtrack(i+1)^2+ytrack(i+1)^2+ztrack(i+1)^2));
                    
                    if (radi<3.5)
                        
                        
                        break;
                    end
                    if (radi > 2*(sqrt(xstart^2+ystart^2+zstart^2))) || (sqrt(ytrack(i).^2) > 45)
                        volume=0;
                        break;
                        
                    end
                end
                
                
                
                for i=1:1000
                    
                    bxtrack=LFM3DInterp(x,y,z,bx,xtrack2(i),ytrack2(i),ztrack2(i));
                    bytrack=LFM3DInterp(x,y,z,by,xtrack2(i),ytrack2(i),ztrack2(i));
                    bztrack=LFM3DInterp(x,y,z,bz,xtrack2(i),ytrack2(i),ztrack2(i));
                    
                    btrack=sqrt(bxtrack^2+bytrack^2+bztrack^2);
                    
                    xtrack2(i+1)=bxtrack/btrack*(-dl)+xtrack2(i);
                    ytrack2(i+1)=bytrack/btrack*(-dl)+ytrack2(i);
                    ztrack2(i+1)=bztrack/btrack*(-dl)+ztrack2(i);
                    
                    volume=dl/btrack+volume;
                    
                    radi=sqrt((xtrack2(i+1)^2+ytrack2(i+1)^2+ztrack2(i+1)^2));
                    
                    if (radi<3.5)
                        
                        break;
                    end
                    
                    if (radi > 2*(sqrt(xstart^2+ystart^2+zstart^2))) || (sqrt(ytrack2(i).^2) > 45)
                        volume=0;
                        break;
                    end
                    
                    
                end
                
                
                volumetotal(j,k)=volume;
                pressuretotal(j,k)=press;
                entrop(j,k)=press*volume.^(1.6666);
            else
                entrop(j,k)=0;
                volumetotal(j,k)=0;
                pressuretotal(j,k)=0;
                bzfailures(j,k) = 1;
            end
            
        end
    end
    
    volumetotal(1:10,1)
    axd = 'end volume'
    pressuretotal(1:10,1)
    axd = 'end pressure'
    entrop(1:10,1)
    axd = 'end entrop'
    
    %entropytime(jj-start+1,:,:)=entrop(:,:);
    %volumetime(jj-start+1,:,:)=volumetotal(:,:);
    %pressuretime(jj-start+1,:,:)=pressuretotal(:,:);
end

save('/thayerfs/research/lfm/entropy/b10v400/b10v400data.mat', 'entropytime', 'volumetime', 'pressuretime');
