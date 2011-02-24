clear all; clc;
positions=load('/thayerfs/research/lfm/entropy/positions.txt'); % load positions

x=positions(:,1);
y=positions(:,2);
z=positions(:,3);

x1=reshape(x,53,50,65);
y1=reshape(y,53,50,65);
z1=reshape(z,53,50,65);



NI=1:53; % first 10 i-shells 
NJ=1:50; % ~ night side
NK=1:65; % k=1:32, northern hemisphere

x=x1(NI,NJ,NK);
y=y1(NI,NJ,NK);
z=z1(NI,NJ,NK);
count=1;
dl=0.3;
start=300;
%finish=600;
finish=300;



for jj=start:finish
  jj
   filename=['/thayerfs/research/lfm/entropy/b10v400/',int2str(jj)]
   
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
   
   bztest=LFM3DInterp(x,y,z,bz,xstart,ystart,zstart)
      end
  end
end