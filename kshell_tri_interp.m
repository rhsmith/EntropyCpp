function d=kshell_tri_interp(x,y,z,data,x0,y1,z1,KK)

% transfer a 3D K-shell to 2D plane, first KK shell
% (x,y,z) to (p,q)
[NI,NJ,NK]=size(x);

p=zeros(NI,NJ);
q=zeros(NI,NJ);
data1=zeros(NI,NJ);
d=0;

for i=1:NI
    for j=1:NJ
        p(i,j)=x(i,j,KK);
        q(i,j)=sqrt(y(i,j,KK)^2+z(i,j,KK)^2);
        p1=x0;
        q1=sqrt(y1^2+z1^2);
        data1(i,j)=data(i,j,KK);
    end
end


% find out which triangle is (p1.q1) in
% search through (i,j) pairs, each cell is divided to two trangles:
% 1 (i,j) (i+1,j),(i,j+1)
% 2 (i+1,j+1) (i+1,j),(i,j+1)
%
for i=1:NI-1
    for j=1:NJ-1
        s1=[p(i,j)-p1 q(i,j)-q1];
        s2=[p(i+1,j)-p1 q(i+1,j)-q1];                 
        s3=[p(i+1,j+1)-p1 q(i+1,j+1)-q1];  
        s4=[p(i,j+1)-p1 q(i,j+1)-q1]; 

%         if (norm(s1)==0)
%             d=data1(i,j);
%             break;
%         end
%         if (norm(s2)==0)
%             d=data1(i+1,j);
%             break;
%         end        
%         if (norm(s3)==0)
%             d=data1(i+1,j+1);
%             break;
%         end     
%         if (norm(s4)==0)
%             d=data1(i,j+1);
%             break;
%         end           
        
        % triangle 1, ANG(12)+ANG(24)+ang(41)=2*pi
        theta12=acos((s1(1)*s2(1)+s1(2)*s2(2))/sqrt((s1(1)^2+s1(2)^2)*(s2(1)^2+s2(2)^2)));
        theta24=acos((s2(1)*s4(1)+s2(2)*s4(2))/sqrt((s2(1)^2+s2(2)^2)*(s4(1)^2+s4(2)^2)));        
        theta41=acos((s4(1)*s1(1)+s4(2)*s1(2))/sqrt((s4(1)^2+s4(2)^2)*(s1(1)^2+s1(2)^2)));           
               
        if (abs(theta12+theta24+theta41-2*pi)<0.001)

            xx1=p(i,j);
            yy1=q(i,j);
            ff1=data1(i,j);
            xx2=p(i+1,j);
            yy2=q(i+1,j);
            ff2=data1(i+1,j);
            xx3=p(i,j+1);
            yy3=q(i,j+1);  
            ff3=data1(i,j+1);
           
            break;
        end
        
        % triangle 2, ANG(23)+ANG(34)+ang(42)=2*pi
        theta23=acos((s2(1)*s3(1)+s2(2)*s3(2))/sqrt((s2(1)^2+s2(2)^2)*(s3(1)^2+s3(2)^2)));
        theta34=acos((s3(1)*s4(1)+s3(2)*s4(2))/sqrt((s3(1)^2+s3(2)^2)*(s4(1)^2+s4(2)^2)));        
        theta42=acos((s4(1)*s2(1)+s4(2)*s2(2))/sqrt((s4(1)^2+s4(2)^2)*(s2(1)^2+s2(2)^2)));           
        
        if (abs(theta23+theta34+theta42-2*pi)<0.001)

            xx1=p(i+1,j+1);
            yy1=q(i+1,j+1);
            ff1=data1(i+1,j+1);
            xx2=p(i+1,j);
            yy2=q(i+1,j);
            ff2=data1(i+1,j);
            xx3=p(i,j+1);
            yy3=q(i,j+1);            
            ff3=data1(i,j+1);
            

            break;
        end        
        
    end
end

%             figure
%             pcolor(p,q,data1);hold on
%             plot(xx1,yy1,'ro');hold on
%             plot(xx2,yy2,'ro');            
%             plot(xx3,yy3,'ro');
%             plot(p1,q1,'b*');hold off   
%             xlabel('p'),ylabel('q');
%             xlim([p1-10 p1+10]);ylim([q1-10 q1+10]);
            
arr1=[xx1 yy1 1]';
arr2=[xx2 yy2 1]';
arr3=[xx3 yy3 1]';
arr= [p1 q1 1]';
d=(ff1*det([arr arr2 arr3]')-ff2*det([arr arr1 arr3]')+ff3*det([arr arr1 arr2]'))/det([arr1 arr2 arr3]');
%g=griddata(p,q,data1,p1,q1);
%g=d;
