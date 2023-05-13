
function [g11 g12 g13 g21 g22 g23 g31 g32 g33]=GR1(h,zx1,zy1,zz1,zx2,zy2,zz2,nn1,nn2,nn3)
ps=0.25;
xs=sqrt((1-2*ps)/(2*(1-ps)));
mu=1;
q=h;
k=q/xs;
x1(1)=zx1;x1(2)=zy1;x1(3)=zz1;x2(1)=zx2;x2(2)=zy2;x2(3)=zz2;
 n(1)=nn1;n(2)=nn2;n(3)=nn3;
 r=sqrt((x1(1)-x2(1))^2+(x1(2)-x2(2))^2+(x1(3)-x2(3))^2);
 
for ii=1:3
   gm(ii)=(x1(ii)-x2(ii))/r;
end
A1(1)=0;A1(2)=0;A1(3)=-i;
A2(1)=-i*xs;A2(2)=i*(2*xs^3-xs);A2(3)=0;
B1(1)=4;B1(2)=-2;B1(3)=-3;
B2(1)=-4*xs^2-1;B2(2)=4*xs^2-1;B2(3)=2*xs^2;
C1(1)=-i*12;C1(2)=i*6;C1(3)=i*6;
C2(1)=i*12*xs;C2(2)=-i*6*xs;C2(3)=-i*6*xs;
D1(1)=-12;D1(2)=6;D1(3)=6;
D2(1)=12;D2(2)=-6;D2(3)=-6;
g=(k*r*A1+B1+C1/(k*r)+D1/(k*r)^2)*exp(-i*k*r)+(k*r*A2+B2+C2/(k*r)+D2/(k*r)^2)*exp(-i*q*r);
f1=(xs^2)*(1-i*2/(q*r)-2/(q*r)^2)*exp(-i*q*r)+(i*2/(k*r)+2/(k*r)^2)*exp(-i*k*r);
f2=(xs^2)*(i/(q*r)+1/(q*r)^2)*exp(-i*q*r)+(1-i/(k*r)-1/(k*r)^2)*exp(-i*k*r);
dt=zeros(3,3);
for ii=1:3
    dt(ii,ii)=1;
    for jj=1:3
GRT(ii,jj)=((g(1)-g(2)-2*g(3))*gm(ii)*gm(jj)*(gm(1)*n(1)+gm(2)*n(2)+gm(3)*n(3))+g(3)*gm(ii)*n(jj)+g(2)*gm(jj)*n(ii)...
    +g(3)*(gm(1)*n(1)+gm(2)*n(2)+gm(3)*n(3))*dt(ii,jj))/(4*pi*r^2);

    end
end
g11=GRT(1,1);g12=GRT(1,2);g13=GRT(1,3);
g21=GRT(2,1);g22=GRT(2,2);g23=GRT(2,3);
g31=GRT(3,1);g32=GRT(3,2);g33=GRT(3,3);






