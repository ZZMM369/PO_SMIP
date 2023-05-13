clc
clear all;
ps=0.25;
gm=sqrt((1-2*ps)/(2*(1-ps)));
mu=1;
lm=2*ps/(1-2*ps)*mu;
h=0.25*pi;
k=h/gm;
r1=1;
dz1=0.01;
NR=16;
NH=15;
NT=NR*NH+2;
NT1=NR*(NH-1)+2;
aa=0;NH2=(NH+1)/2;
for nh=1:NH+1
    ddz=(r1-NH2*dz1)/(NH2*(NH2-1)/2);
    if nh<=(NH+1)/2
        dh(nh)=dz1+ddz*(nh-1);
    else
        dh(nh)=dh(NH-nh+2);
    end
    aa=aa+dh(nh);
    ht(nh)=-r1+aa;
end
for nh=1:NH
    for nr=1:NR+1
        x1((nh-1)*(NR+1)+nr+1)=sqrt(1-ht(nh)^2)*cos(2*pi/NR*(nr-1));
        y1((nh-1)*(NR+1)+nr+1)=sqrt(1-ht(nh)^2)*sin(2*pi/NR*(nr-1));
        z1((nh-1)*(NR+1)+nr+1)=ht(nh);
    end
end
for nh=1:NH-1
    for nr=1:NR
        IDP=(nh-1)*(NR+1)+nr+1;
        IDS=(nh-1)*NR+nr+1;
        xm(IDS)=(x1(IDP)+x1(IDP+1)+x1(IDP+NR+1)+x1(IDP+2+NR))/4;
        ym(IDS)=(y1(IDP)+y1(IDP+1)+y1(IDP+NR+1)+y1(IDP+2+NR))/4;
        zm(IDS)=(ht(nh)+ht(nh+1))/2;
        dx1(IDS)=x1(IDP+NR+1)-x1(IDP);dy1(IDS)=y1(IDP+NR+1)-y1(IDP);
        dx2(IDS)=x1(IDP+2+NR)-x1(IDP);dy2(IDS)=y1(IDP+2+NR)-y1(IDP);
        dx3(IDS)=x1(IDP+1)-x1(IDP);dy3(IDS)=y1(IDP+1)-y1(IDP);
        D(IDS)=ht(nh+1)-ht(nh);
        nnx(IDS)=dy3(IDS)*D(IDS);nny(IDS)=-dx3(IDS)*D(IDS);nnz(IDS)=-dx1(IDS)*dy3(IDS)+dx3(IDS)*dy1(IDS);
        rr(IDS)=sqrt(nnx(IDS)^2+nny(IDS)^2+nnz(IDS)^2);
        nx(IDS)=nnx(IDS)/rr(IDS);
        ny(IDS)=nny(IDS)/rr(IDS);
        nz(IDS)=nnz(IDS)/rr(IDS);
        L1(IDS)=sqrt(dx1(IDS)^2+dy1(IDS)^2+D(IDS)^2);
        L2(IDS)=sqrt((dx2(IDS)-dx1(IDS))^2+(dy2(IDS)-dy1(IDS))^2);
        L3(IDS)=sqrt((dx3(IDS)-dx2(IDS))^2+(dy3(IDS)-dy2(IDS))^2+D(IDS)^2);
        L4(IDS)=sqrt(dx3(IDS)^2+dy3(IDS)^2);
        L5(IDS)=sqrt(dx2(IDS)^2+dy2(IDS)^2+D(IDS)^2);
        C1(IDS)=(L1(IDS)+L2(IDS)+L5(IDS))/2;C2(IDS)=(L3(IDS)+L4(IDS)+L5(IDS))/2;
        S(IDS)=sqrt(C1(IDS)*(C1(IDS)-L1(IDS))*(C1(IDS)-L2(IDS))*(C1(IDS)-L5(IDS)))+sqrt(C2(IDS)*(C2(IDS)-L3(IDS))*(C2(IDS)-L4(IDS))*(C2(IDS)-L5(IDS)));
        R(IDS)=sqrt(S(IDS)/pi);
    end
end
xm(1)=0;ym(1)=0;zm(1)=-1;nx(1)=0;ny(1)=0;nz(1)=-1;
R(1)=sqrt(r1^2-(r1-dz1)^2);S(1)=pi*R(1)^2;
xm(NT1)=0;ym(NT1)=0;zm(NT1)=1;nx(NT1)=0;ny(NT1)=0;nz(NT1)=1;
R(NT1)=sqrt(r1^2-(r1-dz1)^2);S(NT1)=pi*R(NT1)^2;
ylzz=-lm*h^2*exp(i*h*xm);
ylyy=-lm*h^2*exp(i*h*xm);
ylxx=-(lm+2*mu)*h^2*exp(i*h*xm);
ylxz=0;
ylxy=0;
ylyz=0;
ylxn=ylxx.*nx+ylxy.*ny+ylxz.*nz;
ylyn=ylxy.*nx+ylyy.*ny+ylyz.*nz;
ylzn=ylxz.*nx+ylyz.*ny+ylzz.*nz;
ylff=zeros(3*NT1,1);
for nn=1:NT1
    ylff(3*nn-2)=ylxn(nn);ylff(3*nn-1)=ylyn(nn);ylff(3*nn)=ylzn(nn);
end
NN=3*NT1;kyl=zeros(NN,NN);
for n=1:NT1
    for m=1:NT1
        if n==m
            kyl(3*n-2,3*m-2)=-1/2;kyl(3*n-2,3*m-1)=0;kyl(3*n-2,3*m)=0;kyl(3*n-1,3*m-2)=0;kyl(3*n-1,3*m-1)=-1/2;
            kyl(3*n-1,3*m)=0;kyl(3*n,3*m-2)=0;kyl(3*n,3*m-1)=0;kyl(3*n,3*m)=-1/2;
        else
            [kyll(3*n-2,3*m-2) kyll(3*n-2,3*m-1) kyll(3*n-2,3*m) kyll(3*n-1,3*m-2) kyll(3*n-1,3*m-1)...
                kyll(3*n-1,3*m) kyll(3*n,3*m-2) kyll(3*n,3*m-1) kyll(3*n,3*m)]=...
                GR1(h,xm(n),ym(n),zm(n),xm(m),ym(m),zm(m),nx(n),ny(n),nz(n));
            kyl(3*n-2,3*m-2)=S(m)*kyll(3*n-2,3*m-2);kyl(3*n-2,3*m-1)=S(m)*kyll(3*n-2,3*m-1);kyl(3*n-2,3*m)=S(m)*kyll(3*n-2,3*m);
            kyl(3*n-1,3*m-2)=S(m)*kyll(3*n-1,3*m-2);kyl(3*n-1,3*m-1)=S(m)*kyll(3*n-1,3*m-1);kyl(3*n-1,3*m)=S(m)*kyll(3*n-1,3*m);
            kyl(3*n,3*m-2)=S(m)*kyll(3*n,3*m-2);kyl(3*n,3*m-1)=S(m)*kyll(3*n,3*m-1);kyl(3*n,3*m)=S(m)*kyll(3*n,3*m);
            
        end
    end
end
cof=kyl\(-ylff);
ux=zeros(NT1,1);uy=zeros(NT1,1);uz=zeros(NT1,1);
fux=zeros(NT1,1);fuy=zeros(NT1,1);fuz=zeros(NT1,1);
uxt=zeros(NT1,1);uyt=zeros(NT1,1);uzt=zeros(NT1,1);urt=zeros(NT1,1);
for nn=1:NT1
    for mm=1:NT1
        if nn==mm
            FF1=gm^2*(R(mm)/2-i*h*R(mm)^2/3-h*R(mm)^3/9+i*h^3*R(mm)^4/24)+(3*R(mm)/2-2*i*k*R(mm)^2/3-2*k^2*R(mm)^3/9+i*k^3*R(mm)^4/24);
            FF2=gm^2*(R(mm)/2-i*h^3*R(mm)^4/24)+(-R(mm)/2+i*k^3*R(mm)^4/24);
            kwy(3*nn-2,3*mm-2)=(FF1+FF2*nx(mm)*nx(mm))/(4*mu);
            kwy(3*nn-2,3*mm-1)=FF2*nx(mm)*ny(mm)/(4*mu);
            kwy(3*nn-2,3*mm)=FF2*nx(mm)*nz(mm)/(4*mu);
            kwy(3*nn-1,3*mm-2)=FF2*ny(mm)*nx(mm)/(4*mu);
            kwy(3*nn-1,3*mm-1)=(FF1+FF2*ny(mm)*ny(mm))/(4*mu);
            kwy(3*nn-1,3*mm)=FF2*ny(mm)*nz(mm)/(4*mu);
            kwy(3*nn,3*mm-2)=FF2*nz(mm)*nx(mm)/(4*mu);
            kwy(3*nn,3*mm-1)=FF2*nz(mm)*ny(mm)/(4*mu);
            kwy(3*nn,3*mm)=(FF1+FF2*nz(mm)*nz(mm))/(4*mu);
        else
            [kwyy(3*nn-2,3*mm-2) kwyy(3*nn-2,3*mm-1) kwyy(3*nn-2,3*mm) kwyy(3*nn-1,3*mm-2) kwyy(3*nn-1,3*mm-1) kwyy(3*nn-1,3*mm) kwyy(3*nn,3*mm-2) kwyy(3*nn,3*mm-1) kwyy(3*nn,3*mm)]=...
                GR2(h,xm(nn),ym(nn),zm(nn),xm(mm),ym(mm),zm(mm));
            kwy(3*nn-2,3*mm-2)=S(mm)*kwyy(3*nn-2,3*mm-2);kwy(3*nn-2,3*mm-1)=S(mm)*kwyy(3*nn-2,3*mm-1);kwy(3*nn-2,3*mm)=S(mm)*kwyy(3*nn-2,3*mm);
            kwy(3*nn-1,3*mm-2)=S(mm)*kwyy(3*nn-1,3*mm-2);kwy(3*nn-1,3*mm-1)=S(mm)*kwyy(3*nn-1,3*mm-1);kwy(3*nn-1,3*mm)=S(mm)*kwyy(3*nn-1,3*mm);
            kwy(3*nn,3*mm-2)=S(mm)*kwyy(3*nn,3*mm-2);kwy(3*nn,3*mm-1)=S(mm)*kwyy(3*nn,3*mm-1);kwy(3*nn,3*mm)=S(mm)*kwyy(3*nn,3*mm);
        end
    end
end
for n=1:NT1
    for m=1:NT1
        ux(n)=ux(n)+cof(3*m-2)*kwy(3*n-2,3*m-2)+cof(3*m-1)*kwy(3*n-2,3*m-1)+cof(3*m)*kwy(3*n-2,3*m);
        uy(n)=uy(n)+cof(3*m-2)*kwy(3*n-1,3*m-2)+cof(3*m-1)*kwy(3*n-1,3*m-1)+cof(3*m)*kwy(3*n-1,3*m);
        uz(n)=uz(n)+cof(3*m-2)*kwy(3*n,3*m-2)+cof(3*m-1)*kwy(3*n,3*m-1)+cof(3*m)*kwy(3*n,3*m);
        
    end
    
    fuy(n)=0;
    fux(n)=i*h*exp(i*h*(xm(n)));
    fuz(n)=0;
    
    uxt(n)=((fux(n)+ux(n)))/abs(h);
    uyt(n)=((fuy(n)+uy(n)))/abs(h);
    uzt(n)=((fuz(n)+uz(n)))/abs(h);
    urt(n)=uxt(n)*nx(n)+uyt(n)*ny(n)+uzt(n)*nz(n);
end
plot(abs(urt),'r')