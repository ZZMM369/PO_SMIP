function [aaamaxa,wD] = Duhamel(time,gd,ti,Damp)
m=1;kesi=Damp;
ug=gd';
n1=length(ug);
dt=time(2)-time(1);
wn=2*pi/ti;
wD=wn*sqrt(1-kesi^2);
u0=0;
v0=0;
aa0=0;
n2=n1+300;
u(1,1)=u0;
v(1,1)=v0;
aa(1,1)=aa0;
P(1,1)=-m*ug(1,1);
aaa(1,1)=aa(1,1)+ug(1,1);
T2(1,1)=0;
for i=2:n2
    T2(1,i)=(i-1)*dt;
    if i<=n1
        P(1,i)=-m*ug(1,i);
    else
        ug(1,i)=0;
        P(1,i)=-m*ug(1,i);
    end
    for j=2:i
        ut(1,j)=P(1,j)/(m*wD)*dt*exp(-kesi*wn*(T2(1,i)-T2(1,j)))*sin(wD*(T2(1,i)-T2(1,j)));
    end
    u(1,i)=sum(ut);
end
for i=1:n2-1
    T3(1,i)=(i-1)*dt;
    v(1,i)=(u(1,i+1)-u(1,i))/dt;
end
for i=1:n2-2
    T4(1,i)=(i-1)*dt;
    aa(1,i)=(v(1,i+1)-v(1,i))/dt;
    aaa(1,i)=aa(1,i)+ug(1,i);
end
aaamaxa=max(abs(aaa));
end