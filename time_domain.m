function [time,utt6py,peer_peaky,Periody,timey,accelerationy] = time_domain(frequency_spectrum_value)
path='C:\Users...\'; 
Filesname = dir(strcat(path,'*.xls'));
Length = length(Filesname ); 
utt6py=[];
peer_peaky=[];
Periody=[];
timey=[];
accelerationy=[];
for i = 1:Length 
    data_num = xlsread(strcat(path,Filesname (i).name));
    Period=data_num(:,3);
    peer_peak = data_num(:,4);
    time=data_num(:,1);
    acceleration =  data_num(:,2);
    length0=size(acceleration,1);
    maxacceleration=max(abs(acceleration));
    TFxs=0.2/maxacceleration;
    acceleration1=TFxs*acceleration;
    yangben=300;
    IT=0.0005*(1:1:yangben);
    IT1=IT';
    PY=frequency_spectrum_value';
    sh1=PY(1:yangben,100);
    sh0=[sh1];sh=sh0';
    fttp=zeros(length0,1);
    a1=1*acceleration1;a1(length0+1:length0+round(1/3*length0),1)=0;
    lenl=length(a1);
    t0=0.01;
    for n3=1:yangben
        fttp(n3,1)=0;
        for n2=1:lenl
            ta=t0*n2;
            fttp(n3,1)=fttp(n3,1)+a1(n2)*exp(-1i*2*pi*n3*n2/lenl);
        end
        fttp1(n3,1)=fttp(n3,1)/lenl;
    end
    for jj=1:yangben
        utt2(jj,1)=fttp1(jj,1)/(-(2*pi*jj*0.00125)^2);
        utt3py(:,jj)=sh(:,jj)*utt2(jj,1)*(-(2*pi*jj*0.00125)^2);
    end
    for n2=1:length0
        for nn=1:1
            utt4py(n2,nn)=0;
            for jj=1:yangben
                ta=t0*n2;
                tx(n2)=ta;
                utt4py(n2,nn)=1*(utt4py(n2,nn)+utt3py(nn,jj)*exp(1i*2*pi*jj*n2/lenl));
            end
        end
    end
    utt5py=(1/2*pi)*real(utt4py)*2;
    utt6py=[utt6py,utt5py];
    peer_peaky=[peer_peaky,peer_peak];
    Periody=[Periody,Period];
    timey=[timey,time];
    accelerationy=[accelerationy,acceleration1];
end