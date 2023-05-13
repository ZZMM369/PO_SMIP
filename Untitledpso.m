clc;clear
train111=getfield(load('input'),'input');
label=getfield(load('output'),'output');
input=train111;
[minput,ninput]=size(input);
[moutput,noutput]=size(label);
s=zeros(minput,ninput+noutput);
s(1:minput,1:ninput)=input;
s(1:minput,ninput+1:ninput+noutput)=(label(1:minput,1:noutput));
s22=s.';
NN=2000;
s2=s22(1:end,1:NN);
F=0.90;
Ltrain=floor(F*NN);
Ltest=NN-Ltrain;
input_train=s2(1:ninput,1:Ltrain);
output_train=s2(ninput+1:ninput+noutput,1:Ltrain);
input_test=s2(1:ninput,Ltrain+1:NN);
output_test=s2(ninput+1:ninput+noutput,Ltrain+1:NN);
N = size(input_test,2);
inputnum=13;
hidnum=26;
hidnum2=22;
outputnum=1;
D=inputnum*hidnum+hidnum+hidnum*hidnum2+hidnum2*outputnum+outputnum;
[inputn,inputps]=mapminmax(input_train,0,1);
[outputn,outputps]=mapminmax(output_train);
inputn_test=mapminmax('apply',input_test,inputps);
Net=newff(inputn,outputn,[hidnum,hidnum2],{'tansig','tansig','purelin'},'trainlm');
c1=1;
c2=1.49445;
maxgen=1200;
sizepop=1200;
popmax=5;popmin=-5;
Vmax=1;Vmin=-1;
wmax=0.90;
wmin=0.30;
for i=1:sizepop
    pop(i,:)=2*rands(1,2);
    V(i,:)=rand(1,2);
    fitness(i)=fun(pop(i,:));
end
[bestfitness bestindex]=min(fitness);
[worstfitness,worstindex]=max(fitness);
zbest=pop(bestindex,:);
gbest=pop;%
fitnessgbest=fitness;
fitnesszbest=bestfitness;
 for i=1:maxgen
    for j=1:sizepop
        w1=wmax-(wmax-wmin)*(fitness(j)-worstfitness)/(bestfitness-worstfitness);
        V(j,:)=V(j,:)+c1*rand*(gbest(j,:)-pop(j,:))+c2*rand*(zbest-pop(j,:));
        V(j,find(V(j,:)>Vmax))=Vmax;
        V(j,find(V(j,:)<Vmax))=Vmin;
        pop(j,:)=pop(j,:)+0.5*V(j,:);  
        pop(j,find(pop(j,:)>popmax))=popmax;
        pop(j,find(pop(j,:)<popmin))=popmin;
        fitness(j)=fun(pop(j,:));
    end
    for j=1:sizepop
        if fitness(j)<fitnessgbest(j)
            gbest(j,:)=pop(j,:);
            fitnessgbest(j)=fitness(j);
        end
        if fitness(j)<fitnesszbest
            zbest=pop(j,:);
            fitnesszbest=fitness(j);
        end
    end
    yy(i)=fitnesszbest;
 end
D=inputnum*hidnum+hidnum+hidnum*hidnum2+hidnum2+hidnum2*outputnum+outputnum;
w1=pop(1:inputnum*hidnum);
b1=pop(inputnum*hidnum+1:inputnum*hidnum+hidnum);
w2=pop(inputnum*hidnum+hidnum+1:inputnum*hidnum+hidnum+hidnum*hidnum2);
b2=pop(inputnum*hidnum+hidnum+hidnum*hidnum2+1:inputnum*hidnum+hidnum+hidnum*hidnum2+hidnum2);
w3=pop(inputnum*hidnum+hidnum+hidnum*hidnum2+hidnum2+1:inputnum*hidnum+hidnum+hidnum*hidnum2+hidnum2+hidnum2*outputnum);
b3=pop(inputnum*hidnum+hidnum+hidnum*hidnum2+hidnum2+hidnum2*outputnum+1:inputnum*hidnum+hidnum+hidnum*hidnum2+hidnum2+hidnum2*outputnum+outputnum);
Net.iw{1,1}=reshape(w1,hidnum,inputnum);
Net.lw{2,1}=reshape(w2,hidnum2,hidnum);
Net.lw{3,2}=reshape(w3,outputnum,hidnum2);
Net.b{1}=reshape(b1,hidnum,1);
Net.b{2}=reshape(b2,hidnum2,1);
Net.b{3}=reshape(b3,outputnum,1);
Net.trainParam.epochs =100;
Net.trainParam.goal=0.0009;
Net.trainParam.lr=0.001;
Net.trainParam.max_fail=500;
[Net,per2]=train(Net,inputn,outputn);
bn=sim(Net,inputn);
test_simu0=mapminmax('reverse',bn,outputps);
an=sim(Net,inputn_test);
test_simu=mapminmax('reverse',an,outputps);
% save('Net');
% save('bp.mat','Net');





    