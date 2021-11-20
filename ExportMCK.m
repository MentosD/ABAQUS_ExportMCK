clear;
cd('C:\Users\BJUT\Desktop\MCK_to_Matlab');          %Abaqus文件目录
filename = 'Job-1-1_MASS2.mtx';
[data1,data2,data3,data4,date5]=textread(filename,'%n%n%n%n%n','delimiter', ',');
data=[data1 data2 data3 data4 date5];
OnePointDOF=3;                                      %一个节点3个自由度
M = zeros(length(data1),length(data1));
K = zeros(length(data1),length(data1));
C = zeros(length(data1),length(data1));
for n=1:length(data1)
X=(data(n,1)-1)*OnePointDOF+(data(n,2));
Y=(data(n,3)-1)*OnePointDOF+(data(n,4));
M(X,Y)=data(n,5);
M(Y,X)=data(n,5);
end
filename = 'Job-1-1_STIF2.mtx';
[data1,data2,data3,data4,date5]=textread(filename,'%n%n%n%n%n','delimiter', ',');
data=[data1 data2 data3 data4 date5];
for n=1:length(data1)
X=(data(n,1)-1)*OnePointDOF+(data(n,2));
Y=(data(n,3)-1)*OnePointDOF+(data(n,4));
K(X,Y)=data(n,5);
K(Y,X)=data(n,5);
end
% 模态分析
A = pinv(M)* K;
[V,D] = eig(A);
la = diag(D);
ww = sqrt(la);
w = sort(ww);
% Rayleigh阻尼
w1 = w(1);
w2 = w(2);
w1 = 3.81;
w2 = 3.94;
kxis=0.05;
a0=2*kxis/(w(2)+w(1))*w(1)*w(2);
a1=2*kxis/(w(2)+w(1));
C=a0*M+a1*K;