
clear  all
close all
clc
format compact;

%-------relax algorithm估计正弦成分的幅度和频率=========

%=====================
fs=1024;  
t=0:1/fs:1-1/fs;  
%r=2;  xn=10*cos(2*pi*50*t)+20*cos(2*pi*150*t)+randn(size(t));  
r=1;  xn=10*exp(1j*2*pi*50*t)+20*exp(1j*2*pi*150*t)+15*exp(1j*2*pi*160*t)+25*exp(1j*2*pi*200*t)+3*randn(size(t));
%======================= 
N=fs;  
KK=4;  %the number of sinusoidals

w=zeros(1,KK);    %频率的估计
b=zeros(1,KK);   %幅度的估计

y(1,:)=xn;  
yyk=fftshift(fft(y(1,:),N));

[ywk(1),m]=max(abs(yyk));   %类似周期图法计算第一个正弦信号成分
w(1)=-pi+(m-1)*2*pi/N;      %这个地方似乎不应该减1，角频率的估计

figure();
plot(abs(yyk));
[m w(1)/2/pi*fs]

b(1)=y(1,:)*exp(-1j*w(1)*(0:N-1)')/N;  %幅度的估计


%==============================
for K=2:KK  
    
    c=0; cast=10;  Th=0.01;  g=[];   
 
  while  (cast>=Th)
      
    sum=zeros(1,N);
    for ii=1:N
        for i=1:K-1
            sum(ii)=sum(ii)+b(i)*exp(1j*w(i)*(ii-1));   %假设有K个成分，计算除第K个成分外的其他成分的和
        end
    end
    y(K,:)=xn-sum;    %计算第K个成分
    yyk=fftshift(fft(y(K,:),N));
    [~,m]=max(abs(yyk));
    w(K)=-pi+(m-1)*2*pi/N;   %对第K个成分频率的估计
    b(K)=y(K,:)*exp(-1j*w(K)*(0:N-1)')/N;   %对第K个成分幅度的估计
    
    %=======重新估计1->K-1的成分=========
    for k1=1:K-1
        K2=1:K;   K2(k1)=[];    %删除第k1个元素
        sum=zeros(1,N);
        for ii=1:N
            for i=1:length(K2)
                sum(ii)=sum(ii)+b(K2(i))*exp(1j*w(K2(i))*(ii-1));
            end
        end
        y(k1,:)=xn-sum;
        yyk=fftshift(fft(y(k1,:),N));
        [~,m]=max(abs(yyk));
        w(k1)=-pi+(m-1)*2*pi/N;
        b(k1)=y(k1,:)*exp(-1j*w(k1)*(0:N-1)')/N; 
    end
    
    %======迭代更新，直到收敛===========
    sum=zeros(1,N);
    for ii=1:N
        for i=1:K
            sum(ii)=sum(ii)+b(i)*exp(1j*w(i)*(ii-1));  
        end
    end
    v=xn-sum;  g=[g v*v'];
    if c>0
        cast=g(end)-g(end-1);
    end
    c=c+1;
  end
  c
  g
    
    
end  

freq=w/2/pi*N;  
ample=r*abs(b);  %如果只有实部的话，幅度的估计值要乘以2，否则不乘

figure();
stem(freq,ample);   
title('RELAX');  
xlabel('frequency/Hz');ylabel('amplitude'); 
grid on
%===================
[freq; ample]
 