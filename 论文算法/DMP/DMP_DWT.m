clc
clear

K=64;
M=128;
N=256;
X_init=imread('lena.bmp');
X_init=double(X_init);
[N,repeats]=size(X_init);

W=DWT(N);

X_dwt=full(W*sparse(X_init)*W');
% X_dwt1=X_dwt;
% X_dwt=zeros(N,N);
% X_dwt(1:K,1:K)=X_dwt1(1:K,1:K);
% X_recover=full(W'*X_dwt*W);

X_dwt=im2col(X_dwt,[16,16],'distinct');
X_dwt=X_dwt';


PHI=randn(M,N);

Y=PHI*X_dwt;

X_cs=zeros(N,repeats);

for t=1:repeats
    t
    
    PSI=[PHI,Y(:,t)]'*[PHI,Y(:,t)];
    [~,x_cs,~,~]=DMP(PSI,256,64,eps);

    X_cs(:,t)=x_cs;
end

X_cs=X_cs';
X_cs=col2im(X_cs,[16,16],[256,256],'distinct');

X_recover=full(W'*sparse(X_cs)*W);

%  误差(PSNR)
errorx=sum(sum(abs(X_recover-X_init).^2));        %  MSE误差
psnr=10*log10(255*255/(errorx/N/repeats));   %  PSNR

%原始图像
figure(1);
% subplot(1,2,1);
% imshow(uint8(X_init));
% title('原始图像');
% 
% subplot(1,2,2);
imshow(uint8(X_recover));
title(['恢复的图像  信噪比：',num2str(psnr)]);


%  变换图像
figure(2);
imshow(uint8(X_dwt));
title('小波变换后的图像');
