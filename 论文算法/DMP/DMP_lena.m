clear all;

X=imread('lena.bmp');
X1=double(X);
d=dctmtx(16);
dct=@(x)d*x*d';
B=blkproc(X1,[16,16],dct);                                                    %按每个16x16分块大小处理进行DCT变换
mask=[1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0;....
      1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0;....
      1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0;...
      1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0;...
      1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0;....
      1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0;...
      1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0;...
      1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0;...
      1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0;
      1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
      1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
      0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
      0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
      0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
      0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
      0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
      ];
   
   
   
% B2=blkproc(B,[16,16],@(x)mask.*x);   %每个分块保留左上部分系数，其他置0；这样获得稀疏信号及稀疏度K为36；
C=im2col(B,[16 16],'distinct');%每个分块排成一列，每列长度256
M=64;
N=256;
% 测量矩阵可选高斯随机矩阵或哈达玛矩阵
%高斯测量矩阵测量矩阵
phi=randn(M,N);
%哈达玛测量矩阵
% a=hadamard(N);
% q=randperm(N);
% phi=a(q(1:M),:);
%测量
Y=phi*C;
G=zeros(N,N);  %  先建立恢复矩阵的空矩阵
% OMP算法重构
t0=clock;%程序计时开始
for i=1:N  %  列循环
    i
    PSI=[phi,Y(:,i)]'*[phi,Y(:,i)];
    [~,rec,~,~]=DMP(PSI,256,32,eps);
%     rec=cs_omp(Y(:,i),phi,36);%每一列依次用OMP算法重构
   G(:,i)=rec;
end
G2=col2im(G,[16 16],[256 256],'distinct');
inv=@(x)d'*x*d;
X2=blkproc(G2,[16 16],inv);%图像分块DCT逆变换
X2=uint8(X2);
figure(1);
subplot(121);imshow(X)
xlabel('原始图形');
subplot(122);imshow(X2)
xlabel('重构图像');
%为了便于计算误差，减小数据类型对数值影响，统一转换成double类型，且norm函数对uint8 数据无效
 X=double(X);
 X2=double(X2);
disp('显示误差')
% 相对误差
%Relative_error=sqrt(sum(sum(abs(X2-X).^2)))/sqrt(sum(sum(abs(X).^2)))
 Relative_error=norm(X2-X)/norm(X) 

 %  MSE误差 也可用MSE=norm(X2-X).^2/（256*256）
MSE=sum(sum(abs(X2-X).^2))/(256*256)   
 %  PSNR
PSNR=10*log10(255*255/MSE) 
t1=clock;
t=etime(t1,t0)
%程序计时

