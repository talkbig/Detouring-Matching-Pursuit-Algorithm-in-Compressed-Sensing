close all
clear


% load common_date-randn-N=256-M=128-k_width=64-repeats=1000
% load common_date-digit-N=256-M=128-k_width=64-repeats=1000
load common_date-zeros_ones-N=256-M=128-k_width=64-repeats=1000

M=128;
N=256;
tol=eps;

k_width=30;
repeats=100;

dmp_support_err=zeros(repeats,k_width);
dmp_x_err=zeros(repeats,k_width);
dmp_res=zeros(repeats,k_width);
dmp_sf=zeros(repeats,k_width);
dmp_time=zeros(repeats,k_width);

for K=1:k_width
    for count=1:repeats
        PHI=common_PHI(:,:,count);       
        xtrue=common_xtrue(:,K,count);

        xtrue=0.08.*randn(N,1)+xtrue;
        y=PHI*xtrue;
%         y=common_y(:,K,count);
%         y=y+norm(y)/M.*randn(M,1);
        tic;
        PSI=[PHI,y]'*[PHI,y];
        [support,x,res,sf]=DMP(PSI,N,K,tol); 
        time=toc;        
        frame=zeros(N,1);
        rank=common_support(1:K,K,count);
        frame(rank)=ones(K,1);
        frame(support)=zeros(K,1);
        dmp_support_err(count,K)=sum(frame);
        dmp_x_err(count,K)=sum((xtrue-x).^2);
        dmp_res(count,K)=res;
        dmp_sf(count,K)=sf;
        dmp_time(count,K)=time;     
    end
    K
end

dmp_success=sum(dmp_support_err==0,1);
dmp_time_mean=sum(dmp_time,1)/repeats;

figure(1)
subplot(2,1,1);
hold on
plot(dmp_success,'r');
xlabel('稀疏度 K');
ylabel('准确恢复次数');
legend('dmp');
grid on;
title(['M=',num2str(M),' N=',num2str(N),' repeats=',num2str(repeats)]);
subplot(2,1,2);
hold on
plot(dmp_time_mean,'r');
xlabel('稀疏度 K');
ylabel('平均重构运算时间');
legend('dmp');
grid on;
title(['M=',num2str(M),' N=',num2str(N),' repeats=',num2str(repeats)]);

clearvars common_PHI common_xtrue common_y
eval(['save data_repeats-',num2str(repeats),'_',date,'-',num2str(hour(now)),'-',num2str(minute(now))]);