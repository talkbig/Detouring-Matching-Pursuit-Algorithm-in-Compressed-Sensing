close all
clear

M=128;
N=256;
tol=eps;
PHI=randn(M,N);
PHI=PHI./repmat(sqrt(sum(PHI.^2)),M,1);    

k_width=30;
repeats=1000;

omp_support_err=zeros(repeats,k_width);
omp_x_err=zeros(repeats,k_width);
omp_time=zeros(repeats,k_width);

for K=1:k_width
    for count=1:repeats
        rng((K-1)*repeats+count,'v5normal');    
        rank=randperm(N);
        rank=rank(1:K);
        xtrue=zeros(N,1);
%         xtrue(rank)=randn(K,1);
        xtrue(rank)=(1:K)';
        xtrue=xtrue+0.2.*rand(N,1);
        y=PHI*xtrue;
        tic;
        [support,x]=OMP(PHI,y,K);     
        time=toc;         
        frame=zeros(N,1);
        frame(rank)=ones(K,1);
        frame(support)=zeros(K,1);
        omp_support_err(count,K)=sum(frame);
        omp_x_err(count,K)=sum((xtrue-x).^2)/N;
        omp_time(count,K)=time;     
    end
    K
end

omp_success=sum(omp_support_err==0,1);
omp_time_mean=sum(omp_time,1)/repeats;

figure(1)
subplot(2,1,1);
hold on
plot(omp_success,'r');
xlabel('稀疏度 K');
ylabel('准确恢复次数');
legend('OMP');
title(['M=',num2str(M),' N=',num2str(N),' repeats=',num2str(repeats)]);
grid on;
subplot(2,1,2);
hold on
plot(omp_time_mean,'r');
xlabel('稀疏度 K');
ylabel('平均重构运算时间');
legend('OMP');
title(['M=',num2str(M),' N=',num2str(N),' repeats=',num2str(repeats)]);
grid on

eval(['save data_repeats-',num2str(repeats),'_',date,'-',num2str(hour(now)),'-',num2str(minute(now))]);