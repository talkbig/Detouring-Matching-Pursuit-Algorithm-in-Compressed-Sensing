close all
clear

M=128;
N=256;
tol=eps;

k_width=30;
repeats=1000;

sp_support_err=zeros(repeats,k_width);
sp_x_err=zeros(repeats,k_width);
sp_res_norm=zeros(repeats,k_width);
sp_time=zeros(repeats,k_width);

for K=1:k_width
    K
    for count=1:repeats
        rng((K-1)*repeats+count,'v5normal');
        PHI=randn(M,N);
        PHI=PHI./repmat(sqrt(sum(PHI.^2)),M,1);        
        rank=randperm(N);
        rank=rank(1:K);
        xtrue=zeros(N,1);
%         xtrue(rank)=randn(K,1);
%         xtrue(rank)=ones(K,1);
        xtrue(rank)=(1:K)';
        xtrue=xtrue+0.2.*rand(N,1);
        y=PHI*xtrue;
%         y_ref=norm(y)/sqrt(M)
%         y=y+0.2*y_ref*randn(M,1);
        tic;
        [support,x,res_norm]=SP(PHI,y,K);        
        time=toc;        
        frame=zeros(N,1);
        frame(rank)=ones(K,1);
        frame(support)=zeros(K,1);
        sp_support_err(count,K)=sum(frame);
        sp_x_err(count,K)=sum((xtrue-x).^2)/N;
        sp_res_norm(count,K)=res_norm;
        sp_time(count,K)=time;     
    end
end

sp_success=sum(sp_support_err==0,1);
sp_time_mean=sum(sp_time,1)/repeats;

figure(1)
subplot(2,1,1);
hold on
plot(sp_success,'r');
xlabel('稀疏度 K');
ylabel('准确恢复次数');
legend('sp');
grid on;
title(['M=',num2str(M),' N=',num2str(N),' repeats=',num2str(repeats)]);
subplot(2,1,2);
hold on
plot(sp_time_mean,'r');
xlabel('稀疏度 K');
ylabel('平均重构运算时间');
legend('sp');
grid on;
title(['M=',num2str(M),' N=',num2str(N),' repeats=',num2str(repeats)]);

eval(['save data_repeats-',num2str(repeats),'_',date,'-',num2str(hour(now)),'-',num2str(minute(now))]);