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