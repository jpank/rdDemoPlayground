close all

figure(1)
hold on
grid on
title('Reconstrucion time CS algortithms vs. input signal SNR')
plot(mean(aiht_time_vec), '-bx')
plot(mean(cvx_time_vec), '-rx')
plot(mean(spgl1_time_vec), '-gx')
legend('AIHT', 'CVX', 'SPGL1')
set(gca,'xtick',1:size(recSNR_aihtVec,2));
set(gca,'xticklabel',{'3', '4', '5', '6', '7'})
xlabel('Initial SNR [dB]')
ylabel('Time in [s]')




figure(2)
title('Compressed Sensing Reconstruction SNR vs. Input Signal SNR')
plot(mean(recSNR_aihtVec), '-bx')
hold on
plot(mean(recSNR_cvxVec), '-rx')
plot(mean(recSNR_spgl1Vec), '-gx')
plot(mean(inSNRvec), '-magentax')
legend('AIHT', 'CVX', 'SPGL1', 'input')
ylabel('Output SNR [dB]')
xlabel('Initial SNR [dB]')
%set(gca,'xtick',2:7)
% set(gca,'xlim',[3,7])
set(gca,'xtick',1:size(recSNR_aihtVec,2));
%  set(gca,'xticklabel',{'apple','banana'});
set(gca,'xticklabel',{'3', '4', '5', '6', '7'})
grid on

tm_aiht = mean(mean(aiht_time_vec));
tm_cvx = mean(mean(cvx_time_vec));
tm_spgl1 = mean(mean(spgl1_time_vec));

fprintf('--------------------------------------\n');
fprintf('|         Reconstruction times       |\n');
fprintf('|  AIHT [s] |   CVX [s]  | SPGL1 [s] |\n');
fprintf('--------------------------------------\n');
fprintf('|   %3.2f    |    %3.2f   |   %3.2f   |\n', tm_aiht, tm_cvx, tm_spgl1 );
fprintf('--------------------------------------\n');

fprintf('------------------------------------------------\n');
fprintf('|         Reconstrucion SNR comparison         |\n');
fprintf('|  input  |  AIHT [s] |   CVX [s]  | SPGL1 [s] |\n');
fprintf('------------------------------------------------\n');
fprintf('|  %3.2f  |  %3.2f    |    %3.2f   |   %3.2f   |\n', [mean(inSNRvec); mean(recSNR_aihtVec); mean(recSNR_cvxVec); mean(recSNR_spgl1Vec)]);
fprintf('-----------------------------------------------\n');



