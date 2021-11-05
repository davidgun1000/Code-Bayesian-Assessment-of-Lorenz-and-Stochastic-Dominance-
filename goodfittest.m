%goodness of fit test to obtain table 2
%please uncomment for specific model.

load('data_dominance_paper.mat');
load('mixgam3_2008_pred_whole.mat');
x = income_2008_perm;
m = 5000;
table = tabulate(x);
n = size(table(:,1),1);
table(:,3) = table(:,2)./n;
table(:,4) = cumsum(table(:,3));
cum_pred = table(:,4); 
x = table(:,1);

%%
%for gamma2 mixture
for j = 1:m
   cum_dist = W1(j,1)*gamcdf(x,V1(j,1),M1(j,1)/V1(j,1)) + ...
       W2(j,1)*gamcdf(x,V2(j,1),M2(j,1)/V2(j,1)); 
   log_lik_vector(j,1) = sum(log(cum_dist));
   cum_dist_matrix(j,:) = (cum_dist)';     
end
cum_dist_est = mean(cum_dist_matrix);
cum_dist_est = cum_dist_est';
RMSE = sqrt(mean((cum_dist_est - cum_pred).^2))
MAE = mean(abs(cum_dist_est - cum_pred))
d = max(abs(cum_dist_est - cum_pred))
%%
%for gamma3 mixture
% for j = 1:m
%     j
%    cum_dist = mean(W1(:,1))*gamcdf(x,mean(V1(:,1)),mean(M1(:,1))/mean(V1(:,1))) + ...
%        mean(W2(:,1))*gamcdf(x,mean(V2(:,1)),mean(M2(:,1))/mean(V2(:,1))) + ...
%        mean(W3(:,1))*gamcdf(x,mean(V3(:,1)),mean(M3(:,1))/mean(V3(:,1)));
%    cum_dist_matrix(j,:) = (cum_dist)';     
% end
% cum_dist_est = mean(cum_dist_matrix);
% cum_dist_est = cum_dist_est';
% 
% RMSE = sqrt(mean((cum_dist_est - cum_pred).^2))
% MAE = mean(abs(cum_dist_est - cum_pred))
% d = max(abs(cum_dist_est - cum_pred))
%%
%for gamma4 mixture
% for j = 1:m
%    cum_dist = W1(j,1)*gamcdf(x,V1(j,1),M1(j,1)/V1(j,1)) + ...
%        W2(j,1)*gamcdf(x,V2(j,1),M2(j,1)/V2(j,1)) + ...
%        W3(j,1)*gamcdf(x,V3(j,1),M3(j,1)/V3(j,1)) + ...
%        W4(j,1)*gamcdf(x,V4(j,1),M4(j,1)/V4(j,1));
%    log_lik_vector(j,1) = sum(log(cum_dist));
%    cum_dist_matrix(j,:) = (cum_dist)';
% end
% cum_dist_est = mean(cum_dist_matrix);
% cum_dist_est = cum_dist_est';
% sum_log_score = sum(log(cum_dist_est));
% sum_log_score = num2str(sum_log_score,'%.4f')
% max_log_lik = max(log_lik_vector);
% BIC = max_log_lik - (11/2)*log(n);
% BIC_score = num2str(BIC,'%.4f')
% RMSE = sqrt(mean((cum_dist_est - cum_pred).^2))
% MAE = mean(abs(cum_dist_est - cum_pred))
% d = max(abs(cum_dist_est - cum_pred))
%%
%for gamma 5 mixtuer
% for j = 1:m
%    cum_dist = W1(j,1)*gamcdf(x,V1(j,1),M1(j,1)/V1(j,1)) + ...
%        W2(j,1)*gamcdf(x,V2(j,1),M2(j,1)/V2(j,1)) + ...
%        W3(j,1)*gamcdf(x,V3(j,1),M3(j,1)/V3(j,1)) + ...
%        W4(j,1)*gamcdf(x,V4(j,1),M4(j,1)/V4(j,1)) + ...
%        W5(j,1)*gamcdf(x,V5(j,1),M5(j,1)/V5(j,1));
%    log_lik_vector(j,1) = sum(log(cum_dist));
%    cum_dist_matrix(j,:) = (cum_dist)';
% end
% cum_dist_est = mean(cum_dist_matrix);
% cum_dist_est = cum_dist_est';
% sum_log_score = sum(log(cum_dist_est));
% sum_log_score = num2str(sum_log_score,'%.4f')
% max_log_lik = max(log_lik_vector);
% BIC = max_log_lik - (11/2)*log(n);
% BIC_score = num2str(BIC,'%.4f')
% RMSE = sqrt(mean((cum_dist_est - cum_pred).^2))
% MAE = mean(abs(cum_dist_est - cum_pred))
% d = max(abs(cum_dist_est - cum_pred))
%%
%for dagum
% b =     210.6517;
% a =      2.3119;
% p =       2.5380;
% cum_dist_est = dagcdf(x,b,a,p);
% RMSE = sqrt(mean((cum_dist_est - cum_pred).^2))
% MAE = mean(abs(cum_dist_est - cum_pred))
% d = max(abs(cum_dist_est - cum_pred))
%for singh
% b = 289.8104;
% a =  3.4725;
% q =  0.6435;
% cum_dist_est = singhcdf(x,b,a,q);
% RMSE = sqrt(mean((cum_dist_est - cum_pred).^2))
% MAE = mean(abs(cum_dist_est - cum_pred))
%d = max(abs(cum_dist_est - cum_pred))
%%
%for gam1
% m =   454.0787;
% v =   2.5294;
% cum_dist_est = gamcdf(x,v,m/v);
% RMSE = sqrt(mean((cum_dist_est - cum_pred).^2))
% MAE = mean(abs(cum_dist_est - cum_pred))
% d = max(abs(cum_dist_est - cum_pred))
%%
