% a code to produce results given in Table 7
% Classical/Frequentist approach for testing Lorenz and Stochastic
% Dominance.

seed = 170;
rng(seed);
load('data_dominance_paper.mat');

%first order dominance
x = income_2002_perm;
y = income_2008_perm;
n_x = size(x,1);
m_y = size(y,1);
pool_data = [x;y];
com_x = [x ones(n_x,1)];
com_x(:,2) = com_x(:,2)./n_x;
com_x(:,3) = cumsum(com_x(:,2));
com_x(n_x,3) = 1;
com_y = [y ones(m_y,1)];
com_y(:,2) = com_y(:,2)./m_y;
com_y(:,3) = cumsum(com_y(:,2));
com_y(m_y,3) = 1;
t = 0.001:0.001:0.999;
quantile_x = quantile(x,t);
quantile_x = quantile_x';

quantile_y = quantile(y,t);
quantile_y = quantile_y';
grid_point = quantile(pool_data,t);
grid_point = grid_point';
s = size(grid_point,1);

for i = 1:s
    
    ind_x = (x <= grid_point(i,1));
    ind_y = (y <= grid_point(i,1));
    F_x(i,1) = mean(ind_x);
    F_y(i,1) = mean(ind_y);
    pg_hat_x(i,1) = mean((grid_point(i,1) - x).*ind_x);
    pg_hat_y(i,1) = mean((grid_point(i,1) - y).*ind_y);
    L_hat_x(i,1) = sum(x.*(x<=quantile_x(i,1)))/sum(x);
    L_hat_y(i,1) = sum(y.*(y<=quantile_y(i,1)))/sum(y);
end

xFSDy_sup = max(F_x - F_y);
xFSDy_stat = sqrt((n_x*m_y)/(n_x+m_y))*(xFSDy_sup);
yFSDx_sup = max(F_y - F_x);
yFSDx_stat = sqrt((n_x*m_y)/(n_x+m_y))*(yFSDx_sup);

xSSDy_sup = max(pg_hat_x - pg_hat_y);
xSSDy_stat = sqrt((n_x*m_y)/(n_x+m_y))*(xSSDy_sup);

ySSDx_sup = max(pg_hat_y - pg_hat_x);
ySSDx_stat = sqrt((n_x*m_y)/(n_x+m_y))*(ySSDx_sup);

xLDy_sup = max(L_hat_y - L_hat_x);
xLDy_stat = sqrt((n_x*m_y)/(n_x+m_y))*(xLDy_sup);

yLDx_sup = max(L_hat_x - L_hat_y);
yLDx_stat = sqrt((n_x*m_y)/(n_x+m_y))*(yLDx_sup);

p_xFSDy = exp(-2*(xFSDy_stat^2));
p_yFSDx = exp(-2*(yFSDx_stat^2));
xSSDy_accept = 0;
xLDy_accept = 0;
ySSDx_accept = 0;
yLDx_accept = 0;
m = 1000;
for j = 1:m
    for i = 1:n_x
        u1 = rand();
        index = find(u1 <= com_x(:,3),1,'first');
        boot_x(i,1) = com_x(index(1,1),1);      
    end

    for i = 1:m_y
        u1 = rand();
        index = find(u1 <= com_y(:,3),1,'first');
        boot_y(i,1) = com_y(index(1,1),1);   
    end
    
    quantile_x_boot = quantile(boot_x,t);
    quantile_x_boot = quantile_x_boot';
    
    quantile_y_boot = quantile(boot_y,t);
    quantile_y_boot = quantile_y_boot';
    
    
    for i = 1:s
        
        ind_x_boot = (boot_x <= grid_point(i,1));
        ind_y_boot = (boot_y <= grid_point(i,1));
        pg_boot_x(i,1) = mean((grid_point(i,1) - boot_x).*ind_x_boot);
        pg_boot_y(i,1) = mean((grid_point(i,1) - boot_y).*ind_y_boot);
        L_boot_x(i,1) = sum(boot_x.*(boot_x <=quantile_x_boot(i,1)))/sum(boot_x);
        L_boot_y(i,1) = sum(boot_y.*(boot_y <=quantile_y_boot(i,1)))/sum(boot_y);
    end
    
    xSSDy_boot_sup(j,1) = max((pg_boot_x - pg_hat_x) - (pg_boot_y - pg_hat_y));
    xSSDy_boot_stat(j,1) = sqrt((n_x*m_y)/(n_x+m_y))*xSSDy_boot_sup(j,1);
    
    if xSSDy_boot_stat(j,1) >=xSSDy_stat
        xSSDy_accept = xSSDy_accept+1;
    end
    
    ySSDx_boot_sup(j,1) = max((pg_boot_y - pg_hat_y) - (pg_boot_x - pg_hat_x));
    ySSDx_boot_stat(j,1) = sqrt((n_x*m_y)/(n_x+m_y))*ySSDx_boot_sup(j,1);
    
    if ySSDx_boot_stat(j,1) >=ySSDx_stat
        ySSDx_accept = ySSDx_accept+1;
    end
    
    xLDy_boot_sup(j,1) = max((L_boot_y - L_hat_y) - (L_boot_x - L_hat_x));
    xLDy_boot_stat(j,1) = sqrt((n_x*m_y)/(n_x+m_y))*(xLDy_boot_sup(j,1));
    
    if xLDy_boot_stat(j,1) >= xLDy_stat
        xLDy_accept = xLDy_accept + 1;
    end
    
    yLDx_boot_sup(j,1) = max((L_boot_x - L_hat_x) - (L_boot_y - L_hat_y));
    yLDx_boot_stat(j,1) = sqrt((n_x*m_y)/(n_x+m_y))*(yLDx_boot_sup(j,1));
    
    if yLDx_boot_stat(j,1) >= yLDx_stat
        yLDx_accept = yLDx_accept + 1;
    end
    
    
    
    
    
end

p_xSSDy = xSSDy_accept/m;
p_ySSDx = ySSDx_accept/m;

p_xLDy = xLDy_accept/m;
p_yLDx = yLDx_accept/m;


save('dominance.mat','p_xFSDy','p_yFSDx','p_xSSDy','p_ySSDx','p_xLDy','p_yLDx');