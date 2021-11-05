% a main program for assessing Lorenz and stochastic dominance
load('gamma4_1999.mat'); % load the draws for mixture of gamma
load('gamma4_2002.mat');
drawsx = [W1_1999 W2_1999 W3_1999 W4_1999 M1_1999 M2_1999 M3_1999 M4_1999 V1_1999 V2_1999 V3_1999 V4_1999];
drawsy = [W1_2002 W2_2002 W3_2002 W4_2002 M1_2002 M2_2002 M3_2002 M4_2002 V1_2002 V2_2002 V3_2002 V4_2002];
%parameter for dist x
W1_x = drawsx(:,1);
W2_x = drawsx(:,2);
W3_x = drawsx(:,3);
W4_x = drawsx(:,4);
M1_x = drawsx(:,5);
M2_x = drawsx(:,6);
M3_x = drawsx(:,7);
M4_x = drawsx(:,8);
V1_x = drawsx(:,9);
V2_x = drawsx(:,10);
V3_x = drawsx(:,11);
V4_x = drawsx(:,12);
%parameter for distribution y
W1_y = drawsy(:,1);
W2_y = drawsy(:,2);
W3_y = drawsy(:,3);
W4_y = drawsy(:,4);
M1_y = drawsy(:,5);
M2_y = drawsy(:,6);
M3_y = drawsy(:,7);
M4_y = drawsy(:,8);
V1_y = drawsy(:,9);
V2_y = drawsy(:,10);
V3_y = drawsy(:,11);
V4_y = drawsy(:,12);
%mean_x = W1_x.*M1_x + W2_x.*M2_x + W3_x.*M3_x;
%mean_y = W1_y.*M1_y + W2_y.*M2_y + W3_y.*M3_y;

m = size(W1_x,1);
% m=10;
p = (0.001:0.001:0.999)';
n = size(p,1);
for i = 1:m
    i
    u1 = rand(200000,1);
    %get the observations for x
    z1_x = (u1<=W1_x(i,1));
    z2_x = (u1>W1_x(i,1)) & (u1<=(W1_x(i,1)+W2_x(i,1)));
    z3_x = (u1>(W1_x(i,1)+W2_x(i,1)) & u1<=(W1_x(i,1)+W2_x(i,1)+W3_x(i,1)));
    z4_x = (u1>(W1_x(i,1)+W2_x(i,1)+W3_x(i,1))) & (u1<=(W1_x(i,1)+W2_x(i,1)+W3_x(i,1)+W4_x(i,1)));
    n1_x = sum(z1_x);
    n2_x = sum(z2_x);
    n3_x = sum(z3_x);
    n4_x = sum(z4_x);
    obs1_x = gamrnd(V1_x(i,1),M1_x(i,1)/V1_x(i,1),n1_x,1);
    obs2_x = gamrnd(V2_x(i,1),M2_x(i,1)/V2_x(i,1),n2_x,1);
    obs3_x = gamrnd(V3_x(i,1),M3_x(i,1)/V3_x(i,1),n3_x,1);
    obs4_x = gamrnd(V4_x(i,1),M4_x(i,1)/V4_x(i,1),n4_x,1);
    obs_x = [obs1_x;obs2_x;obs3_x;obs4_x];
    obs_x = sort(obs_x,1);
    cdf_x = W1_x(i,1)*gamcdf(obs_x,V1_x(i,1),M1_x(i,1)/V1_x(i,1))+...
    W2_x(i,1)*gamcdf(obs_x,V2_x(i,1),M2_x(i,1)/V2_x(i,1))+...
    W3_x(i,1)*gamcdf(obs_x,V3_x(i,1),M3_x(i,1)/V3_x(i,1))+...
    W4_x(i,1)*gamcdf(obs_x,V4_x(i,1),M4_x(i,1)/V4_x(i,1));
    cdf_x(200000,1) = 1;
    table_x = [obs_x cdf_x];
    
    %get the observation for y
    z1_y = (u1<=W1_y(i,1));
    z2_y = (u1>W1_y(i,1)) & (u1<=(W1_y(i,1)+W2_y(i,1)));
    z3_y = (u1>(W1_y(i,1)+W2_y(i,1)) & u1<=(W1_y(i,1)+W2_y(i,1)+W3_y(i,1)));
    z4_y = (u1>(W1_y(i,1)+W2_y(i,1)+W3_y(i,1))) & (u1<=(W1_y(i,1)+W2_y(i,1)+W3_y(i,1)+W4_y(i,1)));
    n1_y = sum(z1_y);
    n2_y = sum(z2_y);
    n3_y = sum(z3_y);
    n4_y = sum(z4_y);
    obs1_y = gamrnd(V1_y(i,1),M1_y(i,1)/V1_y(i,1),n1_y,1);
    obs2_y = gamrnd(V2_y(i,1),M2_y(i,1)/V2_y(i,1),n2_y,1);
    obs3_y = gamrnd(V3_y(i,1),M3_y(i,1)/V3_y(i,1),n3_y,1);
    obs4_y = gamrnd(V4_y(i,1),M4_y(i,1)/V4_y(i,1),n4_y,1);
    obs_y = [obs1_y;obs2_y;obs3_y;obs4_y];
    obs_y = sort(obs_y,1);
    cdf_y = W1_y(i,1)*gamcdf(obs_y,V1_y(i,1),M1_y(i,1)/V1_y(i,1))+...
    W2_y(i,1)*gamcdf(obs_y,V2_y(i,1),M2_y(i,1)/V2_y(i,1))+...
    W3_y(i,1)*gamcdf(obs_y,V3_y(i,1),M3_y(i,1)/V3_y(i,1))+...
    W4_y(i,1)*gamcdf(obs_y,V4_y(i,1),M4_y(i,1)/V4_y(i,1));
    cdf_y(200000,1) = 1;
    table_y = [obs_y cdf_y];
            
    for t = 1:n
       index_x = find(p(t,1)<=table_x(:,2),1,'first');
       if index_x == 1
          quan_x(t,1) = table_x(index_x(1,1),1);
       else
          quan_x(t,1) = unifrnd(table_x(index_x(1,1)-1,1),table_x(index_x(1,1),1));
       end    
         
       index_y = find(p(t,1)<=table_y(:,2),1,'first');
      if index_y == 1
          quan_y(t,1) = table_y(index_y(1,1),1);
      else
          quan_y(t,1) = unifrnd(table_y(index_y(1,1)-1,1),table_y(index_y(1,1),1));
      end   
    end
    quantile_matrix_x(i,:) = (quan_x)';
    quantile_matrix_y(i,:) = (quan_y)';

    GLD_x = (((W1_x(i,1)).*(M1_x(i,1)).*gamcdf(quan_x,((V1_x(i,1))+1),(M1_x(i,1))/(V1_x(i,1)))+...
             (W2_x(i,1)).*(M2_x(i,1)).*gamcdf(quan_x,((V2_x(i,1))+1),(M2_x(i,1))/(V2_x(i,1)))+...
             (W3_x(i,1)).*(M3_x(i,1)).*gamcdf(quan_x,((V3_x(i,1))+1),(M3_x(i,1))/(V3_x(i,1)))+...
             (W4_x(i,1)).*(M4_x(i,1)).*gamcdf(quan_x,((V4_x(i,1))+1),(M4_x(i,1))/(V4_x(i,1))))) - ... 
             (((W1_x(i,1)).*(M1_x(i,1)).*gamcdf(quan_x(1,1),((V1_x(i,1))+1),(M1_x(i,1))/(V1_x(i,1)))+...
             (W2_x(i,1)).*(M2_x(i,1)).*gamcdf(quan_x(1,1),((V2_x(i,1))+1),(M2_x(i,1))/(V2_x(i,1)))+...
             (W3_x(i,1)).*(M3_x(i,1)).*gamcdf(quan_x(1,1),((V3_x(i,1))+1),(M3_x(i,1))/(V3_x(i,1)))+ ...
             (W4_x(i,1)).*(M4_x(i,1)).*gamcdf(quan_x(1,1),((V4_x(i,1))+1),(M4_x(i,1))/(V4_x(i,1)))));
    mean_x(i,1) = (W1_x(i,1)*M1_x(i,1) + W2_x(i,1)*M2_x(i,1) + W3_x(i,1)*M3_x(i,1) + W4_x(i,1)*M4_x(i,1)) - ... 
        (((W1_x(i,1)).*(M1_x(i,1)).*gamcdf(quan_x(1,1),((V1_x(i,1))+1),(M1_x(i,1))/(V1_x(i,1)))+...
             (W2_x(i,1)).*(M2_x(i,1)).*gamcdf(quan_x(1,1),((V2_x(i,1))+1),(M2_x(i,1))/(V2_x(i,1)))+...
             (W3_x(i,1)).*(M3_x(i,1)).*gamcdf(quan_x(1,1),((V3_x(i,1))+1),(M3_x(i,1))/(V3_x(i,1)))+ ...
             (W4_x(i,1)).*(M4_x(i,1)).*gamcdf(quan_x(1,1),((V4_x(i,1))+1),(M4_x(i,1))/(V4_x(i,1)))));
    LD_x = GLD_x./mean_x(i,1);
    GLD_matrix_x(i,:) = (GLD_x)';
    LD_matrix_x(i,:) = (LD_x)';
    
    GLD_y = (((W1_y(i,1)).*(M1_y(i,1)).*gamcdf(quan_y,((V1_y(i,1))+1),(M1_y(i,1))/(V1_y(i,1)))+...
             (W2_y(i,1)).*(M2_y(i,1)).*gamcdf(quan_y,((V2_y(i,1))+1),(M2_y(i,1))/(V2_y(i,1)))+...
             (W3_y(i,1)).*(M3_y(i,1)).*gamcdf(quan_y,((V3_y(i,1))+1),(M3_y(i,1))/(V3_y(i,1)))+ ...
             (W4_y(i,1)).*(M4_y(i,1)).*gamcdf(quan_y,((V4_y(i,1))+1),(M4_y(i,1))/(V4_y(i,1))))) - ...
             (((W1_y(i,1)).*(M1_y(i,1)).*gamcdf(quan_y(1,1),((V1_y(i,1))+1),(M1_y(i,1))/(V1_y(i,1)))+...
             (W2_y(i,1)).*(M2_y(i,1)).*gamcdf(quan_y(1,1),((V2_y(i,1))+1),(M2_y(i,1))/(V2_y(i,1)))+...
             (W3_y(i,1)).*(M3_y(i,1)).*gamcdf(quan_y(1,1),((V3_y(i,1))+1),(M3_y(i,1))/(V3_y(i,1)))+...
         (W4_y(i,1)).*(M4_y(i,1)).*gamcdf(quan_y(1,1),((V4_y(i,1))+1),(M4_y(i,1))/(V4_y(i,1)))));
    
    mean_y(i,1) = (W1_y(i,1)*M1_y(i,1) + W2_y(i,1)*M2_y(i,1) + W3_y(i,1)*M3_y(i,1) + W4_y(i,1)*M4_y(i,1)) - ...
                (((W1_y(i,1)).*(M1_y(i,1)).*gamcdf(quan_y(1,1),((V1_y(i,1))+1),(M1_y(i,1))/(V1_y(i,1)))+...
             (W2_y(i,1)).*(M2_y(i,1)).*gamcdf(quan_y(1,1),((V2_y(i,1))+1),(M2_y(i,1))/(V2_y(i,1)))+...
             (W3_y(i,1)).*(M3_y(i,1)).*gamcdf(quan_y(1,1),((V3_y(i,1))+1),(M3_y(i,1))/(V3_y(i,1)))+...
         (W4_y(i,1)).*(M4_y(i,1)).*gamcdf(quan_y(1,1),((V4_y(i,1))+1),(M4_y(i,1))/(V4_y(i,1)))));
         
    LD_y = GLD_y./mean_y(i,1);
    GLD_matrix_y(i,:) = (GLD_y)';
    LD_matrix_y(i,:) = (LD_y)';    
end
% obtaining the posterior probabilities of dominance, see the paper Lander
% Et al 2019.
for j = 1:1000
    j
    R = (randperm(m))';
    for i = 1:m
         GLD_matrix_y_perm(i,:) = GLD_matrix_y(R(i,1),:);
         LD_matrix_y_perm(i,:) = LD_matrix_y(R(i,1),:);
         quantile_matrix_y_perm(i,:) = quantile_matrix_y(R(i,1),:);
    end
    
    xGLDy = GLD_matrix_x>=GLD_matrix_y_perm;
    yGLDx = GLD_matrix_y_perm>=GLD_matrix_x;
    xGLDy = double(xGLDy);
    yGLDx = double(yGLDx);
    prop_xGLDy(j,:) = mean(xGLDy);
    prop_yGLDx(j,:) = mean(yGLDx);
    overall_xGLDy(j,1) = mean(prod(xGLDy,2));
    overall_yGLDx(j,1) = mean(prod(yGLDx,2));
    overall_xGLDy_20lowest(j,1) = mean(prod(xGLDy(:,1:200),2));
    overall_yGLDx_20lowest(j,1) = mean(prod(yGLDx(:,1:200),2));
    overall_xGLDy_10lowest(j,1) = mean(prod(xGLDy(:,1:100),2));
    overall_yGLDx_10lowest(j,1) = mean(prod(yGLDx(:,1:100),2));
    
    xLDy = LD_matrix_x>=LD_matrix_y_perm;
    yLDx = LD_matrix_y_perm>=LD_matrix_x;
    xLDy = double(xLDy);
    yLDx = double(yLDx);
    prop_xLDy(j,:) = mean(xLDy);
    prop_yLDx(j,:) = mean(yLDx);
    overall_xLDy(j,1) = mean(prod(xLDy,2));
    overall_yLDx(j,1) = mean(prod(yLDx,2));
    overall_xLDy_20lowest(j,1) = mean(prod(xLDy(:,1:200),2));
    overall_yLDx_20lowest(j,1) = mean(prod(yLDx(:,1:200),2));
    overall_xLDy_10lowest(j,1) = mean(prod(xLDy(:,1:100),2));
    overall_yLDx_10lowest(j,1) = mean(prod(yLDx(:,1:100),2));
    
    xFSDy = quantile_matrix_x>=quantile_matrix_y_perm;
    yFSDx = quantile_matrix_y_perm>=quantile_matrix_x;
    xFSDy = double(xFSDy);
    yFSDx = double(yFSDx);
    prop_xFSDy(j,:) = mean(xFSDy);
    prop_yFSDx(j,:) = mean(yFSDx);
    overall_xFSDy(j,1) = mean(prod(xFSDy,2));
    overall_yFSDx(j,1) = mean(prod(yFSDx,2));
    overall_xFSDy_20lowest(j,1) = mean(prod(xFSDy(:,1:200),2));
    overall_yFSDx_20lowest(j,1) = mean(prod(yFSDx(:,1:200),2));
    overall_xFSDy_10lowest(j,1) = mean(prod(xFSDy(:,1:100),2));
    overall_yFSDx_10lowest(j,1) = mean(prod(yFSDx(:,1:100),2));
    
end

mean_prop_xGLDy = mean(prop_xGLDy);
mean_prop_yGLDx = mean(prop_yGLDx);
mean_overall_xGLDy = mean(overall_xGLDy);
max_overall_xGLDy = max(overall_xGLDy);
min_overall_xGLDy = min(overall_xGLDy);
mean_overall_yGLDx = mean(overall_yGLDx);
max_overall_yGLDx = max(overall_yGLDx);
min_overall_yGLDx = min(overall_yGLDx);
mean_overall_xGLDy_20lowest = mean(overall_xGLDy_20lowest);
max_overall_xGLDy_20lowest = max(overall_xGLDy_20lowest);
min_overall_xGLDy_20lowest = min(overall_xGLDy_20lowest);
mean_overall_yGLDx_20lowest = mean(overall_yGLDx_20lowest);
max_overall_yGLDx_20lowest = max(overall_yGLDx_20lowest);
min_overall_yGLDx_20lowest = min(overall_yGLDx_20lowest);
mean_overall_xGLDy_10lowest = mean(overall_xGLDy_10lowest);
max_overall_xGLDy_10lowest = max(overall_xGLDy_10lowest);
min_overall_xGLDy_10lowest = min(overall_xGLDy_10lowest);
mean_overall_yGLDx_10lowest = mean(overall_yGLDx_10lowest);
max_overall_yGLDx_10lowest = max(overall_yGLDx_10lowest);
min_overall_yGLDx_10lowest = min(overall_yGLDx_10lowest);



mean_prop_xLDy = mean(prop_xLDy);
mean_prop_yLDx = mean(prop_yLDx);
mean_overall_xLDy = mean(overall_xLDy);
max_overall_xLDy = max(overall_xLDy);
min_overall_xLDy = min(overall_xLDy);
mean_overall_yLDx = mean(overall_yLDx);
max_overall_yLDx = max(overall_yLDx);
min_overall_yLDx = min(overall_yLDx);
mean_overall_xLDy_20lowest = mean(overall_xLDy_20lowest);
max_overall_xLDy_20lowest = max(overall_xLDy_20lowest);
min_overall_xLDy_20lowest = min(overall_xLDy_20lowest);
mean_overall_yLDx_20lowest = mean(overall_yLDx_20lowest);
max_overall_yLDx_20lowest = max(overall_yLDx_20lowest);
min_overall_yLDx_20lowest = min(overall_yLDx_20lowest);
mean_overall_xLDy_10lowest = mean(overall_xLDy_10lowest);
max_overall_xLDy_10lowest = max(overall_xLDy_10lowest);
min_overall_xLDy_10lowest = min(overall_xLDy_10lowest);
mean_overall_yLDx_10lowest = mean(overall_yLDx_10lowest);
max_overall_yLDx_10lowest = max(overall_yLDx_10lowest);
min_overall_yLDx_10lowest = min(overall_yLDx_10lowest);

mean_prop_xFSDy = mean(prop_xFSDy);
mean_prop_yFSDx = mean(prop_yFSDx);
mean_overall_xFSDy = mean(overall_xFSDy);
max_overall_xFSDy = max(overall_xFSDy);
min_overall_xFSDy = min(overall_xFSDy);
mean_overall_yFSDx = mean(overall_yFSDx);
max_overall_yFSDx = max(overall_yFSDx);
min_overall_yFSDx = min(overall_yFSDx);
mean_overall_xFSDy_20lowest = mean(overall_xFSDy_20lowest);
max_overall_xFSDy_20lowest = max(overall_xFSDy_20lowest);
min_overall_xFSDy_20lowest = min(overall_xFSDy_20lowest);
mean_overall_yFSDx_20lowest = mean(overall_yFSDx_20lowest);
max_overall_yFSDx_20lowest = max(overall_yFSDx_20lowest);
min_overall_yFSDx_20lowest = min(overall_yFSDx_20lowest);
mean_overall_xFSDy_10lowest = mean(overall_xFSDy_10lowest);
max_overall_xFSDy_10lowest = max(overall_xFSDy_10lowest);
min_overall_xFSDy_10lowest = min(overall_xFSDy_10lowest);
mean_overall_yFSDx_10lowest = mean(overall_yFSDx_10lowest);
max_overall_yFSDx_10lowest = max(overall_yFSDx_10lowest);
min_overall_yFSDx_10lowest = min(overall_yFSDx_10lowest);    

save('prop.mat','prop_xGLDy','prop_yGLDx','prop_xLDy','prop_yLDx','prop_xFSDy','prop_yFSDx');
save('overall_mean.mat','mean_overall_xGLDy','mean_overall_yGLDx','mean_overall_xGLDy_20lowest','mean_overall_yGLDx_20lowest',...
    'mean_overall_xGLDy_10lowest','mean_overall_yGLDx_10lowest','mean_overall_xLDy','mean_overall_yLDx','mean_overall_xLDy_20lowest',...
    'mean_overall_yLDx_20lowest','mean_overall_xLDy_10lowest','mean_overall_yLDx_10lowest','mean_overall_xFSDy','mean_overall_yFSDx',...
    'mean_overall_xFSDy_20lowest','mean_overall_yFSDx_20lowest','mean_overall_xFSDy_10lowest','mean_overall_yFSDx_10lowest');
    
save('overall_max.mat','max_overall_xGLDy','max_overall_yGLDx','max_overall_xGLDy_20lowest','max_overall_yGLDx_20lowest',...
    'max_overall_xGLDy_10lowest','max_overall_yGLDx_10lowest','max_overall_xLDy','max_overall_yLDx','max_overall_xLDy_20lowest',...
    'max_overall_yLDx_20lowest','max_overall_xLDy_10lowest','max_overall_yLDx_10lowest','max_overall_xFSDy','max_overall_yFSDx',...
    'max_overall_xFSDy_20lowest','max_overall_yFSDx_20lowest','max_overall_xFSDy_10lowest','max_overall_yFSDx_10lowest');
    
save('overall_min.mat','min_overall_xGLDy','min_overall_yGLDx','min_overall_xGLDy_20lowest','min_overall_yGLDx_20lowest',...
    'min_overall_xGLDy_10lowest','min_overall_yGLDx_10lowest','min_overall_xLDy','min_overall_yLDx','min_overall_xLDy_20lowest',...
    'min_overall_yLDx_20lowest','min_overall_xLDy_10lowest','min_overall_yLDx_10lowest','min_overall_xFSDy','min_overall_yFSDx',...
    'min_overall_xFSDy_20lowest','min_overall_yFSDx_20lowest','min_overall_xFSDy_10lowest','min_overall_yFSDx_10lowest');
save('overall.mat','overall_xGLDy','overall_yGLDx','overall_xGLDy_20lowest','overall_yGLDx_20lowest',...
     'overall_xGLDy_10lowest','overall_yGLDx_10lowest','overall_xLDy','overall_yLDx','overall_xLDy_20lowest','overall_yLDx_20lowest',...
     'overall_xLDy_10lowest','overall_yLDx_10lowest','overall_xFSDy','overall_yFSDx','overall_xFSDy_20lowest','overall_yFSDx_20lowest',...
     'overall_xFSDy_10lowest','overall_yFSDx_10lowest')
 
 


