%%%
%Simple Univariate model using Normal-Extended gamma variables
%%%

%%%
%Simple Univariate model using Normal-Inverse gamma variables
%%%

clear all
%%
%Simulate the data

rng(123)

%Number of groups/sources
J = 2;

%Two components, so we only need one probability of success
p1_simul = [.2 .1];

%Initialise mu's and sigma^2's
sigma2_simul = [.1 .5 .5 1.5];
mu_simul = [-3 0 0 1];

%Simulate data conditionally on mu,sigma2
y_simul = cell(1,J);
%Number of observations in each group
n = [100 100];
n1 = zeros(1,J); %number of elements in the first group

for j = 1:J
    n1(j) = binornd(n(j),p1_simul(j));
    y_simul{j} = [normrnd(mu_simul((j-1)*2 + 1), sqrt(sigma2_simul((j-1)*2 + 1)), 1, n1(j)) normrnd(mu_simul((j-1)*2 + 2), sqrt(sigma2_simul((j-1)*2 + 2)), 1, n(j) - n1(j))];
end

%Plot of simulating distirbutions and data histograms

%Colors for histograms
Color_Hist = zeros(J*2,3);
Color_Hist(1,:) = [106 168 239]/255; %Component 1 is blue
Color_Hist(2,:) = [183 121 245]/255; %Component 12 is purple-ish
Color_Hist(3,:) = Color_Hist(2,:); %Component 21 is shared
Color_Hist(4,:) = [245 60 60]/255; %Component 3 is red

Color_Lines = zeros(J*2,3);
Color_Lines(1,:) = [0 0 1]; %Component 1 is blue
Color_Lines(2,:) = [120 49 168]/255; %Component 12 is purple-ish
Color_Lines(3,:) = Color_Lines(2,:); %Component 21 is shared
Color_Lines(4,:) = [1 0 0]; %Component 3 is red

figure(1)
dx = (max([y_simul{:}])-min([y_simul{:}]))/50;
xx = min([y_simul{:}])-.5:dx:max([y_simul{:}])+.5;
hold on
for j = 1:J
    [f_k] = hist(y_simul{j}(1:n1(j)),xx);
    bar(xx,p1_simul(j)*f_k/sum(f_k)/dx,'FaceColor',Color_Hist((j-1)*2 + 1,:),'facealpha',.5,'edgecolor','none');
    [f_k] = hist(y_simul{j}(n1(j)+1:n(j)),xx);
    bar(xx,(1-p1_simul(j))*f_k/sum(f_k)/dx,'FaceColor',Color_Hist((j-1)*2 + 2,:),'facealpha',.5,'edgecolor','none');
end
xx = min([y_simul{:}])-.5:0.01:max([y_simul{:}])+.5;
for j = 1:J
    yy = p1_simul(j) * normpdf(xx,mu_simul((j-1)*2 + 1),sqrt(sigma2_simul((j-1)*2 + 1))) + (1 - p1_simul(j)) * normpdf(xx,mu_simul((j-1)*2 + 2),sqrt(sigma2_simul((j-1)*2 + 2)));
    n1_xx = sum(xx<=mu_simul((j-1)*2 + 1)+sqrt(sigma2_simul((j-1)*2 + 1)));
    plot(xx(1:n1_xx),yy(1:n1_xx),'Color',Color_Lines((j-1)*2 + 1,:),'LineWidth',3)
    plot(xx(n1_xx:end),yy(n1_xx:end),'Color',Color_Lines((j-1)*2 + 2,:),'LineWidth',3)
end
axis tight
legend('Comp 1','Comp 12','Comp 21','Comp 3')
xlabel('Y')
%title('Simulating densities and simulated data','FontSize',15)
print(1,'-djpeg','data.jpeg')
close(1)

%%

%Initialise configuration and labels etc...

%Take the simulated data
y = y_simul;

%Gibbs
n_burn1 = 1000;
n_burn2 = 24000;
thin = 1;
n_save = 5000;
n_iter = n_burn1 + n_burn2 + thin*n_save;


%Allocation variables, number of clusters and unique values
Z = cell(1,J);
K = ones(1,J);
N = cell(1,J);

%New latent variables u
u = ones(1,J);
%Non-conjugate
s_u = 1*ones(1,J);
u_accept = zeros(1,J);
u_count = zeros(1,J);



%Parameters of NGG processes
kappa_0 = 0.01;

sigma_0 = .01;

%Hyperparameters of P0
%m0
s1 = 1000;
m0 = 0.25;

%k0
k0 = 0.62;

%Hyperparameters of sigma2
a_sigma2 = 2.07;
b_sigma2 = 0.66;

%For simulation of alpha with groups of different sizes
j1 = [0 1]; j2 = [0 1];
index_j = zeros(1,length(j1)*length(j2));
n_aux = prod(n) * ones(1,length(j1)*length(j2));
%This calculation is not very general!
for i1 = 1:length(j1)
    for i2 = 1:length(j2)
        index_j((i1-1)*length(j2) + i2) = j1(i1) + j2(i2);
        if i1 == 1
            n_aux((i1-1)*length(j2) + i2) = n_aux((i1-1)*length(j2) + i2) / n(1);
        end
        if i2 == 1
            n_aux((i1-1)*length(j2) + i2) = n_aux((i1-1)*length(j2) + i2) / n(2);
        end
    end
end

%Latent varuables in restaurants
% Sigma2 = cell(1,J);
% Mu = cell(1,J);

%All in one cluster
for j = 1:J
    %     Sigma2{j} = 1/gamrnd(a_sigma2,1/b_sigma2);
    %     Mu{j} = sqrt(Sigma2{j})*randn;
    N{j} = n(j);
    Z{j} = ones(1,n(j));
end

%Adaptive variance
tau = .234;
g_w = 0.7;

%Outputs
Z_out = cell(n_save,J);
N_out = cell(n_save,J);
K_out = zeros(n_save,J);
U_out = zeros(n_save,J);
Mu_out = cell(n_save,J);
Sigma2_out = cell(n_save,J);

%%

tic
for iter = 1:n_iter
    
    %For element ji, sample the table value, after removing the element from
    %the list
    
    for j = 1:J
        for i = 1:n(j)
            
            N{j}(Z{j}(i)) = N{j}(Z{j}(i)) - 1;
            %If the cluster is left empty the number of clusters is
            %affected, and I need to remove the unique value
            aloneN = 0;
            if N{j}(Z{j}(i)) == 0
                aloneN = 1;
            end
            
            %Weights for each dish
            f_k = -inf*ones(1,K(j)+1);
            %In the first position is the probability of a new table
            scale = b_sigma2/a_sigma2*(1+k0)/k0;
            f_k(1) = log( tpdf((y{j}(i)-m0)/sqrt(scale),2*a_sigma2)/sqrt(scale) );
            
            for k = 1:K(j)
                
                %People eating dish k excluding data point i
                index = find(Z{j}==k);
                index = setdiff(index, i);
                den = y{j}(index);
                
                if ~isempty(den)
                    sizeden = length(den);
                    sumden = sum(den);
                    meanden = sumden/sizeden;
                                        
                    k0_den = k0 + sizeden;
                    m0_den = (k0*m0 + sumden)/k0_den;
                    a_den = a_sigma2 + sizeden/2;
                    b_den = b_sigma2 + .5*(sum((den - meanden).^2) + k0*sizeden*(meanden - m0)^2/k0_den);
                    scale_aux = b_den/a_den*(1+k0_den)/k0_den;
                    
                    f_k(k + 1) = log( tpdf((y{j}(i)-m0_den)/sqrt(scale_aux),2*a_den)/sqrt(scale_aux) );
                end
            end
            
            w = [kappa_0 * (u(j)+1)^sigma_0 N{j} - sigma_0];
            
            %The (N{j} - sigma_0) term could be < 0!
            if aloneN
                w(1 + Z{j}(i)) = 0;
            end
            
            f_k = exp(f_k-max(f_k)) .* w;
            f_k = f_k/sum(f_k);
            
            cumsum_f_k = cumsum(f_k);
            hh = find(rand <= cumsum_f_k,1);
            
            if hh == 1
                
                %                 %Sample from conditional distirbution of latent
                %                 %variables
                %                 new_Sigma2 = 1/gamrnd(a_sigma2,1/b_sigma2);
                %                 new_Mu = sqrt(new_Sigma2)*randn;
                
                if aloneN%if aloneN==1 I only change the value of Eta{i}(Z(i,j))=Tau_new;
                    %Change the number of clusters that now are one less
                    %(The element was alone!)
                    %Change the unique value (update tau_new)
                    N{j}(Z{j}(i)) = 1;
                    
                    %                     Sigma2{j}(Z{j}(i)) = new_Sigma2;
                    %                     Mu{j}(Z{j}(i)) = new_Mu;
                else
                    %Add a new cluster of dimension one, with a new unique
                    %value
                    K(j) = K(j) + 1;
                    Z{j}(i) = K(j);
                    N{j} = [N{j} 1];
                    %                     Sigma2{j} = [Sigma2{j} new_Sigma2];
                    %                     Mu{j} = [Mu{j} new_Mu];
                end
                
            else %I pick one of the other clusters
                h = hh - 1;
                
                if aloneN
                    %I have one cluster less
                    K(j) = K(j) - 1;
                    h_aux = Z{j}(i);
                    Z{j}(i) = h;
                    N{j}(h) = N{j}(h) + 1;
                    %"Eliminate" the empty cluster
                    %                     Sigma2{j}(h_aux) = [];
                    %                     Mu{j}(h_aux) = [];
                    N{j}(h_aux) = [];
                    Z{j}(Z{j} > h_aux) = Z{j}(Z{j} > h_aux) - 1;
                else
                    %Same number of clusters, I am just moving the element Z(i,j)
                    Z{j}(i) = h;
                    N{j}(h)= N{j}(h) + 1;
                end
            end
        end
    end
    
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Acceleration step
    Sigma2 = cell(1,J);
    Mu = cell(1,J);
    for j = 1:J
        for k = 1:K(j)
            %All the people eating dish k in restaurant j
            dish_k = y{j}(Z{j}==k);
            size_k = length(dish_k);
            sum_k = sum(dish_k);
            mean_k = sum_k/size_k;
            k0_aux = k0 + size_k;
            
            Sigma2{j}(k) = 1/gamrnd(a_sigma2 + size_k/2,1/( b_sigma2 + .5*(sum((dish_k - mean_k).^2) + k0*size_k*(mean_k - m0)^2/k0_aux) ));
            Mu{j}(k) = normrnd((k0*m0+sum_k)/k0_aux,sqrt(Sigma2{j}(k)/k0_aux));
        end
    end
    

    
    %Update u_j
    for j = 1:J
        u_new = exp(log(u(j)) + sqrt(s_u(j))*randn);
        
        log_accept = n(j)*(log(u_new) - log(u(j)));
        log_accept = log_accept - kappa_0/sigma_0*((u_new+1)^sigma_0 - (u(j)+1)^sigma_0);
        log_accept = log_accept - (n(j) - K(j)*sigma_0) * (log(u_new+1) - log(u(j)+1));
        
        if (isreal(log_accept) == 0 )
            stop(1);
        end
        
        accept = 1;
        if (isnan(log_accept) == 1)
            accept = 0;
        elseif ( log_accept < 0 )
            accept = exp(log_accept);
        end
        
        u_accept(j) = u_accept(j) + accept;
        u_count(j) = u_count(j) + 1;
        
        if ( rand < accept )
            u(j) = u_new;
        end
        
        if iter > n_burn1
            s_u(j) = exp(log(s_u(j))+iter^(-g_w)*(accept-tau));
        end
    end
    

    
    %Saving the output
    if iter > (n_burn1 + n_burn2) && mod((iter - n_burn1 - n_burn2),thin) == 0
        
        iter_aux = (iter - n_burn1 - n_burn2)/thin;
        
        for j = 1:J
            Z_out{iter_aux,j} = Z{j};
            N_out{iter_aux,j} = N{j};
            Sigma2_out{iter_aux,j} = Sigma2{j};
            Mu_out{iter_aux,j} = Mu{j};
        end
        K_out(iter_aux,:) = K;
        U_out(iter_aux,:) = u;
    end
end
toc
%profile viewer

%%
save('N_Ind_NGG')


%% PLOTS
load('N_Ind_NGG')

%%
%Predictive
dg = (max([y_simul{:}])-min([y_simul{:}]))/100;
grid_pred = (min([y_simul{:}])-3):dg:(max([y_simul{:}])+3);
l_grid = length(grid_pred);
pred = zeros(J+1,l_grid);

for it = 1:n_save
    scale = b_sigma2/a_sigma2*(k0 + 1)/k0;
    %Predictive in source group
    for j = 1:J
        pred(j,:) = pred(j,:) + kappa_0 * (U_out(it,j)+1)^sigma_0/(kappa_0 * (U_out(it,j)+1)^sigma_0 + n(j) - K_out(it,j)*sigma_0) * tpdf((grid_pred-m0)/sqrt(scale),2*a_sigma2)/sqrt(scale);
        for k = 1:K_out(it,j)
            dish_k = y{j}(Z_out{it,j}==k);
            size_k = length(dish_k);
            sum_k = sum(dish_k);
            mean_k = sum_k/size_k;
            k0_k = k0 + size_k;
            m0_k = (k0*m0 + sum_k)/k0_k;
            a_k = a_sigma2 + size_k/2;
            b_k = b_sigma2 + .5*(sum((dish_k - mean_k).^2) + k0*size_k*(mean_k - m0)^2/k0_k);
            scale_k = b_k/a_k*(1+k0_k)/k0_k;
            
            pred(j,:) = pred(j,:) + (N_out{it,j}(k) - sigma_0)/(kappa_0 * (U_out(it,j)+1)^sigma_0 + n(j) - K_out(it,j)*sigma_0) * tpdf((grid_pred-m0_k)/sqrt(scale_k),2*a_k)/sqrt(scale_k);
        end
    end
    %Predictive of new restaurant
    %Only new cluster as DP's are iid here
    pred(J+1,:) = pred(J+1,:) + tpdf((grid_pred-m0)/sqrt(scale),2*a_sigma2)/sqrt(scale);
end

pred = pred/n_save;
sum(pred'*dg)

%Plot of predictive distirbutions and data histograms

%Predictive in groups
figure(1)
dx = (max([y_simul{:}])-min([y_simul{:}]))/25;
xx = min([y_simul{:}])-.5:dx:max([y_simul{:}])+.5;
hold on
for j = 1:J
    [f_k] = hist([y_simul{j}],xx);
    bar(xx,p1_simul(j)*f_k/sum(f_k)/dx,'FaceColor',Color_Hist((j-1)*2 + 1,:),'facealpha',.5,'edgecolor','none');
    [f_k] = hist(y_simul{j}(n1(j)+1:n(j)),xx);
    bar(xx,(1-p1_simul(j))*f_k/sum(f_k)/dx,'FaceColor',Color_Hist((j-1)*2 + 2,:),'facealpha',.5,'edgecolor','none');
end
xx = grid_pred;
plot(xx,pred(1,:),'-','LineWidth',3)
plot(xx,pred(2,:),'-.','LineWidth',3)
legend({'Comp 1','Comp 2 (Group 1)','Comp 2 (Group 2)','Comp 3','Predictive (Group 1)','Predictive (Group 2)'},'FontSize',15)
xlabel('Y', 'Fontsize', 15)
axis tight
title('')
print(1,'-djpeg','N_Ind_NGG_Sigma_Fixed_Pred_groups.jpeg')
print(1,'-depsc','N_Ind_NGG_Sigma_Fixed_Pred_groups.eps')
print(1,'-deps','N_Ind_NGG_Sigma_Fixed_Pred_groups_bw.eps')
close(1)

%Predictive new group
figure(2)
dx = (max([y_simul{:}])-min([y_simul{:}]))/25;
xx = min([y_simul{:}])-.5:dx:max([y_simul{:}])+.5;
hold on
[f_k] = hist([y_simul{:}],xx);
bar(xx,f_k/sum(f_k)/dx,'FaceColor',Color_Hist((j-1)*2 + 1,:),'facealpha',.5,'edgecolor','none');
xx = grid_pred;
plot(xx,pred(J+1,:),'LineWidth',3)
xlabel('Y', 'Fontsize', 15)
legend({'Data','Predictive new Group'}, 'Fontsize', 15)
axis tight
title('')
print(2,'-djpeg','N_Ind_NGG_Sigma_Fixed_Pred_Newgroup.jpeg')
print(2,'-depsc','N_Ind_NGG_Sigma_Fixed_Pred_Newgroup.eps')
print(2,'-deps','N_Ind_NGG_Sigma_Fixed_Pred_Newgroup_bw.eps')
close(2)

%% 
%Posterior total number of dishes
%AND
%Posterior number of clusters in restaurants

figure(3)
%Sum number of dishes in each restaurant: by definition they will be
%different (absolutely continuous P0!)
Kn_out = sum(K_out,2);

%Total number of dishes - Histogram
xx = unique(Kn_out);
ff = hist(Kn_out,xx);
bar(xx,ff/sum(ff));
axis tight
legend({'Number of Clusters (both Groups)'}, 'Fontsize', 15)
title('')
print(3,'-djpeg','N_Ind_NGG_Sigma_Fixed_Kn.jpeg')
print(3,'-depsc','N_Ind_NGG_Sigma_Fixed_Kn.eps')
print(3,'-deps','N_Ind_NGG_Sigma_Fixed_Kn_bw.eps')
close(3)


figure(4)
NumDish = K_out;

eps = [-0.13 0.13];

hold on
for j = 1:J
    xx = unique(NumDish(:,j));
    ff = hist(NumDish(:,j),xx);
    if j == 1
        bar(xx + eps(j),ff/sum(ff), 0.25,'FaceColor',Color_Hist(j,:),'facealpha',.55,'edgecolor','none');
        xx = unique(NumDish);
        XTick = num2str(unique(NumDish));
        set(gca,'xtick', xx)
        set(gca,'xticklabel', XTick)
    else
        bar(xx + eps(j),ff/sum(ff), 0.25,'FaceColor',Color_Hist(4,:),'facealpha',.55,'edgecolor','none');
    end
end
axis tight
title('')
legend({'Num Dishes Group 1', 'Num Dishes Group 2'}, 'Fontsize', 15)
print(4,'-djpeg','N_Ind_NGG_Sigma_Fixed_NumDishes.jpeg')
print(4,'-depsc','N_Ind_NGG_Sigma_Fixed_NumDishes.eps')
print(4,'-deps','N_Ind_NGG_Sigma_Fixed_NumDishes_bw.eps')
close(4)

%%
%LPML

CPO_j = cell(1,J);
CPO = zeros(1, sum(n));
for it = 1:n_save
    for j = 1:J
        aux = normpdf(y{j}, Mu_out{it,j}(Z_out{it,j}), sqrt(Sigma2_out{it,j}(Z_out{it,j})));
        CPO_j{j} = 1 ./ aux;
    end
    CPO = CPO + [CPO_j{:}];
end
LPML = sum(log( n_save ./ CPO ));
