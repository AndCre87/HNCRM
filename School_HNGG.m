%%%
% Univariate model using Normal-Inverse gamma variables (non-conjugate) and HNGG for School data
%%%

clear all
rng(123)

%% School data

Data = importdata('school_data.txt');
school = Data.data;
school_id = school(:,1);
school = school(:,2);

%Number of groups/sources
J = length(unique(school_id));
y = cell(1,J);
n = zeros(1,J);
for j = 1:J
    y{j} = school(school_id == j);
    n(j) = sum(school_id == j);
end


%Gibbs
n_burn1 = 1000;
n_burn2 = 24000;
thin = 1;
n_save = 5000;
n_iter = n_burn1 + n_burn2 + thin*n_save;

%Table-specific labels
T = cell(1,J);
%Dish-specific labels
K = cell(1,J);

%Numerosity in the clusters
%Number of elements in each group
N = cell(1,J);
%Number of tables eating dish k (sum over groups)
M_k = []; %Overall

%Labels
%Factors (one for each data in each group)

%Table-specific unique values
%One for each table in each group (but the same)
%NOT NECESSARY?

%New latent variables u
u = ones(1,J+1);
%Non-conjugate
s_u = 1*ones(1,J+1);
u_accept = zeros(1,J+1);
u_count = zeros(1,J+1);

%Parameters of NGG processes
mm = 1;
ss = sqrt(1);
a_kappa_0 = (mm/ss)^2;
b_kappa_0 = mm/ss^2;
kappa_0 = 1;

% mm = .5;
% ss = sqrt(1/20);
% a_sigma_0 = (1-mm)*mm^2/ss^2-mm;
% b_sigma_0 = a_sigma_0*(1-mm)/mm;
a_sigma_0 = 2;
b_sigma_0 = 18;
sigma_0 = .1;
%Non-conjugate
s_sigma_0 = 1;
sigma_0_accept = 0;
sigma_0_count = 0;

mm = 1;
ss = sqrt(1);
a_kappa = (mm/ss)^2;
b_kappa = mm/ss^2;
kappa = 1;

% mm = .5;
% ss = sqrt(1/20);
% a_sigma = (1-mm)*mm^2/ss^2-mm;
% b_sigma = a_sigma*(1-mm)/mm;
a_sigma = 2;
b_sigma = 18;
sigma = .1;
%Non-conjugate
s_sigma = 1;
sigma_accept = 0;
sigma_count = 0;

%Hyperparameters of P0
%m0
m0 = 50;
s20 = 25;

%Hyperparameters of sigma2
nu0 = 1;
sigma20 = 100;

%All together
for j = 1:J
    T{j} = ones(1,n(j));
    N{j} = n(j);
    K{j} = 1;
end
M_k = J;
Kn = 1;

%Number of tables in each restaurant
M_j = ones(1,J);

%Neal's 8 aglorithm variables
%Number of auxiliary variables
N_aux = 50;
%Initialise new set of auxiliary variables
Mu_aux = m0 + sqrt(s20)*randn(1, N_aux);
Sigma2_aux = 1./gamrnd(nu0/2, 1/(sigma20*nu0/2), 1, N_aux);

%Initialise
Sigma2 = ones(1,Kn);
Mu = zeros(1,Kn);

%Adaptive variance
tau = .234;
g_w = 0.7;

%Outputs
T_out = cell(n_save,J);
K_out = cell(n_save,J);
N_out = cell(n_save,J);
M_j_out = cell(1,n_save);
M_k_out = cell(1,n_save);
Tau_out = cell(n_save,J);
Sigma2_out = cell(1,n_save);
Mu_out = cell(1,n_save);
kappa_0_out = zeros(1,n_save);
sigma_0_out = zeros(1,n_save);
kappa_out = zeros(1,n_save);
sigma_out = zeros(1,n_save);
Kn_out = zeros(1,n_save);
U_out = zeros(n_save,J+1);

%%

tic
%profile on
for iter = 1:n_iter
    
    %For element ji, sample the table value, after removing the element from the list
    for j = 1:J
        for i = 1:n(j)
            
            aux_t = T{j}(i);
            aux_k = K{j}(aux_t);
            
            N{j}(aux_t) = N{j}(aux_t) - 1;
            
            %Empty clusters
            aloneN = ( N{j}(aux_t)==0 );
            if aloneN
                M_k(aux_k) = M_k(aux_k) - 1;
            end
            aloneMk = ( M_k(aux_k)==0 );
            if aloneMk %Part of the re-use algorithm
                aux_ind = randi(N_aux);
                Mu_aux(aux_ind) = Mu(aux_k);
                Sigma2_aux(aux_ind) = Sigma2(aux_k);
            end
            
            f_k = log( normpdf(y{j}(i), [Mu_aux Mu], sqrt([Sigma2_aux Sigma2])) );
            
            f_k = f_k([1:(N_aux+Kn) N_aux+K{j}]);
            w = [kappa_0 * (u(1)+1)^sigma_0 * ones(1,N_aux)/N_aux (M_k - sigma_0)];
            w = [kappa * (u(j+1)+1)^sigma * w/sum(w) (N{j} - sigma)];
            
            %The (N{j} - sigma) term could be < 0!
            if aloneN
                w(N_aux + Kn + aux_t) = 0;
                %Same thing for the (M_k - sigma_0) term!
                if aloneMk
                    w(N_aux + aux_k) = 0;
                end
            end
            
            f_k = exp(f_k-max(f_k)) .* w;
            f_k = f_k/sum(f_k);
            
            cumsum_f_k = cumsum(f_k);
            hh = find(rand <= cumsum_f_k,1);
            
            if hh <= N_aux
                if aloneN
                    %Select new dish for the new table but use same
                    %label (the cluster was empty)
                    N{j}(aux_t) = 1;
                    
                    if aloneMk
                        %Same number of dishes in the global menu
                        %Number of tables serving this dish is now one
                        M_k(aux_k) = 1;
                        Mu(aux_k) = Mu_aux(hh);
                        Sigma2(aux_k) = Sigma2_aux(hh);
                    else
                        Kn = Kn + 1;
                        M_k = [M_k 1];
                        K{j}(aux_t) = Kn;
                        Mu = [Mu Mu_aux(hh)];
                        Sigma2 = [Sigma2 Sigma2_aux(hh)];
                    end
                else
                    Kn = Kn + 1;
                    M_k = [M_k 1];
                    M_j(j) = M_j(j) + 1;
                    T{j}(i) = M_j(j);
                    K{j} = [K{j} Kn];
                    N{j} = [N{j} 1];
                    Mu = [Mu Mu_aux(hh)];
                    Sigma2 = [Sigma2 Sigma2_aux(hh)];
                end
                
                %Restore auxiliary variable
                Mu_aux(hh) = m0 + sqrt(s20)*randn;
                Sigma2_aux(hh) = 1./gamrnd(nu0/2, 1/(sigma20*nu0/2));
                
            elseif hh <= Kn + N_aux
                %New table, old dish
                h_k = hh - N_aux;
                
                M_k(h_k) = M_k(h_k) + 1;
                
                if aloneN
                    K{j}(aux_t) = h_k;
                    N{j}(aux_t) = 1;
                    if aloneMk
                        Kn = Kn - 1;
                        for i_j = 1:J
                            K{i_j}(K{i_j} > aux_k) = K{i_j}(K{i_j} > aux_k) - 1;
                        end
                        M_k(aux_k) = [];
                        Mu(aux_k) = [];
                        Sigma2(aux_k) = [];
                    end
                else
                    M_j(j) = M_j(j) + 1;
                    T{j}(i) = M_j(j);
                    K{j} = [K{j} h_k];
                    N{j} = [N{j} 1];
                end
            else
                %Old table, old dish
                h_t = hh - Kn - N_aux;
                
                T{j}(i) = h_t;
                N{j}(h_t) = N{j}(h_t) + 1;
                if aloneN
                    %One new person at table h_t
                    %No-one is sitting at table aux -> remove
                    N{j}(aux_t) = [];
                    T{j}(T{j} > aux_t) = T{j}(T{j} > aux_t) - 1;
                    K{j}(aux_t) = [];
                    M_j(j) = M_j(j) - 1;
                    if aloneMk
                        Kn = Kn - 1;
                        for i_j = 1:J
                            K{i_j}(K{i_j} > aux_k) = K{i_j}(K{i_j} > aux_k) - 1;
                        end
                        M_k(aux_k) = [];
                        Mu(aux_k) = [];
                        Sigma2(aux_k) = [];
                    end
                end
            end
        end
    end
    
    %Sample the table dish label, removing all the element at the table
    for j = 1:J
        num_tables = length(unique(T{j}));
        for t = 1:num_tables;
            
            %Dish served at this table
            aux_k = K{j}(t);
            
            %Update number of tables serving dish k
            M_k(aux_k) = M_k(aux_k) - 1;
            %It could have been the only table serving dish k
            aloneMk = ( M_k(aux_k)==0 );
            if aloneMk %Part of the re-use algorithm
                aux_ind = randi(N_aux);
                Mu_aux(aux_ind) = Mu(aux_k);
                Sigma2_aux(aux_ind) = Sigma2(aux_k);
            end
            
            %Select the tables serving dish k
            table_i = find(T{j}==t);
            table = y{j}(table_i);
            sizetable = length(table);
            
            f_k = log( normpdf(table(1), [Mu_aux Mu], sqrt([Sigma2_aux Sigma2])) );
            if sizetable>1
                for i_aux = 2:sizetable
                    f_k = f_k + log( normpdf(table(i_aux), [Mu_aux Mu], sqrt([Sigma2_aux Sigma2])) );
                end
            end
            
            w = [kappa_0*(u(1)+1)^sigma_0*ones(1,N_aux)/N_aux (M_k - sigma_0)];
            
            %The (M_k - sigma_0) term could be < 0!
            if aloneMk
                w(N_aux + aux_k) = 0;
            end
            
            f_k = exp(f_k-max(f_k)) .* w;
            f_k = f_k/sum(f_k);
            cumsum_f_k = cumsum(f_k);
            hh = find(rand <= cumsum_f_k,1);
            
            if hh <= N_aux
                %New dish at this table
                
                if aloneMk
                    %Same number of dishes in the global menu
                    %Number of tables serving this dish is increased by one
                    M_k(aux_k) = 1;
                    Mu(aux_k) = Mu_aux(hh);
                    Sigma2(aux_k) = Sigma2_aux(hh);
                else
                    Kn = Kn + 1;
                    M_k = [M_k 1];
                    K{j}(t) = Kn;
                    Mu = [Mu Mu_aux(hh)];
                    Sigma2 = [Sigma2 Sigma2_aux(hh)];
                end
                
                %Restore auxiliary variables
                Mu_aux(hh) = m0 + sqrt(s20)*randn;
                Sigma2_aux(hh) = 1/gamrnd(nu0/2,1/(sigma20*nu0/2));
            else
                %Old dish at this table
                h_k = hh - N_aux;
                
                M_k(h_k) = M_k(h_k) + 1;
                K{j}(t) = h_k;
                if aloneMk
                    Kn = Kn - 1;
                    for i_j = 1:J
                        K{i_j}(K{i_j} > aux_k) = K{i_j}(K{i_j} > aux_k) - 1;
                    end
                    M_k(aux_k) = [];
                    Mu(aux_k) = [];
                    Sigma2(aux_k) = [];
                end
            end
        end
    end
    
    
    %Acceleration step
    for k = 1:Kn
        %All the people eating dish k
        dish_k = [];
        for j = 1:J
            dish_k = [dish_k; y{j}(K{j}(T{j})==k)];
        end
        size_k = length(dish_k);
        sum_k = sum(dish_k);
        
        Sigma2(k) = 1/gamrnd((nu0 + size_k)/2,2/(nu0*sigma20 + sum((dish_k - Mu(k)).^2)));
        
        %Conditionally to Sigma2(k)
        m0_post = (sum_k*s20 + m0*Sigma2(k))/(size_k*s20 + Sigma2(k));
        s20_post = Sigma2(k)*s20/(size_k*s20 + Sigma2(k));
        Mu(k) = normrnd(m0_post,sqrt(s20_post));
    end
    
    
    
    %Update u0
    u_new = exp(log(u(1)) + sqrt(s_u(1))*randn);
    
    log_accept = sum(M_k)*(log(u_new) - log(u(1)));
    log_accept = log_accept - kappa_0/sigma_0*((u_new+1)^sigma_0 - (u(1)+1)^sigma_0);
    log_accept = log_accept - (sum(M_k) - Kn*sigma_0) * (log(u_new+1) - log(u(1)+1));
    
    if (isreal(log_accept) == 0 )
        stop(1);
    end
    
    accept = 1;
    if (isnan(log_accept) == 1)
        accept = 0;
    elseif ( log_accept < 0 )
        accept = exp(log_accept);
    end
    
    u_accept(1) = u_accept(1) + accept;
    u_count(1) = u_count(1) + 1;
    
    if ( rand < accept )
        u(1) = u_new;
    end
    
    if iter > n_burn1
        s_u(1) = exp(log(s_u(1))+iter^(-g_w)*(accept-tau));
    end
    
    
    for j = 1:J
        %Update uj
        u_new = exp(log(u(j+1)) + sqrt(s_u(j+1))*randn);
        
        log_accept = n(j)*(log(u_new) - log(u(j+1)));
        log_accept = log_accept - kappa/sigma*((u_new+1)^sigma - (u(j+1)+1)^sigma);
        log_accept = log_accept - (n(j) - M_j(j)*sigma) * (log(u_new+1) - log(u(j+1)+1));
        
        if (isreal(log_accept) == 0 )
            stop(1);
        end
        
        accept = 1;
        if (isnan(log_accept) == 1)
            accept = 0;
        elseif ( log_accept < 0 )
            accept = exp(log_accept);
        end
        
        u_accept(j+1) = u_accept(j+1) + accept;
        u_count(j+1) = u_count(j+1) + 1;
        
        if ( rand < accept )
            u(j+1) = u_new;
        end
        
        if iter > n_burn1
            s_u(j+1) = exp(log(s_u(j+1))+iter^(-g_w)*(accept-tau));
        end
    end
    
    
    
    %Update kappa's
    
    kappa_0 = gamrnd(a_kappa_0 + Kn, 1/( b_kappa_0 + ((u(1)+1)^sigma_0-1)/sigma_0 ));
    kappa = gamrnd(a_kappa + sum(M_k), 1/( b_kappa + sum(((u(2:end)+1).^sigma-1)/sigma) ));
    
    %Update sigma's
    
    %sigma_0
    sigma_new = exp(log(sigma_0/(1-sigma_0)) + sqrt(s_sigma_0)*randn);
    sigma_new = 1/(1+1/sigma_new);
    
    log_accept = a_sigma_0*(log(sigma_new) - log(sigma_0)) + b_sigma_0*(log(1-sigma_new) - log(1-sigma_0));
    log_accept = log_accept + Kn * (log(u(1)+1)*(sigma_new-sigma_0) + gammaln(1-sigma_0) - gammaln(1-sigma_new));
    log_accept = log_accept - kappa_0*(((u(1)+1)^sigma_new-1)/sigma_new - ((u(1)+1)^sigma_0-1)/sigma_0);
    log_accept = log_accept + sum(gammaln(M_k - sigma_new) - gammaln(M_k - sigma_0));
    
    if (isreal(log_accept) == 0 )
        stop(1);
    end
    
    accept = 1;
    if (isnan(log_accept) == 1)
        accept = 0;
    elseif ( log_accept < 0 )
        accept = exp(log_accept);
    end
    
    sigma_0_accept = sigma_0_accept + accept;
    sigma_0_count = sigma_0_count + 1;
    
    if ( rand < accept )
        sigma_0 = sigma_new;
    end
    
    if iter > n_burn1
        s_sigma_0 = exp(log(s_sigma_0)+iter^(-g_w)*(accept-tau));
    end
    
    
    %sigma
    sigma_new = exp(log(sigma/(1-sigma)) + sqrt(s_sigma)*randn);
    sigma_new = 1/(1+1/sigma_new);
    
    
    log_accept = a_sigma*(log(sigma_new) - log(sigma)) + b_sigma*(log(1-sigma_new) - log(1-sigma));
    log_accept = log_accept + sum(log(u(2:end)+1) .* M_j)*(sigma_new-sigma) + sum(M_k) * (gammaln(1-sigma) - gammaln(1-sigma_new));
    log_accept = log_accept - kappa*sum( ((u(2:end)+1).^sigma_new-1)/sigma_new - ((u(2:end)+1).^sigma-1)/sigma );
    log_accept = log_accept + sum(gammaln([N{:}] - sigma_new) - gammaln([N{:}] - sigma));
    
    if (isreal(log_accept) == 0 )
        stop(1);
    end
    
    accept = 1;
    if (isnan(log_accept) == 1)
        accept = 0;
    elseif ( log_accept < 0 )
        accept = exp(log_accept);
    end
    
    sigma_accept = sigma_accept + accept;
    sigma_count = sigma_count + 1;
    
    if ( rand < accept )
        sigma = sigma_new;
    end
    
    if iter > n_burn1
        s_sigma = exp(log(s_sigma)+iter^(-g_w)*(accept-tau));
    end
    
    
    
    %Saving the output
    if iter > (n_burn1 + n_burn2) && mod((iter - n_burn1 - n_burn2),thin) == 0
        
        iter_aux = (iter - n_burn1 - n_burn2)/thin;
        
        Mu_out{iter_aux} = Mu;
        Sigma2_out{iter_aux} = Sigma2;
        for j = 1:J
            T_out{iter_aux,j} = T{j};
            K_out{iter_aux,j} = K{j};
            N_out{iter_aux,j} = N{j};
        end
        M_j_out{iter_aux} = M_j;
        M_k_out{iter_aux} = M_k;
        kappa_0_out(iter_aux) = kappa_0;
        sigma_0_out(iter_aux) = sigma_0;
        kappa_out(iter_aux) = kappa;
        sigma_out(iter_aux) = sigma;
        Kn_out(iter_aux) = Kn;
        U_out(iter_aux,:) = u;
    end
end
toc
%profile viewer

%%
save('SCHOOL_HNGG')

%% PLOTS
load('SCHOOL_HNGG')

%%
%Compute the Binder's partition

%Posterior inclusion probabilities
%Based on overall dish indices
pihat_II = zeros(sum(n),sum(n));

for g = 1:n_save
    cij_II = zeros(sum(n),sum(n));
    Z = [];
    for j = 1:J
        %Construct big adjacency matrix
        Z = [Z, K_out{g,j}(T_out{g,j})];
    end
    for j1 = 1:sum(n)
        for j2 = 1:j1
            cij_II(j1,j2) = (Z(j1) == Z(j2));
            cij_II(j2,j1) = cij_II(j1,j2);
        end
    end
    pihat_II = pihat_II + cij_II;
end
pihat_II = pihat_II/n_save;

%Funzione di perdita di Lau-Green (Binder)
Binder_II = zeros(1,n_save);

%Binder LF's parameters are equal in this case
k1 = 1;
k2 = 1;
Binder_coeff = k1/(k1+k2);  %quando è uguale a un mezzo equivale a Dahl (Binder con costi unitari)
%però esistono infinite combinazioni che danno un mezzo!!

for g = 1:n_save
    cij_II = zeros(sum(n),sum(n));
    Z = [];
    for j = 1:J
        %Construct big adjacency matrix
        Z = [Z, K_out{g,j}(T_out{g,j})];
    end
    for j1 = 1:sum(n)
        for j2 = 1:j1
            cij_II(j1,j2) = (Z(j1) == Z(j2));
            cij_II(j2,j1) = cij_II(j1,j2);
        end
    end
    Binder_II(g) = sum(sum(tril(pihat_II-Binder_coeff,-1).*tril(cij_II,-1)));
end
%Consider the half MCMC sample already used!
max_g_II = find(Binder_II == max(Binder_II),1);
% max_g_II1 = n_save;

%Binder's partition
sel_part_Binder = [];
for j = 1:J
    %Construct big adjacency matrix
    sel_part_Binder = [sel_part_Binder, K_out{max_g_II,j}(T_out{max_g_II,j})];
end
save('sel_part_Binder_all.mat', 'sel_part_Binder', 'max_g_II')

%%
%Plot data by increasing sample mean, and colours based on Binder's partition

%Compute sample means of data groups
mean_data = zeros(1,J);
for j = 1:J
    mean_data(j) = mean(y{j});
end
%Sort in ascending order
[mean_data_sorted mean_data_sorted_id] = sort(mean_data);

%Plot data in ascending mean order and clustering colours
Color_Clust = varycolor(length(unique(sel_part_Binder)));
%We have 5 clusters
Markers = {'o','s','d','^','*'};
figure(1)
hold on
for j = 1:J
    aux_j =  mean_data_sorted_id(j);
    y_j = y{aux_j};
    
    partition_j = K_out{max_g_II,aux_j}(T_out{max_g_II,aux_j});
    
    plot([j j], [min(y_j) max(y_j)], '-', 'Color', [160 160 160]/253)
    for i = 1:n(aux_j)
        plot(j, y_j(i), Markers{partition_j(i)}, 'Color', Color_Clust(partition_j(i),:), 'MarkerSize', 3, 'MarkerFaceColor', Color_Clust(partition_j(i),:), 'MarkerSize', 3.5)
    end
end
xlabel('Schools (ordered)', 'Fontsize', 15)
ylabel('Math Scores', 'Fontsize', 15)
axis tight
xlim([-1 101])
print(1,'-djpeg','SCHOOL_ClustMeans_Random1.jpeg')
print(1,'-depsc','SCHOOL_ClustMeans_Random1.eps')
print(1,'-deps','SCHOOL_ClustMeans_Random1_bw.eps')
close(1)

%%

%Plot of posterior expected number of dishes in eaxch group
figure(1)
Kj = zeros(1, J);
for j = 1:J
    Kj_aux = zeros(1, n_save);
    for it = 1:n_save
        Kj_aux(it) = length(K_out{it,j});
    end
    Kj(j) = mean(Kj_aux);
end

[n_sort n_sort_id] = sort(n);

n_sort_id = 1:J;
bar(1:J, Kj(n_sort_id), 0.5, 'k')
xlim([0 J+1])
xlabel('School', 'Fontsize', 15)
print(1,'-djpeg','SCHOOL_TablesMean_Random1.jpeg')
print(1,'-depsc','SCHOOL_TablesMean_Random1.eps')
print(1,'-depsc','SCHOOL_TablesMean_Random1_bw.eps')
close(1)

%%

%Density estimation in selected schools

[Kj_sort Kj_sort_id] = sort(Kj);

Kj_max = Kj(Kj_sort_id(end-2:end));
id_max = zeros(1,3);
id_max(1) = find(Kj == Kj_max(1));
id_max(2) = find(Kj == Kj_max(2));
id_max(3) = find(Kj == Kj_max(3));

Kj_min = Kj(Kj_sort_id(1:3));
id_min = zeros(1,3);
id_min(1) = find(Kj == Kj_min(1));
id_min(2) = find(Kj == Kj_min(2));
id_min(3) = find(Kj == Kj_min(3));

sel_schools = [id_min id_max]; %Increasing! Kj([id_min id_max])
n_schools = length(sel_schools);

%Predictive
dg = (max(school)-min(school))/100;
pred_grid = (min(school)-5):dg:(max(school)+5);
l_grid = length(pred_grid);
pred = zeros(n_schools,l_grid);

for it = 1:n_save    
    %Predictive in source group
    for j_aux = 1:n_schools
        j = sel_schools(j_aux);
        
        %Compute weights as in algorithm
        w = [kappa_0_out(it) * (U_out(it,1)+1)^sigma_0_out(it) (M_k_out{it} - sigma_0_out(it))];
        w = [kappa_out(it) * (U_out(it,j+1)+1)^sigma_out(it) * w/sum(w) (N_out{it,j} - sigma_out(it))];
        w = w/sum(w);
        
        n_int = 100; %Approximate integral from P0
        Mu_aux = m0 + sqrt(s20)*randn(1,n_int);
        Sigma2_aux = 1./gamrnd(nu0/2, 2/sigma20/nu0, 1, n_int);
        for i_int = 1:n_int
            pred(j_aux,:) = pred(j_aux,:) + w(1) * normpdf(pred_grid, Mu_aux(i_int), sqrt(Sigma2_aux(i_int)))/n_int;
        end
        
        for i = 1:Kn_out(it)
            pred(j_aux,:) = pred(j_aux,:) + w(1 + i) * normpdf(pred_grid, Mu_out{it}(i), sqrt(Sigma2_out{it}(i)));
        end
        
        for i = 1:M_j_out{it}(j)
            pred(j_aux,:) = pred(j_aux,:) + w(1 + Kn_out(it) + i) * normpdf(pred_grid, Mu_out{it}(K_out{it,j}(i)), sqrt(Sigma2_out{it}(K_out{it,j}(i))));
        end       
    end
end
pred = pred/n_save;
sum(pred'*dg)

%Plot data in ascending mean order and clustering colours
Color_Clust = varycolor(length(unique(sel_part_Binder)));
greys = [225:-40:5]';
Color_Schools = [greys greys greys]/255;
eps = [-0.005:-0.005:-0.03];
Markers = {'o','s','d','^','*'};

figure(1)
hold on
for j_aux = 1:n_schools
    plot(pred_grid,pred(j_aux,:),'-', 'LineWidth',3,'Color',Color_Schools(j_aux,:))
end
legend({['School ', num2str(sel_schools(1))],['School ' num2str(sel_schools(2))],['School ' num2str(sel_schools(3))],['School ' num2str(sel_schools(4))],['School ' num2str(sel_schools(5))],['School ' num2str(sel_schools(6))]}, 'Fontsize', 12.5)

plot([15 90], [0 0], '-k', 'LineWidth',1)

for j_aux = 1:n_schools
    j =  sel_schools(j_aux);
    y_j = y{j};
    
    partition_j = K_out{max_g_II,j}(T_out{max_g_II,j});
    
    plot([15 90], [eps(j_aux) eps(j_aux)], '-', 'LineWidth',1,'Color',Color_Schools(j_aux,:))
    for i = 1:n(j)
        plot(y_j(i), eps(j_aux), Markers{partition_j(i)}, 'Color', Color_Clust(partition_j(i),:), 'MarkerSize', 5, 'MarkerFaceColor', Color_Clust(partition_j(i),:))
    end
end
axis tight
ylim([-0.035 0.055])
xlabel('Math Score', 'Fontsize', 15)
ylabel('Predictive density', 'Fontsize', 15)
aux = get(gca,'ytick');
labels = {[], [], [], num2str(aux(4:end)')};
set(gca, 'ytick',aux,'yticklabel',labels)
print(1,'-djpeg','SCHOOL_Pred_Groups_MinMax_Random1.jpeg')
print(1,'-depsc','SCHOOL_Pred_Groups_MinMax_Random1.eps')
print(1,'-deps','SCHOOL_Pred_Groups_MinMax_Random1_bw.eps')
close(1)

%%
%Predictive new group
figure(1)
%Data in background
dx = (max(school)-min(school))/25;
xx = min(school)-5:dx:max(school)+5;    [f_k] = hist(school,xx);
bar(xx,f_k/sum(f_k)/dx,'FaceColor',[100 100 100]/255,'edgecolor','none');
hold on
plot(pred_grid,pred(n_schools+1,:),'LineWidth',3,'Color',[0 128 255]/255)
axis tight
xlabel('Math Score', 'Fontsize', 15)
print(1,'-djpeg','SCHOOL_Pred_New_Random1.jpeg')
print(1,'-depsc','SCHOOL_Pred_New_Random1.eps')
print(1,'-deps','SCHOOL_Pred_New_Random1_bw.eps')
close(1)
