%%%
%Simple Univariate model using Normal-Inverse gamma variables and HNGG
%%%

clear all
addpath('../../../../../Codes')
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
y = y_simul;

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
kappa_0 = 0.1;

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
kappa = 0.1;

% mm = .5;
% ss = sqrt(1/20);
% a_sigma = (1-mm)*mm^2/ss^2-mm;
% b_sigma = a_sigma*(1-mm)/mm;
a_sigma = 2;
b_sigma = 18;
sigma = .2;
%Non-conjugate
s_sigma = 1;
sigma_accept = 0;
sigma_count = 0;

%Hyperparameters of P0
%m0
m1 = 0;
s1 = 1000;
m0 = 0.25;

%k0
a_k0 = 1;
b_k0 = 1;
k0 = 0.62;


%Hyperparameters of sigma2
a_sigma2 = 2.07;
b_sigma2 = 0.66;


%All together
for j = 1:J
    T{j} = ones(1,n(j));
    N{j} = n(j);
end
K{1} = 1;
K{2} = 1;
M_k = J;
Kn = 1;

%Number of tables in each restaurant
M_j = ones(1,J);

Sigma2 = 0.1*ones(1,Kn);
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
k0_out = zeros(1,n_save);
m0_out = zeros(1,n_save);
%%
% %Predictive
% dg = (max([y_simul{:}])-min([y_simul{:}]))/100;
% pred_grid = (min([y_simul{:}])-2):dg:(max([y_simul{:}])+2);
% l_grid = length(pred_grid);
% pred_mean = zeros(J+1,l_grid);

tic
%profile on
for iter = 1:n_iter
    
    
    if mod(iter,100) == 0
        iter
    end
    
    %For element ji, sample the table value, after removing the element from
    %the list
    
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
            
            
            %Weights for each dish
            f_k = -inf*ones(1,Kn+1);
            
            %In the first position is the probability of a new table
            scale_one = b_sigma2/a_sigma2*(k0+1)/k0;
            f_k(1) = log( tpdf((y{j}(i)-m0)/sqrt(scale_one),2*a_sigma2)/sqrt(scale_one) );
            
            for k = 1:Kn
                
                %People eating dish k excluding data point y_ji
                den = [];
                for j_aux = 1:J
                    index = find(K{j_aux}(T{j_aux})==k);
                    if j_aux == j
                        index = setdiff(index, i);
                    end
                    den = [den y{j_aux}(index)];
                end
                
                if ~isempty(den)
                    sizeden = length(den);
                    sumden = sum(den);
                    meanden = sumden/sizeden;
                    
                    meanboth = (y{j}(i) + sumden)/(1 + sizeden);
                    
                    k0_den = k0 + sizeden;
                    %k0_aux = k0_den + 1;
                    m0_den = (k0*m0 + sumden)/k0_den;
                    %m0_aux = (k0*m0 + sumden + y{j}(i))/k0_aux;
                    a_den = a_sigma2 + sizeden/2;
                    %a_aux = a_den + 1/2;
                    b_den = b_sigma2 + .5*(sum((den - meanden).^2) + k0*sizeden*(meanden - m0)^2/k0_den);
                    %b_aux = b_sigma2 + .5*(sum(([y{j}(i) den] - meanboth).^2) + k0*(sizeden + 1)*(meanboth - m0)^2/k0_aux);
                    scale_aux = b_den/a_den*(k0_den + 1)/k0_den;
                    
                    f_k(k + 1) = log( tpdf((y{j}(i)-m0_den)/sqrt(scale_aux),2*a_den)/sqrt(scale_aux) );
                    %gammaln(a_aux) - gammaln(a_den) + .5*(log(k0_den) - log(k0_aux)) - 1/2*(log(2*pi)) + a_den*log(b_den) - a_aux*log(b_aux);
                end
            end
            
            f_k = f_k([1:Kn+1 K{j}+1]);
            w = [kappa_0 * (u(1)+1)^sigma_0 (M_k-sigma_0)];
            w = [kappa * (u(j+1)+1)^sigma * w/sum(w) (N{j}-sigma)];
            
            %The (N{j} - sigma) term could be < 0!
            if aloneN
                w(Kn + 1 + aux_t) = 0;
                %Same thing for the (M_k - sigma_0) term!
                if aloneMk
                    w(aux_k + 1) = 0;
                end
            end
            
            f_k = exp(f_k-max(f_k)) .* w;
            f_k = f_k/sum(f_k);
            
            cumsum_f_k = cumsum(f_k);
            hh = find(rand <= cumsum_f_k,1);
            
            if hh == 1
                if aloneN
                    %Select new dish for the new table but use same
                    %label (the cluster was empty)
                    N{j}(aux_t) = 1;
                    
                    if aloneMk
                        %Same number of dishes in the global menu
                        %Number of tables serving this dish is now one
                        M_k(aux_k) = 1;
                    else
                        Kn = Kn + 1;
                        M_k = [M_k 1];
                        K{j}(aux_t) = Kn;
                    end
                else
                    Kn = Kn + 1;
                    M_k = [M_k 1];
                    M_j(j) = M_j(j) + 1;
                    T{j}(i) = M_j(j);
                    K{j} = [K{j} Kn];
                    N{j} = [N{j} 1];
                end
            elseif hh <= Kn+1
                %New table, old dish
                h_k = hh - 1;
                
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
                    end
                else
                    M_j(j) = M_j(j) + 1;
                    T{j}(i) = M_j(j);
                    K{j} = [K{j} h_k];
                    N{j} = [N{j} 1];
                end
            else
                %Old table, old dish
                h_t = hh - Kn - 1;
                
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
            
            f_k = -inf*ones(1,Kn+1);
            
            %Select the tables serving dish k
            table_i = find(T{j}==t);
            table = y{j}(table_i);
            sizetable = length(table);
            sumtable = sum(table);
            meantable = sumtable/sizetable;
            
            k0_table = k0 + sizetable;
            m0_table = (k0*m0 + sumtable)/k0_table;
            a_table = a_sigma2 + sizetable/2;
            b_table = b_sigma2 + .5*(sum((table - meantable).^2) + k0*sizetable*(meantable - m0)^2/k0_table);
            
            f_k(1) = gammaln(a_table) - gammaln(a_sigma2) + .5*(log(k0) - log(k0_table)) - sizetable/2*(log(2*pi)) + a_sigma2*log(b_sigma2) - a_table*log(b_table);
            
            for k = 1:Kn
                
                %People eating dish k excluding data points table_ji
                den = [];
                for j_aux = 1:J
                    index = find(K{j_aux}(T{j_aux})==k);
                    if j_aux == j
                        index = setdiff(index, table_i);
                    end
                    den = [den y{j_aux}(index)];
                end
                
                if ~isempty(den)
                    sizeden = length(den);
                    sumden = sum(den);
                    meanden = sumden/sizeden;
                    
                    meanboth = (sumtable + sumden)/(sizetable + sizeden);
                    
                    k0_den = k0 + sizeden;
                    k0_aux = k0_den + sizetable;
                    m0_den = (k0*m0 + sumden)/k0_den;
                    m0_aux = (k0*m0 + sumden + sumtable)/k0_aux;
                    a_den = a_sigma2 + sizeden/2;
                    a_aux = a_den + sizetable/2;
                    b_den = b_sigma2 + .5*(sum((den - meanden).^2) + k0*sizeden*(meanden - m0)^2/k0_den);
                    b_aux = b_sigma2 + .5*(sum(([table den] - meanboth).^2) + k0*(sizetable + sizeden)*(meanboth - m0)^2/k0_aux);
                    
                    f_k(k + 1) = gammaln(a_aux) - gammaln(a_den) + .5*(log(k0_den) - log(k0_aux)) - sizetable/2*(log(2*pi)) + a_den*log(b_den) - a_aux*log(b_aux);
                end
            end
            
            w = [kappa_0*(u(1)+1)^sigma_0 (M_k-sigma_0)];
            
            %The (M_k - sigma_0) term could be < 0!
            if aloneMk
                w(aux_k + 1) = 0;
            end
            
            f_k = exp(f_k-max(f_k)) .* w;
            f_k = f_k/sum(f_k);
            
            cumsum_f_k = cumsum(f_k);
            hh = find(rand <= cumsum_f_k,1);
            
            if  hh == 1
                if aloneMk
                    %Same number of dishes in the global menu
                    %Number of tables serving this dish is increased by one (equal to 1)
                    M_j(aux_k) = 1;
                    M_k(aux_k) = 1;
                else
                    Kn = Kn + 1;
                    M_k = [M_k 1];
                    K{j}(t) = Kn;
                end
            else
                %Old dish at this table
                h_k = hh - 1;
                
                M_k(h_k) = M_k(h_k) + 1;
                K{j}(t) = h_k;
                if aloneMk
                    Kn = Kn - 1;
                    for i_j = 1:J
                        K{i_j}(K{i_j} > aux_k) = K{i_j}(K{i_j} > aux_k) - 1;
                    end
                    M_k(aux_k) = [];
                end
            end
        end
    end
    
    
    %Acceleration step
    Sigma2 = zeros(1,Kn);
    Mu = zeros(1,Kn);
    for k = 1:Kn
        %All the people eating dish k
        dish_k = [];
        for j = 1:J
            dish_k = [dish_k, y{j}(K{j}(T{j})==k)];
        end
        size_k = length(dish_k);
        sum_k = sum(dish_k);
        mean_k = sum_k/size_k;
        k0_aux = k0 + size_k;
        
        Sigma2(k) = 1/gamrnd(a_sigma2 + size_k/2,1/( b_sigma2 + .5*(sum((dish_k - mean_k).^2) + k0*size_k*(mean_k - m0)^2/k0_aux) ));
        Mu(k) = normrnd((k0*m0+sum_k)/k0_aux,sqrt(Sigma2(k)/k0_aux));
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
    
    %     kappa_0 = gamrnd(a_kappa_0 + Kn, 1/( b_kappa_0 + ((u(1)+1)^sigma_0-1)/sigma_0 ));
    % kappa = gamrnd(a_kappa + sum(M_k), 1/( b_kappa + sum(((u(2:end)+1).^sigma-1)/sigma) ));
    
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
        k0_out(iter_aux) = k0;
        m0_out(iter_aux) = m0;
    end
end
toc
%profile viewer

%%
save('N_HNGG')

%% PLOTS
load('N_HNGG')

%%
%Predictive
dg = (max([y_simul{:}])-min([y_simul{:}]))/100;
pred_grid = (min([y_simul{:}])-2):dg:(max([y_simul{:}])+2);
l_grid = length(pred_grid);
pred = zeros(J+1,l_grid);


for it = 1:n_save
    
    %For new group (but point-wise!)'
    nu_one = 2*a_sigma2;
    m0_one = m0_out(it);
    scale_one = b_sigma2/a_sigma2*(k0_out(it) + 1)/k0_out(it);
    
    %Quantities needed for update
    m0_jk = zeros(1,Kn_out(it));
    nu_jk = zeros(1,Kn_out(it));
    scale_jk = zeros(1,Kn_out(it));
    
    for jk = 1:Kn_out(it)
        %All the people eating dish jk
        dish_jk = [];
        for ind_j = 1:J
            dish_jk = [dish_jk, y{ind_j}(K_out{it,ind_j}(T_out{it,ind_j})==jk)];
        end
        
        size_jk = length(dish_jk);
        sum_jk = sum(dish_jk);
        mean_jk = sum_jk/size_jk;
        
        k0_jk = k0_out(it) + size_jk;
        m0_jk(jk) = (k0_out(it)*m0_out(it) + sum_jk)/k0_jk;
        a_jk = a_sigma2 + size_jk/2;
        nu_jk(jk) = 2*a_jk;
        b_jk = (b_sigma2 + .5*(sum((dish_jk - mean_jk).^2) + k0_out(it)*size_jk*(mean_jk - m0_out(it))^2/k0_jk));
        scale_jk(jk) = b_jk/a_jk*(k0_jk + 1)/k0_jk;
    end
    
    %Predictive in source group
    for j = 1:J
        
        %Compute weights as in algorithm
        nu_it = [nu_one nu_jk([1:Kn_out(it) K_out{it,j}])];
        m0_it = [m0_one m0_jk([1:Kn_out(it) K_out{it,j}])];
        scale_it = [scale_one scale_jk([1:Kn_out(it) K_out{it,j}])];
        w = [kappa_0_out(it) * (U_out(it,1)+1)^sigma_0_out(it) (M_k_out{it} - sigma_0_out(it))];
        w = [kappa_out(it) * (U_out(it,j+1)+1)^sigma_out(it) * w/sum(w) (N_out{it,j} - sigma_out(it))];
        w = w/sum(w);
        
        for ii = 1:length(w)
            pred(j,:) = pred(j,:) + w(ii) * tpdf((pred_grid - m0_it(ii))/sqrt(scale_it(ii)),nu_it(ii))/sqrt(scale_it(ii));
        end
    end
    %Predictive of new restaurant (dish level?)
    p_do = (M_k_out{it} - sigma_0_out(it))/(kappa_0_out(it) * (U_out(it,1)+1)^sigma_0_out(it) + sum(M_k_out{it}) - Kn_out(it) * sigma_0_out(it));
    for jk = 1:Kn_out(it)
        pred(J+1,:) = pred(J+1,:) + p_do(jk) * tpdf((pred_grid-m0_jk(jk))/sqrt(scale_jk(jk)),nu_jk(jk))/sqrt(scale_jk(jk));
    end
    p_dn = kappa_0_out(it) * (U_out(it,1)+1)^sigma_0_out(it)/(kappa_0_out(it) * (U_out(it,1)+1)^sigma_0_out(it) + sum(M_k_out{it}) - Kn_out(it) * sigma_0_out(it));
    pred(J+1,:) = pred(J+1,:) + p_dn * tpdf((pred_grid-m0_one)/sqrt(scale_one),2*a_sigma2)/sqrt(scale_one);
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
    [f_k] = hist(y_simul{j}(1:n1(j)),xx);
    bar(xx,p1_simul(j)*f_k/sum(f_k)/dx,'FaceColor',Color_Hist((j-1)*2 + 1,:),'facealpha',.5,'edgecolor','none');
    [f_k] = hist(y_simul{j}(n1(j)+1:n(j)),xx);
    bar(xx,(1-p1_simul(j))*f_k/sum(f_k)/dx,'FaceColor',Color_Hist((j-1)*2 + 2,:),'facealpha',.5,'edgecolor','none');
end
xx = pred_grid;
plot(xx,pred(1,:),'-','LineWidth',3)
plot(xx,pred(2,:),'-.','LineWidth',3)
axis tight
legend({'Comp 1','Comp 2 (Group 1)','Comp 2 (Group 2)','Comp 3','Predictive (Group 1)','Predictive (Group 2)'},'FontSize',15)
xlabel('Y', 'Fontsize', 15)
print(1,'-djpeg','N_HNGG_Pred_groups.jpeg')
print(1,'-depsc','N_HNGG_Pred_groups.eps')
print(1,'-deps','N_HNGG_Pred_groups_bw.eps')
close(1)


%Predictive new group
figure(1)
dx = (max([y_simul{:}])-min([y_simul{:}]))/25;
xx = min([y_simul{:}])-.5:dx:max([y_simul{:}])+.5;
hold on
[f_k] = hist([y_simul{:}],xx);
bar(xx,f_k/sum(f_k)/dx,'FaceColor',Color_Hist((j-1)*2 + 1,:),'facealpha',.5,'edgecolor','none');
xx = pred_grid;
plot(xx,pred(J+1,:),'LineWidth',3)
axis tight
legend({'Data','Predictive new Group'}, 'Fontsize', 15)
xlabel('Y', 'Fontsize', 15)
print(1,'-djpeg','N_HNGG_Pred_new.jpeg')
print(1,'-depsc','N_HNGG_Pred_new.eps')
print(1,'-deps','N_HNGG_Pred_new_bw.eps')
close(1)


%% 
%Posterior total number of dishes
%AND
%Posterior number of clusters in restaurants

%Total number of dishes - Histogram
figure(1)
xx = unique(Kn_out);
ff = hist(Kn_out,xx);
bar(xx,ff/sum(ff));
axis tight
legend({'Number of Clusters (both Groups)'}, 'Fontsize', 15)
print(1,'-djpeg','N_HNGG_Kn_Sigma_Rand.jpeg')
print(1,'-depsc','N_HNGG_Kn_Sigma_Rand.eps')
print(1,'-deps','N_HNGG_Kn_Sigma_Rand_bw.eps')
close(1)

%Num dishes in each restaurant
figure(2)
NumDish = zeros(n_save,J);
for g = 1:n_save
    for j = 1:J
        NumDish(g,j) = length(unique(K_out{g,j}));
    end
end
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
legend({'Number of Clusters (Group 1)', 'Number of Clusters (Group 2)'}, 'Fontsize', 15)
print(2,'-djpeg','N_HNGG_NumDishes_Sigma_Rand.jpeg')
print(2,'-depsc','N_HNGG_NumDishes_Sigma_Rand.eps')
print(2,'-deps','N_HNGG_NumDishes_Sigma_Rand_bw.eps')
close(2)



%Posterior distribution of sigma, sigma0
figure(1)
dx = (max(sigma_out)-min(sigma_out))/12;
xx = min(sigma_out):dx:1.1;
hold on
[f_k] = hist(sigma_out,xx);
bar(xx,f_k/sum(f_k)/dx,'b','facealpha',.5,'edgecolor','none');
ff = ksdensity(sigma_out, xx);
plot(xx,ff,'LineWidth',1.5)
axis tight
xlim([0 1])
legend({'MCMC sample', 'Kernel Density Est.'}, 'Fontsize', 15)
xlabel('\sigma', 'Fontsize', 15)
print(1,'-djpeg','N_HNGG_SigmaRandom.jpeg')
print(1,'-depsc','N_HNGG_SigmaRandom.eps')
print(1,'-deps','N_HNGG_SigmaRandom_bw.eps')
close(1)

figure(1)
dx = (max(sigma_0_out)-min(sigma_0_out))/10;
xx = min(sigma_0_out):dx:1.1;
hold on
[f_k] = hist(sigma_0_out,xx);
bar(xx,f_k/sum(f_k)/dx,'b','facealpha',.5,'edgecolor','none');
ff = ksdensity(sigma_0_out, xx);
plot(xx,ff,'LineWidth',1.5)
axis tight
xlim([0 1])
legend({'MCMC sample', 'Kernel Density Est.'}, 'Fontsize', 15)
xlabel('\sigma_0', 'Fontsize', 15)
print(1,'-djpeg','N_HNGG_Sigma0Random.jpeg')
print(1,'-depsc','N_HNGG_Sigma0Random.eps')
print(1,'-deps','N_HNGG_Sigma0Random_bw.eps')
close(1)

%%
%LPML

CPO_j = cell(1,J);
CPO = zeros(1, sum(n));
for it = 1:n_save
    for j = 1:J
        index = K_out{it,j}(T_out{it,j});
        aux = normpdf(y{j}, Mu_out{it}(index), sqrt(Sigma2_out{it}(index)));
        CPO_j{j} = 1 ./ aux;
    end
    CPO = CPO + [CPO_j{:}];
end
LPML = sum(log( n_save ./ CPO ));
