function [TP_final, iteration, epsilon, terminatedNormally, performance] = Spare_Decomposition(machines, C, Q, mu, p, gamma, justresults)
%

%% check input arguments
if (nargin < 6 || nargin > 7)
    fprintf('Wrong number of parameters... expecting: machines, C, Q, p, mu, gamma\n\n');
    machines = 3;
    C = 30 * ones(1,machines-1);
    Q = 5 * ones(1,machines);
    mu = 1 * ones(1,machines) * .9;
    p = 0.005 * ones(1,machines) * .1;
    gamma = 0.05 * ones(1,machines) * .9;
    justresults = false;
end

if (~checkSize(C, [1 machines-1]) || ~checkSize(Q, [1 machines]) || ~checkSize(mu, [1 machines]) || ~checkSize(p, [1 machines]) || ~checkSize(gamma, [1 machines]))
    C = C'; Q = Q'; mu = mu'; p = p'; gamma = gamma';
    if (~checkSize(C, [1 machines-1]) || ~checkSize(Q, [1 machines]) || ~checkSize(mu, [1 machines]) || ~checkSize(p, [1 machines]) || ~checkSize(gamma, [1 machines]))
        error("Dimensions mismatch!");
    end
end


%% settings
echoLevel = 3;
pauseAfterIteration = false;
sleepAfterIteration = false;
epsilon = 10^(-4);
omega = 0.5;


%% initialisation
terminatedNormally = true;
if (justresults)
    echoLevel = 0;
    pauseAfterIteration = false;
    sleepAfterIteration = false;
end
TP = NaN(1,machines-1);
TP_old = TP;
mu_u = mu(1:machines-1);
p_u = p(1:machines-1);
gamma_u = gamma(1:machines-1);
Q_u = Q(1:machines-1);
mu_d = mu(2:machines);
p_d = p(2:machines);
gamma_d = gamma(2:machines);
Q_d = Q(2:machines);
MaxSparepartsStation = max(Q);
p_up_empty = zeros(machines-1, MaxSparepartsStation);
p_dn_full = zeros(machines-1, MaxSparepartsStation);
for i = 1:machines-1
    [p_up_empty(i,:), p_dn_full(i,:), pUpMaFail(i), pDnMaFail(i), ps(i), pb(i), p001(i), p011(i), p101(i), Sp00A(i), Sp01A(i), pNm1_10(i), pN10(i), pN11(i), SpNA0(i), SpNA1(i), pSum1(i), pSum2(i), TP(i), nb(i), prob]=corr_TML_MSP(mu_u(i),p_u(i),gamma_u(i),mu_d(i),p_d(i),gamma_d(i),Q_u(i),Q_d(i),C(i), MaxSparepartsStation);
end

iteration = 0;
TP_diff = inf;
TP_diff_cyc = NaN;
cycling = false;
data = NaN(1, 1000);


%% decomposition
if (echoLevel >= 1)
    TP_str = strtrim(sprintf('%2.4f ', TP));
    fprintf('Starting decomposition (iteration #%i).\n  TP: \t\t[%s]\n\n', iteration, TP_str);
end

while (TP_diff > epsilon)
    iteration = iteration + 1;
    
    % Forward pass
    for i=2:machines-1
        iterationFP = 0;
        NoConvergence = 1;
        current_omega = omega;
        
        [p_up_empty(i,:), p_dn_full(i,:), pUpMaFail(i), pDnMaFail(i), ps(i), pb(i), p001(i), p011(i), p101(i), Sp00A(i), Sp01A(i), pNm1_10(i), pN10(i), pN11(i), SpNA0(i), SpNA1(i), pSum1(i), pSum2(i), TP(i), nb(i), sol, IS1(i), IS2(i)]=corr_TML_MSP(mu_u(i),p_u(i),gamma_u(i),mu_d(i),p_d(i),gamma_d(i),Q_u(i),Q_d(i),C(i), MaxSparepartsStation);
        
        while NoConvergence
            % Fixed point iteration
            iterationFP = iterationFP + 1;
            local_p_u = p_u(i);
            local_gamma_u = gamma_u(i);
            local_mu_u = mu_u(i);
            local_TP_i = TP(i);
            
            B_i = 0;
            B_u_i = 0;
            B_d_im1 = 0;
            for alpha = 1 : Q(i)
                B_i = B_i + factorial(Q(i))/factorial(Q(i) - alpha) * (gamma(i) / p(i) )^alpha;
                B_u_i = B_u_i + factorial(Q(i))/factorial(Q(i) - alpha) * (local_gamma_u/ local_p_u )^alpha;
                B_d_im1 = B_d_im1 + factorial(Q(i))/factorial(Q(i)-alpha) * (gamma_d(i-1) / p_d(i-1) )^alpha;
                for l = 1 : alpha - 1
                    B_i = B_i + factorial(Q(i) - l) / factorial(Q(i) - alpha) * (gamma(i) / p(i) )^(alpha-l) *( p_up_empty(i-1,l) / pDnMaFail(i-1) + p_dn_full(i,l) / pUpMaFail(i) );
                    B_u_i = B_u_i + factorial(Q(i) - l) / factorial(Q(i) - alpha) * (local_gamma_u/ local_p_u )^(alpha-l) *( p_dn_full(i,l) / pUpMaFail(i) );
                    B_d_im1 = B_d_im1 + factorial(Q(i) - l) / factorial(Q(i) - alpha) * (gamma_d(i-1) / p_d(i-1) )^(alpha-l) *( p_up_empty(i-1,l) / pDnMaFail(i-1) );
                end
            end
            p_u(i) = current_omega * local_p_u + (1 - current_omega) * (p(i) + mu_d(i-1) * min(max(p101(i-1), p101(i-1)/ (TP(i-1)/local_mu_u - pSum1(i))),1 )+ p_u(i-1) * min(max(Sp01A(i-1), Sp01A(i-1) / (TP(i-1)/local_mu_u - pSum1(i))),1 ));
            gamma_u(i) = current_omega * local_gamma_u + (1 - current_omega) * (gamma(i) + ( Q(i-1)/Q(i)*gamma_u(i-1) - gamma(i))* min(max(Sp00A(i-1), Sp00A(i-1)/ (local_p_u / local_gamma_u / Q(i) * ( TP(i-1) /local_mu_u - pSum1(i)))), 1) );
            mu_u(i) = current_omega * local_mu_u + (1 - current_omega) * (((1+B_u_i)/B_u_i)/ ( 1 / TP(i-1) + 1 / mu(i) * (1+B_i)/B_i - 1 / mu_d(i-1) * (1+B_d_im1)/B_d_im1 ));
            
            % Updating the values of the current line
            [p_up_empty(i,:), p_dn_full(i,:), pUpMaFail(i), pDnMaFail(i), ps(i), pb(i), p001(i), p011(i), p101(i), Sp00A(i), Sp01A(i), pNm1_10(i), pN10(i), pN11(i), SpNA0(i), SpNA1(i), pSum1(i), pSum2(i), TP(i), nb(i), sol, IS1(i), IS2(i)]=corr_TML_MSP(mu_u(i),p_u(i),gamma_u(i),mu_d(i),p_d(i),gamma_d(i),Q_u(i),Q_d(i),C(i), MaxSparepartsStation);
            
            if (abs(TP(i)-local_TP_i)/local_TP_i<0.001 || abs(TP(i)-local_TP_i)/local_TP_i<0.01 && iterationFP>20 || abs(TP(i)-local_TP_i)/local_TP_i<0.1 && iterationFP>50)
                NoConvergence=0;
            end
            
            if (iterationFP>100 && NoConvergence ~=0)
                % Convergence by force
                disp("Convergence by force upstream");
                TP(i) = (TP(i) + local_TP_i)/2;
                NoConvergence = 0;
            end
            
            if (abs(TP(i)-local_TP_i)/local_TP_i < 0.1)
                current_omega = 0;
            end
        end
    end
    
    if (echoLevel >= 3)
        TP_str = strtrim(sprintf('%2.4f ', TP));
        mu_u_str = strtrim(sprintf('%2.4f ', mu_u));
        p_u_str = strtrim(sprintf('%2.4f ', p_u));
        gamma_u_str = strtrim(sprintf('%2.4f ', gamma_u));
        fprintf('Forward pass completed (iteration #%i)\n  TP: \t\t[%s]\n  mu_u: \t[%s]\n  p_u: \t\t[%s]\n  gamma_u: \t[%s]\n\n', iteration, TP_str, mu_u_str, p_u_str, gamma_u_str);
    end
    
    % Backward pass
    for i=machines-2:-1:1
        iterationBP = 0;
        NoConvergence = 1;
        current_omega = omega;
        
        [p_up_empty(i,:), p_dn_full(i,:), pUpMaFail(i), pDnMaFail(i), ps(i), pb(i), p001(i), p011(i), p101(i), Sp00A(i), Sp01A(i), pNm1_10(i), pN10(i), pN11(i), SpNA0(i), SpNA1(i), pSum1(i), pSum2(i), TP(i), nb(i), sol, IS1(i), IS2(i)]=corr_TML_MSP(mu_u(i),p_u(i),gamma_u(i),mu_d(i),p_d(i),gamma_d(i),Q_u(i),Q_d(i),C(i), MaxSparepartsStation);
        
        while NoConvergence
            % Fixed point iteration
            iterationBP = iterationBP + 1;
            local_p_d = p_d(i);
            local_gamma_d = gamma_d(i);
            local_mu_d = mu_d(i);
            local_TP_i = TP(i);
            
            B_ip1 = 0;
            B_u_ip1 = 0;
            B_d_i = 0;
            for alpha = 1 : Q(i+1)
                B_ip1 = B_ip1 + factorial(Q(i+1))/factorial(Q(i+1) - alpha) * (gamma(i+1) / p(i+1) )^alpha;
                B_u_ip1 = B_u_ip1 + factorial(Q(i+1))/factorial(Q(i+1) - alpha) * (gamma_u(i+1) / p_u(i+1) )^alpha;
                B_d_i = B_d_i + factorial(Q(i+1))/factorial(Q(i+1) - alpha) * (local_gamma_d / local_p_d )^alpha;
                for l = 1 : alpha - 1
                    B_ip1 = B_ip1 + factorial(Q(i+1) - l) / factorial(Q(i+1) - alpha) * (gamma(i+1) / p(i+1) )^(alpha-l) * ( p_up_empty(i,l) / pDnMaFail(i) + p_dn_full(i+1,l) / pUpMaFail(i+1) );
                    B_u_ip1 = B_u_ip1 + factorial(Q(i+1) - l) / factorial(Q(i+1) - alpha) * (gamma_u(i+1) / p_u(i+1) )^(alpha-l) * ( p_dn_full(i+1,l) / pUpMaFail(i+1) );
                    B_d_i = B_d_i + factorial(Q(i+1) - l) / factorial(Q(i+1) - alpha) * (local_gamma_d / local_p_d )^(alpha-l) * ( p_up_empty(i,l) / pDnMaFail(i) );
                end
            end
            gamma_d(i) = current_omega * local_gamma_d + (1 - current_omega) * (gamma(i+1) + (Q(i+2)/Q(i+1)*gamma_d(i+1)-gamma(i+1))* min(max(SpNA0(i+1), SpNA0(i+1)/ (local_p_d / local_gamma_d / Q(i+1) * ( TP(i+1)/local_mu_d - pSum2(i)))), 1) );
            p_d(i) = current_omega * local_p_d + (1 - current_omega) * (p(i+1) + mu_u(i+1) * min(max(pNm1_10(i+1), pNm1_10(i+1) / (TP(i+1) / local_mu_d - pSum2(i))), 1)+ p_d(i+1) * min(max(SpNA1(i+1), SpNA1(i+1) / (TP(i+1) / local_mu_d - pSum2(i))), 1));
            mu_d(i) = current_omega * local_mu_d + (1 - current_omega) * (((1+B_d_i)/B_d_i) /( 1 / TP(i+1) + 1 / mu(i+1) * (1+B_ip1)/B_ip1 - 1 / mu_u(i+1) * (1+B_u_ip1)/B_u_ip1));
            
            % Updating the values of the current line
            [p_up_empty(i,:), p_dn_full(i,:), pUpMaFail(i), pDnMaFail(i), ps(i), pb(i), p001(i), p011(i), p101(i), Sp00A(i), Sp01A(i), pNm1_10(i), pN10(i), pN11(i), SpNA0(i), SpNA1(i), pSum1(i), pSum2(i), TP(i), nb(i), sol, IS1(i), IS2(i)]=corr_TML_MSP(mu_u(i),p_u(i),gamma_u(i),mu_d(i),p_d(i),gamma_d(i),Q_u(i),Q_d(i),C(i), MaxSparepartsStation);
            
            if (abs(TP(i)-local_TP_i)/local_TP_i<0.001 || abs(TP(i)-local_TP_i)/local_TP_i<0.01 && iterationBP>20 || abs(TP(i)-local_TP_i)/local_TP_i<0.1 && iterationBP>50)
                NoConvergence = 0;
            end
            
            if (iterationBP>2000 && NoConvergence ~=0)
                % Convergence by force
                disp("Convergence by force downstream");
                TP(i) = (TP(i) + local_TP_i)/2;
                NoConvergence = 0;
            end
            
            if (abs(TP(i)-local_TP_i)/local_TP_i < 0.1)
                current_omega = 0;
            end
        end
    end
    
    if (echoLevel >= 3)
        TP_str = strtrim(sprintf('%2.4f ', TP));
        mu_d_str = strtrim(sprintf('%2.4f ', mu_d));
        p_d_str = strtrim(sprintf('%2.4f ', p_d));
        gamma_d_str = strtrim(sprintf('%2.4f ', gamma_d));
        fprintf('Backward pass completed (iteration #%i)\n  TP: \t\t[%s]\n  mu_d: \t[%s]\n  p_d: \t\t[%s]\n  gamma_d: \t[%s]\n\n', iteration, TP_str, mu_d_str, p_d_str, gamma_d_str);
    end
    
    TP_diff_old = TP_diff;
    TP_diff = abs(TP(1) - TP(end));
    
    data(iteration) = TP_diff;
    
    TP_diff_vec = TP_old - TP;
    if (max(abs(TP_diff_vec)) < 10^(-6))
        fprintf('Convergence to wrong value detected. Stopping decomposition.\n\n');
        terminatedNormally = false;
        epsilon = TP_diff;
    end
    
    if (false && TP_diff > TP_diff_old)
        fprintf('Convergence in wrong direction. Stopping decomposition.\n\n');
        terminatedNormally = false;
        epsilon = TP_diff_old;
        TP = TP_old;
        nb = nb_old;
        IS1 = IS1_old;
        IS2 = IS2_old;
        break;
    end
    
    TP_old = TP;
    nb_old = nb;
    IS1_old = IS1;
    IS2_old = IS2;
    
    %check for cycling
    TP_diff_cyc(iteration) = round(abs(TP(1) - TP(end)), 8);
    for step = 3:6
        lastelement = max(size(TP_diff_cyc));
        if (lastelement-2*step+1 > 0)
            if (all(TP_diff_cyc(lastelement-2*step+1:lastelement-step) == TP_diff_cyc(lastelement-step+1:lastelement)))
                cycling = true;
                TP_diff_cyc = NaN;
            end
        end
    end
    
    if (cycling)
        fprintf('Cycling detected (no further convergence). Adjusting epsilon from %0.6f to %0.6f\n\n', epsilon, epsilon * 10);
        cycling = false;
        epsilon = epsilon * 10;
        if (epsilon >= 0.1 * mu(1))
            terminatedNormally = false;
        end
    end
    
    if (TP_diff == TP_diff_old)
        fprintf('No change in TP (no further convergence). Adjusting epsilon from %0.6f to %0.6f\n\n', epsilon, epsilon * 10);
        epsilon = epsilon * 10;
        if (epsilon >= 0.1 * mu(1))
            terminatedNormally = false;
        end
    end
    
    if (mod(iteration, 100) == 0)
        fprintf('No convergence after iteration #%i. Adjusting epsilon from %0.6f to %0.6f\n\n', iteration, epsilon, epsilon * 10);
        epsilon = epsilon * 10;
        if (epsilon >= 0.1 * mu(1))
            terminatedNormally = false;
        end
    end
    
    if (iteration > 1000)
        fprintf('No convergence at all.\n\n');
        terminatedNormally = false;
    end
    
    if (echoLevel >= 2)
        fprintf("Decomposition iteration #%i completed. Current difference in TP is %2.6f (old was %2.6f).\n\n", iteration, TP_diff, TP_diff_old);
    end
    
    if (pauseAfterIteration), pause; end
    if (sleepAfterIteration), pause(0.25); end
end

%% prepare result
TP_final = TP(end);
if (echoLevel >= 1)
    fprintf("Decomposition converged after iteration #%i. The final TP is %2.4f.\n\n", iteration, TP_final);
end

performance.inventory_buffer = nb;
performance.inventory_buffer_total = sum(nb);
performance.inventory_spares = [IS1 IS2(end)];
performance.inventory_spares_total = sum([IS1 IS2(end)]);


%% function definition

    function [p_empty, p_full, pUpMaFail, pDnMaFail, ps, pb, p001, p011, p101, Sp00A, Sp01A, pNm1_10, pN10, pN11, SpNA0, SpNA1, pSum1, pSum2, PR, nb, sol, IS1, IS2]=corr_TML_MSP(mu1, p1, gam1, mu2, p2, gam2, Q1, Q2, Clocal, MaxSparepartsStation)
        if (mu1 < 0 || p1 < 0 || gam1 < 0 || mu2 < 0 ||  p2 < 0 || gam2 < 0)
            disp("Negative rates during TML!");
        end
        
        %************************************************************
        % Function to map system states to Scilab matrix rows/columns
        % StN means StateNumber
        function [sn]=StN(n, alpha1, alpha2) %, Q1, Q2, C)
            % Extended buffer space including the workspaces at the machines
            N = Clocal + 2;
            
            if n == 0
                if (alpha2 == 0)
                    disp([n,"n="]);
                    disp([alpha1,"alpha1="]);
                    disp([alpha2,"alpha2="]);
                    disp([Q1,"Q1="]);
                    disp([Q2,"Q2="]);
                    disp([Clocal,"C="]);
                    error("Lower boundary error!")
                end
                sn = Q2 * alpha1 + alpha2;
            elseif n < N
                sn = n * (Q1+1)*(Q2+1) + (Q2+1)*alpha1 + (alpha2+1) - (Q1+1);
            else
                if (alpha1 == 0)
                    disp([n,"n="]);
                    disp([alpha1,"alpha1="]);
                    disp([alpha2,"alpha2="]);
                    disp([Q1,"Q1="]);
                    disp([Q2,"Q2="]);
                    disp([Clocal,"C="]);
                    error("Upper boundary error!")
                end
                sn = N*(Q1+1)*(Q2+1) + (Q2+1)*(alpha1-1) + (alpha2+1) - (Q1+1);
            end
        end
        % ***************************************
        
        % Extended buffer space including the workspaces at the machines
        N = Clocal + 2;
        
        NumberOfStates = (N+1)*(Q1+1)*(Q2+1) - (Q1+1) - (Q2+1);
        Qmat = zeros(NumberOfStates, NumberOfStates);
        
        % Now we build up the transition rate matrix
        % For each state we consider its potential next states and
        % the respective transition rate, i.e, we consider * leaving * states
        
        for n=0:N
            for alpha1=0:Q1
                for alpha2=0:Q2
                    if (((n>0) || (alpha2>=1)) && ((n<N) || (alpha1>=1)))
                        
                        % Processing and failures
                        if ((alpha1 > 0) && (n < N))
                            % Finishing a part at the first machine
                            Qmat( StN( n, alpha1, alpha2), StN( n+1, alpha1, alpha2) ) = mu1;
                            % Failure and replacement at the first machine
                            Qmat( StN( n, alpha1, alpha2), StN( n, alpha1-1, alpha2) ) = p1;
                        end
                        
                        if ((alpha2 > 0) && (n > 0))
                            % Finishing a part at the second machine
                            Qmat( StN( n, alpha1, alpha2), StN( n-1, alpha1, alpha2) ) = mu2;
                            % Failure and replacement at the second machine
                            Qmat( StN( n, alpha1, alpha2), StN( n, alpha1, alpha2-1) ) = p2;
                        end
                        % Arrivals of new spare parts
                        if (alpha1 < Q1)
                            % Spare part arrival at the first machine
                            Qmat( StN( n, alpha1, alpha2), StN( n, alpha1+1, alpha2) ) = gam1*(Q1-alpha1);
                        end
                        if (alpha2 < Q2)
                            % Spare part arrival at the second machine
                            Qmat( StN( n, alpha1, alpha2), StN( n, alpha1, alpha2+1) ) = gam2*(Q2-alpha2);
                        end
                    end
                    
                end
            end
        end
        
        % now determine elements of the main diagonal
        for j = 1:NumberOfStates
            for k = 1:NumberOfStates
                if j~=k
                    Qmat(j,j) = Qmat(j,j) - Qmat(j,k);
                end
            end
        end
        
        Qmod = Qmat;
        
        for j=1:NumberOfStates
            Qmod(j,NumberOfStates) = 1;
        end
        
        % modified vector of zeros for the right hand side
        nmod = [zeros(1,NumberOfStates-1),1];
        
        prob = nmod / Qmod;
        
        psum=0;
        
        for s=1:NumberOfStates
            psum=psum + prob(s);
        end
        
        if (abs(psum - 1) > 0.00001)
            error("Probabilities do not sum to 1!")
        end
        
        e1org=0;
        e2org=0;
        psorg=0;
        pborg=0;
        
        psum=0;
        
        pbv = zeros(1,Q1);
        psv = zeros(1,Q2);
        
        for n=0:N
            for alpha1=0:Q1
                for alpha2=0:Q2
                    if (((n>0) || (alpha2>=1)) && ((n<N) || (alpha1>=1)))
                        psum=psum + prob(StN(n,alpha1,alpha2));
                        %mprintf('State (%i,%i,%i) p =%20.20f\n',n,alpha1,alpha2,prob(StN(n,alpha1,alpha2)));
                        if ((alpha1>=1) && (n<N))
                            e1org = e1org + prob(StN(n,alpha1,alpha2));
                        end
                        if ((alpha2>=1) && (n>0))
                            e2org = e2org + prob(StN(n,alpha1,alpha2));
                        end
                        if (n==N)
                            pborg = pborg + prob(StN(n,alpha1,alpha2));
                            pbv(alpha1) = pbv(alpha1) + prob(StN(n,alpha1,alpha2));
                        end
                        if (n==0)
                            psorg = psorg + prob(StN(n,alpha1,alpha2));
                            psv(alpha2) = psv(alpha2) + prob(StN(n,alpha1,alpha2));
                        end
                    end
                end
            end
        end
        
        if (abs(psum - 1) > 0.00001)
            error("Probabilities do not sum to 1!")
        end
        
        nbar=0;
        for n=1:N
            for alpha1=0:Q1
                for alpha2=0:Q2
                    if (n == N && alpha1 == 0), continue; end
                    nbar = nbar + n * prob(StN(n,alpha1,alpha2));
                end
            end
        end
        % First machine cannot fail if it is blocked!
        for alpha1 = 1 : Q1
            for alpha2 = 0 : Q2
                nbar = nbar + N * prob(StN( N ,alpha1,alpha2));
            end
        end
        
        % production rate
        PR1 = 0;
        % first machine must be up
        for alpha1 = 1 : Q1
            % if n=0  second machine cannot be down
            for alpha2 = 1 : Q2
                PR1 = PR1 + mu1 * prob(StN( 0, alpha1, alpha2));
            end
        end
        
        for n=1:N-1
            for alpha1 = 1 : Q1
                for alpha2 = 0 : Q2
                    PR1 = PR1 + mu1 * prob(StN(n,alpha1,alpha2));
                end
            end
        end
        
        % production rate
        PR2 = 0;
        % second machine must be up
        for n=1:N-1
            for alpha1 = 0 : Q1
                for alpha2 = 1 : Q2
                    PR2 = PR2 + mu2 * prob(StN(n,alpha1,alpha2));
                end
            end
        end
        
        % if n=N  second macine cannot be down
        for alpha1 = 1 : Q1
            for alpha2 = 1 : Q2
                PR2 = PR2 + mu2 * prob(StN( N , alpha1, alpha2));
            end
        end
        
        nb = nbar;
        PR = PR1;
        
        Sp00A = 0;
        Sp01A = 0;
        p001= prob(StN(0,0,1));
        p011= prob(StN(0,1,1));
        p101 = prob(StN(1,0,1));
        
        pNm1_10 = prob(StN(N-1,1,0));
        pN10= prob(StN(N,1,0));
        pN11= prob(StN(N,1,1));
        
        SpNA0 = 0;
        SpNA1 = 0;
        
        for alpha1=1:Q1
            SpNA0 = SpNA0 + prob(StN(N,alpha1,0));
            SpNA1 = SpNA1 + prob(StN(N,alpha1,1));
        end
        
        for alpha2=1:Q2
            Sp00A = Sp00A + prob(StN(0,0,alpha2));
            Sp01A = Sp01A + prob(StN(0,1,alpha2));
        end
        
        pSum1=0;
        for n = 0:N-1
            for alpha1 = 2:Q1
                for alpha2 = 0:Q2
                    if (n>0 || alpha2>0)
                        pSum1 = pSum1 + prob(StN(n,alpha1,alpha2));
                    end
                end
            end
        end
        
        pSum2=0;
        for n = 1:N
            for alpha1 = 0:Q1
                for alpha2 = 2:Q2
                    if (n<N || alpha1>0)
                        pSum2 = pSum2 + prob(StN(n,alpha1,alpha2));
                    end
                end
            end
        end
        
        p_full = zeros(1,MaxSparepartsStation);
        p_empty= zeros(1,MaxSparepartsStation);
        
        for alpha1 = 1:Q1
            p_full(alpha1) = pbv(alpha1);
        end
        
        for alpha2 = 1:Q2
            p_empty(alpha2) = psv(alpha2);
        end
        
        pb = sum(pbv);
        ps = sum(psv);
        
        pUpMaFail=0;
        for n = 0:N-1
            for alpha2 = 0:Q2
                if (n>0 || alpha2>0)
                    pUpMaFail = pUpMaFail + prob(StN(n,0,alpha2));
                end
            end
        end
        
        pDnMaFail=0;
        for n = 1:N
            for alpha1 = 0:Q1
                if (n<N || alpha1>0)
                    pDnMaFail = pDnMaFail + prob(StN(n,alpha1,0));
                end
            end
        end
        
        sol = zeros(N + 1, Q1 + 1, Q2 + 1);
        for n = 0:N
            for m1 = 0:Q1
                for m2 = 0:Q2
                    if ~((n == 0 && m2 == 0) || (n == N && m1 == 0))
                        f = StN(n, m1, m2);
                        sol(n+1, m1+1, m2+1) = prob(f);
                    end
                end
            end
        end
        
        % average spare part inventory
        IS1 = 0;
        for m1 = 1:Q1
            tmp = 0;
            for n = 0:N
                for m2 = 0:Q2
                    if (n == 0 && m2 == 0), continue; end
                    tmp = tmp + prob(StN(n, m1, m2));
                end
            end
            IS1 = IS1 + (m1 - 1) * tmp;
        end
        
        IS2 = 0;
        for m2 = 1:Q2
            tmp = 0;
            for n = 0:N
                for m1 = 0:Q1
                    if (n == N && m1 == 0), continue; end
                    tmp = tmp + prob(StN(n, m1, m2));
                end
            end
            IS2 = IS2 + (m2 - 1) * tmp;
        end
        
    end


%checkSize  checks if a numeric array has the desired size
%   checkSize(A, desired_dims_sizes) returns:
%       if A is an array and has the desired dimension sizes, nothing
%       else, an error
%
%   Programmed with Matlab's 'size' function in mind
%       i.e. checkSize( A, size(A) ) should pass the test
%   Inputs:
%       A - input to test whether it is a numeric array of desired size
%       desired_dims_sizes - 1 by d array with the sizes for each dimension
%           if this param is not given, A just needs to be a numeric array
%           if this param is empty, A must be empty
%           otherwise A's size must be equal to this param
%           An Inf value at any place in this param allows
%               for any size in that dimension
%
%   Example usages:
%       checkSize( A, [Inf 2] ) - i.e., no restrictions on dimension 1 and size 2 for dimension 2
%       checkSize( A, [Inf Inf] ) - for an MxN array
%       checkSize( A, [3 4 5]) - for a 3x4x5 array
%       checkSize( A, []) - for an empty array
%       Note, for instance, that checkSize( A, [Inf Inf] ) will fail if A has more than two dimensions

% By Paulo Abelha (p.abelha@abdn.ac.uk) 2016
% Thanks to Jan Simon for pointing out simplifications and bugs
% Changed by Florian Sachs to cover boolean return and our situational needs (originally: checkNumericArraySize)

    function [result] = checkSize( A, dim_sizes )
        result = false;
        if ~exist('A','var')
            return;
        end
        if ~isnumeric(A)
            return;
        end
        if ~exist('dim_sizes','var')
            return;
        end
        if ~isnumeric(dim_sizes)
            return;
        end
        if isempty(dim_sizes)
            return;
        end
        if isempty(A)
            return;
        end
        if size(dim_sizes,1) ~= 1 || size(dim_sizes,2) < 2
            return;
        end
        if ndims(A) ~= size(dim_sizes,2)
            return;
        end
        if ~all((size(A) == dim_sizes | isinf(dim_sizes)))
            return;
        end
        result = true;
    end

end