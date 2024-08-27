function [TP, performance] = Spare_2M(C, Q1, Q2, mu1, mu2, p1, p2, gamma1, gamma2)
%

%% check input arguments
if (nargin ~= 9)
    error('Wrong number of parameters... expecting: C, Q1, Q2, mu1, mu2, p1, p2, gamma1, gamma2');
end


%% initialisation
N = C + 2;
numberOfStatesTotal = (N + 1) * (Q1 + 1) * (Q2 + 1);
numberOfStates = (N + 1) * (Q1 + 1) * (Q2 + 1) - (Q1 + 1) - (Q2 + 1);
Q = zeros(numberOfStatesTotal, numberOfStatesTotal);


%% get generator matrix
for n = 0:N
    for m1 = 0:Q1
        for m2 = 0:Q2
            %calculate state numbers
            f = calculateStateNumber(n, m1, m2, Q1, Q2, C);
            %fill line of Q
            if ((n == N && m1 == 0) || (n == 0 && m2 == 0))
                %prepare Q for impossible states
                Q(f,f) = NaN;
            else
                fnp = calculateStateNumber(n+1, m1, m2, Q1, Q2, C);
                fnm = calculateStateNumber(n-1, m1, m2, Q1, Q2, C);
                f1p = calculateStateNumber(n, m1+1, m2, Q1, Q2, C);
                f1m = calculateStateNumber(n, m1-1, m2, Q1, Q2, C);
                f2p = calculateStateNumber(n, m1, m2+1, Q1, Q2, C);
                f2m = calculateStateNumber(n, m1, m2-1, Q1, Q2, C);
                if (fnm ~= 0 && m1 ~= 0),    Q(f,fnm) = mu1;                    end
                if (fnp ~= 0 && m2 ~= 0),    Q(f,fnp) = mu2;                    end
                if (f1p ~= 0 && n ~= N),     Q(f,f1p) = p1;                     end
                if (f2p ~= 0 && n ~= 0),     Q(f,f2p) = p2;                     end
                if (f1m ~= 0),               Q(f,f1m) = (Q1 - m1 + 1) * gamma1; end
                if (f2m ~= 0),               Q(f,f2m) = (Q2 - m2 + 1) * gamma2; end
            end
        end
    end
end


%% solve

%remove rows and columns of impossible states
possible = ~isnan(diag(Q));
Q = Q(possible,possible);
impossibleStateNumbers = 1:numberOfStatesTotal;
impossibleStateNumbers = impossibleStateNumbers(~possible);

%add main diagonal elements
Q = Q + diag(-sum(Q));

%prepare right-hand side
tnew = zeros(numberOfStates, 1);

%add normalization condition (new: overwrite one balance equation; system is overdetermined)
Qnew = Q;
Qnew(numberOfStates, :) = ones(1,numberOfStates)';
tnew(numberOfStates, 1) = 1;

%solve
r = Qnew \ tnew;


%% transform solution from vector to matrix for easier access
if (sum(abs(r)) - 1 > 10^(-10))
    fprintf("Sum: %2.4f\n  C: %i\n  Q1: %i\n  Q2: %i\n  mu1: %2.4f\n  mu2: %2.4f\n  p1: %2.4f\n  p2: %2.4f\n  gamma1: %2.4f\n  gamma2: %2.4f\n\n", sum(abs(r)), C, Q1, Q2, mu1, mu2, p1, p2, gamma1, gamma2);
    error("Probabilities do not sum up to 1!");
end

sol = zeros(N + 1, Q1 + 1, Q2 + 1);
minus = 0;
for n = 0:N
    for m1 = 0:Q1
        for m2 = 0:Q2
            f = calculateStateNumber(n, m1, m2, Q1, Q2, C);
            if (ismember(f,impossibleStateNumbers))
                minus = minus + 1;
                sol(n+1, m1+1, m2+1) = 0;
            else
                sol(n+1, m1+1, m2+1) = r(f - minus);
            end
        end
    end
end

TP = (sum(sol(:,:,2:end), 'all') - sum(sol(1,:,2:end), 'all')) * mu2;

performance.processing = [sum(sol(1:(end-1),2:end,:), 'all') sum(sol(2:end,:,2:end), 'all')];
performance.availability = [sum(sol(:,2:end,:), 'all') sum(sol(:,:,2:end), 'all')];
performance.p_block = sum(sol(end,:,:), 'all');
performance.p_starv = sum(sol(1,:,:), 'all');


%% function definition
    function [ f ] = calculateStateNumber( n, m1, m2, Q1, Q2, C)
        if (nargin == 5)
            C = inf;
        end
        if (n < 0 || m1 < 0 || m2 < 0 || n > (C + 2) || m1 > Q1 || m2 > Q2)
            f = 0;
        else
            f = 1 + m2 + (Q2 + 1) * m1 + (Q2 + 1) * (Q1 + 1) * n;
        end
    end

end