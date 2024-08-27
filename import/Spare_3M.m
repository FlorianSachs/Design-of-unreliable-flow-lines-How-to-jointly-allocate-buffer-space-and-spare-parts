function [ TP, performance ] = Spare_3M( C1, C2, Q1, Q2, Q3, mu1, mu2, mu3, p1, p2, p3, gamma1, gamma2, gamma3 )
%

%% check input arguments
if (nargin ~= 14)
    error('Wrong number of parameters... expecting: C1, C2, Q1, Q2, Q3, mu1, mu2, mu3, p1, p2, p3, gamma1, gamma2, gamma3');
end


%% initialisation
echo = 0;
S1 = Q1 - 1;
S2 = Q2 - 1;
S3 = Q3 - 1;
numberOfStates = (C1 + 3) * (C2 + 3) * (S1 + 2) * (S2 + 2) * (S3 + 2);
Q = sparse(numberOfStates, numberOfStates);
i = 1;
%varbytes = @(x)getfield(whos('x'),'bytes');
%fprintf("The matrix stores %2.2f MB (%2.2f GB).\n\n", varbytes(Q)/1024^2, varbytes(Q)/1024^3);


%% calculation
tic();
for n1 = 0:1:(C1+2)
    for n2 = 0:1:(C2+2)
        for m1 = 0:1:(S1+1)
            for m2 = 0:1:(S2+1)
                for m3 = 0:1:(S3+1)
                    %calculate state numbers
                    f = calculateStateNumber(n1, n2, m1, m2, m3, S1, S2, S3, C1, C2);
                    fn_1p_2m = calculateStateNumber(n1+1, n2-1, m1, m2, m3, S1, S2, S3, C1, C2);
                    fn1m = calculateStateNumber(n1-1, n2, m1, m2, m3, S1, S2, S3, C1, C2);
                    fn2p = calculateStateNumber(n1, n2+1, m1, m2, m3, S1, S2, S3, C1, C2);
                    f1p = calculateStateNumber(n1, n2, m1+1, m2, m3, S1, S2, S3, C1, C2);
                    f1m = calculateStateNumber(n1, n2, m1-1, m2, m3, S1, S2, S3, C1, C2);
                    f2p = calculateStateNumber(n1, n2, m1, m2+1, m3, S1, S2, S3, C1, C2);
                    f2m = calculateStateNumber(n1, n2, m1, m2-1, m3, S1, S2, S3, C1, C2);
                    f3p = calculateStateNumber(n1, n2, m1, m2, m3+1, S1, S2, S3, C1, C2);
                    f3m = calculateStateNumber(n1, n2, m1, m2, m3-1, S1, S2, S3, C1, C2);
                    %fill line of Q
                    if ((n1 == C1 + 2 && m1 == S1 + 1) || (n1 == 0 && m2 == S2 + 1) || (n2 == C2 + 2 && m2 == S1 + 1) || (n2 == 0 && m3 == S3 + 1))
                        %prepare Q for impossible states
                        Q(f,f) = NaN;
                    else
                        if (fn1m ~= 0 && m1 ~= S1 + 1); Q(f,fn1m) = mu1; end
                        if (fn2p ~= 0 && m3 ~= S3 + 1); Q(f,fn2p) = mu3; end
                        if (fn_1p_2m ~= 0 && m2 ~= S2 + 1); Q(f,fn_1p_2m) = mu2; end
                        if (f1p ~= 0); Q(f,f1p) = (m1 + 1) * gamma1; end
                        if (f1m ~= 0 && n1 ~= C1+2); Q(f,f1m) = p1; end
                        if (f2p ~= 0); Q(f,f2p) = (m2 + 1) * gamma2; end
                        if (f2m ~= 0 && n1 ~= 0 && n2 ~= C2+2); Q(f,f2m) = p2; end
                        if (f3p ~= 0); Q(f,f3p) = (m3 + 1) * gamma3; end
                        if (f3m ~= 0 && n2 ~= 0); Q(f,f3m) = p3; end
                    end
                    i = i + 1;
                end
            end
        end
    end
end

%prepare Q for impossible states
for f = 1:numberOfStates
    if (isnan(Q(f,f)))
        Q(f,f) = 1;
    else
        Q(f,f) = -sum(Q(:,f));
    end
end

%add normalization condition
Q(numberOfStates + 1, :) = ones(1,numberOfStates)';
%disp(Q);


%% solve
t = zeros(numberOfStates, 1);
t(numberOfStates + 1, 1) = 1;
%r = linsolve(Q,t);
S = sparse(Q);
%S = gpuArray(Q);
r = S\t;
%r = gather(r);
r = abs(r);
sum(r);
time = toc();
%fprintf("The complete matrix stores %2.2f MB (%2.2f GB).\n\n", varbytes(Q)/1024^2, varbytes(Q)/1024^3);


%% transform solution from vector to matrix for easier access
sol = zeros(C1+3, C2+3, S1 + 2, S2 + 2);
for n1 = 0:1:(C1+2)
    for n2 = 0:1:(C2+2)
        for m1 = (S1+1):-1:0%for m1 = 0:1:(S1+1)
            for m2 = (S2+1):-1:0%for m2 = 0:1:(S2+1)
                for m3 = (S3+1):-1:0%for m3 = 0:1:(S3+1)
                    sol(n1+1, n2+1, m1+1, m2+1, m3+1) = r(calculateStateNumber(n1, n2, m1, m2, m3, S1, S2, S3, C1, C2));
                end
            end
        end
    end
end


%% calculate output

% probabilities and throughput
p_down = [sum(sol(:,:,S1+2,:,:),'all') sum(sol(:,:,:,S2+2,:),'all') sum(sol(:,:,:,:,S3+2),'all')];
p_block = [sum(sol(C1+3,:,:,:,:),'all') sum(sol(:,C2+3,:,:,:),'all') 0];
p_starv = [0 sum(sol(1,:,:,:,:),'all') sum(sol(:,1,:,:,:),'all')];
p_block_and_starv = sum(sol(1,C2+3,:,:,:),'all');
TP = mu3 * (1-p_down(3)-p_starv(3));

% average extended buffer levels
bufferlevel1 = 0;
for n1 = 1:C1+2
    tmp = 0;
    for n2 = 0:C2+2
        for m1 = 0:S1+1
            for m2 = 0:S2+1
                for m3 = 0:S3+1
                    tmp = tmp + sol(n1+1,n2+1,m1+1,m2+1,m3+1);
                end
            end
        end
    end
    bufferlevel1 = bufferlevel1 + n1 * tmp;
end

bufferlevel2 = 0;
for n2 = 1:C2+2
    tmp = 0;
    for n1 = 0:C1+2
        for m1 = 0:S1+1
            for m2 = 0:S2+1
                for m3 = 0:S3+1
                    tmp = tmp + sol(n1+1,n2+1,m1+1,m2+1,m3+1);
                end
            end
        end
    end
    bufferlevel2 = bufferlevel2 + n2 * tmp;
end

% average spare part inventories
IS1 = 0;
for m1 = 0:S1
    tmp = 0;
    for n1 = 0:C1+2
        for n2 = 0:C2+2
            for m2 = 0:S2+1
                for m3 = 0:S3+1
                    tmp = tmp + sol(n1+1,n2+1,m1+1,m2+1,m3+1);
                end
            end
        end
    end
    IS1 = IS1 + (S1-m1) * tmp;
end

IS2 = 0;
for m2 = 0:S2
    tmp = 0;
    for n1 = 0:C1+2
        for n2 = 0:C2+2
            for m1 = 0:S1+1
                for m3 = 0:S3+1
                    tmp = tmp + sol(n1+1,n2+1,m1+1,m2+1,m3+1);
                end
            end
        end
    end
    IS2 = IS2 + (S2-m2) * tmp;
end

IS3 = 0;
for m3 = 0:S3
    tmp = 0;
    for n1 = 0:C1+2
        for n2 = 0:C2+2
            for m1 = 0:S1+1
                for m2 = 0:S2+1
                    tmp = tmp + sol(n1+1,n2+1,m1+1,m2+1,m3+1);
                end
            end
        end
    end
    IS3 = IS3 + (S3-m3) * tmp;
end

% store information in object
performance.inventory_buffer = [bufferlevel1 bufferlevel2];
performance.inventory_buffer_total = bufferlevel1 + bufferlevel2;
performance.inventory_spares = [IS1 IS2 IS3];
performance.inventory_spares_total = IS1 + IS2 + IS3;
performance.p_block = p_block;
performance.p_starv = p_starv;
performance.p_down = p_down;
performance.p_block_and_starv = p_block_and_starv;
performance.time = time;


%% function for state number
    function [ f ] = calculateStateNumber( n1, n2, s1, s2, s3, S1, S2, S3, N1, N2)
        if (nargin == 9)
            N2 = inf;
        end
        if (n1 < 0 || n2 < 0 || s1 < 0 || s2 < 0 || s3 < 0 || n1 > (N1 + 2) || n2 > (N2 + 2) || s1 > (S1 + 1) || s2 > (S2 + 1) || s3 > (S3 + 1))
            f = 0;
        else
            f = 1 + s3 + (S3 + 2) * s2 + (S3 + 2) * (S2 + 2) * s1 + (S3 + 2) * (S2 + 2) * (S1 + 2) * n1 + (S3 + 2) * (S2 + 2) * (S1 + 2) * (N1 + 3) * n2;
        end
    end

end