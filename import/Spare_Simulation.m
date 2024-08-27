function [TP_final, halfwidth, performance] = Spare_Simulation(machines, C, Q, mu, p, gamma)
%

%% check input arguments
if (nargin ~= 6)
    fprintf('Wrong number of parameters... expecting: machines, C, Q, p, mu, gamma\n\n');
    machines = 5;
    C = 5 * ones(1,machines-1);
    Q = 2 * ones(1,machines);
    mu = 1 * ones(1,machines);
    p = 0.005 * ones(1,machines);
    gamma = 0.05 * ones(1,machines);
end

if (~checkSize(C, [1 machines-1]) || ~checkSize(Q, [1 machines]) || ~checkSize(mu, [1 machines]) || ~checkSize(p, [1 machines]) || ~checkSize(gamma, [1 machines]))
    C = C'; Q = Q'; mu = mu'; p = p'; gamma = gamma';
    if (~checkSize(C, [1 machines-1]) || ~checkSize(Q, [1 machines]) || ~checkSize(mu, [1 machines]) || ~checkSize(p, [1 machines]) || ~checkSize(gamma, [1 machines]))
        error("Dimensions mismatch!");
    end
end


%% settings
echo = false;
warmup = 1000;
runtime = 100000;
minruns = 10;


%% initialisation
prod = 1:machines;
fail = machines+1:2*machines;
delivery = 2*machines+1:3*machines;
halfwidth = Inf;

pDown = zeros(1000,machines);
pStarv = zeros(1000,machines);
pBlock = zeros(1000,machines);
pBlockStarv = zeros(1000,machines);
rho = zeros(1000,1);
rho_new = zeros(1000,1);

k = 0;


%% simulation
if (echo)
    fprintf('Starting simulation.\n\n');
end

stepNumber = 1;
while (k < minruns || isnan(halfwidth) || halfwidth > 0.01)
    %parfor i = 1:stepNumber
    for i = 1:stepNumber
        n = zeros(1, machines-1);
        m = Q;
        warmupReady = 0;
        time = 0;
        event = 0;
        eventP = 0;
        eventF = 0;
        eventD = 0;
        output = 0;
        timesafe = 0;
        timeDown = zeros(1, machines);
        timeDownBuffer = zeros(1, machines);
        timeBlock = zeros(1, machines);
        timeBlockBuffer = zeros(1, machines);
        timeStarv = zeros(1, machines);
        timeStarvBuffer = zeros(1, machines);
        timeBlockStarv = zeros(1, machines);
        timeBlockStarvBuffer = zeros(1, machines);
        while ((time - warmup) < runtime || warmupReady == 0)
            if (time >= warmup && warmupReady == 0)
                warmupReady = 1;
                event = 0;
                eventP = 0;
                eventF = 0;
                eventD = 0;
                output = 0;
                timeDown = zeros(1, machines);
                timeBlock = zeros(1, machines);
                timeStarv = zeros(1, machines);
                timeBlockStarv = zeros(1, machines);
                timesafe = time;
            end
            
            starvAll = [0 n == 0];
            blockAll = [n == (C+2) 0];
            muC = (m > 0 & ~starvAll & ~blockAll) .* mu;
            pC = (m > 0 & ~starvAll & ~blockAll) .* p;
            gammaC = (Q - m) .* gamma;
            rate = sum(muC) + sum(gammaC) + sum(pC);
            prob = cumsum([muC pC gammaC]);
            random = rand * rate;
            timePos = find((random < prob)~=0, 1, 'first');
            timeVal = exprnd(1 / rate);
            timePosC = mod(timePos-1, machines) + 1;
            
            time = time + timeVal;
            event = event + 1;
            
            if (ismember(timePos, prod))
                if (timePosC > 1)
                    n(timePosC-1) = n(timePosC-1) - 1;
                    if (n(timePosC - 1) == 0)
                        timeStarvBuffer(timePosC) = time;
                        if (timePosC < machines)
                            if (n(timePosC) == C(timePosC) + 2)
                                timeBlockStarvBuffer(timePosC) = time;
                            end
                        end
                    end
                    if (n(timePosC-1) == C(timePosC-1) + 1)
                        timeBlock(timePosC-1) = timeBlock(timePosC-1) + time - timeBlockBuffer(timePosC-1);
                        timeBlockBuffer(timePosC-1) = 0;
                        if (timeBlockStarvBuffer(timePosC-1) > 0)
                            timeBlockStarv(timePosC-1) = timeBlockStarv(timePosC-1) + time - timeBlockStarvBuffer(timePosC-1);
                            timeBlockStarvBuffer(timePosC-1) = 0;
                        end
                    end
                end
                
                if (timePosC < machines)
                    n(timePosC) = n(timePosC) + 1;
                    if (n(timePosC) == 1)
                        timeStarv(timePosC+1) = timeStarv(timePosC+1) + time - timeStarvBuffer(timePosC+1);
                        timeStarvBuffer(timePosC+1) = 0;
                        if (timeBlockStarvBuffer(timePosC+1) > 0)
                            timeBlockStarv(timePosC+1) = timeBlockStarv(timePosC+1) + time - timeBlockStarvBuffer(timePosC+1);
                            timeBlockStarvBuffer(timePosC+1) = 0;
                        end
                    end
                    if (n(timePosC) == C(timePosC) + 2)
                        timeBlockBuffer(timePosC) = time;
                        if (timePosC > 1)
                            if (n(timePosC - 1) == 0)
                                timeBlockStarvBuffer(timePosC) = time;
                            end
                        end
                    end
                end
                
                if (timePosC == machines); output = output + 1; end
                eventP = eventP + 1;
            elseif (ismember(timePos, fail))
                m(timePosC) = m(timePosC) - 1;
                eventF = eventF + 1;
                if (m(timePosC) == 0)
                    timeDownBuffer(timePosC) = time;
                end
            elseif (ismember(timePos, delivery))
                m(timePosC) = m(timePosC) + 1;
                eventD = eventD + 1;
                if (m(timePosC) == 1)
                    timeDown(timePosC) = timeDown(timePosC) + time - timeDownBuffer(timePosC);
                    timeDownBuffer(timePosC) = 0;
                end
            end
        end
        
        %i = i + 1
        pDown(k+i,:) = timeDown / (time - timesafe);
        pBlock(k+i,:) = timeBlock / (time - timesafe);
        pStarv(k+i,:) = timeStarv / (time - timesafe);
        pBlockStarv(k+i,:) = timeBlockStarv / (time - timesafe);
        tmp = timeDown / (time - timesafe) + timeStarv / (time - timesafe);
        rho(k+i) = (1-tmp(end)) * mu(end);
        rho_new(k+i) = output / (time - timesafe);
    end
    k = k + stepNumber;
    %rho(1:k) = rho_new(1:k);
    rho_avg = sum(rho(1:k))/k;
    samplevar = sum((rho(1:k)-rho_avg).^2)/(k-1);
    zalpha = tinv(1-0.05/2, k-1);
    halfwidth = zalpha*sqrt(samplevar/k);
    if (echo)
        fprintf('%2i: \t rho: %2.4f \t halfwidth: %2.4f \t [%2.4f, %2.4f]\n', k, rho_avg, halfwidth, rho_avg - halfwidth, rho_avg + halfwidth);
        %fprintf('Events: %i \t productions: %i \t failures: %i \t deliveries: %i\n', event, eventP, eventF, eventD)
    end
end


%% prepare result
format longG
if (echo)
    disp(' ')
    %fprintf('Events: %i \t productions: %i \t failures: %i \t deliveries: %i\n', event, eventP, eventF, eventD)
end

pBlock_avg = sum(pBlock(1:k,:))/k;
pStarv_avg = sum(pStarv(1:k,:))/k;
pDown_avg = sum(pDown(1:k,:))/k;
pBlockStarv_avg = sum(pBlockStarv(1:k,:))/k;
rho_avg = sum(rho(1:k))/k;
TP_final = rho_avg(end);

performance.inventory_buffer = NaN(1, machines - 1);
performance.inventory_spares = NaN(1, machines);
performance.p_block = pBlock_avg;
performance.p_starv = pStarv_avg;
performance.p_down = pDown_avg;
performance.p_block_and_starv = pBlockStarv_avg;


%% function definition

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