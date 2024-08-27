function [isoAV, isoTP, start_spares] = get_availability_characteristics(machines, Q, mu, p, gamma, rhoNeeded)
%get_availability_characteristics Calculates some characteristics of a flow line
%   Calculates isolated availabilities
    start_spares = 1 * ones(1, machines);
    isoAV = NaN(1, machines);
    isoTP = NaN(1, machines);
    for ii = 1:machines
        for Q_local = 1:Q(ii)
            numerator = 0; denominator = 0;
            for jj = 0:Q_local
                if (jj < Q_local)
                    numerator = numerator + factorial(Q_local-1) / factorial(jj) * gamma(ii)^(Q_local - 1 - jj) * p(ii)^jj;
                end
                denominator = denominator + factorial(Q_local) / factorial(jj) * gamma(ii)^(Q_local - jj) * p(ii)^jj;
            end
            isoAV(ii) = (Q_local * gamma(ii) * numerator) / denominator;
            isoTP(ii) = isoAV(ii) * mu(ii);
            if (isoTP(ii) > rhoNeeded)
                start_spares(ii) = Q_local;
                break;
            end
        end
    end
end

