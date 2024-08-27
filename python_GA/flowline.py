import numpy as np
import types
from math import factorial

class flowline:
    def __init__(self):
        self.machines = np.nan
        self.C = np.nan
        self.C_min = np.nan
        self.C_max = np.nan
        self.Q = np.nan
        self.Q_min = np.nan
        self.Q_max = np.nan
        self.mu = np.nan
        self.p = np.nan
        self.gamma = np.nan
        self.costs_buffers = np.nan
        self.costs_spares = np.nan
        self.sol_2M = np.nan
        self.characteristics = types.SimpleNamespace()
        self.characteristics_2M = types.SimpleNamespace()
        self.characteristics_dec = types.SimpleNamespace()

    def fill(self, machines, C_scalar, Q_scalar, mu_scalar, p_scalar, gamma_scalar, costs_buffers_scalar, costs_spares_scalar):
        if machines <= 1 or not np.int32(machines):
            print('Number of machines is not correct.')
            raise

        self.machines = machines
        self.C = C_scalar * np.ones(machines - 1, dtype=int)
        self.Q = Q_scalar * np.ones(machines, dtype=int)

        self.mu = mu_scalar * np.ones(machines, dtype=float)
        self.p = p_scalar * np.ones(machines, dtype=float)
        self.gamma = gamma_scalar * np.ones(machines, dtype=float)

        self.costs_buffers = costs_buffers_scalar * np.ones(machines - 1, dtype=float)
        self.costs_spares = costs_spares_scalar * np.ones(machines, dtype=float)

        self.check()
    
    def set_limits(self, C_min, C_max, Q_min, Q_max):
        self.C_min = C_min
        self.C_max = C_max
        self.Q_min = Q_min
        self.Q_max = Q_max
    
    def fill_tml(self, C, Q1, Q2, mu1, mu2, p1, p2, gamma1, gamma2, costs_buffer=1, costs_spares1=1, costs_spares2=1):
        self.machines = 2
        self.C = np.array([C])
        self.Q = np.array([Q1, Q2])
        self.mu = np.array([mu1, mu2])
        self.p = np.array([p1, p2])
        self.gamma = np.array([gamma1, gamma2])
        self.costs_buffers = np.array([costs_buffer])
        self.costs_spares = np.array([costs_spares1, costs_spares2])

    def check(self):
        if self.machines <= 1 or not np.int32(self.machines):
            print('Number of machines is not correct.')
            raise
        if np.any(self.C < 1) or np.any(self.C.round() != self.C) or np.ndim(self.C) != 1 or np.size(self.C) != self.machines - 1:
            print('Buffer capacities are not correct.')
            raise
        if np.any(self.Q < 1) or np.any(self.Q.round() != self.Q) or np.ndim(self.Q) != 1 or np.size(self.Q) != self.machines:
            print('Component numbers are not correct.')
            raise
        
        if np.ndim(self.mu) != 1 or np.size(self.mu) != self.machines: #np.any(self.mu < 0) or 
            print('Processing rates are not correct.')
            raise
        if np.ndim(self.p) != 1 or np.size(self.p) != self.machines: #np.any(self.p < 0) or 
            print('Failure rates are not correct.')
            raise
        if np.ndim(self.gamma) != 1 or np.size(self.gamma) != self.machines: #np.any(self.gamma < 0) or 
            print('Replenishment rates are not correct.')
            raise
        
        if np.any(self.costs_buffers < 0) or np.ndim(self.costs_buffers) != 1 or np.size(self.costs_buffers) != self.machines - 1:
            print('Buffer costs are not correct.')
            raise
        if np.any(self.costs_spares < 0) or np.ndim(self.costs_spares) != 1 or np.size(self.costs_spares) != self.machines:
            print('Component costs are not correct.')
            raise
        
        return True

    def solve(self):
        self.check()
        
        if self.machines == 2:
            self.solve_2M()
            self.characteristics = self.characteristics_2M
        elif self.machines >= 3:
            self.solve_dec()
            self.characteristics = self.characteristics_dec
        
    def solve_2M(self):
        self.check()

        if self.machines != 2:
            print('Number of machines must be 2.')
            raise

        sol = self.spare_2M()
        self.sol_2M = sol

        characteristics = self.get_characteristics_2M(sol)
        self.characteristics_2M = characteristics

    def solve_dec(self):
        self.check()

        if self.machines < 3:
            print('Number of machines must be larger than or equal to 3.')
            raise
        
        characteristics = self.spare_decomposition()
        self.characteristics_dec = characteristics

    def spare_decomposition(self):
        self.check()
        justresults = True

        machines = self.machines
        C = self.C
        Q = self.Q
        mu = self.mu
        p = self.p
        gamma = self.gamma

        ## settings
        echo_level = 1
        pause_after_iteration = False
        sleep_after_iteration = False
        epsilon = 10 ** (-4)
        omega = 0.5

        ## initialisation
        terminated_normally = True
        if justresults:
            echo_level = 0
            pause_after_iteration = False
            sleep_after_iteration = False

        TP = np.array([np.nan] * (machines - 1))
        TP_old = np.copy(TP)

        Q_u = self.Q[:-1].tolist()
        mu_u = self.mu[:-1].tolist()
        p_u = self.p[:-1].tolist()
        gamma_u = self.gamma[:-1].tolist()

        Q_d = self.Q[1:].tolist()
        mu_d = self.mu[1:].tolist()
        p_d = self.p[1:].tolist()
        gamma_d = self.gamma[1:].tolist()

        values = [np.nan] * (machines - 1)
        tml = flowline()

        for i in range(0, machines - 1):
            tml.fill_tml(C[i], Q_u[i], Q_d[i], mu_u[i], mu_d[i], p_u[i], p_d[i], gamma_u[i], gamma_d[i])
            tml.solve_2M()
            values[i] = tml.characteristics_2M
            TP[i] = tml.characteristics_2M.TP

        iteration = 0
        TP_diff = np.infty
        TP_diff_cyc = []
        cycling = False
        data = []

        ## decomposition
        if echo_level >= 1:
            TP_str = ', '.join(map('{:.4f}'.format, TP))
            print(f'Starting decomposition (iteration #{iteration}).\n  TP: \t\t[{TP_str}]\n\n')

        while TP_diff > epsilon and terminated_normally:
            iteration += 1

            # Forward pass
            for i in range(1, machines - 1):
                iteration_FP = 0
                no_convergence = 1
                current_omega = omega

                tml.fill_tml(C[i], Q_u[i], Q_d[i], mu_u[i], mu_d[i], p_u[i], p_d[i], gamma_u[i], gamma_d[i])
                tml.solve_2M()
                values[i] = tml.characteristics_2M
                TP[i] = tml.characteristics_2M.TP

                while no_convergence:
                    # Fixed point iteration
                    iteration_FP += 1
                    local_p_u = p_u[i]
                    local_gamma_u = gamma_u[i]
                    local_mu_u = mu_u[i]
                    local_TP_i = TP[i]

                    B_i = 0
                    B_u_i = 0
                    B_d_im1 = 0
                    for alpha in range(0, Q[i]):
                        B_i = B_i + factorial(Q[i]) / factorial(Q[i] - (alpha+1)) * (gamma[i] / p[i]) ** (alpha+1)
                        B_u_i = B_u_i + factorial(Q[i]) / factorial(Q[i] - (alpha+1)) * (local_gamma_u / local_p_u) ** (alpha+1)
                        B_d_im1 = B_d_im1 + factorial(Q[i]) / factorial(Q[i] - (alpha+1)) * (gamma_d[i-1] / p_d[i-1]) ** (alpha+1)
                        for l in range(0, alpha - 1):
                            B_i = B_i + factorial(Q[i] - (l+1)) / factorial(Q[i] - (alpha+1)) * (gamma[i] / p[i]) ** ((alpha+1) - (l+1)) * (values[i-1].p_starving[l] / values[i-1].p_dn_fail + values[i].p_blocking[l] / values[i].p_up_fail)
                            B_u_i = B_u_i + factorial(Q[i] - (l+1)) / factorial(Q[i] - (alpha+1)) * (local_gamma_u / local_p_u) ** ((alpha+1) - (l+1)) * (values[i].p_blocking[l] / values[i].p_up_fail)
                            B_d_im1 = B_d_im1 + factorial(Q[i] - (l+1)) / factorial(Q[i] - (alpha+1)) * (gamma_d[i-1] / p_d[i-1]) ** ((alpha+1) - (l+1)) * (values[i-1].p_starving[l] / values[i-1].p_dn_fail)

                    p_u[i] = current_omega * local_p_u + (1 - current_omega) * (p[i] + mu_d[i-1] * min(max(values[i-1].p_101, values[i].p_101 / (TP[i - 1] / local_mu_u - values[i].p_sum1)), 1) + p_u[i-1] * min(max(values[i-1].p_01A, values[i-1].p_01A / (TP[i - 1] / local_mu_u - values[i].p_sum1)), 1))
                    gamma_u[i] = current_omega * local_gamma_u + (1 - current_omega) * (gamma[i] + (Q[i - 1] / Q[i] * gamma_u[i-1] - gamma[i]) * min(max(values[i-1].p_00A, values[i-1].p_00A / (local_p_u / local_gamma_u / Q[i] * (TP[i-1] / local_mu_u - values[i].p_sum1))), 1))
                    mu_u[i] = current_omega * local_mu_u + (1 - current_omega) * (((1 + B_u_i) / B_u_i) / (1 / TP[i-1] + 1 / mu[i] * (1 + B_i) / B_i - 1 / mu_d[i-1] * (1 + B_d_im1) / B_d_im1))

                    # Updating the values of the current line
                    tml.fill_tml(C[i], Q_u[i], Q_d[i], mu_u[i], mu_d[i], p_u[i], p_d[i], gamma_u[i], gamma_d[i])
                    tml.solve_2M()
                    values[i] = tml.characteristics_2M
                    TP[i] = tml.characteristics_2M.TP

                    if np.sum(tml.sol_2M < 0) != 0:
                        print('Problem with solution!')


                    if (abs(TP[i] - local_TP_i) / local_TP_i < 0.001 or
                            abs(TP[i] - local_TP_i) / local_TP_i < 0.01 and iteration_FP > 20 or
                            abs(TP[i] - local_TP_i) / local_TP_i < 0.1 and iteration_FP > 50):
                        no_convergence = 0

                    if iteration_FP > 100 and no_convergence !=0:
                        # Convergence by force
                        print('Convergence by force upstream')
                        TP[i] = (TP[i] + local_TP_i) / 2
                        no_convergence = 0

                    if abs(TP[i] - local_TP_i) / local_TP_i < 0.1:
                        current_omega = 0

            if echo_level >= 3:
                TP_str = ', '.join(map('{:.4f}'.format, TP))
                mu_u_str = ', '.join(map('{:.4f}'.format, mu_u))
                p_u_str = ', '.join(map('{:.4f}'.format, p_u))
                gamma_u_str = ', '.join(map('{:.4f}'.format, gamma_u))
                print(f'Forward pass completed (iteration #{iteration}).\n TP: \t\t[{TP_str}]\n '
                      f'mu_u: \t\t[{mu_u_str}]\n p_u: \t\t[{p_u_str}]\n gamma_u: \t[{gamma_u_str}]\n\n')

            # Backward pass
            for i in range(machines-3, -1, -1):
                iteration_BP = 0
                no_convergence = 1
                current_omega = omega

                tml.fill_tml(C[i], Q_u[i], Q_d[i], mu_u[i], mu_d[i], p_u[i], p_d[i], gamma_u[i], gamma_d[i])
                tml.solve_2M()
                values[i] = tml.characteristics_2M
                TP[i] = tml.characteristics_2M.TP

                while no_convergence:
                    # Fixed point iteration
                    iteration_BP = iteration_BP + 1
                    local_p_d = p_d[i]
                    local_gamma_d = gamma_d[i]
                    local_mu_d = mu_d[i]
                    local_TP_i = TP[i]

                    B_ip1 = 0
                    B_u_ip1 = 0
                    B_d_i = 0
                    for alpha in range(0, Q[i + 1]):
                        B_ip1 = B_ip1 + factorial(Q[i+1]) / factorial(Q[i+1] - (alpha+1)) * (gamma[i+1] / p[i+1]) ** (alpha+1)
                        B_u_ip1 = B_u_ip1 + factorial(Q[i+1]) / factorial(Q[i+1] - (alpha+1)) * (gamma_u[i+1] / p_u[i+1]) ** (alpha+1)
                        B_d_i = B_d_i + factorial(Q[i+1]) / factorial(Q[i+1] - (alpha+1)) * (local_gamma_d / local_p_d) ** (alpha+1)
                        for l in range(0, alpha - 1):
                            B_ip1 = B_ip1 + factorial(Q[i+1] - (l+1)) / factorial(Q[i+1] - (alpha+1)) * (gamma[i+1] / p[i+1]) ** ((alpha+1) - (l+1)) * (values[i].p_starving[l] / values[i].p_dn_fail + values[i+1].p_blocking[l] / values[i+1].p_up_fail)
                            B_u_ip1 = B_u_ip1 + factorial(Q[i+1] - (l+1)) / factorial(Q[i+1] - (alpha+1)) * (gamma_u[i+1] / p_u[i+1]) ** ((alpha+1) - (l+1)) * (values[i+1].p_blocking[l] / values[i+1].p_up_fail)
                            B_d_i = B_d_i + factorial(Q[i+1] - (l+1)) / factorial(Q[i+1] - (alpha+1)) * (local_gamma_d / local_p_d) ** ((alpha+1) - (l+1)) * (values[i].p_starving[l] / values[i].p_dn_fail)

                    gamma_d[i] = current_omega * local_gamma_d + (1 - current_omega) * (gamma[i+1] + (Q[i + 2] / Q[i+1] * gamma_d[i+1] - gamma[i+1]) * min(max(values[i+1].p_NA0, values[i+1].p_NA0 / (local_p_d / local_gamma_d / Q[i+1] * (TP[i+1] / local_mu_d - values[i].p_sum2))), 1))
                    p_d[i] = current_omega * local_p_d + (1 - current_omega) * (
                    p[i+1] + mu_u[i+1] * min(max(values[i+1].p_Nx10, values[i+1].p_Nx10 / (TP[i+1] / local_mu_d - values[i].p_sum2)), 1) + p_d[i+1] * min(max(values[i+1].p_NA1, values[i+1].p_NA1 / (TP[i+1] / local_mu_d - values[i].p_sum2)), 1))
                    mu_d[i] = current_omega * local_mu_d + (1 - current_omega) * (((1 + B_d_i) / B_d_i) / (1 / TP[i+1] + 1 / mu[i+1] * (1 + B_ip1) / B_ip1 - 1 / mu_u[i+1] * (1 + B_u_ip1) / B_u_ip1))

                    # Updating the values of the current line
                    tml.fill_tml(C[i], Q_u[i], Q_d[i], mu_u[i], mu_d[i], p_u[i], p_d[i], gamma_u[i], gamma_d[i])
                    tml.solve_2M()
                    values[i] = tml.characteristics_2M
                    TP[i] = tml.characteristics_2M.TP

                    if abs(TP[i] - local_TP_i) / local_TP_i < 0.001 or (abs(TP[i] - local_TP_i) / local_TP_i < 0.01 and iteration_BP > 20) or (abs(TP[i] - local_TP_i) / local_TP_i < 0.1 and iteration_BP > 50):
                        no_convergence = 0

                    if iteration_BP > 2000 and no_convergence !=0:
                        # Convergence by force
                        print('Convergence by force downstream')
                        TP[i] = (TP[i] + local_TP_i) / 2
                        no_convergence = 0

                    if abs(TP[i] - local_TP_i) / local_TP_i < 0.1:
                        current_omega = 0

            if echo_level >= 3:
                TP_str = ', '.join(map('{:.4f}'.format, TP))
                mu_d_str = ', '.join(map('{:.4f}'.format, mu_d))
                p_d_str = ', '.join(map('{:.4f}'.format, p_d))
                gamma_d_str = ', '.join(map('{:.4f}'.format, gamma_d))
                print(f'Backward pass completed (iteration #{iteration}).\n TP: \t\t[{TP_str}]\n '
                      f'mu_d: \t\t[{mu_d_str}]\n p_d: \t\t[{p_d_str}]\n gamma_d: \t[{gamma_d_str}]\n\n')

            TP_diff_old = TP_diff
            TP_diff = abs(TP[0] - TP[-1])
            data.append(TP_diff)

            TP_diff_vec = TP_old - TP
            if max(abs(TP_diff_vec)) < epsilon / 100:
                print('Convergence to wrong value detected. Stopping decomposition.\n\n')
                terminated_normally = False
                epsilon = TP_diff

            TP_old = np.copy(TP)

            # check for cycling
            TP_diff_cyc.append(round(abs(TP[0] - TP[-1]), 8))

            for step in range(2,6):
                lastelement = len(TP_diff_cyc)
                if lastelement - 2 * step + 1 > 0:
                    if np.all(TP_diff_cyc[lastelement - 2 * step:lastelement - step - 1] == TP_diff_cyc[lastelement - step: lastelement - 1]):
                        cycling = True
                        TP_diff_cyc = []

            if cycling:
                print(f'Cycling detected (no further convergence). Adjusting epsilon from {epsilon:.6f} to {epsilon * 10:.6f}\n\n')
                cycling = False
                epsilon *= 10
                if epsilon >= 0.1 * mu[0]:
                    terminated_normally = False

            if TP_diff == TP_diff_old:
                print(f'No change in TP (no further convergence). Adjusting epsilon from {epsilon:.6f} to {epsilon * 10:.6f}\n\n')
                epsilon *= 10
                if epsilon >= 0.1 * mu[0]:
                    terminated_normally = False

            if iteration % 50 == 0:
                print(f'No convergence after iteration {iteration}. Adjusting epsilon from {epsilon:.6f} to {epsilon * 10:.6f}\n\n')
                epsilon *= 10
                if epsilon >= 0.1 * mu[1]:
                    terminated_normally = False

            if iteration > 200:
                print('No convergence at all.\n\n')
                terminated_normally = False

            if echo_level >= 2:
                print(f'Decomposition iteration #{iteration} completed. Current difference in TP is {TP_diff:.6f} (old was {TP_diff_old:.6f}).\n\n')

        ## prepare
        TP_final = TP[-1]
        wip = []
        inventory_spares1 = []
        inventory_spares2 = []
        for i in range(0, machines - 1):
            wip.append(values[i].wip)
            inventory_spares1.append(values[i].inventory_spares[0])
            inventory_spares2.append(values[i].inventory_spares[1])

        result = types.SimpleNamespace()
        result.TP = TP_final
        result.iteration = iteration
        result.epsilon = epsilon
        result.terminated_normally = terminated_normally
        result.wip = wip
        result.wip_total = np.sum(result.wip)
        result.inventory_spares = inventory_spares1
        result.inventory_spares.append(inventory_spares2[-1])
        result.inventory_spares_total = np.sum(result.inventory_spares)

        if echo_level >= 1:
            wip_str = ', '.join(map('{:.2f}'.format, result.wip))
            inventory_spares_str = ', '.join(map('{:.2f}'.format, result.inventory_spares))
            if terminated_normally:
                print(f'Decomposition converged after iteration #{iteration} with epsilon = {epsilon:.6f}.')
            else:
                print(f'Decomposition _terminated_ after iteration #{iteration} with epsilon = {epsilon:.6f}.')
            print(f'The final TP is {TP_final:.4f}, \t WIP: {result.wip_total:.2f}, [{wip_str}], \t inventories: {result.inventory_spares_total:.2f}, [{inventory_spares_str}].\n\n')

        return result
    
    
    ## state numbers
    def calculate_state_number(self, n, m1, m2, Q1, Q2, C=np.infty):
        if n < 0 or m1 < 0 or m2 < 0 or n > (C + 2) or m1 > Q1 or m2 > Q2:
            return -1
        else:
            return m2 + (Q2 + 1) * m1 + (Q2 + 1) * (Q1 + 1) * n
    
    
    ## two-machine line
    def spare_2M(self):
        self.check()

        C = self.C[0]
        Q1 = self.Q[0]
        Q2 = self.Q[1]
        mu1 = self.mu[0]
        mu2 = self.mu[1]
        p1 = self.p[0]
        p2 = self.p[1]
        gamma1 = self.gamma[0]
        gamma2 = self.gamma[1]

        N = C + 2
        number_of_states_total = (N + 1) * (Q1 + 1) * (Q2 + 1)
        number_of_states = (N + 1) * (Q1 + 1) * (Q2 + 1) - (Q1 + 1) - (Q2 + 1)
        Q = np.zeros([number_of_states_total, number_of_states_total])

        ## get generator matrix
        for n in range(0, N + 1):
            for m1 in range(0, Q1 + 1):
                for m2 in range(0, Q2 + 1):
                    # calculate state numbers
                    f = self.calculate_state_number(n, m1, m2, Q1, Q2, C)

                    # fill line of Q
                    if (n == N and m1 == 0) or (n == 0 and m2 == 0):
                        # prepare Q for impossible states
                        Q[f, f] = np.nan
                    else:
                        fnp = self.calculate_state_number(n + 1, m1, m2, Q1, Q2, C)
                        fnm = self.calculate_state_number(n - 1, m1, m2, Q1, Q2, C)
                        f1p = self.calculate_state_number(n, m1 + 1, m2, Q1, Q2, C)
                        f1m = self.calculate_state_number(n, m1 - 1, m2, Q1, Q2, C)
                        f2p = self.calculate_state_number(n, m1, m2 + 1, Q1, Q2, C)
                        f2m = self.calculate_state_number(n, m1, m2 - 1, Q1, Q2, C)
                        if fnm != -1 and m1 != 0:
                            Q[f, fnm] = mu1
                        if fnp != -1 and m2 != 0:
                            Q[f, fnp] = mu2
                        if f1p != -1 and n != N:
                            Q[f, f1p] = p1
                        if f2p != -1 and n != 0:
                            Q[f, f2p] = p2
                        if f1m != -1:
                            Q[f, f1m] = (Q1 - m1 + 1) * gamma1
                        if f2m != -1:
                            Q[f, f2m] = (Q2 - m2 + 1) * gamma2

        ## solve
        Q_mod = Q

        # remove rows and columns of impossible states
        possible = ~np.isnan(np.diag(Q_mod))
        Q_mod = Q_mod[possible, :]
        Q_mod = Q_mod[:, possible]
        impossible_state_numbers = np.arange(0, number_of_states_total)
        impossible_state_numbers = impossible_state_numbers[~possible]

        # add main diagonal elements
        Q_mod = Q_mod + np.diag(-np.sum(Q_mod, axis=0))

        # prepare right-hand side
        t = np.zeros(number_of_states)

        # add normalization condition (new: overwrite one balance equation; system is overdetermined)
        Q_mod[number_of_states - 1, :] = np.ones(number_of_states)
        t[number_of_states - 1] = 1
        # Q_mod[number_of_states, :] = np.ones(number_of_states)
        # t[number_of_states] = 1

        # solve
        r = np.linalg.solve(Q_mod, t)

        ## transform solution from vector to matrix for easier access
        if np.sum(np.abs(r)) - 1 > 10 ** (-10):
            print('Probabilities do not sum up to 1!')
            raise

        sol = np.zeros([N + 1, Q1 + 1, Q2 + 1])
        minus = 0
        for n in range(0, N + 1):
            for m1 in range(0, Q1 + 1):
                for m2 in range(0, Q2 + 1):
                    f = self.calculate_state_number(n, m1, m2, Q1, Q2, C)
                    if f in impossible_state_numbers:
                        minus = minus + 1
                        sol[n, m1, m2] = 0
                    else:
                        sol[n, m1, m2] = r[f - minus]

        return sol
    
    
    ## characteristics of two-machine lines
    def get_characteristics_2M(self, sol):
        C = self.C[0]
        Q1 = self.Q[0]
        Q2 = self.Q[1]
        mu2 = self.mu[1]
        N = C + 2

        characteristics = types.SimpleNamespace()
        characteristics.terminated_normally = True
        characteristics.p_blocking = np.sum(sol[-1, :, :], axis=1)[1:]
        characteristics.p_starving = np.sum(sol[0, :, :], axis=0)[1:]
        characteristics.p_up_fail = np.sum(sol[:, 0, :])
        characteristics.p_dn_fail = np.sum(sol[:, :, 0])
        characteristics.p_001 = sol[0, 0, 1]
        characteristics.p_011 = sol[0, 1, 1]
        characteristics.p_101 = sol[1, 0, 1]
        characteristics.p_Nx10 = sol[-2, 1, 0]
        characteristics.p_N10 = sol[-1, 1, 0]
        characteristics.p_N11 = sol[-1, 1, 1]

        characteristics.p_00A = np.sum(sol[0, 0, :])
        characteristics.p_01A = np.sum(sol[0, 1, :])
        characteristics.p_NA0 = np.sum(sol[-1, :, 0])
        characteristics.p_NA1 = np.sum(sol[-1, :, 1])

        characteristics.p_sum1 = np.sum(sol[1:-1, 2:, :]) + np.sum(sol[:-1, 2:, 1:]) - np.sum(sol[1:-1, 2:, 1:])
        characteristics.p_sum2 = np.sum(sol[1:-1, :, 2:]) + np.sum(sol[1:, 1:, 2:]) - np.sum(sol[1:-1, 1:, 2:])

        characteristics.TP = (np.sum(sol[:, :, 1:]) - np.sum(sol[0, :, 1:])) * mu2

        characteristics.p_processing = [np.sum(sol[0:-1, 1:, :]), np.sum(sol[1:, :, 1:])]
        characteristics.p_available = [np.sum(sol[:, 1:, :]), np.sum(sol[:, :, 1:])]
        characteristics.p_blocking_total = np.sum(sol[-1, :, :])
        characteristics.p_starving_total = np.sum(sol[0, :, :])

        range1 = np.arange(Q1 + 1) - 1
        range1[0] = 0
        range2 = np.arange(Q2 + 1) - 1
        range2[0] = 0
        characteristics.wip = np.sum(np.sum(np.sum(sol, axis=1), axis=1) * np.arange(N + 1))

        characteristics.inventory_spares = np.array([
            np.sum(np.sum(np.sum(sol, axis=0), axis=1) * range1),
            np.sum(np.sum(np.sum(sol, axis=0), axis=0) * range2),
        ])

        return characteristics
