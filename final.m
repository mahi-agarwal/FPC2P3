% FPC 2 Part 3: 
% Supersonic Civilian Transport Nozzle Optimization Code
% Team Nozzle Be There — Diego Delgado, Faz Zaidi, Mahi Agarwal
% Application: Supersonic Civilian Transport

clear all; close all; clc;

% Defining gas properties for air
gamma = 1.4;
R     = 287;

% Chosen operating Conditions for cruise
P_0_cruise   = 110000;   % stagnation pressure (Pa)
T_0_cruise   = 1000;     % stagnation temperature (K)
P_inf_cruise = 7505;     % ambient at 18 km (Pa)

% Net Pressure Ratio 
NPR = P_0_cruise / P_inf_cruise;

% Finding Mach number for ideal perfect expansion
Me_perfect = sqrt(2/(gamma-1) * (NPR^((gamma-1)/gamma) - 1));
fprintf('Cruise: P_0=%.0f kPa, T_0=%.0f K, P_inf=%.1f kPa\n', ...
        P_0_cruise/1e3, T_0_cruise, P_inf_cruise/1e3);
fprintf('NPR = %.1f, Perfectly-expanded Me = %.2f\n', NPR, Me_perfect);

% Optimization weight
w_pressure = 1.5;

% We normalized our dimensions using the throat radius
r_t = 1;

% Conical geometry with parameters half-angle theta (deg) 
% and length L
theta_vec = linspace(3, 15, 40);
L_vec     = linspace(1.0, 4, 40);

% Bell geometry with parameters exit radius and 
% exit half-angle, with a fixed L = 5*r_t
re_vec      = linspace(1.1, 1.5, 40);
theta_e_vec = linspace(2, 10, 40);

% Conical geometry optimization using function
results_conical = optimize_nozzle('conical', theta_vec, L_vec, ...
    P_0_cruise, T_0_cruise, P_inf_cruise, gamma, R, r_t, w_pressure);

% Bell geometry optimization using function
results_bell = optimize_nozzle('bell', re_vec, theta_e_vec, ...
    P_0_cruise, T_0_cruise, P_inf_cruise, gamma, R, r_t, w_pressure);

% Generating figures using the function we defined
plot_results(results_conical, results_bell, ...
             P_0_cruise, T_0_cruise, P_inf_cruise, gamma, R, r_t);

function qoi = compute_qois(sol, P_inf, P_0, gamma, R)

% This function is used to solve for our quanitites of interest for this
% application
%
% Inputs:
%       sol   - output from solve_nozzle
%       P_inf - ambient pressure (Pa)
%       P_0   - stagnation pressure (Pa)
%       gamma - ratio of specific heats (1.4 for air)
%       R     - specific gas constant (J/kg.K)
%
% Output: struct data type for QoIs:
%       C_F       - thrust coefficient
%       V_e       - exit velocity (m/s)
%       Pe_Pinf   - exit pressure ratio P_e / P_inf
%       regime    - flow regime string
%       M_e       - exit Mach number
%       T_e       - exit temperature (K)

% In Matlab, nargin means "number of arguments in"
% We use this to ensure a default setup in case the number of arguments
% entered are not as expected

    if nargin < 4, gamma = 1.4; end
    if nargin < 5, R = 287; end

    A_t = sol.A(1);
    A_e = sol.A(end);
    P_e = sol.P(end);
    V_e = sol.V(end);
    M_e = sol.M(end);
    T_e = sol.T(end);
    rho_e = sol.rho(end);

    % Thrust coefficient (formula is also defined in our slides): 
    % C_F = (m_dot * V_e + (P_e - P_inf) * A_e) / (P_0 * A_t)
    % 
    % Mass flow rate at throat (choked): 
    % m_dot = rho_e * V_e * A_e

    % Then, we can find the necessary values
    m_dot = rho_e * V_e * A_e;
    thrust = m_dot * V_e + (P_e - P_inf) * A_e;
    C_F = thrust / (P_0 * A_t);

    % Exit pressure ratio
    Pe_Pinf = P_e / P_inf;

    % Output using the qoi struct
    qoi.C_F     = C_F;
    qoi.V_e     = V_e;
    qoi.Pe_Pinf = Pe_Pinf;
    qoi.regime  = sol.regime;
    qoi.M_e     = M_e;
    qoi.T_e     = T_e;
end

function [T_ratio, P_ratio, rho_ratio, A_ratio] = isentropic_flow(M, gamma)
% This function is used to compute stagnation ratio values 
% using isentropic relations for this flow
%
%   Inputs:
%       M     - Mach number
%       gamma - 1.4 for air
%
%   Outputs:
%       T_ratio   - T/T0 
%       P_ratio   - P/P0 
%       rho_ratio - rho/rho0
%       A_ratio   - A/A*

    if nargin < 2
        gamma = 1.4;
    end

    % Isentropic relations
    % To simplify, we are computing the recurring factor separately
    factor = 1 + (gamma - 1) / 2 .* M.^2;

    % Then, we are using the isentropic relations definitions to find the
    % ratios
    T_ratio   = 1 ./ factor;
    P_ratio   = factor .^ (-gamma / (gamma - 1));
    rho_ratio = factor .^ (-1 / (gamma - 1));

    % Defined Area-Mach relation: 
    % A/A* = (1/M) * [(2/(gamma+1)) * (1 + (gamma-1)/2 * M^2)]^((gamma+1)/(2*(gamma-1)))
    % Again, separately finding the exponent simplifies things
    exponent = (gamma + 1) / (2 * (gamma - 1));
    % Finally, solving for the Area ratio
    A_ratio  = (1 ./ M) .* ((2 / (gamma + 1)) .* factor) .^ exponent;
end

function M = mach_from_area(A_ratio, gamma, branch)
% This function helps us employ the Area-Mach relation discussed in class
% to find the area ratio which will help us optimize for both geometry
% scenarios
%
%   Inputs:
%       A_ratio - A/A* (needs to be >= 1)
%       gamma   - 1.4 for air
%       branch  - 'sub' for subsonic flow and 'sup' for supersonic flow
%
%   Output:
%       M - Mach number

    if nargin < 2 || isempty(gamma)
        gamma = 1.4;
    end
    if nargin < 3 || isempty(branch)
        branch = 'sup';
    end

    M = zeros(size(A_ratio));

    for i = 1:numel(A_ratio)
        AR = A_ratio(i);

        if AR < 1 - 1e-10
            error('A/A* = %.4f is less than 1. Not physical.', AR);
        end

        if abs(AR - 1) < 1e-10
            M(i) = 1.0;
            continue;
        end

        % Define residual: f(M) = A/A*(M) - target A/A*
        exponent = (gamma + 1) / (2 * (gamma - 1));
        f = @(m) (1/m) * ((2/(gamma+1)) * (1 + (gamma-1)/2 * m^2))^exponent - AR;

        if strcmpi(branch, 'sub')
            % Subsonic root: M in (0, 1)
            M(i) = fzero(f, [1e-6, 1 - 1e-10]);
        else
            % Supersonic root: M in (1, ~large)
            % Upper bound estimate: for large M, A/A* ~ M, so M_max ~ 2*AR
            M_upper = max(2, 2 * AR);
            M(i) = fzero(f, [1 + 1e-10, M_upper]);
        end
    end
end

function [M2, P2_P1, T2_T1, rho2_rho1, P02_P01] = normal_shock(M1, gamma)
% This function helps us find the properties across a normal shock
%
%   Inputs:
%       M1    - upstream Mach number (for supersonic i.e. >1)
%       gamma - 1.4 for air
%
%   Outputs:
%       M2       - downstream Mach number
%       P2_P1    - static pressure ratio across shock
%       T2_T1    - static temperature ratio across shock
%       rho2_rho1 - density ratio across shock
%       P02_P01  - stagnation pressure ratio (< 1)

    if nargin < 2
        gamma = 1.4;
    end

    g = gamma;

    % Here, we use known normal shock relations
    M2 = sqrt((1 + (g-1)/2 .* M1.^2) ./ (g .* M1.^2 - (g-1)/2));

    P2_P1 = 1 + 2*g/(g+1) .* (M1.^2 - 1);

    T2_T1 = P2_P1 .* (2 + (g-1) .* M1.^2) ./ ((g+1) .* M1.^2);

    rho2_rho1 = (g+1) .* M1.^2 ./ (2 + (g-1) .* M1.^2);

    % Then, we find the stagnation pressure ratio
    P02_P01 = ((((g+1).*M1.^2) ./ ((g-1).*M1.^2 + 2)).^(g/(g-1))) ...
              .* ((2*g.*M1.^2 - (g-1)) ./ (g+1)).^(-1/(g-1));
end


function [x, r, A] = nozzle_geometry(geom_type, params, r_t, N)
% This function helps us define the nozzle geometries we have chosen for
% this project -- conical and bell

    if nargin < 3 || isempty(r_t), r_t = 1; end
    if nargin < 4 || isempty(N), N = 200; end

    % We used switch case for this(given our previous knowledge from other 
    % coding languages such as Java and Python)

    % The switch case choses between the conical or bell geometry
    switch lower(geom_type)
        case 'conical'
            % The conical geometry takes in theta and length as its
            % parameters
            % Then, we use the function defined for this geometry, which
            % can also be found on our slides
            theta_deg = params(1);
            L = params(2);
            x = linspace(0, L, N);
            r = r_t + x .* tan(deg2rad(theta_deg));

        case 'bell'
            % For this case, the parameters are the exit radius, exit
            % theta, and length 
            r_e            = params(1);
            theta_exit_deg = params(2);
            L              = params(3);
            theta_exit = deg2rad(theta_exit_deg);
            % Then, we solve for the coefficients in the parabolic bell
            % curve equation that we set for the geometry
            a = (L * tan(theta_exit) - (r_e - r_t)) / L^2;
            b = tan(theta_exit) - 2 * a * L;
            x = linspace(0, L, N);
            r = r_t + a .* x.^2 + b .* x;

        % This is for a "default" or exception if the case inputs are
        % neither conical nor bell
        otherwise
            error('Unknown geometry type: %s', geom_type);
    end
    % Using radius to solve for area at that point
    A = pi .* r.^2;
end

function sol = solve_nozzle(x, A, P_b, P_0, T_0, gamma, R)
% SOLVE_NOZZLE Solve quasi-1D steady-state flow through a C-D nozzle.
%
%   sol = solve_nozzle(x, A, P_b, P_0, T_0, gamma, R)
%
%   Inputs:
%       x     - axial coordinate vector (diverging section, x=0 at throat)
%       A     - area distribution A(x), same size as x
%       P_b   - back pressure (Pa)
%       P_0   - stagnation pressure (Pa)
%       T_0   - stagnation temperature (K)
%       gamma - ratio of specific heats (default 1.4)
%       R     - specific gas constant (default 287 J/kg/K)
%
%   Output: struct sol with fields:
%       M, P, T, rho, V  - flow properties along x
%       regime            - 'subsonic', 'shock_in_nozzle', 'overexpanded',
%                           'perfectly_expanded', 'underexpanded'
%       shock_x           - shock location (NaN if no shock inside nozzle)
%       shock_idx         - index of shock location

    if nargin < 6 || isempty(gamma), gamma = 1.4; end
    if nargin < 7 || isempty(R), R = 287; end

    N = length(x);
    A_t = A(1);          % throat area (x=0)
    A_e = A(end);        % exit area
    A_star = A_t;        % sonic throat
    AR = A ./ A_star;    % area ratio distribution
    AR_e = A_e / A_star; % exit area ratio

    Pb_P0 = P_b / P_0;

    % Computing critical pressure ratios to classify the flows

    % Case 1: Fully subsonic
    M_sub_exit = mach_from_area(AR_e, gamma, 'sub');
    [~, P_sub_exit] = isentropic_flow(M_sub_exit, gamma);
    % P_sub_exit is P_e/P_0

    % Case 2: Normal shock at exit plane
    M_sup_exit = mach_from_area(AR_e, gamma, 'sup');
    [~, P_sup_exit_isen] = isentropic_flow(M_sup_exit, gamma);
    [M_post_shock, P2_P1, ~, ~, ~] = normal_shock(M_sup_exit, gamma);
    P_shock_at_exit = P_sup_exit_isen * P2_P1;
    % P_shock_at_exit is P_e/P_0

    % Case 3: Fully supersonic and isentropic
    % P_sup_exit_isen is P_e/P_0

    % Classifying the regime based on back pressure

    if Pb_P0 >= P_sub_exit
        % Back pressure is too high: fully subsonic in diverging section
        % (or subsonic with some deceleration)
        sol = solve_subsonic(x, AR, P_0, T_0, gamma, R);
        sol.regime = 'subsonic';
        sol.shock_x = NaN;
        sol.shock_idx = NaN;

    elseif Pb_P0 >= P_shock_at_exit
        % Normal shock inside the diverging section
        sol = solve_with_shock(x, AR, Pb_P0, P_0, T_0, gamma, R);
        sol.regime = 'shock_in_nozzle';

    elseif Pb_P0 >= P_sup_exit_isen
        % Overexpanded: supersonic throughout, P_e < P_b
        sol = solve_supersonic(x, AR, P_0, T_0, gamma, R);
        sol.regime = 'overexpanded';
        sol.shock_x = NaN;
        sol.shock_idx = NaN;

    elseif abs(Pb_P0 - P_sup_exit_isen) < 1e-6
        % Perfectly expanded -- what we want for this application
        sol = solve_supersonic(x, AR, P_0, T_0, gamma, R);
        sol.regime = 'perfectly_expanded';
        sol.shock_x = NaN;
        sol.shock_idx = NaN;

    else
        % Underexpanded: P_e > P_b
        sol = solve_supersonic(x, AR, P_0, T_0, gamma, R);
        sol.regime = 'underexpanded';
        sol.shock_x = NaN;
        sol.shock_idx = NaN;
    end

    % Returning outputs
    sol.x = x;
    sol.A = A;
end


function sol = solve_subsonic(x, AR, P_0, T_0, gamma, R)
% Fully subsonic flow
    N = length(x);
    M = zeros(1, N);
    for i = 1:N
        M(i) = mach_from_area(AR(i), gamma, 'sub');
    end
    [T_ratio, P_ratio, rho_ratio, ~] = isentropic_flow(M, gamma);
    sol.M   = M;
    sol.P   = P_ratio .* P_0;
    sol.T   = T_ratio .* T_0;
    sol.rho = sol.P ./ (R .* sol.T);
    sol.V   = M .* sqrt(gamma .* R .* sol.T);
end


function sol = solve_supersonic(x, AR, P_0, T_0, gamma, R)
% Fully supersonic and isentropic flow
    N = length(x);
    M = zeros(1, N);
    M(1) = 1.0; % sonic at throat so M = 1
    for i = 2:N
        M(i) = mach_from_area(AR(i), gamma, 'sup');
    end
    [T_ratio, P_ratio, rho_ratio, ~] = isentropic_flow(M, gamma);
    sol.M   = M;
    sol.P   = P_ratio .* P_0;
    sol.T   = T_ratio .* T_0;
    sol.rho = sol.P ./ (R .* sol.T);
    sol.V   = M .* sqrt(gamma .* R .* sol.T);
end


function sol = solve_with_shock(x, AR, Pb_P0, P_0, T_0, gamma, R)
% This function is used to solve for flow properties in a nozzle with the
% shock present in the diverging section

    N = length(x);
    AR_e = AR(end);

    % Search for shock location by trying each station in the diverging section
    % At station i, supersonic Mach is M1. After shock, M2 and P02/P01
    % Post-shock flow decelerates subsonically to exit with new A/A*
    % Find where exit pressure matches P_b for perfect expansion output

    % Building the supersonic Mach distribution first
    M_sup = zeros(1, N);
    M_sup(1) = 1.0;
    for i = 2:N
        M_sup(i) = mach_from_area(AR(i), gamma, 'sup');
    end

    % Compute exit pressure for shock at each station
    P_exit_ratio = zeros(1, N);
    for i = 2:N-1
        M1 = M_sup(i);
        [M2, ~, ~, ~, P02_P01] = normal_shock(M1, gamma);

        % New A* after shock (A*_new = A_shock / A_ratio(M2))
        [~, ~, ~, A_ratio_M2] = isentropic_flow(M2, gamma);
        A_star_new = AR(i) * (pi * 1^2) / A_ratio_M2; % in terms of A_t
        % Actually, we need A/A*_new at the exit
        AR_exit_new = AR(end) * (pi * 1^2) / A_star_new;

        % Subsonic Mach at exit with new A*
        if AR_exit_new >= 1
            M_exit_sub = mach_from_area(AR_exit_new, gamma, 'sub');
            [~, P_exit_isen] = isentropic_flow(M_exit_sub, gamma);
            % P_exit / P_0 = (P_exit/P02) * (P02/P01) * (P01/P_0)
            % P01 = P_0 (isentropic upstream of shock)
            P_exit_ratio(i) = P_exit_isen * P02_P01;
        else
            P_exit_ratio(i) = NaN;
        end
    end

    % Find shock location by interpolation: where P_exit_ratio = Pb_P0
    valid = ~isnan(P_exit_ratio(2:end-1));
    idx_range = 2:N-1;
    idx_valid = idx_range(valid);

    if isempty(idx_valid)
        % Fallback: no valid shock location found, treat as supersonic
        sol = solve_supersonic(x, AR, P_0, T_0, gamma, R);
        sol.shock_x = NaN;
        sol.shock_idx = NaN;
        return;
    end

    % Interpolate to find exact shock index
    P_valid = P_exit_ratio(idx_valid);
    % P_exit_ratio increases as shock moves toward throat (higher M1 = stronger shock)
    % Find where P_valid crosses Pb_P0
    shock_idx = idx_valid(1); % default
    for k = 1:length(P_valid)-1
        if (P_valid(k) - Pb_P0) * (P_valid(k+1) - Pb_P0) <= 0
            % Linear interpolation for better accuracy
            frac = (Pb_P0 - P_valid(k)) / (P_valid(k+1) - P_valid(k));
            shock_idx = idx_valid(k) + round(frac);
            break;
        end
    end
    shock_idx = max(2, min(shock_idx, N-1));

    % Build flow field
    M1_shock = M_sup(shock_idx);
    [M2_shock, ~, ~, ~, P02_P01_shock] = normal_shock(M1_shock, gamma);

    % Upstream of shock: supersonic isentropic
    M = zeros(1, N);
    M(1) = 1.0;
    for i = 2:shock_idx
        M(i) = M_sup(i);
    end

    % New A* after shock
    [~, ~, ~, A_ratio_M2] = isentropic_flow(M2_shock, gamma);
    A_star_new_ratio = AR(shock_idx) / A_ratio_M2; % A*_new / A_t

    % Downstream of shock: subsonic isentropic with new A*
    M(shock_idx) = M2_shock; % post-shock Mach (overwrites pre-shock)
    for i = shock_idx+1:N
        AR_new = AR(i) / A_star_new_ratio;
        if AR_new >= 1
            M(i) = mach_from_area(AR_new, gamma, 'sub');
        else
            M(i) = mach_from_area(max(AR_new, 1.0), gamma, 'sub');
        end
    end

    % Compute thermodynamic properties
    % Upstream of shock: use P_0
    % Downstream of shock: use P_02 = P_0 * P02_P01
    P_02 = P_0 * P02_P01_shock;

    sol.M = M;
    sol.P = zeros(1, N);
    sol.T = zeros(1, N);

    % Upstream of shock (isentropic from P_0, T_0)
    [T_r_up, P_r_up, ~, ~] = isentropic_flow(M(1:shock_idx-1), gamma);
    sol.P(1:shock_idx-1) = P_r_up .* P_0;
    sol.T(1:shock_idx-1) = T_r_up .* T_0;

    % Downstream of shock (isentropic from P_02, T_0 -- T_0 unchanged across shock)
    [T_r_dn, P_r_dn, ~, ~] = isentropic_flow(M(shock_idx:N), gamma);
    sol.P(shock_idx:N) = P_r_dn .* P_02;
    sol.T(shock_idx:N) = T_r_dn .* T_0;

    sol.rho = sol.P ./ (R .* sol.T);
    sol.V   = M .* sqrt(gamma .* R .* sol.T);
    sol.shock_x   = x(shock_idx);
    sol.shock_idx = shock_idx;
end


function results = optimize_nozzle(geom_type, param1_vec, param2_vec, ...
    P_0, T_0, P_inf, gamma, R, r_t, w)
% Our main optimization function to check for performance of the nozzle at
% different input conditions
%   J = C_F * lambda - w * (Pe/Pinf - 1)^2

    N1 = length(param1_vec);
    N2 = length(param2_vec);
    
    L_bell = 5 * r_t;

    CF       = NaN(N1, N2);
    Pe_Pinf  = NaN(N1, N2);
    Me       = NaN(N1, N2);
    Ve       = NaN(N1, N2);
    Te       = NaN(N1, N2);
    AR       = NaN(N1, N2);
    J        = NaN(N1, N2);
    regime   = strings(N1, N2);

    % We used a nested for-loop to iterate through 

    for i = 1:N1
        for j = 1:N2
            try
                if strcmpi(geom_type, 'conical')
                    params = [param1_vec(i), param2_vec(j)];
                else
                    params = [param1_vec(i), param2_vec(j), L_bell];
                end

                [x, r, A] = nozzle_geometry(geom_type, params, r_t);

                if A(end) <= A(1), continue; end
                if any(r <= 0), continue; end
                if any(diff(A) < -1e-10 * A(1)), continue; end
                if A(end)/A(1) > 10, continue; end

                AR(i,j) = A(end)/A(1);

                sol = solve_nozzle(x, A, P_inf, P_0, T_0, gamma, R);
                qoi = compute_qois(sol, P_inf, P_0, gamma, R);

                % Divergence correction
                if strcmpi(geom_type, 'conical')
                    th = atan((r(end)-r(1))/(x(end)-x(1)));
                    lam = (1 + cos(th))/2;
                else
                    lam = (1 + cos(deg2rad(param2_vec(j))))/2;
                end

                CF(i,j)      = qoi.C_F * lam;
                Pe_Pinf(i,j) = qoi.Pe_Pinf;
                Me(i,j)      = qoi.M_e;
                Ve(i,j)      = qoi.V_e;
                Te(i,j)      = qoi.T_e;
                regime(i,j)  = string(qoi.regime);
                J(i,j)       = CF(i,j) - w * (Pe_Pinf(i,j) - 1)^2;
            catch
                continue;
            end
        end
        if mod(i,10)==0, fprintf('  %d/%d\n', i, N1); end
    end

    % Grid optimum
    Js = J; Js(isnan(Js)) = -Inf;
    [~, idx] = max(Js(:));
    [oi, oj] = ind2sub([N1,N2], idx);
    fprintf('Grid opt: [%.4f, %.4f], CF=%.5f, Pe/Pinf=%.4f\n', ...
        param1_vec(oi), param2_vec(oj), CF(oi,oj), Pe_Pinf(oi,oj));

    % Local refinement with bounds
    p0 = [param1_vec(oi), param2_vec(oj)];
    opts = optimset('Display','off','TolX',1e-5,'TolFun',1e-6);
    obj = @(p) nozzle_obj(p, geom_type, P_0, T_0, P_inf, gamma, R, r_t, w, L_bell, ...
                          param1_vec, param2_vec);
    [p_ref, neg_J] = fminsearch(obj, p0, opts);

    % Evaluate refined point
    if strcmpi(geom_type,'conical')
        pr = [p_ref(1), p_ref(2)];
    else
        pr = [p_ref(1), p_ref(2), L_bell];
    end
    [xr,rr,Ar] = nozzle_geometry(geom_type, pr, r_t);
    sr = solve_nozzle(xr, Ar, P_inf, P_0, T_0, gamma, R);
    qr = compute_qois(sr, P_inf, P_0, gamma, R);
    if strcmpi(geom_type,'conical')
        th = atan((rr(end)-rr(1))/(xr(end)-xr(1)));
        lr = (1+cos(th))/2;
    else
        lr = (1+cos(deg2rad(p_ref(2))))/2;
    end

    fprintf('Refined:  [%.4f, %.4f], CF=%.5f, Pe/Pinf=%.4f\n', ...
        p_ref, qr.C_F*lr, qr.Pe_Pinf);

    % Pack results
    results.param1_vec  = param1_vec;
    results.param2_vec  = param2_vec;
    results.CF_cruise   = CF;
    results.Pe_Pinf_cruise = Pe_Pinf;
    results.Me_cruise   = Me;
    results.Ve_cruise   = Ve;
    results.Te_cruise   = Te;
    results.AR_cruise   = AR;
    results.J_cruise    = J;
    results.regime_cruise = regime;
    results.geom_type   = geom_type;
    results.opt_idx     = [oi, oj];
    results.opt_params  = [param1_vec(oi), param2_vec(oj)];
    results.opt_CF      = CF(oi,oj);
    results.opt_Pe_Pinf = Pe_Pinf(oi,oj);
    results.opt_Me      = Me(oi,oj);
    results.ref_params  = p_ref;
    results.ref_CF      = qr.C_F * lr;
    results.ref_Pe_Pinf = qr.Pe_Pinf;
    results.ref_Me      = qr.M_e;
    results.ref_Ve      = qr.V_e;
    results.ref_Te      = qr.T_e;
    results.ref_J       = -neg_J;
    results.ref_regime  = qr.regime;
end


function neg_J = nozzle_obj(p, geom_type, P_0, T_0, P_inf, gamma, R, r_t, w, L_bell, ...
                            p1_vec, p2_vec)
    try
        % Bound within grid range
        if strcmpi(geom_type, 'conical')
            if p(1) < min(p1_vec) || p(1) > max(p1_vec) || ...
               p(2) < min(p2_vec) || p(2) > max(p2_vec)
                neg_J = 1e6; return;
            end
            params = [p(1), p(2)];
        else
            if p(1) < min(p1_vec) || p(1) > max(p1_vec) || ...
               p(2) < min(p2_vec) || p(2) > max(p2_vec)
                neg_J = 1e6; return;
            end
            params = [p(1), p(2), L_bell];
        end

        [x, r, A] = nozzle_geometry(geom_type, params, r_t);
        if A(end)<=A(1) || any(r<=0) || any(diff(A)<-1e-10*A(1))
            neg_J = 1e6; return;
        end

        sol = solve_nozzle(x, A, P_inf, P_0, T_0, gamma, R);
        qoi = compute_qois(sol, P_inf, P_0, gamma, R);

        if strcmpi(geom_type,'conical')
            th = atan((r(end)-r(1))/(x(end)-x(1)));
            lam = (1+cos(th))/2;
        else
            lam = (1+cos(deg2rad(p(2))))/2;
        end

        neg_J = -(qoi.C_F * lam - w * (qoi.Pe_Pinf - 1)^2);
    catch
        neg_J = 1e6;
    end
end

function plot_results(results_conical, results_bell, ...
                      P_0_cruise, T_0_cruise, P_inf_cruise, ...
                      gamma, R, r_t)

    geom_names   = {'Conical', 'Bell'};
    all_results  = {results_conical, results_bell};
    param_labels = {{'Half-angle \theta (deg)', 'Length L / r_t'}, ...
                    {'Exit radius r_e / r_t',   'Exit half-angle \theta_{exit} (deg)'}};
    line_colors  = {[0 0.447 0.741], [0.850 0.325 0.098]};

    r_inlet = 2.0 * r_t;
    L_conv  = 1.5 * r_t;

    % Figures 1-2: C_F contour maps
    for g = 1:2
        res = all_results{g};
        [P1, P2] = meshgrid(res.param1_vec, res.param2_vec);

        regime_num = NaN(size(res.regime_cruise));
        regime_num(res.regime_cruise == "underexpanded")     = 1;
        regime_num(res.regime_cruise == "perfectly_expanded") = 2;
        regime_num(res.regime_cruise == "overexpanded")       = 3;
        regime_num(res.regime_cruise == "shock_in_nozzle")    = 4;
        regime_num(res.regime_cruise == "subsonic")           = 5;

        figure('Color','w','Position',[50+600*(g-1) 500 560 450]);
        contourf(P1', P2', res.CF_cruise, 20, 'LineStyle','none');
        hold on;
        if any(~isnan(regime_num(:)))
            [~,hc] = contour(P1', P2', regime_num, ...
                'LineColor',[0.3 0.3 0.3],'LineWidth',1,'LineStyle',':');
            set(hc,'HandleVisibility','off');
        end
        plot(res.opt_params(1), res.opt_params(2), 'rp', ...
             'MarkerSize',16,'MarkerFaceColor','r','LineWidth',1.5,...
             'DisplayName','Grid optimum');
        hold off;
        cb = colorbar; ylabel(cb,'C_F','FontName','Arial','FontSize',11);
        colormap(gca, parula);
        xlabel(param_labels{g}{1},'FontName','Arial','FontSize',12);
        ylabel(param_labels{g}{2},'FontName','Arial','FontSize',12);
        title(sprintf('%s Nozzle — Thrust Coefficient C_F', geom_names{g}), ...
              'FontName','Arial','FontSize',13);
        set(gca,'FontName','Arial');
    end

    % Figures 3-4: Pe/Pinf contour maps
    for g = 1:2
        res = all_results{g};
        [P1, P2] = meshgrid(res.param1_vec, res.param2_vec);

        figure('Color','w','Position',[50+600*(g-1) 50 560 450]);
        contourf(P1', P2', res.Pe_Pinf_cruise, 20, 'LineStyle','none');
        hold on;
        if any(res.Pe_Pinf_cruise(:)<1) && any(res.Pe_Pinf_cruise(:)>1)
            [~,hc] = contour(P1', P2', res.Pe_Pinf_cruise, [1 1], ...
                'LineColor','k','LineWidth',2.5);
            set(hc,'HandleVisibility','off');
        end
        plot(res.opt_params(1), res.opt_params(2), 'rp', ...
             'MarkerSize',16,'MarkerFaceColor','r','LineWidth',1.5,...
             'HandleVisibility','off');
        plot(res.ref_params(1), res.ref_params(2), 'gd', ...
             'MarkerSize',12,'MarkerFaceColor','g','LineWidth',1.5,...
             'HandleVisibility','off');
        hold off;
        cb = colorbar; ylabel(cb,'P_e / P_\infty','FontName','Arial','FontSize',11);
        try clim([0 3]); catch, caxis([0 3]); end
        colormap(gca, turbo);
        xlabel(param_labels{g}{1},'FontName','Arial','FontSize',12);
        ylabel(param_labels{g}{2},'FontName','Arial','FontSize',12);
        title(sprintf('%s Nozzle — Pressure Ratio P_e / P_\\infty', geom_names{g}), ...
              'FontName','Arial','FontSize',13);
        set(gca,'FontName','Arial');
    end

    % Figures 5-6: P(x), T(x), rho(x) profiles
    for g = 1:2
        res = all_results{g};
        if strcmpi(res.geom_type,'conical')
            params = [res.ref_params(1), res.ref_params(2)];
        else
            params = [res.ref_params(1), res.ref_params(2), 5*r_t];
        end

        [x,~,A] = nozzle_geometry(res.geom_type, params, r_t);
        sol = solve_nozzle(x, A, P_inf_cruise, P_0_cruise, T_0_cruise, gamma, R);

        figure('Color','w','Position',[100 100 750 650]);

        subplot(3,1,1);
        plot(x/r_t, sol.P/P_0_cruise, '-', 'Color',line_colors{g},'LineWidth',2);
        ylabel('P / P_0','FontName','Arial','FontSize',12);
        title(sprintf('%s Nozzle — Flow Properties at Cruise', geom_names{g}), ...
              'FontName','Arial','FontSize',13);
        grid on; set(gca,'FontName','Arial');

        subplot(3,1,2);
        plot(x/r_t, sol.T, '-', 'Color',[0.8 0.1 0.1],'LineWidth',2);
        ylabel('T (K)','FontName','Arial','FontSize',12);
        yline(1573,'k--','Inconel 718 limit','LineWidth',1.5,...
              'FontName','Arial','LabelHorizontalAlignment','left');
        grid on; set(gca,'FontName','Arial');

        subplot(3,1,3);
        plot(x/r_t, sol.rho, '-', 'Color',[0.2 0.6 0.2],'LineWidth',2);
        xlabel('x / r_t','FontName','Arial','FontSize',12);
        ylabel('\rho (kg/m^3)','FontName','Arial','FontSize',12);
        grid on; set(gca,'FontName','Arial');
    end

    % Figures 7-8: Individual nozzle geometry (full, with converging)
    x_conv = linspace(-L_conv, 0, 50);
    r_conv = (r_inlet+r_t)/2 + (r_inlet-r_t)/2 .* cos(pi*(x_conv+L_conv)/L_conv);

    for g = 1:2
        res = all_results{g};
        if strcmpi(res.geom_type,'conical')
            params = [res.ref_params(1), res.ref_params(2)];
        else
            params = [res.ref_params(1), res.ref_params(2), 5*r_t];
        end
        [xd, rd, ~] = nozzle_geometry(res.geom_type, params, r_t);

        xf = [x_conv, xd(2:end)];
        rf = [r_conv, rd(2:end)];

        figure('Color','w','Position',[50+600*(g-1) 200 700 400]);
        plot(xf/r_t, rf/r_t, '-','Color',line_colors{g},'LineWidth',2.5);
        hold on;
        plot(xf/r_t,-rf/r_t, '-','Color',line_colors{g},'LineWidth',2.5);
        plot(0, 1,'ko','MarkerSize',7,'MarkerFaceColor','k');
        plot(0,-1,'ko','MarkerSize',7,'MarkerFaceColor','k');
        xline(0,'k:','Throat','LineWidth',1,'FontName','Arial',...
              'LabelVerticalAlignment','top');
        hold off;
        xlabel('x / r_t','FontName','Arial','FontSize',12);
        ylabel('r / r_t','FontName','Arial','FontSize',12);
        title(sprintf('%s Nozzle Geometry', geom_names{g}),...
              'FontName','Arial','FontSize',13);
        grid on; axis equal; set(gca,'FontName','Arial');
    end

    % Figure 9: Bar chart comparison
    figure('Color','w','Position',[100 100 650 500]);
    metrics = {'C_F','M_e','P_e/P_\infty','A_e/A_t'};
    vals = [results_conical.ref_CF,      results_bell.ref_CF; ...
            results_conical.ref_Me,      results_bell.ref_Me; ...
            results_conical.ref_Pe_Pinf, results_bell.ref_Pe_Pinf; ...
            results_conical.AR_cruise(results_conical.opt_idx(1),results_conical.opt_idx(2)), ...
            results_bell.AR_cruise(results_bell.opt_idx(1),results_bell.opt_idx(2))];
    b = bar(vals);
    b(1).FaceColor = line_colors{1};
    b(2).FaceColor = line_colors{2};
    set(gca,'XTickLabel',metrics,'FontName','Arial','FontSize',11);
    ylabel('Value','FontName','Arial','FontSize',12);
    title('Conical vs Bell — Performance Comparison','FontName','Arial','FontSize',14);
    legend('Conical','Bell','Location','northwest','FontName','Arial');
    grid on;

    % Print summary
    fprintf('\n=== Results Summary ===\n');
    for g = 1:2
        res = all_results{g};
        fprintf('%s: CF=%.4f, Me=%.2f, Pe/Pinf=%.4f, Ve=%.1f m/s, Te=%.1f K\n', ...
            geom_names{g}, res.ref_CF, res.ref_Me, res.ref_Pe_Pinf, res.ref_Ve, res.ref_Te);
    end
end
