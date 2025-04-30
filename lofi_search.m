% Main file to perform a single optimisation for a randomly generated
% initial guess at low-fidelity. The user enters mission-defining constants
% and bounds on variables, and must supply ephemeris data through SPICE.
% Part of the Tool for Initial Low-Thrust Design (TILTD).
% Copyright 2022 Darcey Graham

function output = lofi_search(varargin)
assert( ~isempty(varargin), 'Not enough input arguments. Please specify an input file');


% Set some default values before loading inputs.
MBH_noLoops = 1;
rng_seed = 1;
use_parallel = false;   
displayFmincon = 'none';

% Load inputs
for i = 1:length(varargin)
    inputPath = strsplit(varargin{i},{'/','\'});
    fileName = inputPath(end);
    parentPath = inputPath(1:end-1);
    if size(parentPath)>0
        addpath(fullfile(parentPath{:}))
    end
    disp(['Reading input ''', fileName{:}, '''']);
    eval(fileName{:});
end

disp(["MBH: ", MBH_noLoops]);
disp(["Parallel: ", use_parallel]);
disp(["Seed: ", rng_seed]);

rootDir = pwd;

% User: add paths to your MICE library for SPICE
% 'fullfile' allows contructing paths OS independently.
addpath('Functions', fullfile(mice_path{:}, 'src', 'mice'), fullfile(mice_path{:}, 'lib'), rootDir);

% Load the generic kernel(s) for ephemeris data
cspice_furnsh( { fullfile(kernel_path{:} ) } );

if use_parallel
    pc = parcluster('local'); % Get cluster settings for local run.
    tmpdir = getenv('TMPDIR');
    num_cpus = str2num(getenv('SLURM_CPUS_PER_TASK'))/2;
    if isfolder(tmpdir)
        pc.JobStorageLocation = tmpdir; % Set to use system TMPDIR (faster run).
    end
    if isinteger(num_cpus) % If inside slurm script, use SLURM_CPUS_PER_TASK.
                       % Otherwise matlab decides.
        parpool(pc, num_cpus)
    else
        parpool(pc);
    end    
end

rng(rng_seed);

% Keep track of some indices for MBH
ind_vinfi = zeros(3, N_flybys);
ind_vinff = zeros(3, N_flybys);
dtIndex = zeros(1, Np);

% Control vector. Need one in the x,y,z directions for each impulse in each
% phase that permits thrust. Shouldn't need to use anything besides bounds
% of -1 and 1.
u_min = (-1)*ones(N_thrust*3, N);
u_max = ones(N_thrust*3, N);

%% No more user inputs needed, can skip to the end to name the save files
% Ephemeris data

% Min/max wider time range for gathering ephemeris data, in TU
t_ephem_max_et = cspice_str2et(t_ephem_max)/TU;
t_ephem_min_et = cspice_str2et(t_ephem_min)/TU;

% Scaled time over which to get data, in TU
et = linspace(t_ephem_min_et, t_ephem_max_et,N_ephem)*TU;

statePerts = zeros(6*N_perts, N_ephem);

for i = 1:N_perts
    solarPertBool = strcmp(pertNames{i}, 'Sun');        % This is 1 if there IS solar perturbation, 0 if there is NO solar pert
end

if solarPertBool == 1
    for i = 1:N_perts
        if i == solarPertCheck
            [statePerts((1:6)+(solarPertCheck-1)*6, :), ~] = cspice_spkezr(centralBodyName, et, 'J2000', 'LT+S', 'SOLAR SYSTEM BARYCENTER');
            statePerts((1:6)+(solarPertCheck-1)*6, :) = -statePerts((1:6)+(solarPertCheck-1)*6, :); % Take negative to flip direction
        else
            [statePerts((1:6)+(i-1)*6, :), ~] = cspice_spkezr(pertNames{i}, et, 'J2000', 'LT+S', pertRelativeTo);
        end
        statePerts((1:3)+(i-1)*6, :) = statePerts((1:3)+(i-1)*6, :)/DU;
        statePerts((4:6)+(i-1)*6, :) = statePerts((4:6)+(i-1)*6, :)/VU;
    end
else
    for i = 1:N_perts
        [statePerts((1:6)+(i-1)*6, :), ~] = cspice_spkezr(pertNames{i}, et, 'J2000', 'LT+S', pertRelativeTo);
        statePerts((1:3)+(i-1)*6, :) = statePerts((1:3)+(i-1)*6, :)/DU;
        statePerts((4:6)+(i-1)*6, :) = statePerts((4:6)+(i-1)*6, :)/VU;
    end
end

stateLong = zeros(1, 6*N_perts*N_ephem);
for i = 1:N_perts
    for j = 1:6
        stateLong((j+6*(i-1)-1)*N_ephem + 1 : (j+6*(i-1))*N_ephem) = statePerts(j+(i-1)*6,:);        
    end    
end

%% Initial guess for decision variables generated from normal distribution between bounds

% For initial phase: initial epoch, initial RA DEC, initial vinf relative to Jupiter
% For each phase: tof, mf, vinfi and vinff for the gravity assist at the
% end of the trajectory, control vector u
% For final phase: vinfi but no vinff because not doing a GA

% Initial phase
t0_min_et = cspice_str2et(t0_min)/TU;
t0_max_et = cspice_str2et(t0_max)/TU;
t0 = (t0_max_et - t0_min_et)*rand + t0_min_et;

% If start on SOI
if startBody == 1
    RAi = (RAi_max - RAi_min)*rand + RAi_min;
    DECi = (DECi_max - DECi_min)*rand + DECi_min;
    vi = [vi_bound*(2*rand - 1), vi_bound*(2*rand - 1), vi_bound*(2*rand - 1)];
    
% Else if start in orbit
elseif startBody == 2
    ai = (ai_max - ai_min)*rand + ai_min;
    AoPi = (AoPi_max - AoPi_min)*rand + AoPi_min;                                % Argument of periapsis omega in radians
    ei = (ei_max - ei_min)*rand + ei_min;                                                 % Eccentricity
    inci = (inci_max - inci_min)*rand + inci_min;                                       % Inclination in radians
    RAANi = (RAANi_max - RAANi_min)*rand + RAANi_min;                     % Right ascension of ascending node Omega in radians
    TAi = (TAi_max - TAi_min)*rand + TAi_min;

    % Else if start at perturbing body
else
    viRel = [viRel_bound*(2*rand - 1), viRel_bound*(2*rand - 1), viRel_bound*(2*rand - 1)];
end

% Final phase
if endBody == 1
    RAf = (RAf_max - RAf_min)*rand + RAf_min;
    DECf = (DECf_max - DECf_min)*rand + DECf_min;
    vf = [vf_bound*(2*rand - 1), vf_bound*(2*rand - 1), vf_bound*(2*rand - 1)];
elseif endBody == 2
    af = (af_max - af_min)*rand + af_min;
    AoPf = (AoPf_max - AoPf_min)*rand + AoPf_min;                                % Argument of periapsis omega in radians
    ef = (ef_max - ef_min)*rand + ef_min;                                                 % Eccentricity
    incf = (incf_max - incf_min)*rand + incf_min;                                       % Inclination in radians
    RAANf = (RAANf_max - RAANf_min)*rand + RAANf_min;                     % Right ascension of ascending node Omega in radians
    TAf = (TAf_max - TAf_min)*rand + TAf_min;
else
    vfRel = [vfRel_bound*(2*rand - 1), vfRel_bound*(2*rand - 1), vfRel_bound*(2*rand - 1)];
end

% Every phase
dt_all = zeros(1, Np);
vinfi_all = zeros(3,N_flybys);
vinff_all = zeros(3,N_flybys);
r_free = zeros(3,(Np - 1 - N_flybys));
v_free = zeros(3,(Np - 1 - N_flybys));
% Thrust phases
mf_all = zeros(1, N_thrust);
u_all = zeros(3*N_thrust, N);

for i = 1 : Np-1
    dt_all(i) = (dt_max(i) - dt_min(i))*rand + dt_min(i);
end

% For flyby phases
for i = 1 : N_flybys
    vinfi_all(:,i) = [vinf_bound(i)*(2*rand - 1); vinf_bound(i)*(2*rand - 1); vinf_bound(i)*(2*rand - 1)];
    vinff_all(:,i) = [vinf_bound(i)*(2*rand - 1); vinf_bound(i)*(2*rand - 1); vinf_bound(i)*(2*rand - 1)];
end

% For non-flyby phases
for i = 1 : (Np-1-N_flybys)
    r_free(:,i) = [(x_free_max - x_free_min)*rand + x_free_min; (y_free_max - y_free_min)*rand + y_free_min; ...
        (z_free_max - z_free_min)*rand + z_free_min];
    v_free(:,i) = [(vx_free_max - vx_free_min)*rand + vx_free_min; (vy_free_max - vy_free_min)*rand + vy_free_min; ...
        (vz_free_max - vz_free_min)*rand + vz_free_min];
end

% Thrust arcs
for i = 1:N_thrust
    mf_all(i) = (mf_max(i) - mf_min(i))*rand + mf_min(i);
    u_all((1:3)+(i-1)*3,:) = 2*rand(3, N) - 1;
end

% Stuff that all phases have
dt_all(Np) = (dt_max(N_thrust) - dt_min(N_thrust))*rand + dt_min(N_thrust); % THIS DIFFERS IN MBH 'dt_all(Np) = (dt_max(i) - dt_min(i))*rand + dt_min(i);'
% dt_all(Np) = (dt_max(Np) - dt_min(N_thrust))*rand + dt_min(N_thrust); % THIS DIFFERS IN MBH 'dt_all(Np) = (dt_max(i) - dt_min(i))*rand + dt_min(i);'

%% Convert units

% Velocity unit from distance/time
VU = DU/TU;

m0 = m0/MU;

mu_central = mu_central/((DU^3)/(TU^2));        % Sun gravitational parameter from km^3/s^2 to DU_S1^3/TU_S1^2 
mu_perts = mu_perts/((DU^3)/(TU^2));        % Sun gravitational parameter from km^3/s^2 to DU_S1^3/TU_S1^2 

% Put flyby radii in order of the GAs in the trajectory
r_flybys = zeros(1,N_flybys);
if N_flybys > 0
    for i = 1:N_flybys
        r_flybys(i) = r_perts(flybySequence(i));
    end
else
    r_flybys = 0;
end

% Convert constants used into correct units. Moon radius is already in DU
% so don't need to convert that
Tmax = Tmax/((MU*DU)/(TU^2)); % Max thrust from kN to MU_S1*DU_S1/TU_S1
Isp = Isp/TU;                                                % Specific impulse from s to TU_S1
g0 = g0/(DU/(TU^2));                           % Acceleration due to gravity at Earth's surface from km/s^s to DU_S1/TU_S1^2
r_central = r_central/DU;
r_flybys = r_flybys/DU;
SOI = SOI/DU;
h_mins = h_mins/DU;
h_maxs = h_maxs/DU;

%% Optimise phase by phase
A = [];               % A.x <= b
b = [];
Aeq = [];           % Aeq.x = beq
beq = [];
options = optimoptions('fmincon','Algorithm','sqp','ConstraintTolerance', NLP_feas_tol, 'OptimalityTolerance', ...
    NLP_tol, 'StepTolerance', NLP_steptol, 'MaxFunctionEvaluations', NLP_iter_max,'Display', displayFmincon, 'UseParallel', use_parallel);

NpCurrent = 1;

consts = [Np, N, N_perts, N_ephem, N_thrust, N_flybys, NpCurrent, mu_central, r_central, SOI, m0, N_act, D, Tmax, ...
    Isp, g0, DU, TU, MU, whichThrust, flybySequence, startBody, endBody, whichFlyby, mu_perts, r_flybys, h_mins, h_maxs, et, stateLong];

noDecision = 6 + Np + N_thrust + 3*Np*N*N_thrust + 6*(Np-1);
x_store = zeros(Np, noDecision);
optim_store = zeros(Np, noDecision);
viol_store = zeros(Np);
lb_store = zeros(Np, noDecision);
ub_store = zeros(Np, noDecision);
objInd = zeros(1,Np);
phaseSizes = zeros(1,Np);

fprintf("Initial Optimisation.... Phase(%u/%u)\n", 1, Np);

% First phase if starts on SOI
if startBody == 1
    x = [t0, RAi, DECi, vi, dt_all(1)];
    lb = [t0_min_et, RAi_min, DECi_min, -vi_bound, -vi_bound, -vi_bound, dt_min(1)];
    ub = [t0_max_et, RAi_max, DECi_max, vi_bound, vi_bound, vi_bound, dt_max(1)];
    indLastDt = 7;
    
% If start in orbit
elseif startBody == 2
    x = [t0, ai, AoPi, ei, inci, RAANi, TAi, dt_all(1)];
    lb = [t0_min_et, ai_min, AoPi_min, ei_min, inci_min, RAANi_min, TAi_min, dt_min(1)];
    ub = [t0_max_et, ai_max, AoPi_max, ei_max, inci_max, RAANi_max, TAi_max, dt_max(1)];
    indLastDt = 8;
    
% If start at perturbing body, ignored SOI
else
    x = [t0, viRel, dt_all(1)];
    lb = [t0_min_et, -viRel_bound, -viRel_bound, -viRel_bound, dt_min(1)];
    ub = [t0_max_et, viRel_bound, viRel_bound, viRel_bound, dt_max(1)];
    indLastDt = 5;
end

% If first phase ends in flyby
if whichFlyby(1) == 1
    indNow = length(x);
    x = [x, vinfi_all(:,1).', vinff_all(:,1).'];
    lb = [lb, -vinf_bound(1), -vinf_bound(1), -vinf_bound(1), -vinf_bound(1), -vinf_bound(1), -vinf_bound(1)];
    ub = [ub vinf_bound(1), vinf_bound(1), vinf_bound(1), vinf_bound(1), vinf_bound(1), vinf_bound(1)];
    if MBH_noLoops>1
        % Need index of velocity for MBH
        ind_vinfi(:,1) = indNow+1:indNow+3;
        ind_vinff(:,1) = indNow+4:indNow+6;
    end
% Else the first phase ends in a point in free space
else
    x = [x, r_free(:,1).', v_free(:,1).'];
    lb = [lb, x_free_min(1), y_free_min(1), z_free_min(1), vx_free_min(1), vy_free_min(1), vz_free_min(1)];
    ub = [ub, x_free_max(1), y_free_max(1), z_free_max(1), vx_free_max(1), vy_free_max(1), vz_free_max(1)];
end

thrustInd = 1;
sizeX = length(x);
dtIndex(1) = indLastDt;

% If first phase is a thrust arc
if whichThrust(1) == 1
    x = [x, mf_all(1), u_all(1,:), u_all(2,:), u_all(3,:)];
    lb = [lb, mf_min(1), u_min(1,:), u_min(2,:), u_min(3,:)];
    ub = [ub, mf_max(1), u_max(1,:), u_max(2,:), u_max(3,:)];
    thrustInd = 2;
end

% Supply index of last final mass, or tell optimiser to optimise the
% minimum flight time if no thrust arcs
if whichThrust(1) == 1

    indLastMf = sizeX + 1;
    objInd(1) = indLastMf;
    [optimised,~,~,output] = fmincon(@(x)obj_lofiSF(x,indLastMf),x,A,b,Aeq,beq,lb,ub, @(x)con_lofiSF(x,consts), options);
    phaseSizes(1) = length(optimised);
else
    objInd(1) = indLastDt;
    [optimised,~,~,output] = fmincon(@(x)obj_lofiSF_coast(x,indLastDt),x,A,b,Aeq,beq,lb,ub, @(x)con_lofiSF(x,consts), options);
    phaseSizes(1) = length(optimised);
end

%% Intermediate phases

NpCurrent = 2;
for i = 2:Np-1
    fprintf("Initial Optimisation.... Phase(%u/%u)\n", i, Np);

    sizeOptim = length(x);
    consts(7) = NpCurrent;
    x = [optimised, dt_all(i)];
    lb = [lb, dt_min(i)];
    ub = [ub, dt_max(i)];
    indLastDt = sizeOptim+1;
    dtIndex(i) = indLastDt;
    % If this phase ends in a flyby
    if any(i == whichFlyby)
        x = [x, vinfi_all(:,i).', vinff_all(:,i).'];
        lb = [lb, -vinf_bound(i), -vinf_bound(i), -vinf_bound(i), -vinf_bound(i), -vinf_bound(i), -vinf_bound(i)];
        ub = [ub, vinf_bound(i), vinf_bound(i), vinf_bound(i), vinf_bound(i), vinf_bound(i), vinf_bound(i)];

    % Else, this phase ends at a point in free space
    else
        x = [x, r_free(:,i).', v_free(:,i).'];
        lb = [lb, x_free_min(i), y_free_min(i), z_free_min(i), vx_free_min(i), vy_free_min(i), vz_free_min(i)];
        ub = [ub, x_free_max(i), y_free_max(i), z_free_max(i), vx_free_max(i), vy_free_max(i), vz_free_max(i)];
    end

    % If phase is a thrust arc
    if any(i==whichThrust)
        x = [x, mf_all(thrustInd), u_all(1+(thrustInd-1)*3,:), u_all(2+(thrustInd-1)*3,:), u_all(3+(thrustInd-1)*3,:)];
        lb = [lb, mf_min(thrustInd), u_min(1+(thrustInd-1)*3,:), u_min(2+(thrustInd-1)*3,:), u_min(3+(thrustInd-1)*3,:)];
        ub = [ub, mf_max(thrustInd), u_max(1+(thrustInd-1)*3,:), u_max(2+(thrustInd-1)*3,:), u_max(3+(thrustInd-1)*3,:)];
        thrustInd = thrustInd+1;
        
        indLastMf = sizeOptim + 8;
        objInd(i) = indLastMf;
        [optimised,~,~,output] = fmincon(@(x)obj_lofiSF(x,indLastMf),x,A,b,Aeq,beq,lb,ub, @(x)con_lofiSF(x,consts), options);
        phaseSizes(i) = length(optimised);
    else
        % If there has been a thrust arc in the past, optimise for its
        % final mass
        if any(i <= whichThrust)
            objInd(i) = indLastMf;
            [optimised,~,~,output] = fmincon(@(x)obj_lofiSF(x,indLastMf),x,A,b,Aeq,beq,lb,ub, @(x)con_lofiSF(x,consts), options);
            phaseSizes(i) = length(optimised);
        else
            indLast = sizeOptim + 1;
            objInd(i) = indLastDt;
            [optimised,~,~,output] = fmincon(@(x)obj_lofiSF_coast(x,indLastDt),x,A,b,Aeq,beq,lb,ub, @(x)con_lofiSF(x,consts), options);
            phaseSizes(i) = length(optimised);
        end
    end
    
    NpCurrent = NpCurrent+1;
end

%% Final phase

fprintf("Initial Optimisation.... Phase(%u/%u)\n", Np, Np);
consts(7) = Np;

sizeOptim = length(x);
indLastDt = sizeOptim+1;

% If ends on SOI
if endBody == 1
    x = [optimised, dt_all(end), RAf, DECf, vf];
    lb = [lb, dt_min(end), RAf_min, DECf_min, -vf_bound,-vf_bound,-vf_bound];
    ub = [ub, dt_max(end), RAf_max, DECf_max, vf_bound,vf_bound,vf_bound];

    % Or, if trajectory ends in orbit around central body
elseif endBody == 2
    x = [optimised, dt_all(end), af, AoPf, ef, incf, RAANf, TAf];
    lb = [lb, dt_min(end), af_min, AoPf_min, ef_min, incf_min, RAANf_min, TAf_min];
    ub = [ub, dt_max(end), af_max, AoPf_max, ef_max, incf_max, RAANf_max, TAf_max];

    % Or, if trajectory ends at location of perturbing body
else
    x = [optimised, dt_all(end), vfRel];
    lb = [lb, dt_min(end), -vfRel_bound,-vfRel_bound,-vfRel_bound];
    ub = [ub, dt_max(end), vfRel_bound,vfRel_bound,vfRel_bound];
end

dtIndex(end) = indLastDt;

% If phase is a thrust arc
if whichThrust(end) == Np
    x = [x, mf_all(end), u_all(1+(N_thrust-1)*3,:), u_all(2+(N_thrust-1)*3,:), u_all(3+(N_thrust-1)*3,:)];
    lb = [lb, mf_min(end), u_min(1+(N_thrust-1)*3,:), u_min(2+(N_thrust-1)*3,:), u_min(3+(N_thrust-1)*3,:)];
    ub = [ub, mf_max(end), u_max(1+(N_thrust-1)*3,:), u_max(2+(N_thrust-1)*3,:), u_max(3+(N_thrust-1)*3,:)];
    
    indLastMf = sizeOptim + 5;
    objInd(Np) = indLastMf;
    [optimised,m_optim, exitflag,output] = fmincon(@(x)obj_lofiSF(x,indLastMf),x,A,b,Aeq,beq,lb,ub, @(x)con_lofiSF(x,consts), options);
    phaseSizes(Np) = length(optimised);
else
    % If there has been a thrust arc in the past, optimise for its
    % final mass
    if N_thrust ~= 0
        objInd(Np) = indLastMf;
        [optimised,m_optim, exitflag,output] = fmincon(@(x)obj_lofiSF(x,indLastMf),x,A,b,Aeq,beq,lb,ub, @(x)con_lofiSF(x,consts), options);
        phaseSizes(Np) = length(optimised);
        else
        indLast = sizeOptim + 1;
        objInd(Np) = indLastDt;
        [optimised,m_optim, exitflag,output] = fmincon(@(x)obj_lofiSF_coast(x,indLastDt),x,A,b,Aeq,beq,lb,ub, @(x)con_lofiSF(x,consts), options);
        phaseSizes(Np) = length(optimised);
    end
end

fprintf("Initial Optimisation.... Done\n");

optim_archive = zeros(MBH_noLoops, length(optimised)); % For optimised decision variables at each iteration
m_archive = zeros(MBH_noLoops, 1);                    % For final mass value with each iteration
violation_archive = zeros(MBH_noLoops, 1);          % For the violation indicators with each iteration
times = zeros(MBH_noLoops, 1);

% Save to archive if the optimised trajectory is feasible
optim_archive(1,:) = optimised;         % Save all optimised decision variables
m_archive(1) = m_optim;               % Save final mass from that optimised run
violation_archive(1) = output.constrviolation;
best = optimised;

if MBH_noLoops>1
    disp('Starting Monotonic Basin Hopping');
    best = optim_archive(1,:);    %% MBH setup

    % First phase
    if startBody == 1
        sigmas = [s_t0, s_angles, s_angles, s_v, s_v, s_v];
        
        % Or, if trajectory starts in orbit around central body
    elseif startBody == 2
        sigmas = [s_t0, s_a, s_angles, s_e, s_angles, s_angles, s_angles];
        
        % Or, if trajectory starts at location of perturbing body
    else
        sigmas = [s_t0, s_v, s_v, s_v];
    end

    % If first phase ends in flyby
    if whichFlyby(1) == 1
        sigmas = [sigmas, s_dt, s_v,s_v,s_v, s_v,s_v,s_v];
    else
        sigmas = [sigmas, s_dt, s_r, s_r, s_r, s_v, s_v, s_v];
    end

    if whichThrust(1) == 1
        sigmas = [sigmas, s_m, s_u*ones(1,3*N)];
    end

    % Intermediate phases
    for i = 2:Np-1
        sigmas = [sigmas, s_dt];
        
        % If ends in flyby
        if any(i == whichFlyby)
            sigmas = [sigmas, s_v,s_v,s_v, s_v,s_v,s_v];
        else
            sigmas = [sigmas, s_r, s_r, s_r, s_v, s_v, s_v];
        end

        % If phase is a thrust arc
        if any(i==whichThrust)
            sigmas = [sigmas, s_m, s_u*ones(1,N*3)];
        end
    end

    % Final phase
    if endBody == 1
        sigmas = [sigmas, s_dt, s_angles, s_angles, s_v,s_v,s_v];

        % Or, if trajectory ends in orbit around central body
    elseif endBody == 2
        sigmas = [sigmas, s_dt, s_a, s_angles, s_e, s_angles, s_angles, s_angles];
        
        % Or, if trajectory ends at location of perturbing body
    else
        sigmas = [sigmas, s_dt, s_v,s_v,s_v];
    end

    if whichThrust(end) == Np
        sigmas = [sigmas, s_m, s_u*ones(1,3*N)];
    end

    % times(1) = toc;

    % For resolving, know where optimising over dt and where mf

%    [minViolation, iMinViolation] = min(nonzeros(violation_archive));
    minViolation = violation_archive(1);
    
    % Disable fmincon parallelisation here to avoid nesting parpool.
    options.UseParallel=false;

    % Define partial function to simplify loop.
    f = @(k, best)basinhop(k, best, sigmas, MBH_tail, MBH_theta, rho_hop, t0Hop, dtHop, dtIndex, lb, ub, ind_vinfi, ind_vinff, s_v, N_flybys, phaseSizes, objInd, whichThrust, consts, options, A, b, Aeq, beq, Np);

    %% MBH loop
    if use_parallel
        % Delete lockfile if exists.
        jobStorageLocation = pc.JobStorageLocation;

        if exist(fullfile(jobStorageLocation,'.lock'), 'file')==2
            delete(fullfile(jobStorageLocation, '.lock'));
        end
        asyncsave(jobStorageLocation, 'best.mat', struct('best', best, 'minViolation', minViolation));
        parfor k = 2:MBH_noLoops-1
            s = load(fullfile(jobStorageLocation, 'best.mat'), 'best');
            best = s.best;
            
            [optim_archive(k, :), m_archive(k), violation_archive(k)] = f(k, best);
    
            % Re-check file, incase has been updated since. 
            s = load(fullfile(jobStorageLocation, 'best.mat'), 'minViolation');
            minViolation = s.minViolation;

            if violation_archive(k) < minViolation
                asyncsave(jobStorageLocation, 'best.mat', struct('best', optim_archive(k, :), 'minViolation', violation_archive(k)));
                fprintf("New best value found, %5.5G vs %5.5G\n", violation_archive(k), minViolation);
            else
                fprintf("Best value unchanged, %5.5G vs %5.5G\n", violation_archive(k), minViolation);
            end
        end
        % Serial run to make sure final run of workers arnt left out.
        [optim_archive(MBH_noLoops, :), m_archive(MBH_noLoops), violation_archive(MBH_noLoops)] = f(MBH_noLoops, best);
    else
        for k = 2:MBH_noLoops
            [optim_archive(k, :), m_archive(k), violation_archive(k)] = f(k, best);
            if violation_archive(k) < minViolation
                fprintf("New best value found, %5.5G vs %5.5G\n", violation_archive(k), minViolation);
                best = optim_archive(k, :);
                minViolation = violation_archive(k);
            else
                fprintf("Best value unchanged, %5.5G vs %5.5G\n", violation_archive(k), minViolation);
            end
        end
    end
end

output = struct('optim_archive', optim_archive, 'm_archive', m_archive, 'violation_archive', violation_archive, ...
    'consts', consts, 'lb', lb, 'ub', ub);

end