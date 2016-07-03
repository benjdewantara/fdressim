%% Setup and simulate the problem
% Set up a grid of 5 x 1 x 1 blocks, where each block has dimensions of
% 1000 x 1000 x 75 feet.
dims = [5, 1, 1];
physDims = dims.*[1000, 1000, 75]*ft;
G = cartGrid(dims, physDims);
% Add geometry information to grid
G = computeGeometry(G);
% Create rock with constant permeability and porosity
rock = makeRock(G, 15*milli*darcy, 0.18);

% Load modules for solving compressible flow
mrstModule add ad-core ad-blackoil ad-props

% Set up three phase fluid - we will only use the second fluid phase (oil),
% set to a viscosity of 10cP.
fluid = initSimpleADIFluid('mu', [1, 10, 1]*centi*poise);
% Define fluid compressibility and reference pressure, for a linear inverse
% formation volume factor, i.e. b = 1/B, which is the convention in MRST.
c_l = 3.5e-6/psia();
pRef = 6000*psia();
fluid.bO = @(p)(1 + c_l.*(p-pRef));

% Add a rate controlled producer at cell number 4
W = [];
W = addWell(W, G, rock, 4, 'Type', 'rate', 'val', -150*stb/day, 'comp_i', [0, 1]);

% Plot the grid and the well
figure(1); clf
plotGrid(G);
plotWell(G, W);
view(-20, 20);
axis tight

% Set initial state to 6000 PSI and completely oil filled
state0 = initResSol(G, pRef, [0, 1]);

% Total simulation time of about 1 year
T = 360*day;
% Divide into time-steps of length 10 days
dt = 10*day;
% Create a vector of all timesteps
dT = repmat(dt, ceil(T/dt), 1);
% Combine well and timesteps into schedule
schedule = simpleSchedule(dT, 'W', W);
% Set up two-phase flow model for compressible problem
model = TwoPhaseOilWaterModel(G, rock, fluid);
%% Simulate the problem
[ws, states] = simulateScheduleAD(state0, model, schedule);
%% Output to screen the different cell pressures
for i = 1:numel(states)
    fprintf('%d days ', sum(schedule.step.val(1:i))/day);
    fprintf('%4.2f ', states{i}.pressure/psia);
    fprintf('\n')
end
%% Plot well bottom hole pressure
figure(1), clf
time = cumsum(dT);
bhp = getWellOutput(ws, 'bhp');
plot(time/day, bhp/psia)
xlabel('Time [days]')
ylabel('Bottom hole pressure [psi]')
%% Plot pressure development
x = G.cells.centroids(:, 1);
for i = 1:numel(states)
    figure(1); clf;
    plot(x/ft, states{i}.pressure/psia);
    ylim([4000, 6000]);
    title(formatTimeRange(time(i)))
    ylabel('Reservoir pressure [psi]')
    xlabel('X [feet]');
    pause(0.1)
end
