close all
clear all
clc

set(groot,'defaultLineMarkerSize',20);
set(groot,'defaultLineLineWidth',1.3);

addpath('package_2023')
addpath("Lab 4/F160_GrowthCurve")
addpath("Lab 4/F160_Phase/")
addpath("Lab 4/F160_VelProf/")

load('/scratch/baconnet/courses/FSG2221-stability/Lab/package_2023/LinTheory_F160_2023.mat')

%% Constants
Patm = 100790.00; % [Pa]
T    = 19.86;     % [degrees Celsius]
[rho, nu, mu] = suthlaw(Patm, T);

%% Figure 2 - Scaled velocity profiles + compute delta_1

% Read velocity profile data and construct arrays
fileID = fopen("F_160_evolution_Stats.txt");
read_data = textscan(fileID, "%f %f %f %f %f %f %f %f %f", 'HeaderLines', 1);

X = read_data{2};
X_positions = unique(X);

% We have two x positions, X = 180 and X = 480, therefore two profiles
% (Note: data{3} corresponds to Y values and data{5} U velocity values)
Y_180  = read_data{3}(X == X(1)) ; Y_480  = read_data{3}(X == X(end));
Uy_180 = read_data{5}(X == X(1)) ; Uy_480 = read_data{5}(X == X(end));

% Get amplitude at branch I and II for later
A_branch_I  = read_data{7}(X == X(1))  ; A_branch_I  = flip(A_branch_I); 
A_branch_II = read_data{7}(X == X(end)); A_branch_II = flip(A_branch_II);

% Get phases at branch I and II for later
phi_branch_I  = read_data{8}(X == X(1))  ; phi_branch_I  = flip(phi_branch_I); 
phi_branch_II = read_data{8}(X == X(end)); phi_branch_II = flip(phi_branch_II);

% Flip the vectors from descending to ascending order
Y_180 = flip(Y_180);
Y_480 = flip(Y_480);
Uy_180 = flip(Uy_180);
Uy_480 = flip(Uy_480);

Y_data = [Y_180 Y_480];
Uy_data = [Uy_180 Uy_480];

Ue = zeros(1,2);
delta_1 = zeros(1,2);
ny = zeros(1,2);
x_zero = zeros(1,2);
Ywall = zeros(1,2);

labels = ["Branch I", "Branch II"];

f2 = figure('Name','Scaled Velocity Plots');
f2.Position(3:4) = [1200,1200];

% Loop through each branch at a time
for i = 1:2

    Y_raw = Y_data(:,i); 
    Uy_raw = Uy_data(:,i);

    % Compute wall coordinate 
    [Ywall(i), ny(i)] = TS_LAB_JF_P1(Y_raw, Uy_raw);

    % Adjust the Y coordinates and add 0 boundary condition
    Y  = [0; Y_raw - Ywall(i)];
    Uy = [0; Uy_raw];

    % Compute normalization variables Ue and delta_1 (displacement
    % thickness)    
    Ue(i) = mean(Uy(end-5:end));
    delta_1(i) = trapz(Y, 1 - Uy/Ue(i));

    % Compute the estimated x_0 position (virtual origin)
    x_zero(i) = X_positions(i)*0.001 - (delta_1(i)*0.001/1.7208)^2 * Ue(i)/nu;
    x_zero(i) = x_zero(i)*1000; % Convert to mm

    fprintf(labels(i))
    fprintf("\n-------------------------------------\n")
    fprintf("Ywall            : %f mm\n", Ywall(i))
    fprintf("ny               : %d\n", ny(i))
    fprintf("Ue               : %f m/s\n", Ue(i))
    fprintf("delta_1          : %f mm\n", delta_1(i))
    fprintf("Estimated x0     : %f mm\n", x_zero(i))
    fprintf("Absolute distance: %f mm\n", X_positions(i))
    fprintf("Distance from LE : %f mm\n", X_positions(i) - x_zero(i))

    fprintf("\n")
    figure(f2)
    hold on
    plot(Uy./Ue(i), Y./delta_1(i), "--+")

end

figure(f2);
title("Scaled velocity profiles")
legend(labels(1), labels(2))
xlabel("u/Ue")
ylabel("y/\delta_1")
set(gca,'FontSize',20)
grid("on")
set(findall(gcf,'type','text'),'FontSize',20)

%% Figure 1 - Wave Growth Rate
% Read velocity profile data and construct arrays
fileID = fopen("F160_GC_Stats.txt");
read_data = textscan(fileID, "%f %f %f %f %f %f %f %f %f", 'HeaderLines', 1);

Uinf = 5.941282; % As report in F160_GC_Coefficients.txt
X = read_data{2};

Amplitudes = read_data{7};
N = log(Amplitudes/min(Amplitudes));

delta = sqrt((X - mean(x_zero))/1000 * nu/Uinf);
Re_delta = Uinf*delta/nu;

% Plot
f1 = figure('Name','Wave Growth Rate');
f1.Position(3:4) = [1200,1200];
plot(Re_delta, N, 'o')
hold on
plot(Red, Nfac)

legend("Experimental", "OS Theory")
title("Wave Growth Rate")
xlabel("Re_{\delta}")
ylabel("N = ln(A/A_0)")
set(gca,'FontSize',20)
grid("on")
set(findall(gcf,'type','text'),'FontSize',20)

%% Figures 3-a Measured wall normal amplitude and phase at branch I

A_branches = [A_branch_I A_branch_II];
phi_branches = [phi_branch_I phi_branch_II];
uB = [uBI uBII];
phiB = [phiBI phiBII];

for i=1:2

    Y = [0 ; Y_data(:,i) - Ywall(i)];

    [f,~,~,~,eta_fs] = FS_solver_JF(ny(i));
    b = eta_fs(end) - f(end);
    
    % Plot
    % f3a = figure('Name','Wall-Normal Amplitudes');
    % f3a.Position(3:4) = [1200,1200];
    % plot([0; A_branches(:,i)]./max(A_branches(:,i)), Y./delta_1(i), '--o')
    % hold on
    % plot(uB(:,i)/max(uB(:,i)), eta/b)
    % legend("Experimental (amplitude)", "Linear Theory (amplitude)")
    % title(labels(i))
    % ylabel("y/\delta_1")
    % xlabel("Normalized Amplitude")
    % set(gca,'FontSize',20)
    
    f3aphase = figure('Name','Wall-Normal Phases');
    f3aphase.Position(3:4) = [1200,1200];
    plot([0; phi_branches(:,i)/180], Y./delta_1(i), "--o")
    hold on
    plot(phiB(:,i)/3.1415926, eta./b)
    
    legend("Experimental (phase)", "Linear Theory (phase)")
    title("Wall-Normal Phase")
    ylabel("y/\delta_1")
    title(labels(i))
    xlabel("Phase")
    set(gca,'FontSize',20)
    grid("on")
    set(findall(gcf,'type','text'),'FontSize',20)

end

%% Figure 4 - Phase distribution along streamwise direction
fileID = fopen("F160_Phase_Stats.txt");
read_data = textscan(fileID, "%f %f %f %f %f %f %f %f", 'HeaderLines', 1);

phase = unwrap(read_data{8}*3.1415926/180);
X = read_data{2} - mean(x_zero);

% Linear regression
coefs = polyfit(X/1000, phase, 1);

Uinf = 5.950682; % From F160_Phase_Coefficients.txt
alpha_r = coefs(1);
lambda = 2*3.1415926/alpha_r;
freq_TS = Uinf^2 * 160 / (2*3.1415926 * nu * 1e6);

fprintf("Wave number       : %f m^(-1)\n", alpha_r)
fprintf("Wavelength        : %f m\n", lambda)
fprintf("Frequency         : %f Hz\n", freq_TS)
fprintf("Phase speed c     : %f m/s\n", lambda * freq_TS )
fprintf("Phase speed c/Uinf: %f \n", lambda * freq_TS / Uinf)

f4 = figure('Name','Phase along streamwise direction');
f4.Position(3:4) = [1200,1200];
plot(X, phase, "+")
hold on
plot(X, X*coefs(1)/1000 + coefs(2))

legend("Experimental", "Linear Regression")
ylabel("\phi (rad)")
xlabel("Distance from leading edge (mm)")
set(gca,'FontSize',20)
grid("on")
set(findall(gcf,'type','text'),'FontSize',20)
