close all
clear all
clc

set(groot,'defaultLineMarkerSize',20);
set(groot,'defaultLineLineWidth',1.3);

addpath('package_2023')
addpath("Lab 4/F160_GrowthCurve")
addpath("Lab 4/F160_Phase/")
addpath("Lab 4/F160_VelProf/")

%% Constants
Patm = 100790.00; % [Pa]
T    = 19.86;     % [degrees Celsius]
[rho, nu, mu] = suthlaw(Patm, T);

%% Read velocity profile data and construct arrays
fileID = fopen("F_160_evolution_Stats.txt");
read_data = textscan(fileID, "%f %f %f %f %f %f %f %f %f", 'HeaderLines', 1);

X = read_data{2};
X_positions = unique(X);

% We have two x positions, X = 180 and X = 480, therefore two profiles
% (Note: data{3} corresponds to Y values and data{5} U velocity values)
Y_180  = read_data{3}(X == X(1)) ; Y_480  = read_data{3}(X == X(end));
Uy_180 = read_data{5}(X == X(1)) ; Uy_480 = read_data{5}(X == X(end));

% Flip the vectors from descending to ascending order
Y_180 = flip(Y_180);
Y_480 = flip(Y_480);
Uy_180 = flip(Uy_180);
Uy_480 = flip(Uy_480);

%% Figure 2 - Scaled velocity profiles + compute delta_1
Y_data = [Y_180 Y_480];
Uy_data = [Uy_180 Uy_480];

Ue = zeros(1,2);
delta_1 = zeros(1,2);
ny = zeros(1,2);
x_zero = zeros(1,2);

labels = ["Branch I (180 mm)", "Branch II (480 mm)"];

f2 = figure('Name','Scaled Velocity Plots');
f2.Position(3:4) = [1200,1200];

for i = 1:2

    Y_raw = Y_data(:,i); 
    Uy_raw = Uy_data(:,i);

    % Compute wall coordinate 
    [Ywall, ny(i)] = TS_LAB_JF_P1(Y_raw, Uy_raw);

    % Adjust the Y coordinates and add 0 boundary condition
    Y  = [0; Y_raw - Ywall];
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
    fprintf("Ywall        : %f mm\n", Ywall)
    fprintf("ny           : %d\n", ny(i))
    fprintf("Ue           : %f m/s\n", Ue(i))
    fprintf("delta_1      : %f mm\n", delta_1(i))
    fprintf("Estimated x_0: %f mm\n", x_zero(i)*1000)

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

%% Read velocity profile data and construct arrays
fileID = fopen("F160_GC_Stats.txt");
read_data = textscan(fileID, "%f %f %f %f %f %f %f %f %f", 'HeaderLines', 1);

load('/scratch/baconnet/courses/FSG2221-stability/Lab/package_2023/LinTheory_F160_2023.mat')

Uinf = 5.941282; % As report in F160_GC_Coefficients.txt
X = read_data{2};

Amplitudes = read_data{7};
N = log(Amplitudes/cBI);

delta = sqrt((X - mean(x_zero))/1000 * nu/Uinf);
Re_delta = Uinf*delta/nu;

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


