% This script collects the statistics and builds the plot for errors in
% Control System in case of the "invisible" request signals from the
% position of the receiver. The reference aircraft position errors are
% varied

addpath ..\src
addpath additional_functions

%% initialisation

TOSAVE = 0; %save intermediate results and plot as a figure or not
iterations = 1e1;
N = 1;

SNR = [15 10 5];
t_lms_req_SNR = [1.2*1e-8 2.2*1e-8 3.6*1e-8]; %LMS of TOA errors for each SNR value for the request signals (are known a priori)
t_lms_ans_SNR = [5.9*1e-9 1.3*1e-8 2.2*1e-8]; %LMS of TOA errors for each SNR value for the answer signals (are known a priori)

Ts = 5e-3; %request signals period, s
Ta = 10; %antenna rotation period, s
v = 1/Ta; %antenna rotation speed, rotations/s
ssr_parameters = struct('PRI', Ts, 'speed', v);
v_lms = 0.01; %LMS of SSR antenna rotation speed errors
Ts_lms = 1e-8; %LMS of errors in the request signals period
Ts_lms = sqrt(Ts_lms^2*12); % uniform distribution correction
ssr_errors = struct('PRI_error', Ts_lms, 'speed_error', v_lms);

% SSR, receiver and reference aircraft positions
receiver = [0 0 0]; %receiver position, [x, y, z], m
ssr = [30e3 0 0]; %SSR position, [x, y, z], m
reference_aircraft = [20e3 10e3 10e3]; %reference aircraft position, [x, y, z], m

% aircraft positions
L_real = 100e3; %all concidered "path difference": the difference between "SSR-Aircraft-Receiver" and "SSR-Receiver" distances, m
h_aircraft = 10e3; %all concidered height of aircraft, m
phi_real = [0 : 2*pi/9 : 2*pi]; %all concidered "phi" values (phi is the angle between "SSR-Aircraft" and "SSR-Receiver")

Pssr = 2000; %the minimum power of request signals in Watts
Rtransp = 125; %the minimum power of answer signals in Watts

reference_aircraft_LMS = [5 10 100]; % LMS of all considered reference aircraft position errors

if TOSAVE
    % prefix for file names (if it is necessary to save anything)
    datestring = datetime('now','Format','yyyy-MM-dd_HH-mm');
    
    if ~exist('results', 'dir')
        mkdir('results')
    end
    cd results
    mkdir(char(datestring))
    cd ../
    
    prefix = ['results\' char(datestring) '\'];

    % save all initial values if it is necessary
    save([prefix 'init.mat'], 'iterations', 'N', 'SNR', 't_lms_req_SNR', ....
            't_lms_ans_SNR', 'ssr_errors', 'ssr_parameters', 'reference_aircraft_LMS');
end

%% preparations before the main cycle

% We concider all phi_real for each L_real. In this case all aircraft
% positions lay on the ellipse. The parameters of this ellipse are
% calculated below

r_receiver_ssr = get_distance(receiver, ssr); %the distance between the receiver and the SSR
a = (L_real - r_receiver_ssr)/2 + (r_receiver_ssr/2);
b = sqrt((L_real/2).^2 - (r_receiver_ssr/2).^2);

% So we can get the aircraft coordinates
x_aircraft = a.*cos(phi_real);
y_aircraft = b.*sin(phi_real);
x_aircraft = x_aircraft + ssr(1)/2;

% Structure for TOA errors
t_errors = struct('requests_on_receiver', 0, 'requests_on_aircraft', 0, 'answers', 0);

%% collect statistics and find errors for the Control System in different cases

disp(['----- Start main cycle -----']);
for current_ref_aircraft_LMS = reference_aircraft_LMS
    disp(' ')
    disp('--------------------------------------------------------')
    disp(['----- Control System, reference aircraft LMS is ' num2str(current_ref_aircraft_LMS) ' -----']);
    disp('--------------------------------------------------------')
    all_errors_lms = cell(1, length(SNR));
    all_errors_mean = cell(1, length(SNR));
    for isnr = 1 : length(SNR)
        disp(' ')
        disp(['----- SNR is ' num2str(SNR(isnr)) ' -----']);
        r_error_lms = zeros(1, length(x_aircraft));
        r_error_mean = zeros(1, length(x_aircraft));
        for iPhi = 1 : length(x_aircraft)
            disp(['----- phi = ' num2str(phi_real(iPhi)) ' radians -----']);
            
            aircraft = [x_aircraft(iPhi) y_aircraft(iPhi) h_aircraft(1)];
            b = get_distance(receiver, ssr);

            error_R = []; % initialize vector which contains errors in aircraft position determination (the length is unknown a priori)

            % Calculation of all TOA errors according to distances, types of 
            % signals, signals' power and SNR
            hssr_rec_sq = SNR(isnr);
            hssr_rec_sq_times = 10^(hssr_rec_sq/10);
            hssr_aircraft_sq_times = hssr_rec_sq_times*ssr(1)^2/get_distance(ssr,aircraft)^2;
            hp_receiver_sq_times = (hssr_rec_sq_times*Rtransp/Pssr)*ssr(1)^2/get_distance(receiver,aircraft)^2;
            hssr_receiver_ssr_sq_times = hssr_rec_sq_times*ssr(1)^2/get_distance(ssr,receiver)^2;
            t_errors.requests_on_receiver = t_lms_req_SNR(isnr)*hssr_rec_sq_times/hssr_receiver_ssr_sq_times;
            t_errors.requests_on_aircraft = t_lms_req_SNR(isnr)*hssr_rec_sq_times/hssr_aircraft_sq_times;
            t_errors.answers = 2*t_lms_ans_SNR(isnr)*hssr_rec_sq_times/hp_receiver_sq_times;
            
            for i = 1 : iterations
                if ~rem(i,100)
                    disp(['iteration No ' num2str(i)]);
                    %disp(datetime)
                end
                [L, phi, ok] = L_phi_calc_invisible( ssr_parameters, ssr_errors, t_errors, receiver, ssr, aircraft, reference_aircraft, false, N, current_ref_aircraft_LMS );
                [R1, R2] = R1R2_function(L, get_distance(ssr, receiver), phi, aircraft(3));
                if isempty(R1)
                    continue;
                end
                [x, y, ok] = get_coordinates_from_distances(sqrt(R2^2 - aircraft(3)^2), sqrt(R1^2 - aircraft(3)^2), receiver(1:2), ssr(1:2), aircraft(1:2));
                if ((~ok) || isempty(x) || isempty(y))
                    continue;
                end

                error_R = [error_R; sqrt((aircraft(1) - x)^2 + (aircraft(2) - y)^2)];
            end
            r_error_lms(iPhi) = sqrt(mean((error_R - mean(error_R)).^2));
            r_error_mean(iPhi) = mean(error_R);
        end
        all_errors_lms{isnr} = r_error_lms;
        all_errors_mean{isnr} = r_error_mean;
    end
    if TOSAVE
        save([prefix 'all_errors_ref_aircraft_lms_' num2str(current_ref_aircraft_LMS) '.mat'], 'all_errors_lms', 'r_error_mean');
    end
end

%% plot

set(0,'DefaultAxesFontSize',11,'DefaultAxesFontName','Times New Roman');
figure();
plot(phi_real, all_errors_lms{1}, 'ko-');
hold on
plot(phi_real, all_errors_lms{2}, 'kx-');
plot(phi_real, all_errors_lms{3}, 'k*-');
grid on
xlabel('\gamma, radians');
ylabel('\sigma_{\deltaR}, m');
title(['LMS of aircraft position determination errors, \sigma_{\delta_{Tci}} = ' num2str(ssr_errors.PRI_error) ', \sigma_{\delta_{vi}} = ' num2str(ssr_errors.speed_error)]);
legend(['d^2_{SSR} = ' num2str(SNR(1)) ' dB'], ['d^2_{SSR} = ' num2str(SNR(2)) ' dB'], ['d^2_{SSR} = ' num2str(SNR(3)) ' dB'])

% Correct axis
highest_point = max([max(all_errors_lms{1}), max(all_errors_lms{2}), max(all_errors_lms{3})]);
axis([0,2*pi,0,highest_point*1.05])