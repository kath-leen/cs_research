% This script collects the statistics and builds the plot for errors in
% bistatic SSR and the remote control system for different aircraft
% positions (which is defined by the angle "phi")

addpath ..\src
addpath additional_functions

%% initialisation

TOSAVE = 0; %save intermediate results and plot as a figure or not
iterations = 1e1; %the number of iterations to collect the statistics
N = 5; %the number of SSR antenna rotations which are taken to collect the statistics for the control system

% initial verror values
SNR = [15 10 5]; %signal-to-noise ratio for analysis
t_lms_req_SNR = [1.2*1e-8 2.2*1e-8 3.6*1e-8]; %LMS of TOA errors for each SNR value for the request signals (are known a priori)
t_lms_ans_SNR = [5.9*1e-9 1.3*1e-8 2.2*1e-8]; %LMS of TOA errors for each SNR value for the answer signals (are known a priori)

% SSR parameters
Ts = 5e-3; %request signals period, s
Ta = 10; %antenna rotation period, s
v = 1/Ta; %antenna rotation speed, rotations/s
ssr_parameters = struct('PRI', Ts, 'speed', v);
v_lms = 0.01; %LMS of SSR antenna rotation speed errors
Ts_lms = 1e-8; %LMS of errors in the request signals period
Ts_lms = sqrt(Ts_lms^2*12); % uniform distribution correction
ssr_errors = struct('PRI_error', Ts_lms, 'speed_error', v_lms);

% SSR and receiver positions
receiver = [0 0 0]; %receiver position, [x, y, z], m
ssr = [30e3 0 0]; %SSR position, [x, y, z], m

% aircraft positions
L_real = 100e3; %all concidered "path difference": the difference between "SSR-Aircraft-Receiver" and "SSR-Receiver" distances, m
h_aircraft = 10e3; %all concidered height of aircraft, m
phi_real = [0 : 2*pi/9 : 2*pi]; %all concidered "phi" values (phi is the angle between "SSR-Aircraft" and "SSR-Receiver")

Pssr = 2000; %the minimum power of request signals in Watts
Rtransp = 125; %the minimum power of answer signals in Watts

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
    save([prefix 'init.mat'], 'iterations', 'N', 'SNR', 't_lms_req_SNR', 't_lms_ans_SNR', 'ssr_errors', 'ssr_parameters');
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

%% collect statistics and find errors for monostatic SSR different types of bistatic SSR and the Control System

disp(['----- Start main cycle -----']);
lms_bistatic_ssr = cell(1, length(SNR));
lms_bistatic_ssr_no_Ts = cell(1, length(SNR));
lms_control_system_N = cell(1, length(SNR));
lms_mssr = cell(1, length(SNR));
for isnr = 1 : length(SNR)
    disp(' ')
    disp(['----- SNR is ' num2str(SNR(isnr)) ' -----']);
    
    % Initialize vectors to collect statistics
    lms_BSSR_known_Ts = zeros(1, length(x_aircraft));
    lms_BSSR_calculated_Ts = zeros(1, length(x_aircraft));
    lms_MSSR = zeros(1, length(x_aircraft));
    lms_control_system = zeros(1, length(x_aircraft));
    
    % Varying aricraft position
    for iPhi = 1 : length(x_aircraft)
        disp(['----- phi = ' num2str(phi_real(iPhi)) ' radians -----']);
        aircraft = [x_aircraft(iPhi) y_aircraft(iPhi) h_aircraft(1)];
        b = get_distance(receiver, ssr);
        
        % Initialize vectors which contains errors in aircraft position
        % determination; the length is unknown a priori
        error_R_BSSR_known_Ts = [];
        error_R_BSSR_calculated_Ts = [];
        error_R_MSSR = [];
        error_R_control_system = [];

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
        
        % Collection of the statistics
        for i = 1 : iterations
            if ~rem(i,100)
                disp(['iteration No ' num2str(i)]);
                %disp(datetime)
            end
            
            % Bistatic SSR, Ts is calculated on the receiver
            flags = struct('isSSR', 1, 'isMonostatic', 0, 'isCalculatedTs', 1);
            [L, phi] = Lphi_calc( ssr_parameters, ssr_errors, t_errors, receiver, ssr, aircraft, flags, 1);
            [R1, R2] = R1R2_function(L, get_distance(ssr, receiver), phi, aircraft(3));
            if ~isempty(R1)
                [x, y, ok] = coordR_function(sqrt(R2^2 - aircraft(3)^2), sqrt(R1^2 - aircraft(3)^2), receiver, ssr, aircraft);
                if ((ok) && ~isempty(x) && ~isempty(y))
                    error_R_BSSR_calculated_Ts = [error_R_BSSR_calculated_Ts; sqrt((aircraft(1) - x)^2 + (aircraft(2) - y)^2)];
                end
            end
            
            % Bistatic SSR, Ts is known a priori
            flags = struct('isSSR', 1, 'isMonostatic', 0, 'isCalculatedTs', 0);
            [L, phi] = Lphi_calc( ssr_parameters, ssr_errors, t_errors, receiver, ssr, aircraft, flags, 1 );
            [R1, R2] = R1R2_function(L, get_distance(ssr, receiver), phi, aircraft(3));
            if ~isempty(R1)
                [x, y, ok] = coordR_function(sqrt(R2^2 - aircraft(3)^2), sqrt(R1^2 - aircraft(3)^2), receiver, ssr, aircraft);
                if ((ok) && ~isempty(x) && ~isempty(y))
                    error_R_BSSR_known_Ts = [error_R_BSSR_known_Ts; sqrt((aircraft(1) - x)^2 + (aircraft(2) - y)^2)];
                end
            end
            
            % Monostatic SSR
            flags = struct('isSSR', 1, 'isMonostatic', 1, 'isCalculatedTs', 1);
            [L, phi] = Lphi_calc( ssr_parameters, ssr_errors, t_errors, receiver, ssr, aircraft, flags, 1 );
            [ok, x, y] = MSSR(L, phi, aircraft(3), ssr);
            if ((ok) && ~isempty(x) && ~isempty(y))
                error_R_MSSR = [error_R_MSSR; sqrt((aircraft(1) - x)^2 + (aircraft(2) - y)^2)];
            end
            
            % Control System
            flags = struct('isSSR', 0, 'isMonostatic', 0, 'isCalculatedTs', 0);
            [L, phi] = Lphi_calc( ssr_parameters, ssr_errors, t_errors, receiver, ssr, aircraft, flags, N );
            [R1, R2] = R1R2_function(L, get_distance(ssr, receiver), phi, aircraft(3));
            if ~isempty(R1)
                [x, y, ok] = coordR_function(sqrt(R2^2 - aircraft(3)^2), sqrt(R1^2 - aircraft(3)^2), receiver, ssr, aircraft);
                if ((ok) && ~isempty(x) && ~isempty(y))
                    error_R_control_system = [error_R_control_system; sqrt((aircraft(1) - x)^2 + (aircraft(2) - y)^2)];
                end
            end
            
        end
        
        % Save LMS of errors
        lms_BSSR_calculated_Ts(iPhi) = sqrt((mean((error_R_BSSR_calculated_Ts - mean(error_R_BSSR_calculated_Ts)).^2)));
        lms_BSSR_known_Ts(iPhi) = sqrt((mean((error_R_BSSR_known_Ts - mean(error_R_BSSR_known_Ts)).^2)));
        lms_MSSR(iPhi) = sqrt((mean((error_R_MSSR - mean(error_R_MSSR)).^2)));
        lms_control_system(iPhi) = sqrt((mean((error_R_control_system - mean(error_R_control_system)).^2)));
    end
    
    % Form structures with all necessary results
    lms_bistatic_ssr{isnr} = lms_BSSR_known_Ts;
    lms_bistatic_ssr_no_Ts{isnr} = lms_BSSR_calculated_Ts;
    lms_control_system_N{isnr} = lms_control_system;
    lms_mssr{isnr} = lms_MSSR;
end
if TOSAVE
    % Save these results
    save([prefix 'lms_bistatic_ssr_all.mat'], 'lms_bistatic_ssr');
    save([prefix 'lms_bistatic_ssr_no_Ts.mat'], 'lms_bistatic_ssr_no_Ts');
    save([prefix 'lms_control_system_N.mat'], 'lms_control_system_N');
    save([prefix 'lms_mssr.mat'], 'lms_mssr');
end

%% plot

set(0,'DefaultAxesFontSize',11,'DefaultAxesFontName','Times New Roman');
figure();
markers = ['o' 'x' '*' '+' 'o' 'x' '*']; % the collection of plot markers

% Monostatic SSR and Control System (DAS - distance analysis station)
subplot(3,1,1);
title(['Monostatic SSRS and DAS arccuracies comparison']);
hold on
leg = [];
for i = 1 : length(SNR)
    plot(phi_real, lms_control_system_N{i}, ['-k' markers(i)]);
    leg{i} = ['d^2_{SSR} = ' num2str(SNR(i)) ' dB'];
end
legend(leg);
for i = 1 : length(SNR)
    plot(phi_real, lms_mssr{i}, [':k' markers(i)]);
end
grid on
xlabel('angle \gamma, rad');
ylabel('RAS error standard deviation \sigma_{\deltaR}, m');

% Bistatic SSR (Ts is known a priori) and Control System
subplot(3,1,2);
title(['Bistatic SSRS (with known t_{ri} moments) and DAS arccuracies comparison']);
hold on
leg = [];
for i = 1 : length(SNR)
    plot(phi_real, lms_control_system_N{i}, ['-k' markers(i)]);
    leg{i} = ['d^2_{SSR} = ' num2str(SNR(i)) ' dB'];
end
legend(leg);
for i = 1 : length(SNR)
    plot(phi_real, lms_bistatic_ssr{i}, [':k' markers(i)]);
end
grid on
xlabel('angle \gamma, rad');
ylabel('RAS error standard deviation \sigma_{\deltaR}, m');

% Bistatic SSR (Ts is calculated on the receiver) and Control System
subplot(3,1,3);
title(['Bistatic SSRS (without known t_{ri} moments) and DAS arccuracies comparison']);
hold on
leg = [];
count = 1;
for i = 1 : length(SNR)
    plot(phi_real, lms_control_system_N{i}, ['-k' markers(i)]);
    leg{i} = ['d^2_{SSR} = ' num2str(SNR(i)) ' dB'];
end
legend(leg);
for i = 1 : length(SNR)
    plot(phi_real, lms_bistatic_ssr_no_Ts{i}, [':k' markers(i)]);
end
grid on
xlabel('angle \gamma, rad');
ylabel('RAS error standard deviation \sigma_{\deltaR}, m');

% Correct axis
subplot(3,1,1)
axis([0,2*pi,0,max(lms_control_system_N{length(SNR)})*1.05])
subplot(3,1,2)
axis([0,2*pi,0,max(lms_control_system_N{length(SNR)})*1.05])
subplot(3,1,3)
axis([0,2*pi,0,max(lms_control_system_N{length(SNR)})*1.05])