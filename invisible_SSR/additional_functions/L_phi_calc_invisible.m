function [ L, phi, ok ] = L_phi_calc_invisible( ssr_parameters, ssr_errors, toa_lms, receiver, ssr, aircraft, reference_aircraft, isCalculatedTs, N, reference_aircraft_error )
% output
L = [];
phi = [];
ok = 0;

%% initialisation

addpath ..\src
Ta = 1/ssr_parameters.speed; %rotation period, s
ampl_lms = 1e-5; %approximate lms of signals amplitude
t_processing = 3*1e-6; %approximate signal processing time on aorcrafts, s
t_processing_lms = 1e-10; %approximate LMS of time of signal processing on aircrafts, s
theta = 0.5*pi/180; %the SSR antenna radiation pattern width, rad
beta = 20*pi/180; %start position of SSR antenna, rad (beta is the angle between "SSR-receiver" and "SSR-aircraft")

b = get_distance(receiver, ssr);

%% initial values for the main cycle below

% all lengths are unknown a priori
t_1030 = [];
amp_1030 = [];
t_1090 = [];
amp_1090 = [];
t_1090_ref = [];
amp_1090_ref = [];
all_t_1030_for_SSR = [];

% values that will be changed in cycles
t = 0; %time stamp
angle = -beta; %antenna rotaton angle
    
%% distances initial calculation

r_receiver_ssr = get_distance(ssr, receiver);
r_receiver_ssr_prj = get_distance(ssr(1:2), receiver(1:2)); %projection
r_receiver_aircraft = get_distance(receiver, aircraft);
r_receiver_aircraft_prj = get_distance(receiver(1:2), aircraft(1:2)); %projection
r_aircraft_ssr = get_distance(aircraft, ssr);
r_aircraft_ssr_prj = get_distance(aircraft(1:2), ssr(1:2)); %projection
r_receiver_aircraft_ref = get_distance(receiver, reference_aircraft);
r_receiver_aircraft_ref_prj = get_distance(receiver(1:2), reference_aircraft(1:2)); %projection
r_aircraft_ref_ssr = get_distance(reference_aircraft, ssr);
r_aircraft_ref_ssr_prj = get_distance(reference_aircraft(1:2), ssr(1:2)); %projection

% alpha correction
alpha = acos((r_aircraft_ssr_prj^2 + r_receiver_ssr_prj^2 - r_receiver_aircraft_prj^2) / (2*r_aircraft_ssr_prj*r_receiver_ssr_prj));
if (aircraft(2) < 0)
    alpha = -alpha;
end

% alpha for the reference plane correction
alpha_ref = acos((r_aircraft_ref_ssr_prj^2 + r_receiver_ssr_prj^2 - r_receiver_aircraft_ref_prj^2) / (2*r_aircraft_ref_ssr_prj*r_receiver_ssr_prj));
if (reference_aircraft(2) < 0)
    alpha_ref = -alpha_ref;
end

%% the main cycle - calculates the TOA and the amplitude of answer and request signals on the receiver and the aircraft

constr_1090 = [rem(beta+alpha-theta/2+4*pi, 2*pi), rem(beta+alpha+theta/2+4*pi,2*pi)]; %this is the antenna sector in which the aircraft will receive the request signals and therefore will transmit the answer signals
constr_1090_ref = [rem(beta+alpha_ref-theta/2+4*pi, 2*pi), rem(beta+alpha_ref+theta/2+4*pi, 2*pi)]; %this is the antenna sector for the reference aircraft

while (angle <= (beta + theta/2 + 4*pi)) % the main cycle
    angle_to_compare = rem(angle + 4*pi, 2*pi);
    all_t_1030_for_SSR = [all_t_1030_for_SSR t];
    
    if (angle_to_compare >= (beta - theta/2)) && (angle_to_compare <= (beta + theta/2))
        amp_1030 = [amp_1030 abs(sin((angle_to_compare - beta)/theta + pi/2)) + ampl_lms*randn()];
        t_1030 = [t_1030 t + r_receiver_ssr / physconst('LightSpeed') + randn*toa_lms.requests_on_receiver];%t_error(SNR, '3ATx', reference_signal_req, 0, Fs)];
    elseif (angle_to_compare >= (beta - theta/2 + 2*pi)) && (angle_to_compare <= (beta + theta/2 + 2*pi))
        amp_1030 = [amp_1030 abs(sin((angle - beta - 2*pi)/theta + pi/2)) + ampl_lms*randn()];
        t_1030 = [t_1030 t + r_receiver_ssr / physconst('LightSpeed') + randn*toa_lms.requests_on_receiver];%t_error(SNR, '3ATx', reference_signal_req, 0, Fs)];
    end
    
    if (angle_to_compare >= constr_1090(1)) && (angle_to_compare <= constr_1090(2))
        amp_1090 = [amp_1090 1 + ampl_lms*randn()];
        t_1090 = [t_1090 t + r_receiver_aircraft / physconst('LightSpeed') + t_processing + randn*t_processing_lms + r_aircraft_ssr / physconst('LightSpeed') + randn*toa_lms.requests_on_aircraft + randn*toa_lms.answers];%t_error(SNRreq1, '3ATx', reference_signal_req, 0, Fs) + t_error(SNRans, 'Rx', reference_signal_ans, code, Fs)];
    end
    
    if (angle_to_compare >= constr_1090_ref(1)) && (angle_to_compare <= constr_1090_ref(2))
        amp_1090_ref = [amp_1090_ref 1 + ampl_lms*randn()];
        t_1090_ref = [t_1090_ref t + r_receiver_aircraft_ref / physconst('LightSpeed') + t_processing + randn*t_processing_lms + r_aircraft_ref_ssr / physconst('LightSpeed') + randn*toa_lms.requests_on_aircraft + randn*toa_lms.answers];%t_error(SNRreq1, '3ATx', reference_signal_req, 0, Fs) + t_error(SNRans, 'Rx', reference_signal_ans, code, Fs)];
    end
    
    % t and angle update
    t_last = t;
    t = t + ssr_parameters.PRI + ssr_errors.PRI_error*rand - ssr_errors.PRI_error/2;
    angle = angle + (t - t_last)*(ssr_parameters.speed + ssr_errors.speed_error*randn) * 2*pi;
end

%% the additional cycle for the Control System in case the SSR parameters will be analysed based on the N antenna rotations

angle = angle - 2*pi;
count = 0;

% all lengths are unknown a priori
t_1030_long = [t_1030];
amp_1030_long = [amp_1030];
t_1090_long = [t_1090];
amp_1090_long = [amp_1090];
t_1090_ref_long = [t_1090_ref];
amp_1090_ref_long = [amp_1090_ref];

while (count < (N - 1))
    angle_to_compare = rem(angle + 4*pi, 2*pi);
    if (angle_to_compare >= (beta - theta/2)) && (angle_to_compare <= (beta + theta/2))
        amp_1030_long = [amp_1030_long abs(sin((angle - beta)/theta + pi/2)) + ampl_lms*randn()];
        t_1030_long = [t_1030_long t + r_receiver_ssr / physconst('LightSpeed') + randn*toa_lms.requests_on_receiver];%t_error(SNR, '3ATx', reference_signal_req, 0, Fs)];
    elseif angle_to_compare >= (beta - theta/2 + 2*pi)
        amp_1030_long = [amp_1030_long abs(sin((angle - beta - 2*pi)/theta + pi/2)) + ampl_lms*randn()];
        t_1030_long = [t_1030_long t + r_receiver_ssr / physconst('LightSpeed') + randn*toa_lms.requests_on_receiver];%t_error(SNR, '3ATx', reference_signal_req, 0, Fs)];
    end

    if (angle_to_compare >= constr_1090(1)) && (angle_to_compare <= constr_1090(2))
        amp_1090_long = [amp_1090_long 1 + ampl_lms*randn()];
        t_1090_long = [t_1090_long t + r_receiver_aircraft / physconst('LightSpeed') + t_processing + randn*t_processing_lms + r_aircraft_ssr / physconst('LightSpeed') + randn*toa_lms.requests_on_aircraft + randn*toa_lms.answers];%t_error(SNRreq1, '3ATx', reference_signal_req, 0, Fs) + t_error(SNRans, 'Rx', reference_signal_ans, code, Fs)];
    end

    if (angle_to_compare >= constr_1090_ref(1)) && (angle_to_compare <= constr_1090_ref(2))
        amp_1090_ref_long = [amp_1090_ref_long 1 + ampl_lms*randn()];
        t_1090_ref_long = [t_1090_ref_long t + r_receiver_aircraft_ref / physconst('LightSpeed') + t_processing + randn*t_processing_lms + r_aircraft_ref_ssr / physconst('LightSpeed') + randn*toa_lms.requests_on_aircraft + randn*toa_lms.answers];%t_error(SNRreq1, '3ATx', reference_signal_req, 0, Fs) + t_error(SNRans, 'Rx', reference_signal_ans, code, Fs)];
    end

    t_last = t;
    t = t + ssr_parameters.PRI + ssr_errors.PRI_error*rand - ssr_errors.PRI_error/2;
    angle = angle + (t - t_last)*(ssr_parameters.speed + ssr_errors.speed_error*randn) * 2*pi;
    if (angle > (beta + theta/2 + 2*pi))
        angle = angle - 2*pi;
        count = count + 1;
    end
end

%% Additional check in case of any possible errors

if isempty(t_1090) || isempty(t_1030) || (length(t_1030) == 1) || (min(t_1030) > max(t_1090)) || isempty(t_1090_ref)
    return;
end

%% the addition of the reference aircraft position error

reference_aircraft(1) = reference_aircraft(1) + randn()*reference_aircraft_error;
reference_aircraft(2) = reference_aircraft(2) + randn()*reference_aircraft_error;

%% Ts calculations for the Control System with the corresponding flag

delta_t_1090 = [t_1090_long(2:end) - t_1090_long(1:end-1)];
delta_t_1090_ref = [t_1090_ref_long(2:end) - t_1090_ref_long(1:end-1)];
if isCalculatedTs
    [max_meaning, ~] = max(delta_t_1090);
    [min_meaning, ~] = min(delta_t_1090);
    threshold = (max_meaning - min_meaning) / 2;
    [max_meaning, ~] = max(delta_t_1090_ref);
    [min_meaning, ~] = min(delta_t_1090_ref);
    threshold_ref = (max_meaning - min_meaning) / 2;

    Ts_hypothetical = [delta_t_1090(delta_t_1090 < threshold) delta_t_1090_ref(delta_t_1090_ref < threshold_ref)];
    Ts_experimental = mean(Ts_hypothetical);
    
    %Ta_hypothetical = [delta_t_1090(delta_t_1090 >= threshold) delta_t_1090_ref(delta_t_1090_ref >= threshold_ref)];
    %Ta_experimental = mean(Ta_hypothetical);
else
    Ts_experimental = ssr_parameters.PRI;
    %Ta_experimental = Ta;
end

%% find the first packages of the answer signals from the desired and the reference aircrafts

idx = find(delta_t_1090 > Ta/2);
t_1090_first_package = t_1090(1 : idx - 1);
idx = find(delta_t_1090_ref > Ta/2);
t_1090_ref_first_package = t_1090_ref(1 : idx - 1);

%% beta calculation based on all answer signals from the desired and the reference aircrafts
% (beta is the angle between the "SSR - desired aircraft" and the "SSR - reference aircraft")

delta_toa_in_sec = t_1090_first_package(1) - t_1090_ref_first_package(1);
angle = 2*pi * delta_toa_in_sec / Ta;
beta = rem((angle + 2*pi), 2*pi);

%% L calculation based on all answer signals from the desired and the reference aircrafts

if t_1090_first_package(1) < t_1090_ref_first_package(1)
    earlier_package = t_1090_first_package;
    later_package = t_1090_ref_first_package;
else
    earlier_package = t_1090_ref_first_package;
    later_package = t_1090_first_package;
end

while earlier_package(end) < later_package(end)
    earlier_package = [earlier_package earlier_package(end) + Ts_experimental];
end
deltas = [];
for i = 1 : length(later_package)
    deltas = [deltas min(abs(earlier_package - later_package(i)))];
end
delta_l = mean(deltas) * physconst('LightSpeed');
needed_L = delta_l + get_distance(reference_aircraft, ssr) + get_distance(reference_aircraft, receiver);
L = needed_L - get_distance(receiver, ssr);

%% phi calculation

% all coordinates will be temporary rotate to the alpha angle because it is
% more convenient
ssr_init = ssr; % the initial coordinates
receiver_init = receiver; % the initial coordinates
reference_aircraft_init = reference_aircraft; % the initial coordinates
aircraft_init = aircraft; % the initial coordinates

% also the coordinates center will be moved to the SSR position
receiver = [receiver - ssr];
reference_aircraft = [reference_aircraft - ssr];
aircraft = [aircraft - ssr];
ssr = [0 0 0];

alpha = atan(abs(receiver(2))/abs(receiver(1)));
% alpha correction
if (receiver(1) < 0) && (receiver(2) >= 0)
    alpha = pi - alpha;
elseif (receiver(1) < 0) && (receiver(2) < 0)
    alpha = pi + alpha;
elseif (receiver(1) > 0) && (receiver(2) < 0)
    alpha = -alpha;
end
alpha = 2*pi - alpha;

% coordinates rotation
receiver = [(receiver(1)*cos(alpha) - receiver(2)*sin(alpha)) (receiver(1)*sin(alpha) + receiver(2)*cos(alpha))];
ssr = [(ssr(1)*cos(alpha) - ssr(2)*sin(alpha)) (ssr(1)*sin(alpha) + ssr(2)*cos(alpha))];
reference_aircraft = [(reference_aircraft(1)*cos(alpha) - reference_aircraft(2)*sin(alpha)) (reference_aircraft(1)*sin(alpha) + reference_aircraft(2)*cos(alpha))];
aircraft = [(aircraft(1)*cos(alpha) - aircraft(2)*sin(alpha)) (aircraft(1)*sin(alpha) + aircraft(2)*cos(alpha))];

% the angle between "SSR - reference aircraft" and "SSR - receiver"
angle_reference_aircraft_ssr_receiver = atan(abs(reference_aircraft(2)/reference_aircraft(1)));
if (reference_aircraft(1) < 0) && (reference_aircraft(2) >= 0)
    angle_reference_aircraft_ssr_receiver = pi - angle_reference_aircraft_ssr_receiver;
elseif (reference_aircraft(1) < 0) && (reference_aircraft(2) < 0)
    angle_reference_aircraft_ssr_receiver = pi + angle_reference_aircraft_ssr_receiver;
elseif (reference_aircraft(1) > 0) && (reference_aircraft(2) < 0)
    angle_reference_aircraft_ssr_receiver = -angle_reference_aircraft_ssr_receiver;
end

% phi is the angle consists of the sum of the beta and the previous angle
sum_angle = beta - angle_reference_aircraft_ssr_receiver;
phi = rem(sum_angle, 2*pi);

% return all coordinates to their initial state
ssr = ssr_init;
receiver = receiver_init;
reference_aircraft = reference_aircraft_init;
aircraft = aircraft_init;

%% the return flag is now equal to 1

ok = 1;