function [ L, phi ] = Lphi_calc( ssr_parameters, ssr_errors, toa_lms, receiver, ssr, aircraft, flags, N )
% output
L = [];
phi = [];

%% initialisation

Ta = 1/ssr_parameters.speed; %rotation period, s
ampl_lms = 1e-5; %approximate lms of signals amplitude
t_processing = 3*1e-6; %approximate signal processing time on aorcrafts, s
t_processing_lms = 1e-10; %approximate LMS of time of signal processing on aircrafts, s
theta = 0.5*pi/180; %the SSR antenna radiation pattern width, rad
beta = 20*pi/180; %start position of SSR antenna, rad (beta is the angle between "SSR-receiver" and "SSR-aircraft")

%% initial values for the main cycle below

% all lengths are unknown a priori
t_1030 = [];
amp_1030 = [];
t_1090 = [];
amp_1090 = [];
all_t_1030_for_SSR = [];

% values that will be changed in cycles
t = 0; %time stamp
angle = -beta; %antenna rotaton angle
    
%% distances initial calculation

r_receiver_ssr = get_distance(ssr, receiver);
r_receiver_ssr_prj = get_distance(ssr(1:2), receiver(1:2)); % projection
r_receiver_aircraft = get_distance(receiver, aircraft);
r_receiver_aircraft_prj = get_distance(receiver(1:2), aircraft(1:2)); % projection
r_aircraft_ssr = get_distance(aircraft, ssr);
r_aircraft_ssr_prj = get_distance(aircraft(1:2), ssr(1:2)); % projection

% alpha correction
alpha = acos((r_aircraft_ssr_prj^2 + r_receiver_ssr_prj^2 - r_receiver_aircraft_prj^2) / (2*r_aircraft_ssr_prj*r_receiver_ssr_prj));
if (aircraft(2) < 0)
    alpha = -alpha;
end

%% the main cycle - calculates the TOA and the amplitude of answer and request signals on the receiver and the aircraft

constr_1090 = [beta+alpha-theta/2, beta+alpha+theta/2]; %this is the antenna sector in which the aircraft will receive the request signals and therefore will transmit the answer signals
if (alpha < 0)
    constr_1090 = constr_1090 + 2*pi;
end

while (angle <= (beta + theta/2 + 2*pi)) % the main cycle
    all_t_1030_for_SSR = [all_t_1030_for_SSR t]; %this is all timestamps in which the request signals are transmitted (in is necessary for MSSR)
    
    if (angle >= (beta - theta/2)) && (angle <= (beta + theta/2)) %this is the antenna sector in which the request signals will be received on the receiver
        amp_1030 = [amp_1030 abs(sin((angle - beta)/theta + pi/2)) + ampl_lms*randn()];
        t_1030 = [t_1030 t + r_receiver_ssr / physconst('LightSpeed') + randn*toa_lms.requests_on_receiver];
    elseif (angle >= (beta - theta/2 + 2*pi)) %the additional check for the +2pi chase
        amp_1030 = [amp_1030 abs(sin((angle - beta - 2*pi)/theta + pi/2)) + ampl_lms*randn()];
        t_1030 = [t_1030 t + r_receiver_ssr / physconst('LightSpeed') + randn*toa_lms.requests_on_receiver];%t_error(SNR, '3ATx', reference_signal_req, 0, Fs)];
    end
    
    if (angle >= constr_1090(1)) && (angle <= constr_1090(2)) %in this sector the request signals will be received on the aircraft
        if (flags.isMonostatic)
            amp_1090 = [amp_1090 1 + ampl_lms*randn()];
            t_1090 = [t_1090 t + r_aircraft_ssr / physconst('LightSpeed') + t_processing+ r_aircraft_ssr / physconst('LightSpeed')  + randn*t_processing_lms + randn*toa_lms.requests_on_aircraft + randn*toa_lms.answers];
        else
            amp_1090 = [amp_1090 1 + ampl_lms*randn()];
            t_1090 = [t_1090 t + r_receiver_aircraft / physconst('LightSpeed') + t_processing + randn*t_processing_lms + r_aircraft_ssr / physconst('LightSpeed') + randn*toa_lms.requests_on_aircraft + randn*toa_lms.answers];
        end
    end
    
    % t and angle update
    t_last = t;
    t = t + ssr_parameters.PRI + ssr_errors.PRI_error*rand - ssr_errors.PRI_error/2;
    angle = angle + (t - t_last)*(ssr_parameters.speed + ssr_errors.speed_error*randn) * 2*pi;
end

%% the additional cycle for the Control System in case the SSR parameters will be analysed based on the N antenna rotations

if ~flags.isSSR % if it is the Control System
    angle = angle - 2*pi;
    % in this case we will need the long vectors of TOA and amplitudes of all request signals during the N antenna rotations
    t_1030_long = [t_1030];
    amp_1030_long = [amp_1030];
    
    count = 0; %antenna rotation period number
    while (count < (N - 1))
        if (angle >= (beta - theta/2)) && (angle <= (beta + theta/2))
            amp_1030_long = [amp_1030_long abs(sin((angle - beta)/theta + pi/2)) + ampl_lms*randn()];
            t_1030_long = [t_1030_long t + r_receiver_ssr / physconst('LightSpeed') + randn*toa_lms.requests_on_receiver];%t_error(SNR, '3ATx', reference_signal_req, 0, Fs)];
        elseif (angle >= (beta - theta/2 + 2*pi))
            amp_1030_long = [amp_1030_long abs(sin((angle - beta - 2*pi)/theta + pi/2)) + ampl_lms*randn()];
            t_1030_long = [t_1030_long t + r_receiver_ssr / physconst('LightSpeed') + randn*toa_lms.requests_on_receiver];%t_error(SNR, '3ATx', reference_signal_req, 0, Fs)];
        end

        t_last = t;
        t = t + ssr_parameters.PRI + ssr_errors.PRI_error*rand - ssr_errors.PRI_error/2;
        angle = angle + (t - t_last)*(ssr_parameters.speed + ssr_errors.speed_error*randn) * 2*pi;
        if (angle > (beta + theta/2 + 2*pi))
            angle = angle - 2*pi;
            count = count + 1;
        end
    end
end

%% Additional check in case of any possible errors

if isempty(t_1090) || isempty(t_1030) || (length(t_1030) == 1) || (min(t_1030) > max(t_1090))
    return;
end

%% L and phi calculations based on all answer and request signals for the Monostatic SSR

if flags.isMonostatic
    phi = (max(amp_1090) + min(amp_1090))/2;
    all_deltas = zeros(1, length(t_1090));
    for i = 1 : length(t_1090)
        all_deltas(i) = min(abs(all_t_1030_for_SSR - t_1090(1)))- t_processing;
    end
    L = physconst('LightSpeed')*mean(all_deltas);
    return;
end

%% Ts calculations for the Control System with the corresponding flag

if (flags.isCalculatedTs) && (~flags.isSSR)
    delta_t_1030_all = t_1030_long(2:end) - t_1030_long(1:end-1);
    [max_meaning, ~] = max(delta_t_1030_all);
    [min_meaning, ~] = min(delta_t_1030_all);
    threshold = (max_meaning - min_meaning) / 2;
    delta_t_1030_filtered = delta_t_1030_all(delta_t_1030_all < threshold);
    Ts_experimental = mean(delta_t_1030_filtered);
else
    Ts_experimental = ssr_parameters.PRI;
end

%% L and phi calculations based on all answer and request signals for the Control System and the Bistatic SSR

delta_t_1030 = t_1030(2:end) - t_1030(1:end-1);
[~, max_index] = max(delta_t_1030); % max_index is the end of the first package (because t_1030 contains only two packages)
[~, maxindreq] = max(amp_1030(1 : max_index)); % the index of the "most powerful" impulse in the first request signals package (and also the package center)
[~, maxindreqnext] = max(amp_1030(max_index + 1 : end));
maxindans = round(maxindreq * length(t_1090)/max_index); % we assume that it is the index of the index of the impulse referred to the answer signals package center

% check all negative scenarios
if isempty(maxindans) || isempty(maxindreq) || (maxindans > length(t_1090)) || ...
        (maxindreq > length(t_1030)) || (~maxindans) || (~maxindreq)
    return;
end

% calculate the phi value based on the package centres
if (t_1090(maxindans) - t_1030(maxindreq) < Ta/2)
    phi = ((t_1090(maxindans) - t_1030(maxindreq)) / (t_1030(max_index + 1) - t_1030(max_index))) * 2*pi;
else
    phi = 2*pi - ((t_1030(max_index + maxindreqnext) - t_1090(maxindans)) / (t_1030(max_index + 1) - t_1030(max_index))) * 2*pi;
end

% calculate the L value for the Bistatic SSR
if flags.isSSR
    all_deltas = [];
    for i = 1 : length(t_1090)
        all_deltas = [all_deltas min(abs(all_t_1030_for_SSR + r_receiver_ssr / physconst('LightSpeed') - t_1090(i)))- t_processing];
    end
    L = physconst('LightSpeed')*mean(all_deltas);
    return;
end

% calculate the L value for the Control System
if t_1030(max_index) > t_1090(maxindans) % all first request package TOA are higher that t_1090(maxindans)
    L = physconst('LightSpeed') * (min(abs(t_1030(1 : max_index)-t_1090(maxindans))) - t_processing);
elseif (t_1030(max_index + 1) < t_1090(maxindans)) % all second request package TOA are lower that t_1090(maxindans)
    L = physconst('LightSpeed') * (min(abs(t_1030-t_1090(maxindans))) - t_processing);
elseif ((t_1090(maxindans) - t_1030(max_index)) < (t_1030(max_index + 1)- t_1090(maxindans)))  % t_1090(maxindans) is between two request packages and closer to the first request package than to the second one
    t_1030_first_package = t_1030(1 : max_index); % the final length is unknown a priori
    % interpolate request signals TOA
    while (t_1030_first_package(end) < t_1090(maxindans))
        t_1030_first_package = [t_1030_first_package t_1030_first_package(end) + Ts_experimental];
    end
    L = physconst('LightSpeed') * (t_1090(maxindans) - t_1030_first_package(end - 1) - t_processing);
else   % t_1090(maxindans) is between two request packages and closer to the second request package than to the first one
    t_1030_second_package = t_1030(max_index + 1 : end); % the final length is unknown a priori
    % interpolate request signals TOA
    while (t_1030_second_package(1) > t_1090(maxindans))
        t_1030_second_package = [t_1030_second_package(1) - Ts_experimental t_1030_second_package];
    end
    L = physconst('LightSpeed') * (t_1090(maxindans) - t_1030_second_package(1)- t_processing);
end