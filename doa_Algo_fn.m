function doa_Algo_fn(M,theta_Degree)
%% Simualtion Parameters
ipFreq = 100;
samplingFreq = 100*ipFreq;
samplingPeriod = 1/samplingFreq;
sampleCount = 2*samplingFreq/ipFreq;

time = samplingPeriod * [0:1:sampleCount-1];
snr_db = 10;
snr_linear = 10^(snr_db/10);
noise_var = 1;

complexNoise = sqrt(noise_var/2) * (randn(1, length(time)) + 1i*randn(1, length(time)));
noisePwr = (complexNoise*complexNoise') * 1/length(time);
signalPower = snr_linear*noisePwr;

% number of Signal Sources
d = 2;
%theta_Degree = [-5 5];
theta = pi/180 * theta_Degree;

phi1 = 2*pi*rand(1,sampleCount);
phi2 = 2*pi*rand(1,sampleCount);


% Two orthogonal/uncorrelated sources prototypes
ipWave1 = exp(1i * (2*pi*ipFreq*time + phi1));
E1 = (ipWave1*ipWave1') * 1/length(time);
ipWave1 = sqrt(signalPower/E1) * exp(1i * (2*pi*ipFreq*time + phi1));

ipWave2 = exp(1i * (2*pi*ipFreq*time + phi2));
E2 = (ipWave2*ipWave2') * 1/length(time);
ipWave2 = sqrt(signalPower/E2) * exp(1i * (2*pi*ipFreq*time + phi2));

ipWaveSource = [ipWave1;ipWave2];
%plot(time,ipWave1);
cov = 1/length(time) * (ipWave1 *ipWave2');


% number of ULA Antenna Elements
%M = 8;
spatialLength = 0.5;
beta = 2*pi*spatialLength*sin(theta);
% Construct the array steering matrix for the ULA
steeringMatrix = exp(1i * (0:M-1).' * beta);

% Form the ipWaveData vector as d x M x length(time)
ipWaveData = zeros(d,M,length(time)); 

for i = 1:d
    for j=1:M
        ipWaveData(i,j,:) = ipWaveSource(i,:) .* exp(1i*beta(i)*j) + (sqrt(noise_var/2) * (randn(1, length(time)) + 1i*randn(1, length(time))));
        % signal = reshape(ipWaveData(1,1,:),1,200);
        % plot(time,real(signal));
        % hold on;
    end
end


%% At recevier , the Matrix X will be

X = reshape(sum(ipWaveData,1),M,length(time));

Rxx = 1/length(time) * X*X';

%% Coventional Beam former 

theta_grid_degrees = [-90 : 1 : 90];
theta_grid = theta_grid_degrees * pi/180;

for g = 1:length(theta_grid)
    steeringVector = exp(1i * 2*pi*spatialLength * sin(theta_grid(g)) * (0:M-1)');
    beamformerPower(g) = 10*log(abs(steeringVector' * Rxx * steeringVector) / M^2);
end

%% Capon Beam Former
for g = 1:length(theta_grid)
    steeringVector = exp(1i * 2*pi*spatialLength * sin(theta_grid(g)) * (0:M-1)');
    deno = steeringVector' * inv(Rxx) * steeringVector;
    capon_beamformerPower(g) = 10*log(abs(1/deno));
end


%% Plot Results with Full Details
figure('Color', 'w', 'Position', [100, 100, 800, 500]);
p1 = plot(theta_grid_degrees, beamformerPower, 'b-', 'LineWidth', 1.8);
hold on; grid on; grid minor;
p2 = plot(theta_grid_degrees, capon_beamformerPower, 'm-', 'LineWidth', 1.8);


% Mark True Source Angles
y_limits = [min(beamformerPower)-2, max(beamformerPower)+3];
for k = 1:d
    xline(theta_Degree(k), 'r--', sprintf('True: %d^o', theta_Degree(k)), ...
        'LineWidth', 1.5, 'LabelVerticalAlignment', 'bottom', ...
        'LabelHorizontalAlignment', 'center', 'FontSize', 10, 'Color', [0.8 0 0]);
end

% Automatically detect and mark peaks
[pks, locs] = findpeaks(beamformerPower, theta_grid_degrees, 'MinPeakProminence', 2);
p3 = plot(locs, pks, 'ro', 'MarkerSize', 7, 'MarkerFaceColor', 'r');

% Annotate detected peaks on the plot
for p = 1:length(locs)
    text(locs(p), pks(p) + 1.2, sprintf('Peak: %.1f^o', locs(p)), ...
        'HorizontalAlignment', 'center', 'FontWeight', 'bold', 'FontSize', 9);
end

% Graph Details
title(sprintf('Conventional vs. Capon (MVDR) Beamformer Spectrum (M = %d, SNR = %d dB)', M, snr_db), ...
    'FontSize', 12, 'FontWeight', 'bold');
xlabel('Angle of Arrival \theta (Degrees)', 'FontSize', 11);
ylabel('Spatial Power Spectrum (dB)', 'FontSize', 11);
xlim([-90 90]);
ylim(y_limits);
legend([p1, p2, p3], ...
    {'Conventional (Delay-and-Sum)', 'Capon (MVDR) Spectrum', 'Detected Peaks'}, ...
    'Location', 'northeast');
end