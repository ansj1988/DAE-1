% Antenna gain estimates from noisy measurement data
% Developed by: Artur Nogueira de São José
% Graduate Program in Electrical Engineering - UFMG, Brazil

% Technique: DAE-1
% Current status: under review by a Ph.D. jury

% If you use this code, please cite this paper:

% @article{NDESAOJOSE2020107720,
% title = "Improving antenna gain estimations in non-ideal test sites with auto-tunable filters",
% journal = "Measurement",
% volume = "159",
% pages = "107720",
% year = "2020",
% issn = "0263-2241",
% doi = "https://doi.org/10.1016/j.measurement.2020.107720",
% url = "http://www.sciencedirect.com/science/article/pii/S026322412030258X",
% author = "Artur {N. de São José} and Virginie Deniau and Úrsula {do C. Resende} and Ricardo Adriano"
% }

%% Step 1: loading variables and configuring the tool

clear
close all
clc

% Examples of input signals:
% 1. Directional patch antenna
% 2. Omnidirectional patch antenna

load('directional.mat');
% load('omni.mat');

%% Step 2: calculations

% Signal length
N = length(hA0_1); 

% Average S12 curve per test setup arrangement
M_arrange1 = (hA0_1 + hA0_2 + hA0_3 + hA0_4 + hA0_5 + hA0_6 + hA0_7 + hA0_8 + hA0_9 + hA0_10)./10;
M_arrange2 = (hA45_1 + hA45_2 + hA45_3 + hA45_4 + hA45_5 + hA45_6 + hA45_7 + hA45_8 + hA45_9 + hA45_10)./10;
M_arrange3 = (hB0_1 + hB0_2 + hB0_3 + hB0_4 + hB0_5 + hB0_6 + hB0_7 + hB0_8 + hB0_9 + hB0_10)./10;
M_arrange4 = (hB0_1 + hB0_2 + hB0_3 + hB0_4 + hB0_5 + hB0_6 + hB0_7 + hB0_8 + hB0_9)./9; %omni
M_arrange5 = (hB45_1 + hB45_2 + hB45_3 + hB45_4 + hB45_5 + hB45_6 + hB45_7 + hB45_8 + hB45_9 + hB45_10)./10;

s12_dae = (M_arrange1 + M_arrange2 +  M_arrange3 + M_arrange4 + M_arrange5)./5;

%% Step 3: benchmark analysis 
% Implementing the technique proposed by Froes et al (2019) - see https://ieeexplore.ieee.org/document/8606248 
% We call it DAE-2

% B: max curves, S: min curves, M: mean curves

for k=1:N

    B(k,1) = abs(max([M_arrange1(k); M_arrange2(k); M_arrange4(k); M_arrange5(k)]));
    S(k,1) = abs(min([M_arrange1(k); M_arrange2(k); M_arrange4(k); M_arrange5(k)]));
    M(k,1) = abs(s12_dae(k));
    
    % How close is M from B?
    p(k,1) = 1 - ( B(k,1) - M(k,1) ) / ( B(k,1) - S(k,1) ); 
    
    % How close is M from S?
    q(k,1) = 1 - ( M(k,1) - S(k,1) ) / ( B(k,1) - S(k,1) );
    
    % Weights used to recover the original S12 curve
    w1(k,1) = exp(p(k,1) - 0.5);
    w2(k,1) = exp(q(k,1) - 0.5);

    % Retrieved S12 curve
    s12_froes(k,1) = ( w1(k,1)*( ( B(k,1)+M(k,1) )/2 ) + w2(k,1)*( ( M(k,1)+S(k,1) )/2 ) ) / ( w1(k,1) + w2(k,1) ); 

end

%% Step 4: converting S12 into gain

d = 1;                                                                            % Distance between the antennas (in meters)
lambda = (3e8./freqHz)';                                                          % Wavelengths (in meters)
gain_lab = 10*log10( (4*pi*d./lambda).*abs(s12_lea)' );                           % Measured signal
gain_chamber = 10*log10( (4*pi*d./lambda).*abs(s12_camara)' );                    % Reference gain curve
gain_dae1 = 10*log10( (4*pi*d./lambda).*abs(s12_dae)' );                          % Filtered signal using DAE-1
gain_dae2 = 10*log10( (4*pi*d./lambda).*abs(s12_froes)' );                        % Filtered signal using DAE-2

%% Step 5: graphs

figure
hold on
plot(freqHz,gain_lab,'b')
plot(freqHz,gain_dae1,'y')
plot(freqHz,gain_dae2,'k')
plot(freqHz,gain_chamber,'r')
ylabel('Gain (dB)')
xlabel('Frequency (Hz)')
legend('Original signal','DAE-1','DAE-2','Chamber')
hold off

% A filtering quality measure based on the correlation between the cleaned
% and reference (anechoic chamber) gain curves
correl_nofilter = corr(gain_lab',gain_chamber','Type','Pearson')
correl_dae1 = corr(gain_dae1',gain_chamber','Type','Pearson')
correl_dae2 = corr(gain_dae2',gain_chamber','Type','Pearson')

