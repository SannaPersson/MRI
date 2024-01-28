close all; clear all; clc;

fov = 256e-3; % Field of View in meters
Nx = 96; % Matrix size in x-direction
Ny = Nx; % Matrix size in y-direction (assuming square matrix)
sliceThickness = 3e-3; % Slice thickness in meters
Nslices = 1; % Number of slices
Oversampling = 2; % Oversampling factor
phi = pi/2; % Orientation of the readout

% Set system limits
sys = mr.opts('MaxGrad',22,'GradUnit','mT/m',...
    'MaxSlew',160,'SlewUnit','T/m/s',...
    'rfRingdownTime', 30e-6, 'rfDeadtime', 100e-6, 'adcDeadTime', 10e-6);  
seq = mr.Sequence(sys); % Create a new sequence object
warning('OFF', 'mr:restoreShape'); % Disable warning

% Create fat-sat pulse
B0 = 3; % Main magnetic field strength in Tesla
sat_ppm = -3.45; % Chemical shift for fat saturation in ppm
sat_freq = sat_ppm * 1e-6 * B0 * sys.gamma; % Convert ppm to Hz
rf_fs = mr.makeGaussPulse(110*pi/180, 'system', sys, 'Duration', 8e-3,...
    'bandwidth', abs(sat_freq), 'freqOffset', sat_freq);
gz_fs = mr.makeTrapezoid('z', sys, 'delay', mr.calcDuration(rf_fs), 'Area', 1/1e-4); % spoil up to 0.1mm

% Create 90 degree slice selection pulse and gradient
[rf, gz] = mr.makeSincPulse(pi/2, 'system', sys, 'Duration', 3e-3,...
    'SliceThickness', sliceThickness, 'apodization', 0.5, 'timeBwProduct', 4);

% Define k-space and trajectory parameters
deltaK = 1 / (2*fov); % Deltak
kRadius = round(Nx/2)
TR = 0.05; % Repetition time of RF pulses in seconds (example value)
T = TR; % Duration of one shot
A = deltaK / T; % Scaling factor
omega = 2500     %;1;%5*pi%pi / T; % Angular frequency for sinusoidal motion
numPoints = 150; % Number of points in the trajectory (example value)




% Generate the sinusoidal trajectory for the first quadrant
kx = zeros(1, numPoints);
ky = zeros(1, numPoints);
for tIdx = 1:numPoints
    t = (tIdx - 1) * T / numPoints;
    if t <= T / 2
        kx(tIdx) = A * t * (1 + sin(omega * t));
        ky(tIdx) = A * t * (1 - sin(omega * t));
    else
        kx(tIdx) = A * (t + (T - t) * sin(omega * t));
        ky(tIdx) = A * (t - (T - t) * sin(omega * t));
    end
end
plot(kx, ky);
hold on;
% Repeat trajectory calculation for other quadrants by scaling and translating
% Quadrant 2 (Q2), Quadrant 3 (Q3), Quadrant 4 (Q4)
kxQ1 = kx; kyQ1 = ky;
kxQ2 = -kx; kyQ2 = ky;
kxQ3 = -kx; kyQ3 = -ky;
kxQ4 = kx; kyQ4 = -ky;


% Concatenate all quadrants to form the full trajectory
kx_full = [kxQ1; kxQ2; kxQ3; kxQ4];
ky_full = [kyQ1; kyQ2; kyQ3; kyQ4];
%figure(1)
%h= scatter(kx_full, ky_full, "filled")
%set(h,'markersize',4)

ka = [kx_full; ky_full]
safety_margin=0.94; % we need that  otherwise we just about violate the slew rate due to the rounding errors
gzReph = mr.makeTrapezoid('z',sys,'Area',-gz.area/2);
% Generate the sinusoidal trajectory and gradients for each shot
gradient_waveforms = cell(1, 4); % Store gradients for each shot
adc_events = cell(1, 4); % Store ADC events for each shot
spoiler_waveforms = cell(1, 4); %
for shot = 1:4
    % Calculate k-space trajectory for the current shot
    kx_shot = kx_full(shot, :);
    ky_shot = ky_full(shot, :);

    % Calculate gradients from the trajectory
    [ga, sa] = mr.traj2grad([kx_shot; ky_shot]);

    % limit analysis
    safety_margin=0.94; % we need that  otherwise we just about violate the slew rate due to the rounding errors
    dt_gabs=abs(ga(1,:)+1i*ga(2,:))/(sys.maxGrad*safety_margin)*sys.gradRasterTime;
    dt_sabs=sqrt(abs(sa(1,:)+1i*sa(2,:))/(sys.maxSlew*safety_margin))*sys.gradRasterTime;
    dt_smooth=max([dt_gabs;dt_sabs]);
    
    % apply the lower limit not to lose the trajectory detail
    dt_min=4*sys.gradRasterTime/numPoints; % we want at least 4 points per revolution
    dt_smooth(dt_smooth<dt_min)=dt_min;
    t_smooth=[0 cumsum(dt_smooth,2)];
    kopt_smooth=interp1(t_smooth, [kx_shot; ky_shot]', (0:floor(t_smooth(end)/sys.gradRasterTime))*sys.gradRasterTime)';

    % analyze what we've got
    fprintf('duration orig %d us\n', round(1e6*sys.gradRasterTime*length([kx_shot; ky_shot])));
    fprintf('duration opt %d us\n', round(1e6*sys.gradRasterTime*length(kopt_smooth)));
    size(kopt_smooth)
    plot(kopt_smooth(1,:), kopt_smooth(2,:));
    hold on;
    [gos, sos] = mr.traj2grad(kopt_smooth);
    % Rescale gradients to system limits if necessary
    size(gos)
    gos = gos ./ max(max(abs(gos)) / sys.maxGrad, 2);
    sos = sos ./ max(max(abs(sos)) / sys.maxSlew, 2);
    size(gos)
    gx_shot = gos(1, :);
    gy_shot = gos(2, :);
    size(gx_shot)
    
    % Make sure the gradient waveforms start and end at zero
    gos(1, 1) = 0; % Set first point of X gradient to zero%
    gos(2, 1) = 0; % Set first point of Y gradient to zero
    gos(1, end) = 0; % Set last point of X gradient to zero
    gos(2, end) = 0; % Set last point of Y gradient to zero


    % calculate ADC
    % round-down dwell time to 10 ns
    adcTime = sys.gradRasterTime*size(gos,2);
    size(gos,2);
    % actually it is trickier than that: the (Siemens) interpreter sequence 
    % per default will try to split the trajectory into segments <=1000 samples
    % and every of these segments will have to have duration aligned to the
    % gradient raster time
    adcSamplesPerSegment=700; % you may need to play with this number to fill the entire trajectory
    adcSamplesDesired=kRadius*numPoints;
    adcSegments=round(adcSamplesDesired/adcSamplesPerSegment);
    adcSamples=adcSegments*adcSamplesPerSegment;
    adcDwell=round(adcTime/adcSamples/100e-9)*100e-9; % on Siemens adcDwell needs to be aligned to 100ns (if my memory serves me right)
    
    
%    fs = 25e4
%    adcDwell=round(1/fs/100e-9)*100e-9
%    adcSamples = floor(adcTime/adcDwell)
    
    % Calculate the number of segments based on the samples per segment
%    adcSegments = ceil(adcSamples / adcSamplesPerSegment);
    
    % Update the total number of samples to be an integer multiple of the samples per segment
%    adcSamples = adcSegments * adcSamplesPerSegment;

    adcSegmentDuration=adcSamplesPerSegment*adcDwell; % with the 100 samples above and the 100ns alignment we automatically fullfill the segment alignment requirement
    if mod(adcSegmentDuration, sys.gradRasterTime)>eps 
        error('ADC segmentation model results in incorrect segment duration');
    end
    % update segment count
    adcSegments=floor(adcTime/adcSegmentDuration);
    adcSamples=adcSegments*adcSamplesPerSegment;

    % Define gradients and ADC events for the current shot

    adc_events{shot} = mr.makeAdc(adcSamples, 'Dwell', sys.gradRasterTime, 'system', sys);

 %   gos = [gos gos(:, end)]

    gradient_waveforms{shot}.gx = mr.makeArbitraryGrad('x', gx_shot, 'first', 0, 'last', 0, 'system', sys, 'Delay',mr.calcDuration(gzReph));
    gradient_waveforms{shot}.gy = mr.makeArbitraryGrad('y', gy_shot, 'first', 0, 'last', 0, 'system', sys, 'Delay',mr.calcDuration(gzReph));
    %gradient_waveforms{shot}.gx = mr.makeArbitraryGrad('x', gx_shot, 'Delay',mr.calcDuration(gzReph), 'Sys', sys);
    %gradient_waveforms{shot}.gy = mr.makeArbitraryGrad('y', gy_shot, 'Delay',mr.calcDuration(gzReph), 'Sys', sys);


    % spoilers
    gz_spoil=mr.makeTrapezoid('z',sys,'Area',deltaK*Nx*4);
    gx_spoil=mr.makeExtendedTrapezoid('x','times',[0 mr.calcDuration(gz_spoil)],'amplitudes',[gos(1,end),0]); %todo: make a really good spoiler
    gy_spoil=mr.makeExtendedTrapezoid('y','times',[0 mr.calcDuration(gz_spoil)],'amplitudes',[gos(2,end),0]); %todo: make a really good spoiler

    spoiler_waveforms{shot}.gx_spoil = gx_spoil;
    spoiler_waveforms{shot}.gy_spoil = gy_spoil;

end



% because of the ADC alignment requirements the sampling window possibly
% extends past the end of the trajectory (these points will have to be
% discarded in the reconstruction, which is no problem). However, the
% ramp-down parts and the Z-spoiler now have to be added to the readout
% block otherwise there will be a gap inbetween
% gz_spoil.delay=mr.calcDuration(gx);
% gx_spoil.delay=gz_spoil.delay;
% gy_spoil.delay=gz_spoil.delay;
% gx_combined=mr.addGradients([gx,gx_spoil], lims);
% gy_combined=mr.addGradients([gy,gy_spoil], lims);
% gz_combined=mr.addGradients([gzReph,gz_spoil], lims);



for s = 1:Nslices
   % seq.addBlock(rf_fs, gz_fs); % Fat-saturation block

    % Apply RF pulse with frequency offset for slice selection
    rf.freqOffset = gz.amplitude * sliceThickness * (s - 1 - (Nslices - 1) / 2);


    for shot = 1:4
        seq.addBlock(rf, gz);
        % Apply the gradients and ADC events for the current shot

        seq.addBlock(mr.rotate('z',phi,gzReph,gradient_waveforms{shot}.gx,gradient_waveforms{shot}.gy,adc_events{shot}));
        seq.addBlock(mr.rotate('z',phi,spoiler_waveforms{shot}.gx_spoil,spoiler_waveforms{shot}.gy_spoil,gz_spoil));
    end
end

% check whether the timing of the sequence is correct
[ok, error_report]=seq.checkTiming;

if (ok)
    fprintf('Timing check passed successfully\n');
else
    fprintf('Timing check failed! Error listing follows:\n');
    fprintf([error_report{:}]);
    fprintf('\n');
end





%
seq.setDefinition('FOV', [fov fov sliceThickness]);
seq.setDefinition('Name', 'four_shot');
seq.setDefinition('MaxAdcSegmentLength', adcSamplesPerSegment); % this is important for making the sequence run automatically on siemens scanners without further parameter tweaking

seq.write('four_shot.seq');   % Output sequence for scanner

seq.write('outputs/test.seq');   % Output sequence for scanner
%nce is ready, so let's see what we got 
seq.plot();             % Plot sequence waveforms

%% k-space trajectory calculation
[ktraj_adc, t_adc, ktraj, t_ktraj, t_excitation, t_refocusing] = seq.calculateKspacePP();
seq
% plot k-spaces
figure; plot(t_ktraj, ktraj'); title('k-space components as functions of time'); % plot the entire k-space trajectory
figure; plot(ktraj(1,:),ktraj(2,:),'b'); % a 2D plot
hold;plot(ktraj_adc(1,:),ktraj_adc(2,:),'r.'); title('2D k-space');

% seq.install('siemens');
return

%% calculate PNS and check whether we are still within limits

[pns,tpns]=seq.calcPNS('idea/asc/MP_GPA_K2309_2250V_951A_AS82.asc'); % prisma