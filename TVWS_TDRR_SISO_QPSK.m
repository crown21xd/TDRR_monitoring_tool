
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% configurable parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Corresponds to Section II (System Model) of the paper: Defines the SISO-OFDM system parameters,
% including antenna counts, channel taps, modulation, and simulation settings for evaluating jammer impact.
num_receive_antennas = 1;                    % 8 par.B:BTS  number of receive antennas
num_transmit_antennas = 1;                    % 2 par.U:OPERATOR number of transmit antennas and data streams
num_channel_taps = 1;                        % 4 par.L: number of channel taps (multipath delay spread)
cyclic_prefix_length = 16;                   % par.P: cyclic prefix length (to combat inter-symbol interference)
num_subcarriers = 64;                        % par.N_sc: number of total subcarriers (FFT size for OFDM)
used_subcarrier_indices = [-26:-22, -20:-8, -6:-1, 1:6, 8:20, 22:26] + 32; % par.tonelist: list of used subcarriers (data-bearing tones), shifted by 32 for DC centering
num_used_subcarriers = length(used_subcarrier_indices); % par.num_tones: total number of active subcarriers
num_ofdm_symbols = 50;                       % par.num_ofdm_symbols: number of OFDM symbols sent per channel realization (temporal diversity)
jammer_power_dB = 30;                        % par.rho_dB: jammer strength relative to signal strength (in dB) - high value simulates strong jamming
modulation_scheme = 'QPSK';                  % par.mod: transmit constellation (modulation scheme: QPSK for 2 bits/symbol)
random_jammer_power_flag = 0;                % par.random_jammer_power: if true, jammer power varies randomly per trial; else fixed at par.rho_dB
    
num_trials = 100;                            % par.trials: number of Monte-Carlo trials (different channel realizations for statistical averaging)
snr_dB_list = -5:1:50;                       % par.SNRdB_list: list of SNR values (in dB) to simulate (signal-to-noise ratio sweep)

plot_flag = 1;                              % par.plot: flag to enable/disable plotting of simulation results
print_messages_flag = 1;                     % par.printmsg: flag to enable/disable printing of progress messages during simulation

rng(0);                                     % set random seed for reproducibility of results

%% 

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% simulation start
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Corresponds to Section II (System Model): Modulation constellation definition and bit mapping for signal generation.

% create a unique simulation name based on parameters for logging or saving results
%summarzes the given paramenters
simulation_name = [num2str(num_receive_antennas) 'x' num2str(num_transmit_antennas) '_' modulation_scheme '_L', num2str(num_channel_taps) '_P' num2str(cyclic_prefix_length) ...
  '_Nsc' num2str(num_subcarriers) '_Nofd' num2str(num_ofdm_symbols) '_rho' num2str(jammer_power_dB) ...
  '_' num2str(num_trials) 'Trials']; % par.simName
fprintf(simulation_name)

% define the modulation constellation based on the chosen modulation scheme
switch (modulation_scheme)
  case 'BPSK'  % Binary Phase Shift Keying: 2 symbols, 1 bit each
    constellation_symbols = [ -1 1 ]; % par.symbols
  case 'QPSK'  % Quadrature Phase Shift Keying: 4 symbols, 2 bits each
    constellation_symbols = [ ... % par.symbols
      -1-1i,-1+1i, ...
      +1-1i,+1+1i ];
  case '16QAM'  % 16-Quadrature Amplitude Modulation: 16 symbols, 4 bits each
    constellation_symbols = [ ... % par.symbols
      -3-3i,-3-1i,-3+3i,-3+1i, ...
      -1-3i,-1-1i,-1+3i,-1+1i, ...
      +3-3i,+3-1i,+3+3i,+3+1i, ...
      +1-3i,+1-1i,+1+3i,+1+1i ];
  case '64QAM'  % 64-Quadrature Amplitude Modulation: 64 symbols, 6 bits each
    constellation_symbols = [ ... % par.symbols
      -7-7i,-7-5i,-7-1i,-7-3i,-7+7i,-7+5i,-7+1i,-7+3i, ...
      -5-7i,-5-5i,-5-1i,-5-3i,-5+7i,-5+5i,-5+1i,-5+3i, ...
      -1-7i,-1-5i,-1-1i,-1-3i,-1+7i,-1+5i,-1+1i,-1+3i, ...
      -3-7i,-3-5i,-3-1i,-3-3i,-3+7i,-3+5i,-3+1i,-3+3i, ...
      +7-7i,+7-5i,+7-1i,+7-3i,+7+7i,+7+5i,+7+1i,+7+3i, ...
      +5-7i,+5-5i,+5-1i,+5-3i,+5+7i,+5+5i,+5+1i,+5+3i, ...
      +1-7i,+1-5i,+1-1i,+1-3i,+1+7i,+1+5i,+1+1i,+1+3i, ...
      +3-7i,+3-5i,+3-1i,-3-3i,+3+7i,+3+5i,+3+1i,+3+3i ];
end

% normalize the constellation to have unit average energy for fair comparison
constellation_energy = mean(abs(constellation_symbols).^2);  % par.constellation_energy: compute average symbol energy
constellation_symbols = constellation_symbols / sqrt(constellation_energy);  % normalize symbols
normalized_energy = mean(abs(constellation_symbols).^2);  % par.Es: verify normalized energy is 1

% precompute bit-to-symbol mapping for efficient decoding
bits_per_symbol = log2(length(constellation_symbols));  % par.Q: bits per symbol (e.g., 2 for QPSK)
bit_mappings = de2bi(0:length(constellation_symbols)-1, bits_per_symbol, 'left-msb');  % par.bits: binary representations of symbol indices

% track simulation time
time_elapsed = 0;

% intermediate results
y_jammerless_frequency = zeros(num_receive_antennas,num_subcarriers,num_ofdm_symbols); % y_jl_f
y_jammer_compliant_frequency = zeros(num_receive_antennas,num_subcarriers,num_ofdm_symbols); % y_oj_comp_f
y_jammer_noncompliant_frequency = zeros(num_receive_antennas,num_subcarriers,num_ofdm_symbols); % y_oj_noncomp_f
y_compliant_frequency = zeros(num_receive_antennas,num_subcarriers,num_ofdm_symbols); % y_comp_f
y_noncompliant_frequency = zeros(num_receive_antennas,num_subcarriers,num_ofdm_symbols); % y_noncomp_f


% Initialize accumulators for averaging for received signal
avg_y_jammerless = zeros(num_receive_antennas, num_used_subcarriers, num_ofdm_symbols);
avg_y_jammer_compliant = zeros(num_receive_antennas, num_used_subcarriers, num_ofdm_symbols);
avg_y_jammer_noncompliant = zeros(num_receive_antennas, num_used_subcarriers, num_ofdm_symbols);


% results
jammer_compliant_freq_domain = zeros(num_trials, num_used_subcarriers,num_receive_antennas, num_ofdm_symbols); % Y_oj_comp
jammer_noncompliant_freq_domain = zeros(num_trials, num_used_subcarriers,num_receive_antennas, num_ofdm_symbols); % Y_oj_noncomp

mse_jammerless = zeros(1, length(snr_dB_list)); % MSE_jl
mse_compliant = zeros(1, length(snr_dB_list)); % MSE_comp
mse_noncompliant = zeros(1, length(snr_dB_list)); % MSE_noncomp

results.BER_jammerless = zeros(1, length(snr_dB_list)); % res.BER_jl
results.BER_compliant = zeros(1, length(snr_dB_list)); % res.BER_comp
results.BER_noncompliant = zeros(1, length(snr_dB_list)); % res.BER_noncomp
results.BER_projections = zeros(num_channel_taps, length(snr_dB_list)); % res.BER_projs

results.spatial_interference_distribution = zeros(num_receive_antennas,num_used_subcarriers);

% Initialize variables for per-trial errors to compute success rate of mitigating illegal users/channel taps
errors_jammerless_per_trial = zeros(num_trials, length(snr_dB_list)); % errors_jl_per_trial
errors_compliant_per_trial = zeros(num_trials, length(snr_dB_list)); % errors_comp_per_trial
errors_noncompliant_per_trial = zeros(num_trials, length(snr_dB_list)); % errors_noncomp_per_trial

% generate bit stream
bits = randi([0 1],num_trials,num_ofdm_symbols,num_transmit_antennas,bits_per_symbol,num_used_subcarriers); % bits
indices = zeros(num_trials, num_ofdm_symbols, num_transmit_antennas,num_used_subcarriers); % idx

%% 

% trials loop
% Corresponds to Section III (Jammer Modeling and Interference Analysis): Monte-Carlo simulation loop for channel realizations and jammer signal generation.
tic
for trial_index=1:num_trials

  

  % generate jammer's channel in the time domain
  % Corresponds to Section II (System Model): Rayleigh fading channel generation.
  jammer_channel_time = sqrt(0.5)*randn(num_receive_antennas,num_channel_taps) + 1j*sqrt(0.5)*randn(num_receive_antennas,num_channel_taps);  % j_
  jammer_channel_freq = (1/sqrt(num_subcarriers))*fft(jammer_channel_time,num_subcarriers,2);  % j_f

  % generate the UEs' channel in the time domain
  user_channel_time = sqrt(0.5)*randn(num_receive_antennas,num_transmit_antennas,num_channel_taps) + 1j*sqrt(0.5)*randn(num_receive_antennas,num_transmit_antennas,num_channel_taps);  % H_
  user_channel_freq = (1/sqrt(num_subcarriers))*fft(user_channel_time,num_subcarriers,3); % H_f
  % E[||H_f(:,:,i)||_F^2] = B*U*L/N_sc for all i=1,...,N_sc
  % expected receive signal energy per subcarrier is E_s*B*U*L/N_sc = B*U*L/N_sc

  % set jammer power: either random per trial or fixed
  if random_jammer_power_flag
    jammer_power_linear = num_transmit_antennas*10^(jammer_power_dB*rand()/10);  % par.rho
  else
    jammer_power_linear = num_transmit_antennas*10^(jammer_power_dB/10);  % par.rho
  end



  % loop over SNR values to simulate different noise levels
  for snr_index = 1:length(snr_dB_list)
    
    noise_variance = num_channel_taps*num_transmit_antennas*10^(-snr_dB_list(snr_index)/10);  % par.N0
    
    % loop over OFDM symbols for temporal averaging
    for ofdm_symbol_index=1:num_ofdm_symbols
  
      %%% generate transmit signals
      % Corresponds to Section III (Jammer Modeling): Gener ation of compliant and noncompliant jammer signals, and UE signals.
        
      % generate jammer's transmit signal in frequency domain (OFDM-compliant)
      jammer_signal_freq = sqrt(0.5)*randn(1,num_subcarriers)+ 1j*sqrt(0.5)*randn(1,num_subcarriers);  % generates random bits in frequency domain, Gaussian noise
      jammer_signal_time_no_cp = sqrt(num_subcarriers)*ifft(jammer_signal_freq);% converts frequency to time domain, ifft (Inverse Fast Fourier Transform)
      jammer_signal_time = [jammer_signal_time_no_cp(:,num_subcarriers-cyclic_prefix_length+1:num_subcarriers), jammer_signal_time_no_cp];  % add cp
      jammer_signal_time = sqrt(jammer_power_linear)*jammer_signal_time;  % scale jammer power compared to jammerless
  
      % generate jammer's transmit signal in time domain (OFDM-noncompliant)
      jammer_signal_time_noncomp = sqrt(0.5)*randn(1,num_subcarriers+cyclic_prefix_length)+ 1j*sqrt(0.5)*randn(1,num_subcarriers+cyclic_prefix_length); % generates random noise in time  domain
      jammer_signal_time_noncomp = sqrt(jammer_power_linear)*jammer_signal_time_noncomp;  % scale jammer power compared to jammerless       

      % OFDM-compliant TVWS operator: Generation of legitimate user signal in TVWS using OFDM
      % generate UEs' transmit signal in frequency domain (OFDM-compliant)
       for tone_index=1:num_used_subcarriers
            bits_slice = squeeze(bits(trial_index, ofdm_symbol_index, :, :, tone_index)); % size [2 1]
            bits_slice = bits_slice.'; % transpose to [1 2]
            indices(trial_index, ofdm_symbol_index, 1, tone_index) = bi2de(bits_slice, 'left-msb') + 1;
       end
      symbols_used_tones = constellation_symbols(squeeze(indices(trial_index,ofdm_symbol_index,:,:))); % S_f_used_tones
      symbols_freq = zeros(num_transmit_antennas, num_subcarriers);  % S_f
      % Subcarrier assignment: Assign the modulated symbols to the designated used subcarriers in the frequency domain
      symbols_freq(:,used_subcarrier_indices) = symbols_used_tones; % assign symbols only to used subcarriers
      symbols_time_no_cp = sqrt(num_subcarriers)*ifft(symbols_freq,num_subcarriers,2);  % S_no_cp
      symbols_time = [symbols_time_no_cp(:,num_subcarriers-cyclic_prefix_length+1:num_subcarriers), symbols_time_no_cp];  % S

    
      %%% generate receive signals
      % Threat DB 
      receive_signal_legitimate_TVWSO = zeros(num_receive_antennas, num_subcarriers+cyclic_prefix_length+num_channel_taps-1);% y_jl 8 x 64+16+4-1
      receive_signal_jammer_compliant = zeros(num_receive_antennas, num_subcarriers+cyclic_prefix_length+num_channel_taps-1);     % y_oj_comp
      receive_signal_jammer_noncompliant = zeros(num_receive_antennas, num_subcarriers+cyclic_prefix_length+num_channel_taps-1);  % y_oj_noncomp
      receive_signal_compliant = zeros(num_receive_antennas, num_subcarriers+cyclic_prefix_length+num_channel_taps-1);        % y_comp
      receive_signal_noncompliant = zeros(num_receive_antennas, num_subcarriers+cyclic_prefix_length+num_channel_taps-1);     % y_noncomp
    
      

      % convolve channels with transmit signals to simulate received signals
      for receive_ant_index=1:num_receive_antennas
        for transmit_ant_index=1:num_transmit_antennas %1->2
          receive_signal_legitimate_TVWSO(receive_ant_index,:) = receive_signal_legitimate_TVWSO(receive_ant_index,:) + conv(squeeze(user_channel_time(receive_ant_index,transmit_ant_index,:)),symbols_time(transmit_ant_index,:));  % sum of convolutions for each UE 
                %receive_signal_legitimate_TVWSO(1=>8,83) = 
        end
        receive_signal_jammer_compliant(receive_ant_index,:) = conv(jammer_channel_time(receive_ant_index,:),jammer_signal_time);  % jammer convolution compliant
        receive_signal_jammer_noncompliant(receive_ant_index,:) = conv(jammer_channel_time(receive_ant_index,:),jammer_signal_time_noncomp);  % jammer convolution noncompliant
      end
      receive_signal_compliant = receive_signal_legitimate_TVWSO + receive_signal_jammer_compliant;  % total received signal with compliant jammer
      receive_signal_noncompliant = receive_signal_legitimate_TVWSO + receive_signal_jammer_noncompliant;  % total received signal with noncompliant jammer
        
      


      %%% generate noise in the time domain
      noise_time = sqrt(0.5)*randn(num_receive_antennas,num_subcarriers)+ 1j*sqrt(0.5)*randn(num_receive_antennas,num_subcarriers);  % N
      noise_time = sqrt(noise_variance)*noise_time;  % scale noise by noise variance

      %%% truncate signals in time-domain to remove cyclic prefix and channel reverberation
      receive_signal_legitimate_TVWSO_truncated = receive_signal_legitimate_TVWSO(:,cyclic_prefix_length+1:cyclic_prefix_length+num_subcarriers) + noise_time;  % jammerless signal + noise
      receive_signal_jammer_compliant_truncated = receive_signal_jammer_compliant(:,cyclic_prefix_length+1:cyclic_prefix_length+num_subcarriers); % jammer only, no noise
      receive_signal_jammer_noncompliant_truncated = receive_signal_jammer_noncompliant(:,cyclic_prefix_length+1:cyclic_prefix_length+num_subcarriers); % jammer only, no noise
      receive_signal_compliant_truncated = receive_signal_compliant(:,cyclic_prefix_length+1:cyclic_prefix_length+num_subcarriers) + noise_time;  % full signal + noise
      receive_signal_noncompliant_truncated = receive_signal_noncompliant(:,cyclic_prefix_length+1:cyclic_prefix_length+num_subcarriers) + noise_time;  % full signal + noise

      




      %%% convert to frequency domain
      y_jammerless_frequency(:,:,ofdm_symbol_index) = (1/sqrt(num_subcarriers))*fft(receive_signal_legitimate_TVWSO_truncated, num_subcarriers, 2);  % jammerless frequency domain
      y_jammer_compliant_frequency(:,:,ofdm_symbol_index) = (1/sqrt(num_subcarriers))*fft(receive_signal_jammer_compliant_truncated, num_subcarriers, 2);  % jammer compliant freq domain
      y_jammer_noncompliant_frequency(:,:,ofdm_symbol_index) = (1/sqrt(num_subcarriers))*fft(receive_signal_jammer_noncompliant_truncated, num_subcarriers, 2);  % jammer noncompliant freq domain
      y_compliant_frequency(:,:,ofdm_symbol_index) = (1/sqrt(num_subcarriers))*fft(receive_signal_compliant_truncated, num_subcarriers, 2);  % full compliant freq domain
      y_noncompliant_frequency(:,:,ofdm_symbol_index) = (1/sqrt(num_subcarriers))*fft(receive_signal_noncompliant_truncated, num_subcarriers, 2);  % full noncompliant freq domain


      %%% remove subcarriers that are not used
      y_jammerless_frequency_tones = y_jammerless_frequency(:,used_subcarrier_indices,:);  % select used subcarriers only
      y_jammer_compliant_frequency_tones = y_jammer_compliant_frequency(:,used_subcarrier_indices,:);
      y_jammer_noncompliant_frequency_tones = y_jammer_noncompliant_frequency(:,used_subcarrier_indices,:);
      y_compliant_frequency_tones = y_compliant_frequency(:,used_subcarrier_indices,:);
      y_noncompliant_frequency_tones = y_noncompliant_frequency(:,used_subcarrier_indices,:);


      user_channel_frequency_tones = user_channel_freq(:,:,used_subcarrier_indices);  % channel frequency response for used subcarriers
    
      %%% collect jammer statistics for later analysis
      for tone_index=1:num_used_subcarriers
        jammer_compliant_freq_domain(trial_index,tone_index,:,ofdm_symbol_index) = y_jammer_compliant_frequency_tones(:,tone_index,ofdm_symbol_index);  % store jammer compliant freq domain
        jammer_noncompliant_freq_domain(trial_index,tone_index,:,ofdm_symbol_index) = y_jammer_noncompliant_frequency_tones(:,tone_index,ofdm_symbol_index);  % store jammer noncompliant freq domain   
      end

    end % iterate over OFDM symbols with same channel realization


    for tone_index=1:num_used_subcarriers
  
      %%% extract jammer statistics
      % Corresponds to Section IV (Mitigation): SVD for subspace estimation.
      % SVD: What spatial direction the jammer energy points to (U)? & How strong it is in each direction (S)?
      [U_compliant,S_compliant,~] = svd(squeeze(jammer_compliant_freq_domain(trial_index,tone_index,:,:)));  % singular value decomposition of jammer compliant
      [U_noncompliant,S_noncompliant,~] = svd(squeeze(jammer_noncompliant_freq_domain(trial_index,tone_index,:,:)));  % SVD of jammer noncompliant
      singular_values = diag(S_noncompliant);
      dominant_sv = singular_values(1);  % largest singular value
      results.spatial_interference_distribution(:, tone_index) = results.spatial_interference_distribution(:, tone_index) + dominant_sv;

      % form best rank-(num_receive_antennas-1) projectors: those that remove most of the
      % jammer interference per subcarrier
      % Corresponds to Section IV (Mitigation): Nulling matrices.
      jammer_estimate_compliant = U_compliant(:,1);  % dominant jammer subspace vector compliant
      % Mitigation: Projection matrix to nullify the compliant jammer subspace, mitigating its effect on TVWS operator
      projection_matrix_compliant = eye(num_receive_antennas) - jammer_estimate_compliant*jammer_estimate_compliant';  % projection matrix to null jammer compliant
      jammer_estimate_noncompliant = U_noncompliant(:,1);  % dominant jammer subspace vector noncompliant
      projection_matrix_noncompliant = eye(num_receive_antennas) - jammer_estimate_noncompliant*jammer_estimate_noncompliant';  % projection matrix to null jammer noncompliant

      % projectors = zeros(num_channel_taps,num_receive_antennas,num_receive_antennas);  % initialize multiple projectors for different nulling dimensions
      % for b=1:num_channel_taps
      %   projectors(b,:,:) = eye(num_receive_antennas) - U_noncompliant(:,1:b)*U_noncompliant(:,1:b)';  % projector nulling b jammer subspace vectors
      % end
      projectors = repmat(eye(num_receive_antennas), [num_channel_taps, 1, 1]); %no-op projection in SISO
  
      %%% detect signals using zero-forcing (and, in some cases, jammer mitigation)
      % Corresponds to Section IV (Mitigation): Zero-forcing detection with/without nulling., recover the transmitted symbol    
      % symbols_estimate_jammerless = (sqrt(num_subcarriers)*(user_channel_frequency_tones(:,:,tone_index)))\squeeze(y_jammerless_frequency_tones(:,tone_index,:));  % zero-forcing estimate jammerles
      symbols_estimate_jammerless = y_jammerless_frequency_tones(:,tone_index,:) ./ user_channel_frequency_tones(:,:,tone_index);

      % Mitigated TVWS operator symbol estimation using projection matrix to suppress compliant jammer interference
      symbols_estimate_compliant = (sqrt(num_subcarriers)*(projection_matrix_compliant*squeeze(user_channel_frequency_tones(:,:,tone_index))))\(projection_matrix_compliant*squeeze(y_compliant_frequency_tones(:,tone_index,:)));  % zero-forcing with jammer mitigation compliant
      symbols_estimate_noncompliant = (sqrt(num_subcarriers)*(projection_matrix_noncompliant*squeeze(user_channel_frequency_tones(:,:,tone_index))))\(projection_matrix_noncompliant*squeeze(y_noncompliant_frequency_tones(:,tone_index,:)));  % zero-forcing with jammer mitigation noncompliant
      symbols_estimate_projections = zeros(num_receive_antennas,num_transmit_antennas,num_ofdm_symbols);  % initialize estimates for multiple projectors
      for b=1:num_channel_taps
        symbols_estimate_projections(b,:,:) = (sqrt(num_subcarriers)*(squeeze(projectors(b,:,:))*squeeze(user_channel_frequency_tones(:,:,tone_index))))\(squeeze(projectors(b,:,:))*squeeze(y_noncompliant_frequency_tones(:,tone_index,:)));  % zero-forcing with multiple nulling projectors
      end
        
           % -- compute bit outputs by minimum distance decoding
      % Corresponds to Section IV (Mitigation): Decoding.
      symbols_estimate_jammerless_vector = reshape(symbols_estimate_jammerless, [num_ofdm_symbols*num_transmit_antennas,1]);  % vectorize estimated symbols
      [~,indices_hat_vector] = min(abs(symbols_estimate_jammerless_vector*ones(1,length(constellation_symbols))-ones(num_ofdm_symbols*num_transmit_antennas,1)*constellation_symbols).^2,[],2);  % find closest constellation points
      indices_hat_jammerless = reshape(indices_hat_vector, [num_transmit_antennas, num_ofdm_symbols]);  % reshape indices
      bits_hat_jammerless = bit_mappings(indices_hat_jammerless,:);  % map indices to bits

      symbols_estimate_compliant_vector = reshape(symbols_estimate_compliant, [num_ofdm_symbols*num_transmit_antennas,1]);
      [~,indices_hat_vector] = min(abs(symbols_estimate_compliant_vector*ones(1,length(constellation_symbols))-ones(num_ofdm_symbols*num_transmit_antennas,1)*constellation_symbols).^2,[],2);
      indices_hat_compliant = reshape(indices_hat_vector, [num_transmit_antennas, num_ofdm_symbols]);
      bits_hat_compliant = bit_mappings(indices_hat_compliant,:);

      symbols_estimate_noncompliant_vector = reshape(symbols_estimate_noncompliant, [num_ofdm_symbols*num_transmit_antennas,1]);
      [~,indices_hat_vector] = min(abs(symbols_estimate_noncompliant_vector*ones(1,length(constellation_symbols))-ones(num_ofdm_symbols*num_transmit_antennas,1)*constellation_symbols).^2,[],2);
      indices_hat_noncompliant = reshape(indices_hat_vector, [num_transmit_antennas, num_ofdm_symbols]);
      bits_hat_noncompliant = bit_mappings(indices_hat_noncompliant,:);

      bits_hat_projections = zeros(num_channel_taps,num_transmit_antennas*num_ofdm_symbols,bits_per_symbol);
      for b=1:num_channel_taps
        symbols_estimate_vector = reshape(symbols_estimate_projections(b,:,:), [num_ofdm_symbols*num_transmit_antennas,1]);
        [~,indices_hat_vector] = min(abs(symbols_estimate_vector*ones(1,length(constellation_symbols))-ones(num_ofdm_symbols*num_transmit_antennas,1)*constellation_symbols).^2,[],2);
        indices_hat = reshape(indices_hat_vector, [num_transmit_antennas, num_ofdm_symbols]);
        bits_hat_projections(b,:,:) = bit_mappings(indices_hat,:);
      end
      

      %gian marti BER
      % -- compute error metrics by comparing estimated bits to true bits
      % Corresponds to Section IV (Mitigation): BER calculation.
      bit_tensor = squeeze(bits(trial_index,:,:,:,tone_index));  % extract bits for current trial and subcarrier
      bit_tensor = permute(bit_tensor,[3 2 1]);  % reorder dimensions for comparison
      true_bits = bit_tensor(:,:)';  % flatten to 2D matrix
      bits_hat_jammerless = permute(bits_hat_jammerless, [3 2 1]);
        bits_hat_jammerless = bits_hat_jammerless(:,:)';
        
        bits_hat_compliant = permute(bits_hat_compliant, [3 2 1]);
        bits_hat_compliant = bits_hat_compliant(:,:)';
        
        bits_hat_noncompliant = permute(bits_hat_noncompliant, [3 2 1]);
        bits_hat_noncompliant = bits_hat_noncompliant(:,:)';

      results.BER_jammerless(snr_index) = results.BER_jammerless(snr_index) + sum(sum(true_bits~=bits_hat_jammerless));  % count bit errors jammerless
      results.BER_compliant(snr_index) = results.BER_compliant(snr_index) + sum(sum(true_bits~=bits_hat_compliant));  % count bit errors jammer compliant
      results.BER_noncompliant(snr_index) = results.BER_noncompliant(snr_index) + sum(sum(true_bits~=bits_hat_noncompliant));  % count bit errors jammer noncompliant
      for b=1:num_channel_taps
          bits_hat_proj = permute(squeeze(bits_hat_projections(b,:,:,:)), [3 2 1]);
          bits_hat_projections = bits_hat_proj(:,:)';
        results.BER_projections(b,snr_index) = results.BER_projections(b,snr_index) + sum(sum(true_bits~=squeeze(bits_hat_projections(b,:,:))));  % count bit errors for each nulling dimension
      end

    end
%% 

  end % iterate over SNRs

  % -- keep track of simulation time and print estimated remaining time
  if(print_messages_flag)
    if toc>10
      time = toc;
      time_elapsed = time_elapsed + time;
        fprintf('estimated remaining simulation time: %3.0f min.\n', ...
          time_elapsed*(num_trials/trial_index-1)/60);
      tic
    end
  end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         % 
         if mod(trial_index, 1) == 0 % Update plot every 5 trials

            % Plot frequency domain used subcarriers signals (magnitude) for antenna 2 and last OFDM symbol
            figure(11)
            clf
            antenna_idx = 1;
            ofdm_symbol_idx = num_ofdm_symbols; % last OFDM symbol

            subplot(5,1,1)
            plot(used_subcarrier_indices, real(squeeze(y_jammerless_frequency_tones(antenna_idx,:,ofdm_symbol_idx))))
            ylim([-100 100])
            title('Legitimate TVWS Frequency Tones')
            xlabel('Subcarrier Index')
            ylabel('Magnitude')

            subplot(5,1,2)
            plot(used_subcarrier_indices, real(squeeze(y_jammer_compliant_frequency_tones(antenna_idx,:,ofdm_symbol_idx))))
            ylim([-100 100])
            title('OFDM-compliant Illegitimate TVWS Frequency Tones ')
            xlabel('Subcarrier Index')
            ylabel('Magnitude')

            subplot(5,1,3)
            plot(used_subcarrier_indices, real(squeeze(y_jammer_noncompliant_frequency_tones(antenna_idx,:,ofdm_symbol_idx))))
            ylim([-100 100])
            title('Cyclic prefix-violating Illegitimate TVWS Frequency Tones ')
            xlabel('Subcarrier Index')
            ylabel('Magnitude')

            sgtitle(['Received Signals at Used Subcarriers ( OFDM Symbol ', num2str(ofdm_symbol_idx), ') at Trial ', num2str(trial_index)]) % Add a super title with trial number
            drawnow % Force update of the figure for live view

         end

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % Accumulate signals across trials
    avg_y_jammerless = avg_y_jammerless + y_jammerless_frequency_tones;
    avg_y_jammer_compliant = avg_y_jammer_compliant + y_jammer_compliant_frequency_tones;
    avg_y_jammer_noncompliant = avg_y_jammer_noncompliant + y_jammer_noncompliant_frequency_tones;

end % iterate over MC trials


% normalize results
% Corresponds to Section V (Results): Averaging BER over trials and subcarriers.
results.BER_jammerless = results.BER_jammerless/(num_trials*num_transmit_antennas*num_ofdm_symbols*bits_per_symbol*num_used_subcarriers);
results.BER_compliant = results.BER_compliant/(num_trials*num_transmit_antennas*num_ofdm_symbols*bits_per_symbol*num_used_subcarriers);
results.BER_noncompliant = results.BER_noncompliant/(num_trials*num_transmit_antennas*num_ofdm_symbols*bits_per_symbol*num_used_subcarriers);
results.BER_projections = results.BER_projections/(num_trials*num_transmit_antennas*num_ofdm_symbols*bits_per_symbol*num_used_subcarriers);

results.spatial_interference_distribution = results.spatial_interference_distribution./sum(results.spatial_interference_distribution,1);
results.spatial_interference_distribution_sc_std = std(results.spatial_interference_distribution.');
results.spatial_interference_distribution_sc_avg = mean(results.spatial_interference_distribution.');

avg_y_jammerless = avg_y_jammerless / num_trials;
avg_y_jammer_compliant = avg_y_jammer_compliant / num_trials;
avg_y_jammer_noncompliant = avg_y_jammer_noncompliant / num_trials;


%%% plot results
% Corresponds to Section V (Results): Visualization of BER performance and interference distribution.
if plot_flag
   
  figure(1)
  semilogy(snr_dB_list, results.BER_jammerless, 'LineWidth', 2.5)
  hold off
  grid on
  axis([min(snr_dB_list) max(snr_dB_list) 1e-4 1])
  xlabel('average SNR per receive antenna [dB]','FontSize',12)
  ylabel('uncoded bit error-rate (BER)','FontSize',12)
  legend("Legitimate TVWS",'FontSize',12,'Interpreter','none')
  set(gca,'FontSize',12)
  set(gcf,'position',[10,10,400,300])

  figure(2)
  semilogy(snr_dB_list, results.BER_jammerless, 'LineWidth', 2.5)
  hold on
  semilogy(snr_dB_list, results.BER_compliant, 'LineWidth', 2.5)
  hold off
  grid on
  axis([min(snr_dB_list) max(snr_dB_list) 1e-4 1])
  xlabel('average SNR per receive antenna [dB]','FontSize',12)
  ylabel('uncoded bit error-rate (BER)','FontSize',12)
  legend(["Legitimate TVWS", "OFDM-compliant Illegitimate TVWS"],'FontSize',12,'Interpreter','none')
  set(gca,'FontSize',12)
  set(gcf,'position',[10,10,400,300])

  figure(3)
  semilogy(snr_dB_list, results.BER_jammerless, 'LineWidth', 2.5)
  hold on
  semilogy(snr_dB_list, results.BER_noncompliant, 'LineWidth', 2.5)
  hold off
  grid on
  axis([min(snr_dB_list) max(snr_dB_list) 1e-4 1])
  xlabel('average SNR per receive antenna [dB]','FontSize',12)
  ylabel('uncoded bit error-rate (BER)','FontSize',12)
  legend(["Legitimate TVWS", "Cyclic prefix-violating Illegitimate TVWS"],'FontSize',12,'Interpreter','none')
  set(gca,'FontSize',12)
  set(gcf,'position',[10,10,400,300])

  figure(4)
  bar(1:num_receive_antennas,results.spatial_interference_distribution_sc_avg)
  hold on
  er = errorbar(1:num_receive_antennas,results.spatial_interference_distribution_sc_avg, ...
    2*results.spatial_interference_distribution_sc_std, ...
    2*results.spatial_interference_distribution_sc_std);
  er.Color = [0 0 0];
  er.LineStyle = 'none';
  hold off
  axis([0.25 num_receive_antennas+0.75 0 1])
  grid on
  xlabel('ordered spatial dimension index')
  ylabel('fraction of receive interference')
  set(gcf,'position',[10,10,400,300])

  figure(5)
    clf
    
    % Plot jammerless baseline
    semilogy(snr_dB_list, results.BER_jammerless, 'k--', 'LineWidth', 2.5)
    hold on
    
    % Plot projected cases (CP-violating jammer with nulling)
    colors = lines(num_channel_taps); % MATLAB color map for consistency
    for b = 1:num_channel_taps
    semilogy(snr_dB_list, results.BER_projections(b,:), 'LineWidth', 2.5, 'Color', colors(b,:))
    end
    
    grid on
    axis([min(snr_dB_list) max(snr_dB_list) 1e-4 1])
    xlabel('Average SNR per receive antenna [dB]', 'FontSize', 12)
    ylabel('Uncoded Bit Error Rate (BER)', 'FontSize', 12)
    
    % ----- Legend setup -----
    numbers = cell(1, num_channel_taps + 1);
    numbers{1} = "Legitimate TVWS (Jammerless)";
    for b = 1:num_channel_taps
    numbers{b+1} = "CP-violating Jammer: Null " + num2str(b) + " Dimension(s)";
    end
    legend(numbers, 'FontSize', 11, 'Interpreter', 'none', 'Location', 'southwest')
    
    set(gca, 'FontSize', 12)
    set(gcf, 'Position', [10, 10, 600, 400])
    title('BER vs. SNR for Legitimate and Projected CP-violating Jammer Cases')

    % Plot jammer statistics
    figure(6)
    % Average magnitude over trials and OFDM symbols for compliant jammer
    jammer_compliant_avg = mean(mean(abs(jammer_compliant_freq_domain),4),1); % Average over trials and OFDM symbols
    subplot(2,1,1)
    imagesc(squeeze(jammer_compliant_avg))
    title('Jammer Compliant Frequency Domain Magnitude (Avg over Trials and OFDM Symbols)')
    xlabel('Receive Antenna')
    ylabel('Subcarrier Index')
    colorbar
    jammer_noncompliant_avg = mean(mean(abs(jammer_noncompliant_freq_domain),4),1); % Average over trials and OFDM symbols
    subplot(2,1,2)
    imagesc(squeeze(jammer_noncompliant_avg))
    title('Jammer Noncompliant Frequency Domain Magnitude (Avg over Trials and OFDM Symbols)')
    xlabel('Receive Antenna')
    ylabel('Subcarrier Index')
    colorbar
    set(gcf,'position',[10,10,600,400])
    
    % Plot frequency domain average signals for last OFDM symbol
    antenna_idx = 1;
    ofdm_symbol_idx = num_ofdm_symbols; % last OFDM symbol
    figure(7)
    subplot(3,1,1)
    plot(used_subcarrier_indices, real(squeeze(avg_y_jammerless(antenna_idx,:,ofdm_symbol_idx))))
    ylim([-100 100])
    title('Average Legitimate TVWS')
    xlabel('Subcarrier Index')
    ylabel('Magnitude')
    subplot(3,1,2)
    plot(used_subcarrier_indices, real(squeeze(avg_y_jammer_compliant(antenna_idx,:,ofdm_symbol_idx))))
    ylim([-100 100])
    title('Average OFDM-compliant Illegitimate TVWS Frequency Tones')
    xlabel('Subcarrier Index')
    ylabel('Magnitude')
    subplot(3,1,3)
    plot(used_subcarrier_indices, real(squeeze(avg_y_jammer_noncompliant(antenna_idx,:,ofdm_symbol_idx))))
    ylim([-100 100])
    title('Average Cyclic prefix-violating Illegitimate TVWS Frequency Tones')
    xlabel('Subcarrier Index')
    ylabel('Magnitude')
    % Super title for the figure
    sgtitle(['Average Received Signals at Used Subcarriers (OFDM Symbol ', num2str(ofdm_symbol_idx), ') Across ', num2str(num_trials), ' Trials'])
    
    % Figure 1: Average Legitimate TVWS
    figure(8)
    clf
    plot(used_subcarrier_indices, real(squeeze(avg_y_jammerless(antenna_idx,:,ofdm_symbol_idx))))
    ylim([-100 100])
    %title('Average Legitimate TVWS')
    xlabel('Subcarrier Index')
    ylabel('Magnitude')
    sgtitle(['Average Received Signals (Legitimate TVWS) at Used Subcarriers (OFDM Symbol ', num2str(ofdm_symbol_idx), ') Across ', num2str(num_trials), ' Trials'])
    
    % Figure 2: Average OFDM-compliant Illegitimate TVWS
    figure(9)
    clf
    plot(used_subcarrier_indices, real(squeeze(avg_y_jammer_compliant(antenna_idx,:,ofdm_symbol_idx))))
    ylim([-100 100])
    %title('Average OFDM-compliant Illegitimate TVWS Frequency Tones')
    xlabel('Subcarrier Index')
    ylabel('Magnitude')
    sgtitle(['Average Received Signals (OFDM-compliant Illegitimate TVWS)at Used Subcarriers (OFDM Symbol ', num2str(ofdm_symbol_idx), ') Across ', num2str(num_trials), ' Trials'])
    
    % Figure 3: Average Cyclic prefix-violating Illegitimate TVWS
    figure(10)
    clf
    plot(used_subcarrier_indices, real(squeeze(avg_y_jammer_noncompliant(antenna_idx,:,ofdm_symbol_idx))))
    ylim([-100 100])
    %title('Average Cyclic prefix-violating Illegitimate TVWS Frequency Tones')
    xlabel('Subcarrier Index')
    ylabel('Magnitude')
    sgtitle(['Average Received Signal (Cyclic prefix-violating Illegitimate TVWS) at Used Subcarriers (OFDM Symbol ', num2str(ofdm_symbol_idx), ') Across ', num2str(num_trials), ' Trials'])


current_dir = pwd;
figures_dir = fullfile(current_dir, 'figures'); 
  
  % Create folder if it doesn't exist
  if ~exist(figures_dir, 'dir')
      mkdir(figures_dir);
  end

  % Find all figure handles
  figHandles = findall(0, 'Type', 'figure');

   for i = 1:length(figHandles)
      figNum = figHandles(i).Number;
      filename = fullfile(figures_dir, sprintf('%d.png', figNum));
      fprintf('Saving Figure %d → %s\n', figNum, filename);

      % --- Add margin & save high quality ---
      % Set figure background to white
      set(figHandles(i), 'Color', 'white');

      % Adjust paper size slightly larger than figure
      set(figHandles(i), 'PaperPositionMode', 'auto');
      pos = get(figHandles(i), 'Position');
      margin = 50; % pixels of margin
      set(figHandles(i), 'Position', [pos(1)-margin pos(2)-margin pos(3)+2*margin pos(4)+2*margin]);

      % Save as high-resolution PNG
      print(figHandles(i), filename, '-dpng', '-r300');
  end

  fprintf('All figures saved in folder "%s".\n', figures_dir);
  % 
  % pause(2)  % pause 2 seconds to view the plot
  % close all

end
