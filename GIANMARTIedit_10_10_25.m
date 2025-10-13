% ==========================================================================
% -- Single-Antenna Jammers in MIMO-OFDM Can Resemble Multi-Antenna Jammers
% --------------------------------------------------------------------------
% -- (c) 2023 Gian Marti
% -- e-mail: gimarti@ethz.ch
% --------------------------------------------------------------------------
% -- If you use this simulator or parts of it, then you must cite our paper:
%
% -- Gian Marti and Christoph Studer,
% -- "Single-Antenna Jammers in MIMO-OFDM Can Resemble Multi-Antenna Jammers," 
% -- IEEE Communication Letters, 2023
% ==========================================================================

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% configurable parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%
params.num_receive_antennas = 1;                    % number of receive antennas
params.num_transmit_antennas = 1;                    % number of transmit antennas and data streams
params.num_channel_taps = 1;                    % number of channel taps
params.cyclic_prefix_length = 16;                   % cyclic prefix length
params.num_subcarriers = 64;                % number of total subcarriers
params.used_subcarrier_indices = [-26:-22, -20:-8, -6:-1, 1:6, 8:20, 22:26] + 32; % list of used subcarriers, only meaningful for params.num_subcarriers = 64;
params.num_used_subcarriers = length(params.used_subcarrier_indices);
params.num_ofdm_symbols_per_trial = 50;    % number of OFDM symbols sent per channel realization
params.jammer_strength_dB = 30;              % jammer strength relative to signal strength (in dB)
params.modulation_scheme = 'BPSK';             % transmit constellation
params.random_jammer_power_flag = 0;  % if true, then the jammer-power is uniformly (in dB) sampled from
                              % the range [0, params.jammer_strength_dB] for every trial.
                              % Otherwise, it is deterministically equal params.jammer_strength_dB.

params.num_monte_carlo_trials = 100;             % number of Monte-Carlo trials (different channel realizations)
params.snr_dB_list = -5:1:20;     % list of SNR values (in dB) to be simulated

params.plot_results_flag = 1;                 % whether or not to plot the results
params.print_messages_flag = 1;             % whether or not to print log messages

rng(0);                       % random seed


%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% simulation start
%%%%%%%%%%%%%%%%%%%%%%%%%%%

params.simulation_name = [num2str(params.num_receive_antennas) 'x' num2str(params.num_transmit_antennas) '_' params.modulation_scheme '_L', num2str(params.num_channel_taps) '_P' num2str(params.cyclic_prefix_length) ...
  '_Nsc' num2str(params.num_subcarriers) '_Nofd' num2str(params.num_ofdm_symbols_per_trial) '_rho' num2str(params.jammer_strength_dB) ...
  '_' num2str(params.num_monte_carlo_trials) 'Trials'];

switch (params.modulation_scheme)
  case 'BPSK'
    params.constellation_symbols = [ -1 1 ];
  case 'QPSK'
    params.constellation_symbols = [ ...
      -1-1i,-1+1i, ...
      +1-1i,+1+1i ];
  case '16QAM'
    params.constellation_symbols = [ ...
      -3-3i,-3-1i,-3+3i,-3+1i, ...
      -1-3i,-1-1i,-1+3i,-1+1i, ...
      +3-3i,+3-1i,+3+3i,+3+1i, ...
      +1-3i,+1-1i,+1+3i,+1+1i ];
  case '64QAM'
    params.constellation_symbols = [ ...
      -7-7i,-7-5i,-7-1i,-7-3i,-7+7i,-7+5i,-7+1i,-7+3i, ...
      -5-7i,-5-5i,-5-1i,-5-3i,-5+7i,-5+5i,-5+1i,-5+3i, ...
      -1-7i,-1-5i,-1-1i,-1-3i,-1+7i,-1+5i,-1+1i,-1+3i, ...
      -3-7i,-3-5i,-3-1i,-3-3i,-3+7i,-3+5i,-3+1i,-3+3i, ...
      +7-7i,+7-5i,+7-1i,+7-3i,+7+7i,+7+5i,+7+1i,+7+3i, ...
      +5-7i,+5-5i,+5-1i,+5-3i,+5+7i,+5+5i,+5+1i,+5+3i, ...
      +1-7i,+1-5i,+1-1i,+1-3i,+1+7i,+1+5i,+1+1i,+1+3i, ...
      +3-7i,+3-5i,+3-1i,+3-3i,+3+7i,+3+5i,+3+1i,+3+3i ];
end

% normalize transmission constellation
params.constellation_energy = mean(abs(params.constellation_symbols).^2);
params.constellation_symbols = params.constellation_symbols/sqrt(params.constellation_energy);
params.symbol_energy = mean(abs(params.constellation_symbols).^2); % equal to 1

% precompute bit labels
params.bits_per_symbol = log2(length(params.constellation_symbols)); % number of bits per symbol
params.bit_labels = de2bi(0:length(params.constellation_symbols)-1,params.bits_per_symbol,'left-msb');

% track simulation time
time_elapsed = 0;

% intermediate results
rx_signal_freq_jammerless = zeros(params.num_receive_antennas,params.num_subcarriers,params.num_ofdm_symbols_per_trial);
rx_jammer_only_freq_compliant = zeros(params.num_receive_antennas,params.num_subcarriers,params.num_ofdm_symbols_per_trial);
rx_jammer_only_freq_noncompliant = zeros(params.num_receive_antennas,params.num_subcarriers,params.num_ofdm_symbols_per_trial);
rx_signal_with_jammer_freq_compliant = zeros(params.num_receive_antennas,params.num_subcarriers,params.num_ofdm_symbols_per_trial);
rx_signal_with_jammer_freq_noncompliant = zeros(params.num_receive_antennas,params.num_subcarriers,params.num_ofdm_symbols_per_trial);

% results
jammer_signals_freq_compliant_trials = zeros(params.num_monte_carlo_trials, params.num_used_subcarriers,params.num_receive_antennas, params.num_ofdm_symbols_per_trial);
jammer_signals_freq_noncompliant_trials = zeros(params.num_monte_carlo_trials, params.num_used_subcarriers,params.num_receive_antennas, params.num_ofdm_symbols_per_trial);

MSE_jl = zeros(1, length(params.snr_dB_list));
MSE_comp = zeros(1, length(params.snr_dB_list));
MSE_noncomp = zeros(1, length(params.snr_dB_list));

res.BER_jl = zeros(1, length(params.snr_dB_list));
res.BER_comp = zeros(1, length(params.snr_dB_list));
res.BER_noncomp = zeros(1, length(params.snr_dB_list));
res.BER_projs = zeros(params.num_channel_taps, length(params.snr_dB_list));

res.spatial_interference_distr = zeros(params.num_receive_antennas,params.num_used_subcarriers);

% generate bit stream
transmitted_bits = randi([0 1],params.num_monte_carlo_trials,params.num_ofdm_symbols_per_trial,params.num_transmit_antennas,params.bits_per_symbol,params.num_used_subcarriers);
symbol_indices = zeros(params.num_monte_carlo_trials, params.num_ofdm_symbols_per_trial, params.num_transmit_antennas,params.num_used_subcarriers);


% trials loop
tic
for trial_index=1:params.num_monte_carlo_trials

  % generate jammer's channel in the time domain
  jammer_channel_time = sqrt(0.5)*randn(params.num_receive_antennas,params.num_channel_taps) + 1j*sqrt(0.5)*randn(params.num_receive_antennas,params.num_channel_taps);
  jammer_channel_freq = (1/sqrt(params.num_subcarriers))*fft(jammer_channel_time,params.num_subcarriers,2);


  % generate the UEs' channel in the time domain
  user_channel_time = sqrt(0.5)*randn(params.num_receive_antennas,params.num_transmit_antennas,params.num_channel_taps) + 1j*sqrt(0.5)*randn(params.num_receive_antennas,params.num_transmit_antennas,params.num_channel_taps);
  user_channel_freq = (1/sqrt(params.num_subcarriers))*fft(user_channel_time,params.num_subcarriers,3); % E[||H_f(:,:,i)||_F^2] = B*U*L/N_sc for all i=1,...,N_sc
  % so the expected receive signal energy per subcarrier is E_s*B*U*L/N_sc = B*U*L/N_sc

  if params.random_jammer_power_flag
    jammer_power = params.num_transmit_antennas*10^(params.jammer_strength_dB*rand()/10);
  else
    jammer_power = params.num_transmit_antennas*10^(params.jammer_strength_dB/10);
  end

  for snr_index = 1:length(params.snr_dB_list)

    noise_power = params.num_channel_taps*params.num_transmit_antennas*10^(-params.snr_dB_list(snr_index)/10);

    for symbol_index=1:params.num_ofdm_symbols_per_trial

      %%% generate transmit signals

      % generate jammer's transmit signal in the frequency domain for OFDM-compliant case
      jammer_signal_freq_compliant = sqrt(0.5)*randn(1,params.num_subcarriers)+ 1j*sqrt(0.5)*randn(1,params.num_subcarriers);
      jammer_signal_time_no_cp = sqrt(params.num_subcarriers)*ifft(jammer_signal_freq_compliant); %this does in fact not change the Gaussian statistics
      jammer_signal_time_compliant = [jammer_signal_time_no_cp(:,params.num_subcarriers-params.cyclic_prefix_length+1:params.num_subcarriers), jammer_signal_time_no_cp];
      jammer_signal_time_compliant = sqrt(jammer_power)*jammer_signal_time_compliant;

      % generate the jammer's transmit signal in the time domain for OFDM-noncompliant case
      jammer_signal_time_noncompliant = sqrt(0.5)*randn(1,params.num_subcarriers+params.cyclic_prefix_length)+ 1j*sqrt(0.5)*randn(1,params.num_subcarriers+params.cyclic_prefix_length);
      jammer_signal_time_noncompliant = sqrt(jammer_power)*jammer_signal_time_noncompliant;

      % generate UEs' transmit signal in the frequency domain for OFDM-compliant case
      for subcarrier_index=1:params.num_used_subcarriers
        symbol_indices(trial_index,symbol_index,:,subcarrier_index) = bi2de(transmitted_bits(trial_index,symbol_index,:,:,subcarrier_index),'left-msb')+1;
      end
      user_symbols_freq_used = params.constellation_symbols(squeeze(symbol_indices(trial_index,symbol_index,:,:))); % symbol vectors
      user_signal_freq = zeros(params.num_transmit_antennas, params.num_subcarriers);
      user_signal_freq(:,params.used_subcarrier_indices) = user_symbols_freq_used; % only transmit signal over the used subcarriers
      user_signal_time_no_cp = sqrt(params.num_subcarriers)*ifft(user_signal_freq,params.num_subcarriers,2);
      user_signal_time = [user_signal_time_no_cp(:,params.num_subcarriers-params.cyclic_prefix_length+1:params.num_subcarriers), user_signal_time_no_cp];

    
      %%% generate receive signals
      rx_signal_time_jammerless = zeros(params.num_receive_antennas, params.num_subcarriers+params.cyclic_prefix_length+params.num_channel_taps-1);          % jammerless receive signal
      rx_jammer_only_time_compliant = zeros(params.num_receive_antennas, params.num_subcarriers+params.cyclic_prefix_length+params.num_channel_taps-1);     % only compliant jammer
      rx_jammer_only_time_noncompliant = zeros(params.num_receive_antennas, params.num_subcarriers+params.cyclic_prefix_length+params.num_channel_taps-1);  % only noncompliant jammer
      rx_signal_with_jammer_time_compliant = zeros(params.num_receive_antennas, params.num_subcarriers+params.cyclic_prefix_length+params.num_channel_taps-1);        % full receive signal, jammer compliant
      rx_signal_with_jammer_time_noncompliant = zeros(params.num_receive_antennas, params.num_subcarriers+params.cyclic_prefix_length+params.num_channel_taps-1);     % full receive signal, jammer noncompliant

      for antenna_index=1:params.num_receive_antennas
        for transmit_antenna_index=1:params.num_transmit_antennas
          rx_signal_time_jammerless(antenna_index,:) = rx_signal_time_jammerless(antenna_index,:) + conv(squeeze(user_channel_time(antenna_index,transmit_antenna_index,:)),user_signal_time(transmit_antenna_index,:));
        end
        rx_jammer_only_time_compliant(antenna_index,:) = conv(jammer_channel_time(antenna_index,:),jammer_signal_time_compliant);
        rx_jammer_only_time_noncompliant(antenna_index,:) = conv(jammer_channel_time(antenna_index,:),jammer_signal_time_noncompliant);
      end
      rx_signal_with_jammer_time_compliant = rx_signal_time_jammerless + rx_jammer_only_time_compliant;
      rx_signal_with_jammer_time_noncompliant = rx_signal_time_jammerless + rx_jammer_only_time_noncompliant;




      


      %%% generate noise in the time domain
      noise_signal = sqrt(0.5)*randn(params.num_receive_antennas,params.num_subcarriers)+ 1j*sqrt(0.5)*randn(params.num_receive_antennas,params.num_subcarriers);
      noise_signal = sqrt(noise_power)*noise_signal;

      %%% truncate signals in time-domain to remove cylic prefix and channel reverberation
      rx_signal_truncated_jammerless = rx_signal_time_jammerless(:,params.cyclic_prefix_length+1:params.cyclic_prefix_length+params.num_subcarriers) + noise_signal;
      rx_jammer_only_truncated_compliant = rx_jammer_only_time_compliant(:,params.cyclic_prefix_length+1:params.cyclic_prefix_length+params.num_subcarriers); % no noise
      rx_jammer_only_truncated_noncompliant = rx_jammer_only_time_noncompliant(:,params.cyclic_prefix_length+1:params.cyclic_prefix_length+params.num_subcarriers); % no noise
      rx_signal_with_jammer_truncated_compliant = rx_signal_with_jammer_time_compliant(:,params.cyclic_prefix_length+1:params.cyclic_prefix_length+params.num_subcarriers) + noise_signal;
      rx_signal_with_jammer_truncated_noncompliant = rx_signal_with_jammer_time_noncompliant(:,params.cyclic_prefix_length+1:params.cyclic_prefix_length+params.num_subcarriers) + noise_signal;

      %%% convert to frequency domain
      rx_signal_freq_jammerless(:,:,symbol_index) = (1/sqrt(params.num_subcarriers))*fft(rx_signal_truncated_jammerless, params.num_subcarriers, 2);
      rx_jammer_only_freq_compliant(:,:,symbol_index) = (1/sqrt(params.num_subcarriers))*fft(rx_jammer_only_truncated_compliant, params.num_subcarriers, 2);
      rx_jammer_only_freq_noncompliant(:,:,symbol_index) = (1/sqrt(params.num_subcarriers))*fft(rx_jammer_only_truncated_noncompliant, params.num_subcarriers, 2);
      rx_signal_with_jammer_freq_compliant(:,:,symbol_index) = (1/sqrt(params.num_subcarriers))*fft(rx_signal_with_jammer_truncated_compliant, params.num_subcarriers, 2);
      rx_signal_with_jammer_freq_noncompliant(:,:,symbol_index) = (1/sqrt(params.num_subcarriers))*fft(rx_signal_with_jammer_truncated_noncompliant, params.num_subcarriers, 2);


      %%% remove subcarriers that are not used
      rx_signal_freq_tones_jammerless = rx_signal_freq_jammerless(:,params.used_subcarrier_indices,:);
      rx_jammer_only_freq_tones_compliant = rx_jammer_only_freq_compliant(:,params.used_subcarrier_indices,:);
      rx_jammer_only_freq_tones_noncompliant = rx_jammer_only_freq_noncompliant(:,params.used_subcarrier_indices,:);
      rx_signal_with_jammer_freq_tones_compliant = rx_signal_with_jammer_freq_compliant(:,params.used_subcarrier_indices,:);
      rx_signal_with_jammer_freq_tones_noncompliant = rx_signal_with_jammer_freq_noncompliant(:,params.used_subcarrier_indices,:);

      % save signal for plotting (after removing unused subcarriers)
      if trial_index == 1 && snr_index == length(params.snr_dB_list) && symbol_index == 1
          plot_signal_jammerless = rx_signal_freq_tones_jammerless;
      end

      user_channel_freq_tones = user_channel_freq(:,:,params.used_subcarrier_indices);


      %%% collect jammer statistics for later analysis
      for subcarrier_index=1:params.num_used_subcarriers
        jammer_signals_freq_compliant_trials(trial_index,subcarrier_index,:,symbol_index) = rx_jammer_only_freq_tones_compliant(:,subcarrier_index,symbol_index);
        jammer_signals_freq_noncompliant_trials(trial_index,subcarrier_index,:,symbol_index) = rx_jammer_only_freq_tones_noncompliant(:,subcarrier_index,symbol_index);
      end

  
    end % iterate over OFDM symbols with same channel realization


    for subcarrier_index=1:params.num_used_subcarriers

      %%% extract jammer statistics
      [U_comp,S_comp,~] = svd(squeeze(jammer_signals_freq_compliant_trials(trial_index,subcarrier_index,:,:)));
      [U_noncomp,S_noncomp,~] = svd(squeeze(jammer_signals_freq_noncompliant_trials(trial_index,subcarrier_index,:,:)));
      res.spatial_interference_distr(:,subcarrier_index) = res.spatial_interference_distr(:,subcarrier_index) + diag(S_noncomp); % interference distribution over different spatial dimensions



      % form best rank-(B-1) projectors: those that remove most of the
      % jammer interference per subcarrier
      jammer_est_comp = U_comp(:,1);
      projector_comp = eye(params.num_receive_antennas) - jammer_est_comp*jammer_est_comp';
      jammer_est_noncomp = U_noncomp(:,1);
      projector_noncomp = eye(params.num_receive_antennas) - jammer_est_noncomp*jammer_est_noncomp';

      projectors = zeros(params.num_channel_taps,params.num_receive_antennas,params.num_receive_antennas);
      for tap_index=1:params.num_channel_taps
        projectors(tap_index,:,:) = eye(params.num_receive_antennas) - U_noncomp(:,1:tap_index)*U_noncomp(:,1:tap_index)';
      end

      %%% detect signals using zero-forcing (and, in some cases, jammer mitigation)
      symbol_est_jl = (sqrt(params.num_subcarriers)*(user_channel_freq_tones(:,:,subcarrier_index)))\squeeze(rx_signal_freq_tones_jammerless(:,subcarrier_index,:));
      symbol_est_comp = (sqrt(params.num_subcarriers)*(projector_comp*squeeze(user_channel_freq_tones(:,:,subcarrier_index))))\(projector_comp*squeeze(rx_signal_with_jammer_freq_tones_compliant(:,subcarrier_index,:)));
      symbol_est_noncomp = (sqrt(params.num_subcarriers)*(projector_noncomp*squeeze(user_channel_freq_tones(:,:,subcarrier_index))))\(projector_noncomp*squeeze(rx_signal_with_jammer_freq_tones_noncompliant(:,subcarrier_index,:)));
      symbol_est_projs = zeros(params.num_receive_antennas,params.num_transmit_antennas,params.num_ofdm_symbols_per_trial);
      for tap_index=1:params.num_channel_taps
        symbol_est_projs(tap_index,:,:) = (sqrt(params.num_subcarriers)*(squeeze(projectors(tap_index,:,:))*squeeze(user_channel_freq_tones(:,:,subcarrier_index))))\(squeeze(projectors(tap_index,:,:))*squeeze(rx_signal_with_jammer_freq_tones_noncompliant(:,subcarrier_index,:)));
      end

      % -- compute bit outputs
      symbol_est_jl_vec = reshape(symbol_est_jl, [params.num_ofdm_symbols_per_trial*params.num_transmit_antennas,1]);
      [~,idxhat_vec] = min(abs(symbol_est_jl_vec*ones(1,length(params.constellation_symbols))-ones(params.num_ofdm_symbols_per_trial*params.num_transmit_antennas,1)*params.constellation_symbols).^2,[],2);
      idxhat_jl = reshape(idxhat_vec, [params.num_transmit_antennas, params.num_ofdm_symbols_per_trial]);
      bithat_jl = params.bit_labels(idxhat_jl,:);

      symbol_est_comp_vec = reshape(symbol_est_comp, [params.num_ofdm_symbols_per_trial*params.num_transmit_antennas,1]);
      [~,idxhat_vec] = min(abs(symbol_est_comp_vec*ones(1,length(params.constellation_symbols))-ones(params.num_ofdm_symbols_per_trial*params.num_transmit_antennas,1)*params.constellation_symbols).^2,[],2);
      idxhat_comp = reshape(idxhat_vec, [params.num_transmit_antennas, params.num_ofdm_symbols_per_trial]);
      bithat_comp = params.bit_labels(idxhat_comp,:);

      symbol_est_noncomp_vec = reshape(symbol_est_noncomp, [params.num_ofdm_symbols_per_trial*params.num_transmit_antennas,1]);
      [~,idxhat_vec] = min(abs(symbol_est_noncomp_vec*ones(1,length(params.constellation_symbols))-ones(params.num_ofdm_symbols_per_trial*params.num_transmit_antennas,1)*params.constellation_symbols).^2,[],2);
      idxhat_noncomp = reshape(idxhat_vec, [params.num_transmit_antennas, params.num_ofdm_symbols_per_trial]);
      bithat_noncomp = params.bit_labels(idxhat_noncomp,:);

      bithat_projs = zeros(params.num_channel_taps,params.num_transmit_antennas*params.num_ofdm_symbols_per_trial,params.bits_per_symbol);
      for tap_index=1:params.num_channel_taps
        symbol_est_vec = reshape(symbol_est_projs(tap_index,:,:), [params.num_ofdm_symbols_per_trial*params.num_transmit_antennas,1]);
        [~,idxhat_vec] = min(abs(symbol_est_vec*ones(1,length(params.constellation_symbols))-ones(params.num_ofdm_symbols_per_trial*params.num_transmit_antennas,1)*params.constellation_symbols).^2,[],2);
        idxhat = reshape(idxhat_vec, [params.num_transmit_antennas, params.num_ofdm_symbols_per_trial]);
        bithat_projs(tap_index,:,:) = params.bit_labels(idxhat,:);
      end

      % -- compute error metrics
      bit_tensor = squeeze(transmitted_bits(trial_index,:,:,:,subcarrier_index));
      bit_tensor = permute(bit_tensor,[3 2 1]);
      true_bits = bit_tensor(:,:)';
      res.BER_jl(snr_index) = res.BER_jl(snr_index) + sum(sum(true_bits~=bithat_jl));
      res.BER_comp(snr_index) = res.BER_comp(snr_index) + sum(sum(true_bits~=bithat_comp));
      res.BER_noncomp(snr_index) = res.BER_noncomp(snr_index) + sum(sum(true_bits~=bithat_noncomp));
      for tap_index=1:params.num_channel_taps
        res.BER_projs(tap_index,snr_index) = res.BER_projs(tap_index,snr_index) + sum(sum(true_bits~=squeeze(bithat_projs(tap_index,:,:))));
      end

    end


  end % iterate over SNRs

  % -- keep track of simulation time
  if(params.print_messages_flag)
    if toc>10
      time = toc;
      time_elapsed = time_elapsed + time;
        fprintf('estimated remaining simulation time: %3.0f min.\n', ...
          time_elapsed*(params.num_monte_carlo_trials/trial_index-1)/60);
      tic
    end
  end

end % iterate over MC trials

% normalize results
res.BER_jl = res.BER_jl/(params.num_monte_carlo_trials*params.num_transmit_antennas*params.num_ofdm_symbols_per_trial*params.bits_per_symbol*params.num_used_subcarriers);
res.BER_comp = res.BER_comp/(params.num_monte_carlo_trials*params.num_transmit_antennas*params.num_ofdm_symbols_per_trial*params.bits_per_symbol*params.num_used_subcarriers);
res.BER_noncomp = res.BER_noncomp/(params.num_monte_carlo_trials*params.num_transmit_antennas*params.num_ofdm_symbols_per_trial*params.bits_per_symbol*params.num_used_subcarriers);
res.BER_projs = res.BER_projs/(params.num_monte_carlo_trials*params.num_transmit_antennas*params.num_ofdm_symbols_per_trial*params.bits_per_symbol*params.num_used_subcarriers);

res.spatial_interference_distr = res.spatial_interference_distr./sum(res.spatial_interference_distr,1);
res.spatial_interference_distr_sc_std = std(res.spatial_interference_distr.');
res.spatial_interference_distr_sc_avg = mean(res.spatial_interference_distr.');



%%% plot results
if params.plot_results_flag

  figure(1)
  semilogy(params.snr_dB_list, res.BER_jl, 'LineWidth', 2.5)
  hold on
  semilogy(params.snr_dB_list, res.BER_comp, 'LineWidth', 2.5)
  semilogy(params.snr_dB_list, res.BER_noncomp, 'LineWidth', 2.5)
  hold off
  grid on
  axis([min(params.snr_dB_list) max(params.snr_dB_list) 1e-4 1])  
  xlabel('average SNR per receive antenna [dB]','FontSize',12)
  ylabel('uncoded bit error-rate (BER)','FontSize',12)
  legend(["Jammerless", "OFDM-compliant jammer", "Cyclic prefix-violating jammer"],'FontSize',12,'Interpreter','none')
  set(gca,'FontSize',12)
  set(gcf,'position',[10,10,400,300])

  figure(2)
  bar(1:params.num_receive_antennas,res.spatial_interference_distr_sc_avg)
  hold on
  er = errorbar(1:params.num_receive_antennas,res.spatial_interference_distr_sc_avg, ...
    2*res.spatial_interference_distr_sc_std, ...
    2*res.spatial_interference_distr_sc_std);
  er.Color = [0 0 0];
  er.LineStyle = 'none';
  hold off
  axis([0.25 params.num_receive_antennas+0.75 0 1])
  grid on
  xlabel('ordered spatial dimension index')
  ylabel('fraction of receive interference')
  set(gcf,'position',[10,10,400,300])

  figure(3)
  semilogy(params.snr_dB_list, res.BER_projs(1,:), 'LineWidth', 2.5)
  hold on
  for b=2:params.num_channel_taps
    semilogy(params.snr_dB_list, res.BER_projs(b,:), 'LineWidth', 2.5)
  end
  hold off
  grid on
  axis([min(params.snr_dB_list) max(params.snr_dB_list) 1e-4 1])
  xlabel('average SNR per receive antenna [dB]','FontSize',12)
  ylabel('uncoded bit error-rate (BER)','FontSize',12)
  numbers = {};
  for b = 1:params.num_channel_taps
    numbers{b} = "Null " + num2str(b) + " Dimensions";
  end
  legend(numbers,'FontSize',12,'Interpreter','none')
  set(gca,'FontSize',12)
  set(gcf,'position',[10,10,400,300])

  figure(4)
  plot(1:params.num_used_subcarriers, abs(squeeze(plot_signal_jammerless(1,:,1))), 'LineWidth', 2)
  xlabel('Used Subcarrier Index','FontSize',12)
  ylabel('Signal Magnitude','FontSize',12)
  title('Jammerless Signal After Removing Unused Subcarriers','FontSize',12)
  grid on
  set(gca,'FontSize',12)
  set(gcf,'position',[10,10,400,300])

end
