%% Nerve EEG for Guitar. May 2017

clear all;
close all;
clc;

cd('C:\Users\Thomas\Documents\MATLAB\SMC8tk2\SMC8SEM_PROJECT\0905_NEEG_DATA');

%% Load data

% Load data
load('NEEG090517_dummy.mat');
dummy=y(2:end,4800:end);
load('NEEG090517_faststrum.mat');
fast_strum=y(2:end,4800:end);
load('NEEG090517_LoudSlow2.mat');
slow_loud_strum2=y(2:end,4800:end);
load('NEEG090517_loudslowstrum.mat');
slow_loud_strum=y(2:end,4800:end);
load('NEEG090517_picking.mat');
picking=y(2:end,4800:end);
load('NEEG090517_pickingFastLoud.mat');
fast_loud_picking=y(2:end,4800:end);
load('NEEG090517_pickingQuiteSlow.mat');
slow_picking=y(2:end,4800:end);
load('NEEG090517_slowstrum.mat');
slow_strum=y(2:end,4800:end);

%Clear unnecessary variables from workspace)
clear('y');

%% Reference


%% Remove areas of noise

% Smile Noise
%datamat_T008_OSmile(:,[1:(5*4800) (17*4800):(29*4800) (45*4800):(56*4800) (66*4800):end])=[];


% %% Create noise samples
%
% % Smile Noise
% noise_smile=[datamat_T008_OSmile(:,(17*4800):(29*4800)) datamat_T008_OSmile(:,...
%     (46*4800):(56*4800))];
%
% % Frown Noise
% noise_frown=[datamat_T009_Frown(:,(13*4800):(18*4800)) datamat_T008_OSmile(:,...
%     (29*4800):(37*4800))];

%% Plot All Channels
% datatime=(1:length(fast_strum))/4800;
% figure();
% plot(datatime,(fast_strum([1 3 4 5 6 7 8 9],:))); hold on
% xlabel('Time(s)');
% ylabel('Voltage(mV)');
% title('Raw fast_strum');
% %ylim([-5000 5000]);
% grid minor ;
% clear('datatime')
% 
% %% Apply bandpass and comb filter
% 
% % Clear Unnecessary Variables
% clear('z_mag', 'z_phas');

% Index Wanted Variables
worksp_idx= evalin('base','whos');
[work_L ~]=size(worksp_idx);

% % Design bandpass
d_bp1 = fdesign.bandpass('N,Fp1,Fp2,Ap', 60, 80, 800, .5, 4800);
% %high pass (above resting heart rate), low pass above action potential maximum freq
h_bp1 = design(d_bp1, 'cheby1');

% View filter z_magnitude response
%fvtool(h_bp1);

% Design Comb Filter
d_comb=fdesign.comb('notch','L,BW,GBW,Nsh',96,5,-2,2,4800);
h_comb=design(d_comb);
%fvtool(h_comb);

% Apply Filters
idx_var=0;
for idx_var=0:work_L-1
    idx_var=idx_var+1;
    name_var=(worksp_idx(idx_var).name);
    curr_var=eval(name_var);
    idx_chan=0;
    [curr_var_L ~]=size(curr_var);
    for idx_chan=0:curr_var_L-1;
        idx_chan=idx_chan+1;
        curr_var(idx_chan,:)=filter(h_bp1,(curr_var(idx_chan,:)));
    end
    assignin('base' , (worksp_idx(idx_var).name) , curr_var);
    idx_chan=0;
    [curr_var_L ~]=size(curr_var);
    for idx_chan=0:curr_var_L-1;
        idx_chan=idx_chan+1;
        curr_var(idx_chan,:)=filter(h_comb,curr_var(idx_chan,:));
    end
    assignin('base' , (worksp_idx(idx_var).name) , curr_var);
end
clear( 'curr_var', 'curr_var_L', 'idx_chan',...
    'idx_var', 'name_var', 'work_L', 'worksp_idx',...
    'd_bp1', 'h_bp1','h_comb', 'd_comb');

% %% Cross correlate signal against white noise
%
% % Generate White Noise
%
% wn=wgn(1,length(datamat_T008_OSmile(8,:)),0.1);
%
% % Perform and plot cross correlation
% [acor,lag] = xcorr(wn(1:1000),datamat_T008_OSmile(8,:));
% [~,I] = max(abs(acor));
% lagDiff = lag(I);
% timeDiff = lagDiff/48000;
% figure()
% plot(wn(1:1000)); hold on;
% % idx_np=find((acor/max(acor))<0.9);
% % acor(:,idx_np)=0;
% plot(lag,acor/max(acor))

% %% View z_magnitude and z_phase (per channel)
% 
% % Create cells for z_magnitude and z_phase per channel
% worksp_idx= evalin('base','whos');% Create workspace index matrix
% [work_L ~]=size(worksp_idx);
% idx_var=0;
% for idx_var=0:work_L-1
%     idx_var=idx_var+1;
%     name_var=(worksp_idx(idx_var).name);
%     curr_var=eval(name_var);
%     idx_chan=0;
%     [curr_var_L ~]=size(curr_var);
%     for idx_chan=0:curr_var_L-1;
%         idx_chan=idx_chan+1;
%         y{idx_var}(idx_chan,:)=fft(curr_var(idx_chan,:));
%         z_mag{idx_var}(idx_chan,:)=abs(y{idx_var}(idx_chan,:));
%         z_phas{idx_var}(idx_chan,:)=unwrap(angle(y{idx_var}(idx_chan,:)));
%         % Normalise z_magnitude and z_phase
%         z_mag{idx_var}(idx_chan,:)=z_mag{idx_var}(idx_chan,:)/max(abs(z_mag{idx_var}(idx_chan,:)));
%         z_phas{idx_var}(idx_chan,:)=z_phas{idx_var}(idx_chan,:)/max(abs(z_phas{idx_var}(idx_chan,:)));
%     end
% end
% clear( 'curr_var','curr_var_L', 'idx_chan','y',...
%     'idx_var', 'name_var','work_L', 'worksp_idx');

% % Plot z_magnitude and z_phase
% fs=4800;
% work_L=length(z_mag);
% idx_var=0;
% for idx_var=0:work_L-1
%     idx_var=idx_var+1;
%     n=length(z_mag{idx_var}(1,:))-1;
%     f=0:fs/n:fs;
%     figure(idx_var)
%     idx_chan=0;
%     [curr_var_L ~]=size(z_mag{idx_var});
%     for idx_chan=0:curr_var_L-1;
%         idx_chan=idx_chan+1;
%         subplot(2,1,1);
%         plot(f,z_mag{idx_var}(idx_chan,:)); hold on;
%         xlim([0 2400]);
%         if idx_var==1
%             title('Dummy')
%         elseif idx_var==2
%             title('Fast_loud_picking')
%         elseif idx_var==3
%             title('Fast_strum')
%         elseif idx_var==4
%             title('Picking')
%         elseif idx_var==5
%             title('Slow_loud_strum')
%         elseif idx_var==6
%             title('Slow_loud_strum2')
%         elseif idx_var==7
%             title('Slow_Picking')
%         elseif idx_var==8
%             title('Slow_strum')
%         else
%         end
%         subplot(2,1,2);
%         plot(f,z_phas{idx_var}(idx_chan,:)); hold on;
%         xlim([0 2400]);
%     end
% end
% clear( 'curr_var','curr_var_L', 'idx_chan','y',...
%     'idx_var', 'name_var', 'work_L', 'worksp_idx',...
%     'n', 'fs', 'f');
% 
% %% Plot All Channels (Recommended to do after each stage to monitor signal)
% mains=[(-11000)*(ones(size(datamat_T008_OSmile(1,:))));(11000)*(ones(size(datamat_T008_OSmile(1,:))))];
% datatime=1:length(datamat_T008_OSmile(1,:));
% datatime=datatime/4800;
% figure();
% plot(datatime,datamat_T008_OSmile(:,:)); hold on
% plot(datatime,mains(1:2,:),'k');
% xlabel('Time(s)');
% ylabel('Voltage(mV)');
% title('Bandpass Filtered Multi-Channel T008 Open Smile');
% ylim([-1000 1000]);
% grid minor ;
% clear('dattime','mains')
% %% Plot All Channels (Recommended to do after each stage to monitor signal)
% mains=[(-11000)*(ones(size(datamat_T009_Frown(1,:))));(11000)*(ones(size(datamat_T009_Frown(1,:))))];
% datatime=1:length(datamat_T009_Frown(1,:));
% datatime=datatime/4800;
% figure();
% plot(datatime,datamat_T009_Frown(:,:)); hold on
% plot(datatime,mains(1:2,:),'k');
% xlabel('Time(s)');
% ylabel('Voltage(mV)');
% title('Bandpass Filtered Multi-Channel T009 Frown');
% ylim([-1000 1000]);
% grid minor ;
% clear('dattime','mains')

%% Smooth/resample data
worksp_idx= evalin('base','whos');
[work_L ~]=size(worksp_idx);
idx_var=0;
for idx_var=0:work_L-1
    idx_var=idx_var+1;
    name_var=(worksp_idx(idx_var).name);
    curr_var=eval(name_var);
    idx_chan=0;
    [curr_var_L ~]=size(curr_var);
    for idx_chan=0:curr_var_L-1;
        idx_chan=idx_chan+1;
        curr_var2(idx_chan,:)=interp(curr_var(idx_chan,:),2);
        clear('idx_low','idx_high','idx_sig');
        % idx_high=find((curr_var2(idx_chan,:)) >110);
        %curr_var2(idx_chan,idx_high)=110*ones(size(curr_var2(idx_chan,idx_high)));%Cap Positive Signal
        %idx_low=find((curr_var2(idx_chan,:)) <-110);
        %curr_var2(idx_chan,idx_low)=-110*ones(size(curr_var2(idx_chan,idx_low)));%Cap Positive Signal
    end
    assignin('base' , (worksp_idx(idx_var).name) , curr_var2);
    clear('curr_var2')
end
clear( 'curr_var', 'curr_var_L', 'idx_chan',...
    'idx_var', 'name_var', 'work_L', 'worksp_idx',...
    'chan_L', 'idx_sig', 'idx_spike', 'var_L',...
    'idx_high','idx_low','curr_data');

%% Rectify Signal
%
%  datamat_T009_Frown=abs(datamat_T009_Frown);
%  datamat_T008_OSmile=abs(datamat_T008_OSmile);

%% Measure pulse widths (pw) and their onsets and offsets
worksp_idx= evalin('base','whos');
[work_L ~]=size(worksp_idx);
idx_var=0;
for idx_var=0:work_L-1
    idx_var=idx_var+1;
    name_var=(worksp_idx(idx_var).name);
    curr_var=eval(name_var);
    idx_chan=0;
    [curr_var_L ~]=size(curr_var);
    for idx_chan=0:curr_var_L-1;
        idx_chan=idx_chan+1;
        [zpw{idx_var,idx_chan},zINITCROSS{idx_var,idx_chan},...
            zFINALCROSS{idx_var,idx_chan}]=pulsewidth(curr_var(idx_chan,:),...
            zfs,'StateLevels',[-1 1]);
        %figure()
        %pulsewidth(curr_var(idx_chan,:),zfs,'StateLevels',[-1 1]);
    end
    assignin('base' , (worksp_idx(idx_var).name) , curr_var);
end
clear( 'curr_var', 'curr_var_L', 'idx_chan',...
    'idx_var', 'name_var', 'work_L', 'worksp_idx');
%% Create cell of pulsewidths with indices for onsets and offsets
idx_var=0;
[var_L chan_L]=size(zpw);
for idx_var=0:var_L-1
    idx_var=idx_var+1;
    idx_sig=0;
    for idx_sig=0:chan_L-1;
        idx_sig=idx_sig+1;
        spikecell{idx_var,idx_sig}=[round(zpw{idx_var,idx_sig},...
            3,'significant');round(zINITCROSS{idx_var,idx_sig},6,...
            'significant');round(zFINALCROSS{idx_var,idx_sig},6,'significant')];
    end
end
clear('idx_sig', 'idx_var', 'var_B','chan_L','idx_high',...
    'idx_low','zFINALCROSS', 'zINITCROSS', 'zpw','var_L');

%% Add peak amplitude to cell of spikes
worksp_idx= evalin('base','whos');
[work_L ~]=size(worksp_idx);
idx_var=0;
for idx_var=0:work_L-2
    idx_var=idx_var+1;
    name_var=(worksp_idx(idx_var).name);
    curr_var=eval(name_var);
    idx_chan=0;
    [var_B var_L]=size(spikecell);
    for idx_chan=0:var_L-1
        idx_chan=idx_chan+1;
        curr_data=curr_var(idx_chan,:);
        chan_L=length(spikecell{idx_var,idx_chan});
        idx_sig=0;
        for idx_sig=0:chan_L-1;
            idx_sig=idx_sig+1;
            spikecell{idx_var,idx_chan}(4,idx_sig)=max(abs(curr_data(1,floor(spikecell{idx_var,idx_chan}(2,...
                idx_sig)*zfs):ceil(spikecell{idx_var,idx_chan}(3,idx_sig)*zfs))));
        end
        spikecell{idx_var,idx_chan}(4,:)=spikecell{idx_var,idx_chan}(4,:)/max(spikecell{idx_var,idx_chan}(4,:)); %Normalise amplitude
        %spikecell{idx_var,idx_chan}(4,:)=round(spikecell{idx_var,idx_chan}(4,:)/0.01)*0.01;
    end
end
clear( 'curr_var', 'curr_var_L', 'curr_var_B', 'idx_chan',...
    'idx_var', 'name_var', 'work_L', 'worksp_idx','var_B',...
    'chan_L', 'idx_sig', 'idx_spike', 'var_L','curr_data');
%%  Check distrubution of peaks per amplitude
figure()
histogram(spikecell{1,3}(4,:),'BinWidth',0.001);
hold on;
line([(mean(spikecell{1,3}(4,:))) (mean(spikecell{1,3}(4,:)))],...
    [0 max(histcounts(spikecell{1,3}(4,:),'BinWidth',0.001))],'Color','r');
grid minor;
title('Open Smile Distribution of Peaks per Amplitude')
xlabel('Normalised Amplitude')
ylabel('Occurences')
text((mean(spikecell{1,3}(4,:))),(max(histcounts(spikecell{1,...
    5}(4,:),'BinWidth',0.001))/2),'red = mean amplitude')

figure()
histogram(spikecell{2,5}(4,:),'BinWidth',0.001);
hold on;
line([(mean(spikecell{2,5}(4,:))) (mean(spikecell{2,5}(4,:)))],...
    [0 max(histcounts(spikecell{2,5}(4,:),'BinWidth',0.001))],'Color','r');
grid minor;
title('Frown Distribution of Peaks per Amplitude')
xlabel('Normalised Amplitude')
ylabel('Occurences')
text((mean(spikecell{2,5}(4,:))),(max(histcounts(spikecell{2,...
    1}(4,:),'BinWidth',0.001))/2),'red = mean amplitude')
%% Remove pulses/spikes with negative peaks or zero values
idx_var=0;
[var_L chan_L]=size(spikecell);
for idx_var=0:var_L-1
    idx_var=idx_var+1;
    idx_sig=0;
    for idx_sig=0:chan_L-1;
        idx_sig=idx_sig+1;
        idx_noise=find(spikecell{idx_var,idx_sig}(4,:) <0);
        spikecell{idx_var,idx_sig}(:,idx_noise)=[];
        clear('idx_noise','idx_low')
        idx_noise=find(spikecell{idx_var,idx_sig}(4,:) ==0);
        spikecell{idx_var,idx_sig}(:,idx_noise)=[];
        clear('idx_noise','idx_low')
    end
end
clear('idx_sig', 'idx_var', 'var_L','chan_L','idx_high','idx_low');

%% Averaged histcounts (number of occurances at particular pulsewidths)
clear('zN','zedges')
idx_var=0;
[var_L chan_L]=size(spikecell);
for idx_var=0:var_L-1
    idx_var=idx_var+1;
    idx_sig=0;
    for idx_sig=0:chan_L-1;
        idx_sig=idx_sig+1;
        [zPulse_Occ{idx_var,idx_sig},zPulse_Wid{idx_var,idx_sig}] = histcounts(spikecell{idx_var,idx_sig}(1,:),...
            'BinWidth',0.0002,'BinLimits',[0.00005,0.0050002]);
        zPulse_Occ{idx_var,idx_sig}=zPulse_Occ{idx_var,idx_sig}/max(zPulse_Occ{idx_var,idx_sig});
        zPulse_Wid{idx_var,idx_sig}(:,1)=[];
    end
end
clear('idx_sig', 'idx_var', 'var_L','chan_L');

%% Weighted (by spike amplitude) mean number of occurances per pulsewidth

[cell_L cell_B]=size(spikecell);
idx_var=0;
for idx_var=0:cell_L-1
    idx_var=idx_var+1;
    idx_chan=0;
    [var_B var_L]=size(spikecell);
    for idx_chan=0:var_L-1
        idx_chan=idx_chan+1;
        chan_L=length(zPulse_Wid{idx_var,idx_chan});
        idx_sig=0;
        for idx_sig=0:chan_L-1;
            idx_sig=idx_sig+1;
            if zPulse_Occ{idx_var,idx_chan}(idx_sig)>0;
                idx_pw=find(spikecell{idx_var,idx_chan}(1,:)==zPulse_Wid{idx_var,idx_chan}(idx_sig));
                zPulse_MPeak{idx_var,idx_chan}(:,idx_sig)=(sum(abs(spikecell{idx_var,...
                    idx_chan}(4,idx_pw))))/zPulse_Occ{idx_var,idx_chan}(idx_sig);
            else
                zPulse_MPeak{idx_var,idx_chan}(:,idx_sig)=0;
            end
        end
        zPulse_MPeak{idx_var,idx_chan}=zPulse_MPeak{idx_var,idx_chan}/max(zPulse_MPeak{idx_var,idx_chan});
        zPulse_AmpOcc{idx_var,idx_chan}=(diag(zPulse_MPeak{idx_var,idx_chan}'*...
            zPulse_Occ{idx_var,idx_chan}))/max(diag(zPulse_MPeak{idx_var,...
            idx_chan}'*zPulse_Occ{idx_var,idx_chan}));
    end
end
clear( 'curr_var', 'curr_var_L', 'curr_var_B', 'idx_chan',...
    'idx_var', 'name_var', 'work_L', 'worksp_idx','var_B',...
    'chan_L', 'idx_sig', 'idx_spike', 'var_L','curr_data',...
    'cell_L','cell_B','idx_pw');

%% Plot normalised number of occurences per pulsewidth
figure()
plot((zPulse_Wid{1,3}),(zPulse_Occ{1,3}),'b');
hold on;
plot((zPulse_Wid{2,5}),(zPulse_Occ{2,5}),'r');
grid minor;
xlabel('Pulsewidth of Action Potential (seconds)');
ylabel('Normalised Occurances');
legend('Open Smile','Frown')
title('Normalised Number of Occurances Per Pulsewidth');
%% Plot normalised mean amplitude per pulsewidth
figure()
plot((zPulse_Wid{1,3}),(zPulse_MPeak{1,3}),'b');
hold on;
plot((zPulse_Wid{2,5}),(zPulse_MPeak{2,5}),'r');
grid minor;
xlabel('Pulsewidth of Action Potential (seconds)');
ylabel('Normalised Occurances');
legend('Open Smile','Frown')
title('Normalised Mean Amplitude Per Pulsewidth');
%% Plot weighted curve for occurences*mean amplitude per pulsewidth
figure()
plot((zPulse_Wid{2,5}),(zPulse_AmpOcc{2,5}),'r');
hold on;
plot((zPulse_Wid{2,7}),(zPulse_AmpOcc{2,7}),'r--');
hold on;
plot((zPulse_Wid{1,3}),(zPulse_AmpOcc{1,3}),'b');
hold on;
plot((zPulse_Wid{1,7}),(zPulse_AmpOcc{1,7}),'b--');
grid minor;
xlabel('Pulsewidth of Action Potential (seconds)');
ylabel('Weighted Occurances*Mean Amplitude');
% line([0.001 0.001 ],[0 1],'Color','g');
% line([0.002 0.002 ],[0 1],'Color','k');
legend('Frown:Muscle','Frown:Nerve','Smile:Muscle','Smile:Nerve','location','northwest')%'Minimum Refractory Period','Maximum AP Pulsewidth'
title('Normalised Mean Power per Action Potential Pulse-width');
set(legend,'FontSize',13);
set(findall(gca, 'Type', 'Line'),'LineWidth',2);
%% Use the curve above to determine noise threshold and remove noise
idx_var=0;
[var_L chan_L]=size(spikecell);
for idx_var=0:var_L-1
    idx_var=idx_var+1;
    idx_sig=0;
    for idx_sig=0:chan_L-1;
        idx_sig=idx_sig+1;
        idx_noise=find(spikecell{idx_var,idx_sig}(1,:) <0.95*10.^-3);
        spikecell{idx_var,idx_sig}(:,idx_noise)=[];
        clear('idx_noise')
        idx_noise=find(spikecell{idx_var,idx_sig}(1,:) >5.05*10.^-3);
        spikecell{idx_var,idx_sig}(:,idx_noise)=[];
        clear('idx_noise','idx_low')
        idx_noise=find(spikecell{idx_var,idx_sig}(4,:) <0.1);
        spikecell{idx_var,idx_sig}(:,idx_noise)=[];
    end
end
clear('idx_sig', 'idx_var', 'var_L','chan_L','idx_high','idx_low');
%%  Check distrubution of peaks per amplitude
figure()
histfit(spikecell{2,5}(4,:),(1/0.001));
hold on;
%histogram(spikecell{2,5}(4,:),'BinWidth',0.001);
%hold on;
line([(mean(spikecell{2,5}(4,:))) (mean(spikecell{2,5}(4,:)))],...
    [0 max(histcounts(spikecell{2,5}(4,:),'BinWidth',0.001))],'Color','r');
grid minor;
title('Open Smile Distribution of Peaks per Amplitude')
xlabel('Normalised Amplitude')
ylabel('Occurences')

figure()
histfit(spikecell{1,1}(4,:),(1/0.001));
hold on;
%histogram(spikecell{1,1}(4,:),'BinWidth',0.001);
%hold on;
line([(mean(spikecell{1,1}(4,:))) (mean(spikecell{1,1}(4,:)))],...
    [0 max(histcounts(spikecell{1,1}(4,:),'BinWidth',0.001))],'Color','r');
grid minor;
title('Frown Distribution of Peaks per Amplitude')
xlabel('Normalised Amplitude')
ylabel('Occurences')

%%  Check distrubution of peaks per pulsewidth
figure()
histfit(spikecell{2,1}(1,:),(1/0.001));
hold on;
%histogram(spikecell{2,5}(4,:),'BinWidth',0.001);
%hold on;
line([(mean(spikecell{2,1}(1,:))) (mean(spikecell{2,1}(1,:)))],...
    [0 max(histcounts(spikecell{2,1}(1,:),'BinWidth',0.001))],'Color','r');
grid minor;
title('Open Smile Distribution of Peaks per Pulsewidth')
xlabel('Normalised Amplitude')
ylabel('Occurences')

figure()
histfit(spikecell{1,1}(1,:),(1/0.001));
hold on;
%histogram(spikecell{1,1}(1,:),'BinWidth',0.001);
%hold on;
line([(mean(spikecell{1,1}(1,:))) (mean(spikecell{1,1}(1,:)))],...
    [0 max(histcounts(spikecell{1,1}(1,:),'BinWidth',0.001))],'Color','r');
grid minor;
title('Frown Distribution of Peaks per Pulsewidth')
xlabel('Normalised Amplitude')
ylabel('Occurences')

%% Averaged histcounts (number of occurances at particular pulsewidths)
clear('zN','zedges')
idx_var=0;
[var_L chan_L]=size(spikecell);
for idx_var=0:var_L-1
    idx_var=idx_var+1;
    idx_sig=0;
    for idx_sig=0:chan_L-1;
        idx_sig=idx_sig+1;
        [zPulse_Occ{idx_var,idx_sig},zPulse_Wid{idx_var,idx_sig}] = histcounts(spikecell{idx_var,idx_sig}(1,:),...
            'BinWidth',0.0002,'BinLimits',[0.00005,0.0050002]);
        zPulse_Occ{idx_var,idx_sig}=zPulse_Occ{idx_var,idx_sig}/max(zPulse_Occ{idx_var,idx_sig});
        zPulse_Wid{idx_var,idx_sig}(:,1)=[];
    end
end
clear('idx_sig', 'idx_var', 'var_L','chan_L');

%% Weighted (by spike amplitude) mean number of occurances per pulsewidth

[cell_L cell_B]=size(spikecell);
idx_var=0;
for idx_var=0:cell_L-1
    idx_var=idx_var+1;
    idx_chan=0;
    [var_B var_L]=size(spikecell);
    for idx_chan=0:var_L-1
        idx_chan=idx_chan+1;
        chan_L=length(zPulse_Wid{idx_var,idx_chan});
        idx_sig=0;
        for idx_sig=0:chan_L-1;
            idx_sig=idx_sig+1;
            if zPulse_Occ{idx_var,idx_chan}(idx_sig)>0;
                idx_pw=find(spikecell{idx_var,idx_chan}(1,:)==zPulse_Wid{idx_var,idx_chan}(idx_sig));
                zPulse_MPeak{idx_var,idx_chan}(:,idx_sig)=(sum(abs(spikecell{idx_var,...
                    idx_chan}(4,idx_pw))))/zPulse_Occ{idx_var,idx_chan}(idx_sig);
            else
                zPulse_MPeak{idx_var,idx_chan}(:,idx_sig)=0;
            end
        end
        zPulse_MPeak{idx_var,idx_chan}=zPulse_MPeak{idx_var,idx_chan}/max(zPulse_MPeak{idx_var,idx_chan});
        zPulse_AmpOcc{idx_var,idx_chan}=(diag(zPulse_MPeak{idx_var,idx_chan}'*...
            zPulse_Occ{idx_var,idx_chan}))/max(diag(zPulse_MPeak{idx_var,...
            idx_chan}'*zPulse_Occ{idx_var,idx_chan}));
    end
end
clear( 'curr_var', 'curr_var_L', 'curr_var_B', 'idx_chan',...
    'idx_var', 'name_var', 'work_L', 'worksp_idx','var_B',...
    'chan_L', 'idx_sig', 'idx_spike', 'var_L','curr_data',...
    'cell_L','cell_B','idx_pw');
%% Plot normalised number of occurences per pulsewidth
figure()
plot((zPulse_Wid{1,3}),(zPulse_Occ{1,3}),'b');
hold on;
plot((zPulse_Wid{2,5}),(zPulse_Occ{2,5}),'r');
grid minor;
xlabel('Pulsewidth of Action Potential (seconds)');
ylabel('Normalised Occurances');
legend('Open Smile','Frown')
title('Normalised Number of Occurances Per Pulsewidth');
%% Plot normalised mean amplitude per pulsewidth
figure()
plot((zPulse_Wid{1,3}),(zPulse_MPeak{1,3}),'b');
hold on;
plot((zPulse_Wid{2,5}),(zPulse_MPeak{2,5}),'r');
grid minor;
xlabel('Pulsewidth of Action Potential (seconds)');
ylabel('Normalised Occurances');
legend('Open Smile','Frown')
title('Normalised Mean Amplitude Per Pulsewidth');
%% Plot weighted curve for occurences*mean amplitude per pulsewidth
figure()
plot((zPulse_Wid{1,3}),(zPulse_AmpOcc{1,3}),'b');
hold on;
plot((zPulse_Wid{2,5}),(zPulse_AmpOcc{2,5}),'r');
grid minor;
xlabel('Pulsewidth of Action Potential (seconds)');
ylabel('Weighted Occurances*Mean Amplitude');
line([0.001 0.001 ],[0 1],'Color','g');
line([0.002 0.002 ],[0 1],'Color','k');
legend('Frown','Open Smile','Minimum Refractory Period','Maximum AP Pulsewidth')
title('Normalised Weighted Occurances*Mean Amplitude Per Pulsewidth');
%% Label Action Potentials (less than 2.01ms) and overlapped Motor Action Potentials (greater than 2.01ms)
idx_var=0;
[var_L chan_L]=size(spikecell);
for idx_var=0:var_L-1
    idx_var=idx_var+1;
    idx_sig=0;
    for idx_sig=0:chan_L-1;
        idx_sig=idx_sig+1;
        idx_ap=find(spikecell{idx_var,idx_sig}(1,:) <2.01*10.^-3);
        spikecell{idx_var,idx_sig}(5,idx_ap)=1;
        idx_cmap=find(spikecell{idx_var,idx_sig}(1,:) >2.01*10.^-3);
        spikecell{idx_var,idx_sig}(5,idx_cmap)=2;
        clear('idx_ap','idx_cmap')
    end
end
clear('idx_sig', 'idx_var', 'var_L','chan_L','idx_high','idx_low');

%% Create peri-time-stimulus histograms

%Smile
timewindow=(0.005*4800);
idx_tw=(5*4800);
idx_tw_end=(67*4800); %length(datamat_T008_OSmile)/zfs
for idx_tw=idx_tw:idx_tw_end
    idx_tw=idx_tw+timewindow;
    if (idx_tw+timewindow)<(idx_tw_end)
        clear('idx_p');
        idx_p=find((spikecell{1,3}(2,:))*4800>idx_tw&(spikecell{1,3}(2,...
            :))*4800<idx_tw+timewindow);
        if isempty(idx_p)
            spike_train{1,1}(1,idx_tw:idx_tw+timewindow)=zeros([1 timewindow+1]);
        else
            clear('spike_p')
            spike_p=mean(spikecell{1,3}(5,idx_p));
            spike_train{1,1}(1,idx_tw:idx_tw+timewindow)=spike_p*ones([1 timewindow+1]);
        end
    else
    end
end
idx_tw=(5*4800);
spike_train{1,1}(1,:)=(spike_train{1,1}(1,:))/2;
spike_train{1,1}(:,1:idx_tw)=[];
clear('idx_sig', 'idx_var', 'var_L','chan_L');
%Frown
timewindow=(0.005*4800);
idx_tw=(5*4800);
idx_tw_end=(55*4800);%length(datamat_T009_Frown)/zfs
for idx_tw=idx_tw:idx_tw_end
    idx_tw=idx_tw+timewindow;
    if (idx_tw+timewindow)<(idx_tw_end)
        clear('idx_p');
        idx_p=find((spikecell{2,5}(2,:))*4800>idx_tw&(spikecell{2,5}(2,...
            :))*4800<idx_tw+timewindow);
        if isempty(idx_p)
            spike_train{2,5}(1,idx_tw:idx_tw+timewindow)=zeros([1 timewindow+1]);
        else
            clear('spike_p')
            spike_p=mean(spikecell{2,5}(5,idx_p));
            spike_train{2,5}(1,idx_tw:idx_tw+timewindow)=spike_p*ones([1 timewindow+1]);
        end
    else
    end
end
idx_tw=(5*4800);
spike_train{2,5}(1,:)=(spike_train{2,5}(1,:))/2;
spike_train{2,5}(:,1:idx_tw)=[];
clear('idx_sig', 'idx_var', 'var_L','chan_L');
%% Plot PSTH
figure()
plot((1:length(spike_train{1,1}))/4800,spike_train{1,1}(1,:))
title('Smile Spike Train PSTH')
xlabel('Time')
ylabel('Probability of Action Potential')
figure()
plot((1:length(spike_train{2,5}))/4800,spike_train{2,5}(1,:))
title('Frown Spike Train PSTH')
xlabel('Time')
ylabel('Probability of Action Potential')

%% Xcorr smilefrown
[acor,lag] = xcorr(spike_train{1,1}(:,1000:10950),spike_train{2,5}(:,:));
[~,I] = max(abs(acor));
lagDiff = lag(I);
timeDiff = lagDiff/4800;
figure()
plot(spike_train{1,1}(:,1000:10950)); hold on;
%idx_np=find((acor/max(acor))<0.9);
%acor(:,idx_np)=0;
plot(lag,acor/max(acor))
title('Smile xcorr frown')
%% Xcorr Smilesmile
[acor,lag] = xcorr(spike_train{1,1}(:,1000:10950),spike_train{1,1}(:,:));
[~,I] = max(abs(acor));
lagDiff = lag(I);
timeDiff = lagDiff/4800;
figure()
plot(spike_train{1,1}(:,1000:10950)); hold on;
%idx_np=find((acor/max(acor))<0.9);
%acor(:,idx_np)=0;
plot(lag,acor/max(acor))
title('Smile xcorr smile')
clear('timewindow', 'timeDiff', 'spike_p', 'lagDiff', 'lag',...
    'idx_tw_end', 'idx_p', 'idx_tw', 'idx_np', 'idx_noise', 'I', 'acor')

%%%%%%%%%%%%%%%%
ai1=find((spikecell{1,1}(5,:))==1);
ti1=spikecell{1,1}(2,ai1);
ti1=floor(ti1*zfs);
ti1e=spikecell{1,1}(3,ai1);
ti1e=ceil(ti1e*zfs);
cm=round(spikecell{1,1}(1,ai1)/max(spikecell{1,1}(1,ai1)),2);
cm=[cm;cm;cm];
ind=0;
xi1mx=0;
for ind=0:length(ti1)-1
    ind=ind+1;
    xi1{ind}=slow_loud_strum(1,ti1(ind)-5:ti1e(ind)+5);hold on;
    if max(xi1{ind})>xi1mx;
        xi1mx=max(xi1{ind});
    else
    end
end
clear('ind');
ind=0;
for ind=0:length(ti1)-1
    ind=ind+1;
    plot3(ones(length(xi1{ind}))*spikecell{1,1}(4,ai1(ind)),...
        (max(xi1{ind})/xi1mx)+(1:length(xi1{ind}))/length(xi1{ind}),...
        xi1{ind},'color',cm(:,ind)); hold on;
    pause(0.1)
end
clear('ind');

ai2=find((spikecell{1,1}(5,:))==2);
ti2=spikecell{1,1}(2,ai2);
ti2=floor(ti2*zfs);
ti2e=spikecell{1,1}(3,ai2);
ti2e=ceil(ti2e*zfs);
cm2=round(spikecell{1,1}(1,ai2)/max(spikecell{1,1}(1,ai2)),2);
cm2=[cm2;cm2;cm2];
ind=0;
xi2mx=0;
for ind=0:length(ti2)-1
    ind=ind+1;
    xi2{ind}=datamat_T008_OSmile(1,ti2(ind)-5:ti2e(ind)+5);hold on;
    if max(xi2{ind})>xi2mx;
        xi2mx=max(xi2{ind});
    else
    end
end
clear('ind');
ind=0;
for ind=0:length(ti2)-1500
    ind=ind+1;
    plot3(ones(length(xi2{ind}))*spikecell{1,1}(4,ai2(ind)),...
        (max(xi2{ind})/xi2mx)+(1:length(xi2{ind}))/length(xi2{ind}),...
        xi2{ind},'color',cm2(:,ind)); hold on;
    pause(0.01)
end
clear('ind');
%%
[stft,freqvector,timevector]=spectrogram(xi2{1}/max(xi2{1}),[],[],[],zfs,'yaxis');

[d,a,costs]=nnmf(abs(spectrogram(xi2{1},...
    [],[],[],zfs,'yaxis')),2);
figure();
plot(a(:,1:2));
figure();
plot(d(1,:));


plot(xi2{10}); hold on;
plot((1:25)+19,xi2{1}(:,fliplr(1:25))); hold on;
plot((1:length(xi2{1}(:,(37:end))))+37-length(xi2{1}(:,(37:end))),xi2{1}(:,fliplr(37:end)));

plot(xi2{1}(:,fliplr(37:end))); hold on;
plot(xi2{1}(:,fliplr(1:6))); hold on;
plot(sum([xi2{1}(:,fliplr(37:44)) ; xi2{1}(:,fliplr(1:8))]));

%%


% figure();
%plot((1:length(diff(ti1)))/length(diff(ti1)),diff(ti1)); hold on ;
%plot((1:length(diff(ti2)))/length(diff(ti2)),diff(ti2));
%
% ind=40;
%         figure();
%      spectrogram(xi2{ind},10,8,4,zfs,'yaxis')