%% Muscle plots
%data from electrode pairs red and blue

mus_fr_win=((1/50)*zfs); %muscle max fire rate (50hz) window 
chan_time=(1:(length(slow_loud_strum(1,:))))/zfs ; %Samples to time

muscle1=sum([(slow_loud_strum(6,:));(slow_loud_strum(7,:))]); %muscle 1
muscle2=sum([(slow_loud_strum(8,:));(slow_loud_strum(9,:))]); %muscle 2


muscle1=((muscle1-mean(muscle1))/std(muscle1))/max(abs((muscle1-mean(muscle1))/std(muscle1))); %center, standardise deviation(gain) and normalise
muscle2=((muscle2-mean(muscle2))/std(muscle2))/max(abs((muscle2-mean(muscle2))/std(muscle2)));

muscle1=envelope(muscle1,mus_fr_win,'rms'); %RMS of muscle
muscle2=envelope(muscle2,mus_fr_win,'rms');

muscle_ntens=median([muscle1;-muscle2]); %net tension across muscles
muscle_tone=sum([muscle1;muscle2]); %total muscle tension
muscle_ntens=muscle_ntens/max(muscle_ntens); %Normalise
muscle_tone=muscle_tone/max(muscle_tone);

figure();
plot(chan_time,muscle1,'b','lineWidth',1.1)
hold on;
plot(chan_time,muscle2,'r','lineWidth',1.1)
hold on;
plot(chan_time,muscle_ntens,'k','lineWidth',1.1)
hold on;
plot(chan_time,muscle_tone,'g','lineWidth',1.1)
grid minor

hand_pos=muscle_ntens.*abs(muscle_tone); %hand posisiton=net tension weighted by median muscle tone
hand_pos=((hand_pos-mean(hand_pos))/std(hand_pos))/max(abs((hand_pos-mean(hand_pos))/std(hand_pos)));

hand_vel=(gradient(hand_pos)); %hand velocity
hand_vel=((hand_vel-mean(hand_vel))/std(hand_vel))/max(abs((hand_vel-mean(hand_vel))/std(hand_vel)));

figure();
%plot(chan_time,muscle_ntens,'k','lineWidth',1.1)
%hold on;
plot(chan_time,muscle_tone,'g','lineWidth',0.3)
hold on;
plot(chan_time,hand_vel,'b','lineWidth',0.3)
hold on;
plot(chan_time,hand_pos,'r','lineWidth',0.3)
grid minor
xlim([38.7 39.6]);
ylabel('Normalised scale');
xlabel('Time(seconds)');
hold on;
legend('muscle tone','velocity','position','Nerve1 Rate','Nerve5 Rate');
set(legend,'FontSize',13);


%% Nerve

nerve1=slow_loud_strum(1,:).*hand_vel; %nerve 1
nerve1=((nerve1-mean(nerve1))/std(nerve1))/max(abs((nerve1-mean(nerve1))/std(nerve1)));
nerve1_rate=envelope(nerve1,mus_fr_win,'rms');
nerve1_rate=nerve1_rate/max(nerve1_rate);
hold on;
plot(chan_time,nerve1_rate,'--','lineWidth',1);

nerve5=slow_loud_strum(5,:).*hand_vel; %nerve 5
nerve5=((nerve5-mean(nerve5))/std(nerve5))/max(abs((nerve5-mean(nerve5))/std(nerve5)));
nerve5_rate=envelope(nerve5,mus_fr_win,'rms');
nerve5_rate=nerve5_rate/max(nerve5_rate);

hold on;
plot(chan_time,nerve5_rate,'--','lineWidth',1);

%plot(chan_time,sum([nerve5_rate;nerve1_rate]),'--','lineWidth',1);






%Plot all
figure();
idx_a=0;
for idx_a=0:8
    idx_a=idx_a+1;
plot(chan_time,slow_loud_strum(idx_a,:)+idx_a,'lineWidth',0.5);hold on;
end


hold on;
plot(chan_time,muscle_tone,'g','lineWidth',1.3)
hold on;
plot(chan_time,hand_vel,'b','lineWidth',1)
hold on;
plot(chan_time,hand_pos,'r','lineWidth',1.3)
grid minor





%% Audio 
sls_audio=audioread('DSC_0025.mov')';%slow loud strum
sls_audio=diff([sls_audio(1,:);-sls_audio(2,:)])/2;

%% Perform and plot cross correlation
emg=(abs(hand_vel));
emg=envelope(emg,4800,'rms');
emg=emg-mean(emg);
emg=emg/max(emg);
aud=envelope(abs(sls_audio(1,:)),4800,'rms');
aud=aud/max(aud);
[acor,lag] = xcorr(emg,aud);
[~,I] = max(abs(acor));
lagDiff = lag(I);
timeDiff=lagDiff/48000;
figure()
plot((1:length(envelope(hand_vel,4800,'rms')))/48000,envelope(abs(hand_vel),4800,'rms')); hold on;
%plot((1:length(hand_vel))/(48000),hand_vel);

% idx_np=find((acor/max(acor))<0.9);
% acor(:,idx_np)=0;
plot(lag/48000,acor/max(acor))

plot(timeDiff+(1:length(aud))/48000,aud); hold on;
plot((1:length(emg))/48000,emg); hold on;

figure();
flp_time=(1:length(sls_audio))/zfs;
plot(flp_time,sls_audio(1,:),'b');
hold on;
plot(chan_time-timeDiff,hand_vel,'g')
hold on;
plot(chan_time-timeDiff,hand_pos,'k')
hold on;
plot(chan_time-timeDiff,muscle_ntens,'r')
hold on;
plot(chan_time-timeDiff,(abs(muscle_ntens).*abs(hand_vel)),'r')
grid minor

%% View Audio
zfs=48000;
segmentTime=1*10^-3;
segmentLength=2*round(segmentTime*zfs/2)+1;
window = hann(segmentLength);
nOverlap =(segmentLength+1)/2; 
nDft = 4800;
spectrogram(slow_loud_strum(1,:),window, nOverlap, nDft, zfs,'yaxis')
[stft1, freqVector1, timeVector1] = spectrogram(slow_loud_strum(1,:), ...
    window, nOverlap, nDft, zfs);
powerSpectrum1 = abs(stft1).^2/length(window);
[dictionary1, activationMatrix1] = nnmf(powerSpectrum1, 10);
freqRange = [60, 1200];
dynRangeDb = 60;
figure()
nmfEquationPlot(dictionary1, activationMatrix1, zfs, ...
segmentLength, nOverlap, powerSpectrum1, freqRange, dynRangeDb)



%% Create set of probability vector for both muscles.
tens_muscle1=(envelope(diff([abs(slow_loud_strum(6,:));abs(slow_loud_strum(7,:))]),80,'rms'));
tens_muscle2=(envelope(diff([abs(slow_loud_strum(8,:));abs(slow_loud_strum(9,:))]),80,'rms'));
tens_muscle1=tens_muscle1/max(tens_muscle1);
tens_muscle2=-1*tens_muscle2/max(tens_muscle2);
plot(tens_muscle1); hold on;
plot(tens_muscle2);

hold on;
plot(diff([tens_muscle1;tens_muscle2])*-1);
%% Add phase code probability to spikecell
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
% Create set of probability vector for both muscles.
tens_muscle1{idx_chan,:}=(envelope(diff([abs(curr_var(6,:));abs(curr_var(7,:))]),mus_fr_win,'rms'));
tens_muscle2{idx_chan,:}=(envelope(diff([abs(curr_var(8,:));abs(curr_var(9,:))]),mus_fr_win,'rms'));
tens_muscle1{idx_chan,:}=(tens_muscle1{idx_chan,:})/max(tens_muscle1{idx_chan,:});
tens_muscle2{idx_chan,:}=(-1*tens_muscle2{idx_chan,:})/max(tens_muscle2{idx_chan,:});
%net_tense{idx_chan,:}=
%         idx_sig=0;
%         for idx_sig=0:chan_L-1;
%             idx_sig=idx_sig+1;
%             spikecell{idx_var,idx_chan}(4,idx_sig)=max((curr_data(1,floor(spikecell{idx_var,idx_chan}(2,...
%                 idx_sig)*fs):ceil(spikecell{idx_var,idx_chan}(3,idx_sig)*fs))));
%         end
%        spikecell{idx_var,idx_chan}(7,:)=spikecell{idx_var,idx_chan}(7,:)/max(spikecell{idx_var,idx_chan}(7,:)); %Normalise probability
       %spikecell{idx_var,idx_chan}(5,:)=spikecell{idx_var,idx_chan}(5,:)/max(spikecell{idx_var,idx_chan}(5,:))*280; %Normalise amplitude
        %spikecell{idx_var,idx_chan}(4,:)=round(spikecell{idx_var,idx_chan}(4,:)/0.01)*0.01;
    end
end
clear( 'curr_var', 'curr_var_L', 'curr_var_B', 'idx_chan',...
    'idx_var', 'name_var', 'work_L', 'worksp_idx','var_B',...
    'chan_L', 'idx_sig', 'idx_spike', 'var_L','curr_data');