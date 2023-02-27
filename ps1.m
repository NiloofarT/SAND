clc;
clear all;
close all;

% adding path and loading data
addpath(genpath("C:\Users\niloofar\Documents\GitHub\teaching\sand\labs"))
%% 1
load("nsa2009_1.mat");
%% 1.a
neuron = 33;
spikeList=data(neuron).spks;

for i=1:length(spikeList)
    spike_count_triall(i)= numel(spikeList{i});
    firing_rate_trial(i)= spike_count_triall(i)/2.048;
end
average_firing_rate_trial=mean(firing_rate_trial);

%% 1.b
nboot = 10000;
bootstat = bootstrp(nboot,'mean',firing_rate_trial);
SE=std(bootstat);
% bootstat = bootstrp(nboot,@(spikeList) std(cat(3,spikeList{:}),[],3),spikeList);
%SE = bootstat/sqrt(nboot);
%% 1.c
% trials 1-75
for i=1:75
    spike_count_triall1(i)= numel(spikeList{i});
    firing_rate_trial1(i)= spike_count_triall1(i)/2.048;
end
average_firing_rate_trial1=mean(firing_rate_trial1);

% trials 76-150
for i=76:150
    spike_count_triall2(i-75)= numel(spikeList{i});
    firing_rate_trial2(i-75)= spike_count_triall2(i-75)/2.048;
end
average_firing_rate_trial2=mean(firing_rate_trial2);
%% 1.d
diff= firing_rate_trial2 - firing_rate_trial1;
[h,p] = ttest(diff,0,'Alpha',0.01);
%% 1.e
FF = var(spike_count_triall)/mean(spike_count_triall)
%% 2
clear all;
clc;
addpath(genpath("C:\Users\niloofar\Documents\GitHub\teaching\sand\labs"))
load("data_v1_binned_moving.mat");
%% 2.a
theta = linspace(0,2*pi,100);
b=1;
a=1;
mu=0;
sig=1;
param(1)=b;
param(2)=a;
param(3)=mu;
param(4)=sig;
figure();
plot(theta,Rcgs(param,theta))
hold on
sig=0.5;
param(4)=sig;
plot(theta,Rcgs(param,theta))
hold on
sig=1.5;
param(4)=sig;
plot(theta,Rcgs(param,theta))
title("circular Gaussian curve")
legend("sig = 1","sig = 0.5","sig = 1.5")
%% 2.b
dataset_number = 1;
neuron_number = 1;

S = data{dataset_number}.spikes(neuron_number,:,:,:);
S = squeeze(S); % make data cube from the original 4d data (only look at 1 neuron for now)
% Sum spike counts over time so that we have a single tuning function
spike_counts = sum(S(:,25:80,:),2);   % only take the data from time bins 25-80
spike_counts = squeeze(spike_counts); % make data cube into a matrix [direction x trials]

stim_direction = linspace(0,2*pi-2*pi/16,16);
% fit model
x_matrix = repmat(stim_direction',1,13); % first we need to make sure 'x' and 'y' are the right sizes and lined up
options=[]; % for now just use default options (see help fminsearch for more)
b_Rcgs = fminsearch('RcgsCost',[0.5 0.5 0.1 0.5],options,x_matrix(:),spike_counts(:)); % do the optimization
b_mse = fminsearch('vonMisesCost',[1 0.1 pi],options,x_matrix(:),spike_counts(:)); % do the optimization
% plot
figure()
scatter(x_matrix(:)*180/pi,spike_counts(:),spike_counts(:)*0+60,'filled','MarkerFaceAlpha',0.4)     % plot the data
hold on
x0 = linspace(-10,345,256)*pi/180;
plot(x0*180/pi,Rcgs(b_Rcgs,x0),'LineWidth',1) % plot the prediction
hold on
plot(x0*180/pi,vonMises(b_mse,x0),'LineWidth',1) % plot the prediction
hold off
xlim([min(x0) max(x0)]*180/pi)
box off; set(gca,'TickDir','out')
ylabel('Spike Rate')
xlabel('Stimulus Direction [deg]')
title("Tuning curve Neuron 28 ")
legend("datapoints","Gaussian","vonMises")
%% 2.c  
MSE_Gauss_N1 = RcgsCost(b_Rcgs,x_matrix(:),spike_counts(:))
MSE_vonMises_N1 = vonMisesCost(b_mse,x_matrix(:),spike_counts(:))

%% 3
clear all;
clc;
addpath(genpath("C:\Users\niloofar\Documents\GitHub\teaching\sand\labs"))
load('chb_sample.mat');
%% 3.a
fs=256;
t=0:1/fs:10;
figure();
for i=1:length(data)
subplot(length(data),1,i)
plot(t,data{1,i}(1:fs*10+1))
ylabel("chl" + i)
end
xlabel("t(s)")
%% or
figure();
for i=1:length(data)
plot(t,data{1,i}(1:fs*10+1)+i*250)
hold on
%plot(t,repmat((i*250/2),length(t),1),'k')
grid on;
end
xlabel("t(s)")
title("first 10 second of EEG signal for 23 Channels")
%% 3.b
CE = ones(length(data),length(data));
for ii=1:length(data)
    for jj=1:length(data)
       C = corrcoef(data{1,jj},data{1,ii});
       CE(ii,jj) = C(1,2);
    end
end

%% plot
figure();
% B = rot90( CE )
imagesc(CE)
axis xy
colorbar
cl1 = get(gca,'CLim');
xlabel('channel')
ylabel('channel')
title("Correlation between electrodes")
set(gca,'Ydir','reverse')
%% 3.c
y=data{1,9};
%Spectrogram
figure()
subplot(3,3,[5 6 8 9])
[S,f,t] = spectrogram(y,4096,4096-100,512,256);
imagesc(t,f,log(abs(S)))
axis xy
xlabel('Time (s)'); ylabel('Frequency (Hz)')
ylim([0 50])
set(gca,'TickDir','out'); box off
subplot(3,3,[2 3])
plot(t,sum(abs(S).^2,1),'LineWidth',2)
axis tight; set(gca,'TickDir','out'); box off
subplot(3,3,[4 7])
plot(f,sum(abs(S),2),'LineWidth',2)
axis tight; box off; set(gca,'TickDir','out')
xlim([0 50])
camroll(90);
