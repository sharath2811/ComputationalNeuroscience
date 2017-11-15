% Created by Eugene M. Izhikevich, February 25, 2003
% Excitatory neurons    Inhibitory neurons
function [FR]=izhikevichNetwork_withStim(szPls,S)

Ne=2;                   Ni=1;
re=rand(Ne,1);          ri=rand(Ni,1);      % Adding noise to currents
a=[0.02*ones(Ne,1);     0.02+0.08*ri];      % Time scale of recovery varibale 0.02(RS) 0.1(FS) 
b=[0.2*ones(Ne,1);      0.25-0.05*ri];      % Sensitivity of the recovery variable
c=[-65+15*re.^2;        -65*ones(Ni,1)];    % Voltage reset value
d=[8-6*re.^2;           2*ones(Ni,1)];      % Reset of recovery   
% S=[0.5*rand(Ne+Ni,Ne),  -rand(Ne+Ni,Ni)];   % Connections between neurons
S=[0 0 0;1 0 -1;1 0 0].*20 %each column is a neuron, each row is a connection

v=-65*ones(Ne+Ni,1);                       % Initial values of v
u=b.*v;                                     % Initial values of u
firings=[];                                 % Spike timings
simTime=10000;                               % Simulation time 100000

%% Creating pulse for the network
szPls=10;
pls=szPls.*[zeros(1,0.45*simTime) ones(1,0.1*simTime) zeros(1,0.45*simTime)];

for t=1:simTime;                            % Simulation time in ms
  I=[4*randn(Ne,1);2*randn(Ni,1)];   % Thalamic input (normally distrubuted) + pulse
  I=I+[pls(t);0;0];
  fired=find(v>=30);                        % indices of spikes
  firings=[firings; t+0*fired,fired];       % Time (t) in first column and neuron number that fired in second
  v(fired)=c(fired);                        % Neurons that fired get voltage reset
  u(fired)=u(fired)+d(fired);               % Neurons that fired get recovery variable reset
  I=I+sum(S(:,fired),2);                    % Sum up inputs (I + fired neurons)
  v=v+0.5*(0.04*v.^2+5*v+140-u+I);  % Update voltage, step 0.5 ms
  v=v+0.5*(0.04*v.^2+5*v+140-u+I);  % for numerical stability
  v_track(:,t)=v;
  u=u+a.*(b.*v-u);                          % Update recovery variable 
end;


subplot(411)
plot(firings(:,1),firings(:,2),'.');
title('Raster Plot');
xlabel('Time');
ylabel('Neuron Number');


%% What was the firing rate of neuron 1?
idx=find(firings(:,2)==1);
spkTm1=firings(idx,1)/1000;
isi=mean(diff(spkTm1));

idx=find(firings(:,2)==1);
spkTm2=firings(idx,1)/1000;
isi=mean(diff(spkTm2));

xc=xcorr(spkTm1,spkTm2);
plot([-50:50],xc(:,2),'k.-');
% 
% %% What's the mean firing rate of all neurons?
% nrnNum=Ne+Ni;
% subplot(414);
% for i=1:nrnNum;
%     idx=find(firings(:,2)==i);
%     spkNum(i)=length(idx);
%     spkTm=firings(idx,1)/1000;
%     FR(i)=1/mean(diff(spkTm));
% end;
% plot(FR,'ko')
% xlabel('Neuron number');
% ylabel('Firing rate (Hz)');
% title('Mean firing rate')
% 
% %% What's the average firing rate at 5 msec?
% % First, need to turn spike times into a continuous time-series
% subplot(412);
% spkMtx=zeros(nrnNum,simTime);
% for i=1:nrnNum;
%     k=find(firings(:,2)==i);
%     if ~isempty(k);
%         spkTm=firings(k,1);
%         spkMtx(i,spkTm)=1;
%     end;
% end;
% %% take mean and smooth data here
% plot(1:simTime,mean(spkMtx));
% ylabel('Firing rate (average spike count)');
% xlabel('Time (msec)');
% title('Firing rate over time');
% 
% %% Plot Pulse to network
% subplot(413);
% plot(1:simTime,pls,'r');
% 
