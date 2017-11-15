clear all
close all
%% PFC - motor cortex (MC) network example

% Set up PFC 
Ne_Pfc=400;                 Ni_Pfc=100;
re=rand(Ne_Pfc,1);          ri=rand(Ni_Pfc,1);      % Adding noise to currents
a=[0.02*ones(Ne_Pfc,1);     0.02+0.08*ri];          % Time scale of recovery varibale 0.02(RS) 0.1(FS) 
b=[0.2*ones(Ne_Pfc,1);      0.25-0.05*ri];          % Sensitivity of the recovery variable
c=[-65+15*re.^2;            -65*ones(Ni_Pfc,1)];    % Voltage reset value
d=[8-6*re.^2;               2*ones(Ni_Pfc,1)];      % Reset of recovery   

% Set up MC 
Ne_MC=400;                  Ni_MC=100;
re=rand(Ne_MC,1);           ri=rand(Ni_MC,1);       % Adding noise to currents
a=[a; 0.02*ones(Ne_MC,1);   0.02+0.08*ri];          % Time scale of recovery varibale 0.02(RS) 0.1(FS) 
b=[b; 0.2*ones(Ne_MC,1);    0.25-0.05*ri];          % Sensitivity of the recovery variable
c=[c; -65+15*re.^2;         -65*ones(Ni_MC,1)];     % Voltage reset value
d=[d; 8-6*re.^2;            2*ones(Ni_MC,1)];       % Reset of recovery 

% PFC - PFC Weight matrix
S1=[0.5*rand(Ne_Pfc+Ni_Pfc,Ne_Pfc),  -rand(Ne_Pfc+Ni_Pfc,Ni_Pfc)];

% MC - MC Weight matrix
S4=[0.5*rand(Ne_MC+Ni_MC,Ne_MC),  -rand(Ne_MC+Ni_MC,Ni_MC)];

% MC -> PFC Weight matrix
NumCon_MP=Ne_MC+Ni_MC;
S2=zeros(NumCon_MP);

% PFC -> MC Weight matrix (Send only excitatory-excitatory connections!)
S3=[0.5*rand(Ne_Pfc+Ni_Pfc,Ne_Pfc), zeros(Ne_Pfc+Ni_Pfc,Ni_Pfc)];

% Put them all together for global weight matrix
S=[S1 S2;S3 S4];
S(501,1)=10;

v=-65*ones(Ne_Pfc+Ni_Pfc+Ne_MC+Ni_MC,1);    % Initial values of v
u=b.*v;                                     % Initial values of u
firings=[];                                 % Spike timings
simTime=50000;                              % Simulation time


%% Creating pulse for the network
szPls=2;
pls=szPls.*[zeros(1,0.45*simTime) ones(1,0.1*simTime) zeros(1,0.45*simTime)];

for t=1:simTime;                            % Simulation time in ms
    %% Stim PFC Only
    I=[4*randn(Ne_Pfc,1); 1*randn(Ni_Pfc,1); 4*randn(Ne_MC,1); 1*randn(Ni_MC,1)];
    pulseIn=[repmat(pls(t),Ne_Pfc+Ni_Pfc,1); zeros(Ne_MC+Ni_MC,1)];
    I=I+pulseIn;
    
    %% Compute voltage and FR
    fired=find(v>=30);                        % indices of spikes
    firings=[firings; t+0*fired,fired];       % Time (t) in first column and neuron number that fired in second
    v(fired)=c(fired);                        % Neurons that fired get voltage reset
    u(fired)=u(fired)+d(fired);               % Neurons that fired get recovery variable reset
    I=I+sum(S(:,fired),2);                    % Sum up inputs (I + fired neurons)
    v=v+0.5*(0.04*v.^2+5*v+140-u+I);          % Update voltage, step 0.5 ms
    v=v+0.5*(0.04*v.^2+5*v+140-u+I);          % for numerical stability
    u=u+a.*(b.*v-u);                          % Update recovery variable
end;

figure(1)
plot(firings(:,1),firings(:,2),'.');
title('Raster Plot');
xlabel('Time');
ylabel('Neuron Number');

figure(2)
nrnNum=size(S,1);
subplot(212);
for i=1:nrnNum;
    idx=find(firings(:,2)==i);
    spkNum(i)=length(idx);
    spkTm=firings(idx,1)/1000;
    FR(i)=1/mean(diff(spkTm));
end;
plot(FR,'ko')
xlabel('Neuron number');
ylabel('Firing rate (Hz)');
title('Mean firing rate');

subplot(211);
spkMtx=zeros(nrnNum,simTime);
for i=1:nrnNum;
    k=find(firings(:,2)==i);
    if ~isempty(k);
        spkTm=firings(k,1);
        spkMtx(i,spkTm)=1;
    end;
end;

%% take mean and smooth data here (notice PFC leading MC)
pfcMnFr=mean(spkMtx(1:500,:));
% pfcMnFr=smooth(pfcMnFr,10);
mcMnFr=mean(spkMtx(501:1000,:));
% mcMnFr=smooth(mcMnFr,10);

plot(1:simTime,pfcMnFr,'b');hold on;
plot(1:simTime,mcMnFr,'r');
ylabel('Firing rate (average spike count)');
xlabel('Time (msec)');
title('Firing rate over time');

selExPair=[spkMtx(1,:)' spkMtx(501,:)'];

% %% Assess xcorr between regions
figure();
mnFr=[pfcMnFr' mcMnFr'];
maxlag=50;
xc=xcorr(mnFr,maxlag);
for i=1:size(mnFr,2)^2;
    subplot(size(mnFr,2),size(mnFr,2),i);
    plot([-maxlag:maxlag],xc(:,i),'k.-');
    axis tight
end;

figure();
maxlag=150;
xc=xcorr(selExPair,maxlag);
bar([-maxlag:maxlag],xc(:,2));

