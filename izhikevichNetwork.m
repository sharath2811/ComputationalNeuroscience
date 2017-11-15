% Created by Eugene M. Izhikevich, February 25, 2003
% Excitatory neurons    Inhibitory neurons
clear all; close all

Ne=800;                 Ni=200;
re=rand(Ne,1);          ri=rand(Ni,1);      % Adding noise to currents
a=[0.02*ones(Ne,1);     0.02+0.08*ri];      % Time scale of recovery varibale 0.02(RS) 0.1(FS) 
b=[0.2*ones(Ne,1);      0.25-0.05*ri];      % Sensitivity of the recovery variable
c=[-65+15*re.^2;        -65*ones(Ni,1)];    % Voltage reset value
d=[8-6*re.^2;           2*ones(Ni,1)];      % Reset of recovery   
S=[0.5*rand(Ne+Ni,Ne),  -rand(Ne+Ni,Ni)];   % Connections between neurons

v=-65*ones(Ne+Ni,1);                        % Initial values of v
u=b.*v;                                     % Initial values of u
firings=[];                                 % Spike timings

for t=1:1000                                % Simulation of 1000 ms
  I=[5*randn(Ne,1);2*randn(Ni,1)];          % Thalamic input (normally distrubuted)
  fired=find(v>=30);                        % indices of spikes
  firings=[firings; t+0*fired,fired];       % Time (t) in first column and neuron number that fired in second
  v(fired)=c(fired);                        % Neurons that fired get voltage reset
  u(fired)=u(fired)+d(fired);               % Neurons that fired get recovery variable reset
  I=I+sum(S(:,fired),2);                    % Sum up inputs (I + fired neurons)
  v=v+0.5*(0.04*v.^2+5*v+140-u+I);  % Update voltage, step 0.5 ms
  v=v+0.5*(0.04*v.^2+5*v+140-u+I);  % for numerical stability
  u=u+a.*(b.*v-u);                          % Update recovery variable 
end;
plot(firings(:,1),firings(:,2),'.');


%% What was the firing rate of neuron 1?
idx=find(firings(:,2)==1);
spkTm=firings(idx,1)/1000;
isi=mean(diff(spkTm));

% What's the mean firing rate of all neurons?
figure();
for i=1:1000;
    idx=find(firings(:,2)==i);
    spkTm=firings(idx,1)/1000;
    FR(i)=1/mean(diff(spkTm));
end;
plot(FR,'ko')
xlabel('Neuron number');
ylabel('Firing rate (Hz)');

% What's the average firing rate at 5 msec?

 
%% How could you reduce the overall firing rate of neurons?
