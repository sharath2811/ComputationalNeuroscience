
%% integrate Pulse wilson RS with EPSP

%%  Integration of Wilson model with the Euler method
clear all; clf;

neuronType='RS';

%% Parameters of the model: 1=K,R  2=Ca,T  3=KCa,H  4=Na
switch neuronType
    case 'bursty'; g(1)=26; g(2)=2.25; g(3)=9.5; g(4)=1;  % Bursty Neuron 
    case 'RS';     g(1)=26; g(2)=0.1; g(3)=5; g(4)=1;     % RS 
    case 'FSI';    g(1)=26; g(2)=0.25; g(3)=0; g(4)=1;    % FSI
end;
 
% % Battery volatge based on [Ion]
E(1)=-.95; E(2)=1.20; E(3)=E(1); E(4)=.50;


%% Setting some constants and initial values for EPSP
E_syn=10;
tau_syn=4; 

 
%% Initial values
dt=0.01;
I_ext=0;
V=-1;
x=zeros(1,4);
switch neuronType
    case {'bursty','RS'}; tau(1)=dt./4.2; tau(2)=dt./14; tau(3)=dt./45; tau(4)=1;
    case 'FSI'; tau(1)=dt./1.5; tau(2)=dt./14; tau(3)=dt./45; tau(4)=1;
end;

plsTime=[50 55];

%% Integration
t_rec=0;m=1;
for t=-100:dt:200;                  % time loop 
    if m>2
        if ismember(t,plsTime); g_syn(m-1)=1;end;
        g_syn(m)= g_syn(m-1) - dt/tau_syn * g_syn(m-1);
        I_ext=-g_syn(m) * (V-E_syn);
    else
        g_syn(m)=0;
    end;
    
    % nonlinear conductances
    x0(1)=1.24  +  3.7*V + 3.2*V^2;  % K,R
    x0(2)=4.205 + 11.6*V + 8  *V^2;  % Ca,T
    x0(3)=3*x(2);                    % KCa,H
    x0(4)=17.8  + 47.6*V +33.8*V^2;  % Na
    
    x=x-tau.*(x-x0);                 % integrate conductances
    I=g.*x.*(V-E);                   % integrate currents
    V=V+dt*(I_ext-sum(I));           % integrate voltage
    
    if t>=0;
        t_rec=t_rec+1;
        x_plot(t_rec)=t;
        y_plot(t_rec)=V;
        I_plot(t_rec)=I_ext;
        trackGate(:,t_rec)=x0;
    end
    m=m+1;
end 
%% Detect spikes
[pks,locs]=findpeaks(100*y_plot,'MinPeakHeight',-30);

%% Plotting reults

plot(x_plot,100*y_plot); xlabel('Time (ms)'); ylabel('Membrane potential');hold on;
plot(x_plot(locs),pks,'gs');hold on;
plot(x_plot,I_plot);

% Getting FR
ISI=diff(x_plot(locs));
FR=1/(mean(ISI)/1000);
title([neuronType ': FR=' num2str(FR)]);

%% Cut and paste below to run as function
% [x_plot,y_plot,I_plot,pks,locs]=wilson_euler_2(0.4,'FSI',1);

