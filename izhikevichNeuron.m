%% Izhikevich simple neuron Chapter 8 Dynamical Systems in neuroscience
clear all                    		
a=0.02;                     % time scale of recovery varibale 0.02(RS) 0.1(FS) 
b=0.2;                      % tensitivity of the recovery variable
c=-50;                      % voltage reset value
d=2;                        % reset of recovery    
vPeak=30

T=1000;tau=0.001; 			% time span and time step (ms)
n=round(T/tau);             % number of simulation steps
v(1)=c; 
u(1)=0;

% puluse of input DC current
stepSize=5;
I=[zeros(1,0.25*n), stepSize*ones(1,0.5*n), zeros(1,0.25*n)];  
%I=[zeros(1,0.1*n), 70*ones(1,0.9*n)];


for i=1:n-1;            % Need to add some burn-in time
   v(i+1)=v(i)+tau*(0.04*v(i).^2+5*v(i)+140-u(i)+I(i));
   u(i+1)=u(i)+tau*a.*(b.*v(i)-u(i));  
   if v(i+1)>=vPeak;	% spike is fired!
	v(i+1)=c;			% membrane voltage rest
	u(i+1)=u(i+1)+d;	% recovery variable update
   end;
end;


% %% Plotting tools
subplot(2,1,1);
plot(tau*(1:n), v, 'k'); hold on; 
xlabel('time (ms)');
ylabel('membrane potential (mV)');
subplot(2,1,2);
plot(u);
xlabel('time (ms)');
ylabel('recovery variable, u');
hold on;

