%set up network variables
Ne=800;                  Ni=200;
re=rand(Ne,1);          ri=rand(Ni,1);      % Adding noise to currents
a=[0.02*ones(Ne,1);     0.02+0.08*ri];      % Time scale of recovery varibale 0.02(RS) 0.1(FS) 
b=[0.2*ones(Ne,1);      0.25-0.05*ri];      % Sensitivity of the recovery variable
c=[-65+15*re.^2;        -65*ones(Ni,1)];    % Voltage reset value
d=[8-6*re.^2;           2*ones(Ni,1)];      % Reset of recovery

v=-65*ones(Ne+Ni,1);                        % Initial values of v
u=b.*v;                                     % Initial values of u
firings=[];
weights=zeros(5000,1000);

% set up connection weights
S=[0.5*rand(Ne+Ni,Ne),  -rand(Ne+Ni,Ni)];   % Connections between neurons

for t=1:1000;                            % Simulation time in ms
    I=[5*randn(Ne,1);2*randn(Ni,1)];   % randomly distributed input
%   pls =[ones(200,1);zeros(600,1);ones(50,1);zeros(150,1)];  %randomly generated pulse
%   I=I.*pls;

  fired=find(v>=30); 
  %if ~isempty(fired)
      %perform hebbian learning
      %increase the connection weights by a 0.1
      firings=[firings; t+0*fired,fired];
      if fired>800
      S(fired,:)=S(fired,:)-0.1;
      else
      S(fired,:)=S(fired,:)+0.1;
      end;
  %end;
  
  v(fired)=c(fired);                        % Neurons that fired get voltage reset
  u(fired)=u(fired)+d(fired);               % Neurons that fired get recovery variable reset
  I=I+sum(S(:,fired),2);                    % Sum up inputs (I + fired neurons)
  v=v+0.5*(0.04*v.^2+5*v+140-u+I);  % Update voltage, step 0.5 ms
  v=v+0.5*(0.04*v.^2+5*v+140-u+I);
  u=u+a.*(b.*v-u);                          % Update recovery variable 
end;

maxthresh =96; %set threshhold for connection weights
minthresh =-95.5; %set threshhold for inhibitory neurons
for i=1:1000;
    for j=1:1000;
        if S(i,j)<maxthresh;
            S(i,j)=0; %set every connection weight that falls below the threshold to 0
        end;
    end;
end;

for i=1:1000;
    for j=1:1000;
        if S(i,j)>minthresh && S(i,j)<0; %removing inhibitory connections
            S(i,j)=0; %set every connection weight that falls below the threshold to 0
        end;
    end;
end;

%Only nodes that have strong connection remain is the network