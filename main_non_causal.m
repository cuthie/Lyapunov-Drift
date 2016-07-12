
clc

% The number of SUs
n=2;

% The number of horizon
m=3;

% Sequence of Spectrum Access Decision
a=reshape((dec2bin(0:2^(n*m)-1)-'0')',n,m,2^(n*m));

% Sequence of Spectrum Sensing decision
theta=(dec2bin(0:2^(m)-1)-'0')';

% Parameter for the channel gains between SU and PU
mug=1;

% Parameter for the channel gains between SU and FC
muh=1;

% Parameter for energy being harvested
muH=1;

% Maximum Power limit
P_max=1;

% Spectrum sensing power
p_s=0.1;

% Maximum Battery level
B_max=0.4:0.2:1.2;

% The number of B_max iteration
b_it=length(B_max);

% Number of iteration through H1and H2
nH=1000;

% Sesning time lower bound
tau_l=0.1;

% Time slot length
T=2;

% Interference power limit
Q_avg=10;

% Maximum sum-capacity for initialization
sum_cap_nc=zeros(1,b_it);

% Vector of channel gains for SU and PU
g=exprnd(mug,n,m);

% Vector of channel gains for SU and FC
h=exprnd(muh,n,m);

% Energy harvesting for the SUs
Eng_h=exprnd(muH,n,m);

% Loop for B_max vs throughput plot
for b_lp=1:b_it
    
    % B_max value for each iteration
    b_m=B_max(b_lp);
    
    % Loop for decision to sense
    for l=2^(n*m)
        
        a_mod=a(:,:,l);
        
        % Loop for spectrum sensing result
        for j=1
            
            theta_mod=theta(:,j);
            
            tic
            
            %% Calling the Heuristic subroutine
            [sum_cap_nc(b_lp)]=fh_non_cau_eh_mex...
                (n,m,a_mod,theta_mod,tau_l,T,g,h,P_max,b_m,p_s,Eng_h,Q_avg,nH);
            
            toc
            
            j
            
        end
        
        l
        
    end
    
    b_lp
    
end

figure
hold on
grid on
plot(B_max,sum_cap_nc,'m*-');
hold off
