function[sum_cap_avg] = fh_non_cau_eh(n,m,a,theta,tau_l,T,g,h,P_max,b_m,p_s,Eng_h,Q_avg,nH)

% The throughput optimization problem
% Finite Horizon
% Non-causal CSI and ESI
% Energy harvesting and sharing

%---------------------------------------------------
% Discretized variable initialization
%---------------------------------------------------

% Battery level values
b1=0:0.2:b_m;
b2=0:0.2:b_m;

% Dicrete Power transmission levels
step1=P_max/3;
step2=P_max/3;
p1=0:step1:P_max;
p2=0:step2:P_max;

% The energy sharing discretization level
e12=0:step1:P_max;
e21=0:step2:P_max;

% The energy transfer efficiency
eta=0.8;

% Noise variance
n_var=1;

% Interference variance
inf_var=2;

% Quantized sensing time level for the SU
t=tau_l:0.95:T;

% The total inteference limit
Q=m*Q_avg;

% Discrete bins of interference
q=0:5:Q;

% Probability of Spectrum efficiency
mu=0.8;

% Value Function
V=zeros(m,length(b1),length(b2),length(q));

% The Optimal Transmission Power
P1_star=zeros(m,length(b1),length(b2),length(q));
P2_star=zeros(m,length(b1),length(b2),length(q));

% The optimal energy transfer
e12_star=zeros(m,length(b1),length(b2),length(q));
e21_star=zeros(m,length(b1),length(b2),length(q));

% The optimal sensing time
t_star=zeros(m,length(b1),length(b2),length(q));

% Parameters for the energy harvesting profile
muH1=1;
muH2=1;

% Loop for different horizon
for k=m:-1:1
    
    % Loop for SU1 battery state
    for b1_loop=1:length(b1)
        
        % Loop for SU2 battery state
        for b2_loop=1:length(b2)
            
            % Loop for the interference bins
            for q_loop=1:length(q)
                
                % Initialization of the maximum value
                % function
                V_max=0;
                
                % Initialization of the maximum index
                max_index_p1=1;
                max_index_p2=1;
                max_index_e12=1;
                max_index_e21=1;
                max_index_t=1;
                
                % Loop for SU1 transmission power
                for ind_p1=1:length(p1)
                    
                    % Loop for SU2 transmission power
                    for ind_p2=1:length(p2)
                        
                        % Loop for the sensing time
                        for ind_t=1:length(t)
                            
                            % Loop for energy sharing from
                            % SU1 to SU2
                            for e12_loop=1:length(e12)
                                
                                % Loop for energy sharing
                                % from SU2 to SU1
                                for e21_loop=1:length(e21)
                                    
                                    %Looping conditions
                                    if(e12(e12_loop)~=0 && e21(e21_loop)~=0)
                                        
                                        continue;
                                        
                                    end
                                    
                                    % Checking the energy
                                    % harvesting constraint
                                    % for this horizon SU1
                                    if(a(n-1,k)*((p_s*t(ind_t))+(p1(ind_p1)*(T-t(ind_t))*(1-theta(k))))...
                                            -e12(e12_loop)+(eta*e21(e21_loop))<=b1(b1_loop))
                                        
                                        % Checking the
                                        % energy harvesting
                                        % constaints for
                                        % this horizon SU2
                                        if(a(n,k)*((p_s*t(ind_t))+(p2(ind_p2)*(T-t(ind_t))*(1-theta(k))))...
                                                -e21(e21_loop)+(eta*e12(e12_loop))<=b2(b2_loop))
                                            
                                            % Interference
                                            % term for this
                                            % horizon
                                            
                                            % SU1
                                            inf_1=p1(ind_p1)*g(n-1,k)*a(n-1,k)*(1-theta(k));
                                            
                                            % SU2
                                            inf_2=p2(ind_p2)*g(n,k)*a(n,k)*(1-theta(k));
                                            
                                            % Total
                                            % interference
                                            % term
                                            intf=mu*((T-t(ind_t))/T)*(inf_1+inf_2);
                                            
                                            % Check if the interfernce constraint
                                            % is satisfied or not
                                            if(intf<=q(q_loop))
                                                
                                                % The capacity expression for the present slot
                                                % SU 1
                                                cp1_term=p1(ind_p1)*h(n-1,k)*a(n-1,k)*(1-theta(k));
                                                
                                                % SU 2
                                                cp2_term=p2(ind_p2)*h(n,k)*a(n,k)*(1-theta(k));
                                                
                                                % Present slot sum capacity
                                                % expression
                                                
                                                % When PU
                                                % is
                                                % present
                                                cp_term=(mu/m)*((T-t(ind_t))/T)*(1/log(2))...
                                                    *log(1+(cp1_term+cp2_term)/(n_var+inf_var));
                                                
                                                % PU is
                                                % absent
                                                cp_term_ab=((1-mu)/m)*((T-t(ind_t))/T)*(1/log(2))...
                                                    *log(1+(cp1_term+cp2_term)/(n_var));
                                                
                                                % Total
                                                % sum-capacity
                                                cp_tem=cp_term+cp_term_ab;
                                                
                                                % Initializing the value function
                                                tem=0;
                                                
                                                % If not the last slot
                                                if(k<m)
                                                    
                                                    % q_index for the next slot
                                                    q_next=q(q_loop)-intf;
                                                    [~,closest_index_q]=min(abs(q_next-q));
                                                    
                                                    % Loop for the future
                                                    % value function simulations
                                                    for nxt_ind=1:nH
                                                        
                                                        % Unknown quanties
                                                        H1=exprnd(muH1);
                                                        H2=exprnd(muH2);
                                                        
                                                        % Battery levels if the
                                                        % Control is choosen above
                                                        b1_next=min(b_m,b1(b1_loop)+H1-a(n-1,k)*(p_s*t(ind_t)+...
                                                            (p1(ind_p1)*(T-t(ind_t))*(1-theta(k))))...
                                                            -e12(e12_loop)+(eta*e21(e21_loop)));
                                                        
                                                        b2_next=min(b_m,b2(b2_loop)+H2-a(n,k)*(p_s*t(ind_t)+...
                                                            (p2(ind_p2)*(T-t(ind_t))*(1-theta(k))))...
                                                            -e21(e21_loop)+(eta*e12(e12_loop)));
                                                        
                                                        % Choosing the appropiate battery index for
                                                        % the next slot
                                                        [~,closest_index_b1]=min(abs(b1_next-b1));
                                                        [~,closest_index_b2]=min(abs(b2_next-b2));
                                                        
                                                        % Calualate the value
                                                        % function if the control
                                                        % is chosen above
                                                        
                                                        % Calculate the value
                                                        % function
                                                        tem=tem+V(k+1,closest_index_b1,closest_index_b2,closest_index_q);
                                                        
                                                    end
                                                    
                                                end
                                                
                                                % Modified value function
                                                tem=cp_tem+tem/nH;
                                                
                                                % Test if the previous result is better than
                                                % the previous one
                                                if(tem>V_max)
                                                    
                                                    % In case test passed
                                                    % reset the maximal index
                                                    max_index_p1=ind_p1;
                                                    max_index_p2=ind_p2;
                                                    max_index_e12=e12_loop;
                                                    max_index_e21=e21_loop;
                                                    max_index_t=ind_t;
                                                    
                                                    % Reset the maximal value function
                                                    V_max=tem;
                                                    
                                                end
                                                
                                            end
                                            
                                        end
                                        
                                    end
                                    
                                end
                                
                            end
                            
                        end
                        
                    end
                    
                end
                
                % Save maximal values
                V(k,b1_loop,b2_loop,q_loop)=V_max;
                
                % Save optimal control power input
                P1_star(k,b1_loop,b2_loop,q_loop)=p1(max_index_p1);
                P2_star(k,b1_loop,b2_loop,q_loop)=p2(max_index_p2);
                
                % Save optimal control energy transfered
                e12_star(k,b1_loop,b2_loop,q_loop)=e12(max_index_e12);
                e21_star(k,b1_loop,b2_loop,q_loop)=e21(max_index_e21);
                
                % Save optimal control sensing time
                t_star(k,b1_loop,b2_loop,q_loop)=t(max_index_t);
                
            end
            
        end
        
    end
    
end

% Throughput vector
sum_cap_avg=0;

% Amount of Interference used
q_used=0;

% set initial values for time simulation
% initial battery levels
B1_tem =0.4;
B2_tem =0.4;

% Loop for the Simulation
for k=1:1:m
    
    % Finding the closest index correponding to B,g,h for the look up table
    [~,closest_index_B1]=min(abs(B1_tem-b1));
    [~,closest_index_B2]=min(abs(B2_tem-b2));
    
    % Finding closest index of q in the present slot
    [~,closest_index_Q]=min(abs(Q-q_used-q));
    
    % Get the optiomal values of sensing time from look up table
    t_tem=t_star(k,closest_index_B1,closest_index_B2,closest_index_Q);
    
    % Modified Power limit
    tp1_lim=B1_tem/(T-t_tem);
    tp2_lim=B2_tem/(T-t_tem);
    
    % Get the optimal values of power from the look up table
    p1_tem=min(P1_star(k,closest_index_B1,closest_index_B2,closest_index_Q),tp1_lim);
    p2_tem=min(P2_star(k,closest_index_B1,closest_index_B2,closest_index_Q),tp2_lim);
    
    % Get the optimal values of energy shared from the look up table
    e12_tem=min(e12_star(k,closest_index_B1,closest_index_B2,closest_index_Q),B1_tem-((p_s*t_tem)+(p1_tem*(T-t_tem))));
    e21_tem=min(e21_star(k,closest_index_B1,closest_index_B2,closest_index_Q),B2_tem-((p_s*t_tem)+(p2_tem*(T-t_tem))));
    
    %Interference caused in the present slot
    in1_tem=p1_tem*g(n-1,k)*a(n-1,k)*(1-theta(k));
    in2_tem=p2_tem*g(n,k)*a(n,k)*(1-theta(k));
    in_tem=mu*((T-t_tem)/T)*(in1_tem+in2_tem);
    
    % The amount of intereference used
    q_used=q_used+in_tem;
    
    % Calculate the Throughput
    cp1_tem=p1_tem*h(n-1,k)*a(n-1,k)*(1-theta(k));
    cp2_tem=p2_tem*h(n,k)*a(n,k)*(1-theta(k));
    
    % Throughput when SU is present
    cp1=(mu/m)*((T-t_tem)/T)*(1/log(2))*log(1+(cp1_tem+cp2_tem)/(n_var+inf_var));
    
    % Throughput when SU is absent
    cp2=((1-mu)/m)*((T-t_tem)/T)*(1/log(2))*log(1+(cp1_tem+cp2_tem)/(n_var));
    
    % Total throughput
    cp_tm=cp1+cp2;
    
    % Battery state for the next slot
    B1_tem=min(b_m,B1_tem+Eng_h(n-1,k)-a(n-1,k)*((p_s*t_tem)+(p1_tem*(T-t_tem)*(1-theta(k))))...
        -e12_tem+(eta*e21_tem));
    B2_tem=min(b_m,B2_tem+Eng_h(n,k)-a(n,k)*((p_s*t_tem)+(p2_tem*(T-t_tem)*(1-theta(k))))...
        -e21_tem+(eta*e12_tem));
    
    
    % Optimal value storage
    sum_cap_avg=sum_cap_avg+cp_tm;
    
end

end

