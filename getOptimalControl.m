function [u_opt, t, r_min, costError, costTotal, Fth, T, sol] = ...
    getOptimalControl(R0,relTol,antagonistic,makePlot, inputGuess)
% Returns optimal control vector and time stamps found using Pontryagin's
% minimum principle on a linearized version of each individual tcpa.
% Electrothermal dyanmics are evaluated at T = avg(T_min, T_max), where
% thermomechanical model is only valid in the linear region T_min<T<T_max
% such that each tcpa must produce some minimum force > 0. 
% antagonistic [bool] if true, modifies the initialized TCPAs to put them
%           in a symmetric antagonsitic configruation, rather than the
%           default unidirecitonal configuration]
% makePlot [bool] 
% inputGuess [struct] input previous sol to improve computation cost
% u_opt [V^2] nxlength(t) vector of optimal control inputs for n tcpas
% t [s] 1xlength(t) vector of timestamps for control inputs
% r_min [N] scalar minimum force all n tcpas in parallel produce at T=T_min
% costError [N^2*s] integral of force error squared of the linear model
% costTotal [N^2*s] integral of (force error squared + control squared)
% Fth [N] force components of each TCPA as function of time
% T [ºC] temperature of each TCPA as function of time
% sol [struct] solution structure returned by bvp4c solver

persistent print;
if isempty(print)
    print = true;
end

if print
    disp('> ')
    disp('> getOptimalControl.m')
    print = false;
end

  
% Create linear TCPA state space model and input constraints 
[T_min, T_max] = getT_lin(); % Min and max of lin range of thermomech model
T_lin = (T_max+T_min)/2; % [ºC] temp electrotherm model is linearized about
T_amb = getT_amb();
global tcpa
n = tcpa(1).n; % updates n if antagonistic setup is used
for i = 1:n
    % Linearize electrothermal model
    gain = interp1(tcpa(i).temp, tcpa(i).gain, T_lin);
    tau = interp1(tcpa(i).temp, tcpa(i).tau, T_lin);
    % Find vector max and min inputs based on steady state gains and temp
    % contraints
    u_max(i,1) = (T_max - T_amb)/gain;
    u_min(i,1) = (T_min - T_amb)/gain;
    tcpa(i).u_max = u_max(i,1);
    % Construct state matricies based on electrothermal model
    A(i,i) = -1/tau;
    B(i,i) = gain/tau;
    % Linear approx of thermomechanical model between T_min and T_max
    linRegion = T_min < tcpa(i).temp & tcpa(i).temp < T_max;
    coeff = polyfit(tcpa(i).temp(linRegion) - T_amb, ...
        tcpa(i).thermalForce(linRegion), 1);
    C(1,i) = -coeff(1);
    C0(i,1) = -coeff(2);
    % Modify coefficients of antagonisitic actuators
    if antagonistic
        u_max(i+n,1) = (T_max - T_amb)/gain;
        u_min(i+n,1) = (T_min - T_amb)/gain;
        A(i+n,i+n) = -1/tau;
        B(i+n,i+n) = gain/tau;
        C(1,i+n) = coeff(1);
        C0(i+n,1) = coeff(2); 
    end
end
if antagonistic
    n = 2*n;
end

% Hamiltonian H = 0.5 C^2 x^2 - r C x + 0.5 r^2 + 0.5 u^T R u + p^T(Ax + Bu)
% where u_min < u < u_max, x = T - T_amb, and r(t) = ref(t) + r_min

R = R0/n*diag(u_max.^-2); % n normalized for number of muscles
% reference defined w/ respect to minimum allowable for produced by muscles
% where F = C x + sum(C0)  ->  C x = C (T - T_amb) = C (F - sum(C0)),
% and r = ref(t) + r_min = ref(t) + C(T_min - T_amb) + sum(C0)
r_min = C*(T_min - T_amb)*ones(n,1) + sum(C0); % minium force produced by the muscles

% Convert Hamiltonian to a split bounded value problem, where u^2 term
% ensures u is continuous and helps multiple tcpas share a reference w/o
% constantly saturating (i.e. bang-bang control)
%u = @(p) median([u_min, -inv(R.')*B.'*p, u_max],2);
u = @(p) min(max(u_min*ones(1,width(p)),-inv(R.')*B.'*p),u_max*ones(1,width(p)));
% Note y = [T - T_amb; p], where y is 2nx1 and dT and p are nx1
dydt = @(t,y) [A*y(1:n,:) + B*u(y(n+1:2*n,:));
    -C.'*C*y(1:n,:) + C.'*(ref(t)+r_min-sum(C0)) - A.'*y(n+1:2*n,:)];
boundCond = @(y0, yf) [y0(1:n)-T_min+T_amb; % == 0
    yf(n+1:2*n)]; % == 0
[~, t_max] = ref(0);
ref_max = max(ref([0:t_max]));
%guess = @(t) [(T_min - T_amb)*ones(n,1); zeros(n,1)];
%scale guess temperature with reference signal
if isempty(inputGuess)
    guess = @(t) [((T_max-T_min)*ref(t)/ref_max+(T_min - T_amb))*ones(n,1); zeros(n,1)];
    guessSol = bvpinit(linspace(0,t_max,10), guess);
else guessSol = inputGuess; end

% manually define jacobians to decrease solve time
fjac = @(t,y) [ A,          -B*inv(R)*B*diag(u_min <= -inv(R)*B*y(n+1:2*n) &  -inv(R)*B*y(n+1:2*n) <= u_max);
                -C.'*C,     -A.'];
options = bvpset('FJacobian',fjac,'BCJacobian',@bcjac,'NMax',ceil(2E6/n),'RelTol',relTol,'Vectorized','on');

%options = bvpset('NMax',floor(100000),'RelTol',relTol)
%options = bvpset('RelTol',1e-2);
sol = bvp4c(dydt, boundCond, guessSol,options);

t = sol.x;
T = sol.y(1:n,:) + T_amb;
p = sol.y(n+1:2*n,:);
for i = 1:length(t)
    u_opt(:,i) = u(p(:,i));
    effortSq(i) = u_opt(:,i).'*R*u_opt(:,i);
    errorSq(i) = (C*(T(:,i) - T_amb) + sum(C0) - ref(t(i)) - r_min).^2;
end

costError = 0.5*trapz(t, errorSq);
costTotal = costError + 0.5*trapz(t, effortSq);

% if sum(T<T_min-1,'all')
%     warning(strcat('u_min constraints unsuccessful at imposing T>=T_min;',...
%         'possible convergance issue for x = ',num2str(tcpa(1).x)))
% end

Fth = C.'*ones(1,length(t)).*(T-T_amb)+C0;

if makePlot
figure
subplot(5,1,1)
plot(t,Fth); hold on
plot(t,C*(T-T_amb)+sum(C0),'k')
plot(t,ref(t)+r_min,'k--')
title(strcat('Force vs temp for R = ',num2str(diag(R).')))
subplot(5,1,2)
plot(t,T); title('Temp vs t')
subplot(5,1,3)
plot(t,p); title('p vs t')
subplot(5,1,4)
plot(t,u_opt./(u_max*ones(1,length(t)))); title('u/u_{max} vs t')
subplot(5,1,5)
plot(t, errorSq); hold on
plot(t, effortSq); legend('squared error','squared weighted effort')
end

end


function [dBCdy0,dBCdyf] = bcjac(y0,yf)
n = length(y0)/2;
I = eye(n);
O = zeros(n);
dBCdy0 = [I O; O O];
dBCdyf = [O O; O I];
end