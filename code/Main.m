clear all
rng('default')
K=40; %number of terminals
M=200; %number of APs
Nm=1; %number of antennas/AP
B=20; %Mhz
T=200;
D=1; %in kilometer
tau=K;
[U,S,V11]=svd(randn(tau,tau));%U includes tau orthogonal sequences 
%U=ones(tau,tau);
Hb = 15; % Base station height in m
Hm = 1.65; % Mobile height in m
f = 1900; % Frequency in MHz
aL = (1.1*log10(f)-0.7)*Hm-(1.56*log10(f)-0.8);
L = 46.3+33.9*log10(f)-13.82*log10(Hb)-aL;

power_f=1; %downlink power: 1W
noise_p = 10^((-203.975+10*log10(B*10^6)+9)/10); %noise power
Pd = power_f/noise_p;%nomalized receive SNR
Pp=0.2/noise_p;%pilot power
sigma_shd=8; %in dB

d0=0.01;%km
d1=0.05;%km

save('../data/system_paramater.mat','M','Nm','K','D','L','Pd',...
    'noise_p','Pp','sigma_shd','d0','d1','tau','U','B','T');

xi=0.1; %penalty coefficient


%Power consumption parameters:
myalpha=(1/0.4)*ones(M,1);
P_fix=0;
P_tc=0.2*ones(M,1);
P_bt=0.25*10^(-3)*ones(M,1);
P_0=0.825*ones(M,1);
P_fix_bar=P_fix + Nm*sum(P_tc) + sum(P_0);


% Large-scale fading matrix
BETAA=get_slow_fading('../data/system_paramater.mat');
Phii=zeros(tau,K);
if tau<K
    Phii(:,1:1:tau)=U;
    for k=(tau+1):K
    Point=randi([1,tau]);
    Phii(:,k)=U(:,Point);
    end
else
    Phii=U(:,1:1:K);
end
Phii_cf = Phii; % pilot set of cell-free systems

% Create Gamma matrix
mau=zeros(M,K);
for m=1:M
    b = BETAA(m,:);
    for k=1:K
        mau(m,k)=norm( b.^(1/2).*(Phii_cf(:,k)'*Phii_cf))^2;
    end
end

Gammaa=tau*Pp*(BETAA.^2)./(tau*Pp*mau + 1);
etaa1=1./(Nm*sum(Gammaa,2)); %Consider equal power allocation
etaa = repmat((etaa1),1,K); 

So = ones(K,1); % QoS constraint


mytheta_prev = sqrt(etaa).*sqrt(Gammaa);
t_now=1;
t_prev=1;
z=mytheta_prev;
y_prev=mytheta_prev;
mytheta_now=mytheta_prev;
v_prev = 1.1*mytheta_prev;

maxPGMIterations = 1000; % maximum number of iteraitons for the inner loop (i.e. the APG method)
maxOuterIterations = 2; % maximum number of iterations of the outer loop


mynu = 0.5; % \nu parameter in Algorithm 2
mydelta = 1e-5; % \delta parameter in Algorithm 2

bestobj = zeros(maxPGMIterations*maxOuterIterations,1);
interval = 10;
m = 1;
for iOuter=1:maxOuterIterations % outer loop
    for iPGM=1:maxPGMIterations % Inner loop (i.e. APG method)

        y = mytheta_now+(t_prev/t_now)*(z-mytheta_now)...
            +((t_prev-1)/t_now)*(mytheta_now-mytheta_prev);
        
        % line search for y
        [z,stepsize_y]  = linesearch_y(xi,So,B,tau,T,M,Nm,K,Gammaa,BETAA,Phii_cf,Pd,...
            P_fix_bar,noise_p,myalpha,mynu,mydelta,y,y_prev,z);
        y_prev = y;


        % Line search for theta
        [v,stepsize_theta] = linesearch_theta(xi,So,B,tau,T,M,Nm,K,Gammaa,BETAA,Phii_cf,Pd,...
            P_fix_bar,noise_p,myalpha,mynu,mydelta,mytheta_now,mytheta_prev,v_prev);

        v_prev = v;

        mytheta_prev = mytheta_now; % Update mytheta_prev

        Fz = computeEnergyEfficiency(xi,So,B,tau,T,M,Nm,K,(Gammaa),BETAA,Phii_cf,Pd,z,...
            P_fix_bar,Pd*noise_p*Nm,myalpha);
        Fv = computeEnergyEfficiency(xi,So,B,tau,T,M,Nm,K,(Gammaa),BETAA,Phii_cf,Pd,v,...
            P_fix_bar,Pd*noise_p*Nm,myalpha);

        if(Fz>=Fv)
            mytheta_now = z;
        else
            mytheta_now = v;
        end


        bestobj(m)=max(Fz,Fv); % save best objective

        [~,SumRate] = gradPenalty(tau,T,M,Nm,K,Gammaa,BETAA,Phii_cf,Pd,mytheta_now,So);
        
        SumPenal = SumRate; % Total penalty function
        
        t_prev = t_now;
        t_now = (sqrt(4*t_prev^2+1)+1)/2;
        
        if (iPGM>interval)
            tol = (bestobj(m)-bestobj(m-interval))/bestobj(m-interval);
            if (abs(tol)<1e-3) % if APG converges, break
                break
            end

        end
        m = m+1;
    end % END of APG method
    if(SumPenal<1e-4) % if feasible point achieved, break the outer loop
        bestobj(m:end)=[];
        break
    else % if not, increase penalty parameter
        xi = xi*10;
    end
end
EE_true = 1./(1./bestobj + sum(P_bt)); % Convert into true EE

hold on 
xlabel('Iteration');
ylabel('Total Energy Efficiency (bits/Joule)');
plot(EE_true, 'r')
saveas(gcf,'../results/Algorithm1_convergence.png')


