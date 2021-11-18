%Basic script: Calculate H (connection probability), Q (capture probability) ,C (coverage probability) probabilities analytically and simulate them
%before and after SIC. Diagrams H,Q,C,Csic - distance

clear
%system parameters
P = 14; % constant power, end-devices transmit with, in dBm
P = power(10, P / 10) / 1000; %convert from dBm
NF = 6; %% receiver noise figure in dBm
BW = 125; % bandwidth of a single channel for uplink, in kHz
NdB = -174 + NF + 10 * log10(BW * 1000); % AWGN noise variance, in dBm
N = power(10, NdB / 10) / 1000; %convert from dBm
w = 1; %(dB) the aggregate SIR ratio required (threshold) for a successful uplink transmission, for co-SF collisions
w = power(10, w / 10); %convert from dB to pure number
f = 868.1e6; % carier frequency in Hz
wavelength = 3e8 / f;%% carrier wavelength in m
n = 3; % path loss exponent: n=2.7 for sub-urban & n=4 for urban areas
kappa = (wavelength / (4 * pi))^2;
qo = [-6, -9, -12, -15, -17.5, -20]; % distance dependent sensitivity, in dBm
l = 2 * [0, 0.5, 1, 1.5, 2, 2.5, 3]; %radius for cyclic regions (between GW-ED), in km
p = 0.0033; %duty cycle

R = l(7); %radius
range = R / 6; % l_i - l_i-1 : witdh of annuli
EDs = 1500; %number of EDs
V = pi * R^2; %total area of EDs
lambda = EDs / V; %density of PPP for EDs

%Calculate the number of EDs per SF, Ned and the density of PPP for the
%active EDs in ring i, a_i
Ned = zeros(1, 6);
a = zeros(1, 6);

for i=1:6
    Ned(i) = lambda * pi * ((l(i+1))^2 - (l(i))^2);
    a(i) = 2 * p * Ned(i); 
end

%simulation parameters 
trials = 1e4; %1e5 for official results

%doo is the variable distance of the main node (desirable)
doo = 0:0.2:R;

%initialization
Pr0 = zeros(1, length(doo));
Pr1 = zeros(1, length(doo));
Pr_bigger_than_one = zeros(1, length(doo));
Hdi = zeros(1, length(doo));
Hdi_sim = zeros(1, length(doo));
Qdi = zeros(1, length(doo));
Qdi_sim = zeros(1, length(doo));
Qdk_1 = zeros(1, length(doo));
Qdk_2 = zeros(1, length(doo));
Qdk_1sim = zeros(1, length(doo));
Qdk_2sim = zeros(1, length(doo));
Qsic_sim = zeros(1, length(doo));

for i=1:length(doo)
    
    if doo(i) == 0
        start = 1;
    else
        start = ceil(doo(i)/range);
    end
    inter_signals = Ned(start);
    
    %Probability if there is no interference
    Pr0(i) = exp(-a(start));
    %Probability if there is one interferer 
    Pr1(i) = a(start) * exp(-a(start));
    %Probabilityif there are more than one interferers
    Pr_bigger_than_one(i) = 1 - Pr0(i)-Pr1(i);
    
    %Path loss attenuation function of the current distance
    g_doo = kappa * (doo(i) * 1000) ^ (-n);
    
    %Calculate the theoretical curve for SNR (eq. 3.9)
    qsf = power(10, qo(start) / 10);
    Hdi(i) = exp(-N * qsf / (P * g_doo));
    
    % Calculate the theoretical curve for SIR (eq. 3.12)
    % Qdi = Qdsic
    Qdi(i) = exp(-2 * pi * 2 * p * lambda *integral(@(x) x.*(1-1./(1 + w.*(doo(i)./x).^n)),l(start),l(start+1)));
    Qdi(1) = 1;
   
    % Calculate the theoretical curve for SIR for interferer k*
    %If interfering=1 (eq. 3.36)
    const = 2 * (doo(i)^n) / (l(start+1)^2 - l(start)^2);
    Qdk_1(i) = const * integral(@(x) x./(w.*x.^n+doo(i).^n),l(start),l(start+1));
    Qdk_1(1) = 0;
    %If interfering>1 (eq. 3.39)
    const = 2 / (l(start+1)^2 - l(start)^2);
    Qdk_2(i) = const * integral(@(dk) dk./(w.*(dk./doo(i)).^n+1).*exp(-2*pi*2*p*lambda*integral(@(x) x./(1/w.*(x./dk).^n+1),l(start),l(start+1))) ,l(start),l(start+1));   
    Qdk_2(1) = 0;
    
    %Create the channel gain |hi|^2 for the desirable signal (modelled by an exponential random variable of mean 1).
    channel = exprnd(1,1,trials);
    
    % Computation: P[SNRi>=qsfi] = p[|hi|^2 >= N*qsf/(p*g)] - Simulation of connection probability
    c = length(find(channel>=N*qsf/(P*g_doo)));
    Hdi_sim(i) = c / trials;
    

    %Simulation of capture probability
    %initializations
    counter = zeros(1, trials);
    counter_dom = zeros(1, trials);
    counter_sic = zeros(1, trials);
    counter_sic_only_one = zeros(1, trials);
    for u=1:trials
        %Calculate the number of interfering signals - PPP for active EDs
        %in ring i, with average number a_i
        interfering = poissrnd(a(start), 1) ;
        
        %Create the channel gain |hk|^2 between each ED k (interfering signal) and GW
        h_interfering = exprnd(1, 1, interfering);
        
        %Create random distances in km for the interfering signals
        if start == 1
            dis = l(start) + (l(start+1)-l(start)) * sqrt(rand(1, interfering));
        elseif start > 1
            dis = sqrt(l(start)^2 + (l(start+1)^2-l(start)^2) * rand(1, interfering));
        end

        %Calculate the path loss attenuation function for each interferer
        g = kappa * (dis * 1000).^(-n);   
        
	%Calculate the mean value of I = sum(P * |hk|^2 * g(dk))
        I = P * sum(h_interfering.*g);
       
       	%Computation: P[SIRi>=w] = p[|hi|^2 >= I*w/(P*g)]
        counter(u) = length(find(channel>=w*I/(P*g_doo)));

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %                            SIC                               %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    

	%Simulate Qdk1
        if interfering == 1
            counter_sic_only_one(u) = length(find(h_interfering>=w*(channel*g_doo)./g));
        %Simulate Qdk2
        elseif interfering > 1
	    %Choose random one interferer
            r = randi(interfering);
            dominant = h_interfering(r).*g(r);
            %Total interference - dominant (chosen) interferer
            rest_inter = I - P * dominant;
            %Examine if the interferer exceeds the SIR threshold
            counter_dom(u) = length(find(P*dominant>=w*(rest_inter+P*channel*g_doo)));
            %Computation: P[SIRi_sic>=w] = p[|hi|^2 >= I_sic*w/(P*g)]
	    %Calculate if the desirable signal exceeds the SIR threshold after SIC
            counter_sic(u) = length(find(channel>=w*rest_inter/(g_doo*P)));
            
        end
    end
    
    %SIR computational- for desirable signal
    Qdi_sim(i) = sum(counter) / trials^2; % Computation P[SIRij>=qsfi] = p[|hij|^2 >= w*I/(p*g)]
   
    %SIR threshold for dominant interferer - Computational
    %If interfering=1
    Qdk_1sim(i) = sum(counter_sic_only_one) / trials^2;
    %If interfering>1
    Qdk_2sim(i) = sum(counter_dom) / trials^2;
   
    %Success probability (decoding the desirable signal) after SIC
    Qsic_sim(i) = sum(counter_sic) / trials^2;

end

%Save .dat files for latex diagrams
Qtotal = [transpose(doo), transpose(Qdi_sim+Qdk_1sim+Qdk_2sim.*Qsic_sim)];
save('~/Thesis/Diagrams/Qtotalsim12km.dat','Qtotal','-ascii');
