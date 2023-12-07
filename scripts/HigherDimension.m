%% Regret Minimization and Separation - 1 Item >2 (Higher Number of) iid (Nominal Distribution) Bidders Case 
% We work with a believed-to-be-true distribution 
% Believed-to-be-true distribution is systematically contamined with a worst-case-inspired distribution
% Performance of our proposed mechanism under the contamined distribution
% is compared to benchmark mechanisms second-price auction (SPA) and
% single-sample mechanism (SSM)

clear;

%% Parameters
%Current values generate right-hand side figure of Figure 4 in main text
%and can be adjusted to generate left-hand side figure (I <-- 10)
I = 100; 
vbar = 1; 
mean = (1/2)*vbar; 
var = (1/10)*vbar;
supportmarg = [0:0.001:vbar];
supN = length(supportmarg);

contLevel=[0:0.05:1]; %contamination level -- informally this is the difference between the training and test distributions

%% Worst-case distribution for one bidder when all others have value 0 (non-standardized, i.e., sums up to 1/I)
conpdf = zeros(supN,1);
for vind = 1:supN
    if supportmarg(vind) >= supportmarg(supN)/exp(1) && vind<supN
        conpdf(vind) = 1 + (1/exp(1))*(supportmarg(supN)/supportmarg(vind)^2);
    end
end
conpdf = (1/I)*(conpdf./sum(conpdf))*(1 - 1/exp(1));
conpdf(supN,1) = (1/I)*(1/exp(1));


%% Nominal marginal distribution for one bidder (non-standardized, i.e., sums up to 1/I)
nomProb = zeros(supN,1);
for vind = 1:supN
    nomProb(vind) = mvnpdf(supportmarg(vind),mean,var);
end
nomProb = nomProb./sum(nomProb);

%% Expected Revenue under Nominal Distribution for SPA, SSM, and our proposed mechanism
revNom = zeros(supN,1); %Expected Revenue under Nominal Distribution for Different Reserve Prices
robRevNom = 0;
tempsum = 0;
for rind = 1:supN
    %Expected payment of one bidder
    revOne = supportmarg(rind)*(sum(nomProb(rind:supN)))*(sum(nomProb(1:rind-1)))^(I-1);
    for tempind = rind+1:supN
        revOne = revOne + supportmarg(tempind)*(sum(nomProb(tempind:supN)))*((sum(nomProb(1:tempind)))^(I-1) - (sum(nomProb(1:tempind-1)))^(I-1));
    end
    %Multiply by number of bidders to compute total
    revNom(rind) = I*revOne;
    if supportmarg(rind)>=vbar/exp(1)
        tempsum = tempsum + ((1+log(supportmarg(rind)/vbar))-(1+log(supportmarg(rind-1)/vbar)));
        robRevNom = robRevNom + ((1+log(supportmarg(rind)/vbar))-(1+log(supportmarg(rind-1)/vbar)))*revNom(rind);
    end
end
SPARevNom = revNom(1);

%% Expected Revenue under Worst-Case Like Distribution for SPA, SSM, and our proposed mechanism
revWorst = zeros(supN,1); %Expected Revenue under Worst-Case Like Distribution for Different Reserve Prices
robRevWorst = 0;
for rind = 1:supN
    revWorst(rind) = supportmarg(rind)*(sum(conpdf(rind:supN)))*I;
    if supportmarg(rind)>=vbar/exp(1)
        robRevWorst = robRevWorst + ((1+log(supportmarg(rind)/vbar))-(1+log(supportmarg(rind-1)/vbar)))*revWorst(rind);
    end
end
SPARevWorst = revWorst(1);

%% Expected Revenue under Contamined Distribution for SPA, SSM, and our proposed mechanism
SSMTrueRev = zeros(length(contLevel),1);
for contInd=1:length(contLevel)
    for rind = 1:supN
        if rind==1
            SSMTrueRev(contInd) = SSMTrueRev(contInd) + ((1-contLevel(contInd))*nomProb(rind) + contLevel(contInd)*((I-1)/I))*((1-contLevel(contInd))*revNom(rind) + contLevel(contInd)*revWorst(rind));
        else
            SSMTrueRev(contInd) = SSMTrueRev(contInd) + ((1-contLevel(contInd))*nomProb(rind) + contLevel(contInd)*conpdf(rind))*((1-contLevel(contInd))*revNom(rind) + contLevel(contInd)*revWorst(rind));
        end
    end
    robRev(contInd) = contLevel(contInd)*robRevWorst + (1-contLevel(contInd))*robRevNom;
    SPARev(contInd) = contLevel(contInd)*SPARevWorst + (1-contLevel(contInd))*SPARevNom;
end

%% Plot results
figure;
set(0,'DefaultTextFontname', 'CMU Serif')
set(0,'DefaultAxesFontName', 'CMU Serif')
plot(contLevel, robRev,'LineWidth',3)
hold on;
plot(contLevel, SSMTrueRev,'-*','LineWidth',3)
hold on;
plot(contLevel, SPARev,'--','LineWidth',3)
hold on;
set(findall(gcf,'-property','FontSize'),'FontSize',18);
legend('Minimax Regret Mechanism','Single Sample Mechanism','Second Price Auction','Location','southwest');
%title("$I =$" + I + ", " +"$\overline{v} =$" + vbar + ", " + "$\mu=$" + mean + ", " + "$\sigma^2 = $" + var,'FontSize',18,'interpreter','latex')
xlabel('Contamination Level ($\epsilon$)','FontSize',24,'interpreter','latex');
ylabel('Revenues','FontSize',24,'interpreter','latex');
ylim([0 1]);




