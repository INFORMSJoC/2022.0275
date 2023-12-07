%% Regret Minimization and Separation - 1 Item 2 Bidders Case 
% Analysis of regret distribution particularly QUANTILES -- not expected revenue
% We work with a believed-to-be-true distribution to compute an expected revenue maximizing mechanism
% Believed-to-be-true distribution is systematically contamined with a worst-case-inspired distribution
% Performance of the computed mechanism is compared with the performance of our proposed mechanism under the contamined distribution
% Comparison performed against second price auction (SPA) and single sample mechanism (SSM) as well

clear;

%% Parameters
%To generate the left-hand side figures in Figures 2,3 in main text and Figures 1,2,3 in online appendix, 
%the following parameters need to be adjusted to the values described in the paper
%Current values generate left-hand side figure of Figure 2 in main text
%Adjust here----
v1bar = 1; 
v2bar = 1;
mean1 = (1/2)*v1bar; 
mean2 = (1/2)*v2bar;
mu = [mean1, mean2];
var1 = (1/10)*v1bar;
var2 = (1/10)*v2bar;
corr = 0;
%----until here

cov = [var1,corr;corr,var2]; 
support1 = [0:0.05:v1bar];
support2 = [0:0.05:v2bar];
supN1 = length(support1);
supN2 = length(support2);

tol=0.001; %used for numerical errors
contLevels=[0:0.05:1]; %contamination level -- informally this is the difference between the training and test distributions
percentile = 0.75;

%% Worst-case inspired distribution to be used for contamination
conpdf = zeros(supN1,supN2);
for v1ind = 1:supN1
    for v2ind = 1:supN2
        if v2ind==1
            if support1(v1ind) >= support1(supN1)/exp(1) && v1ind<supN1
                conpdf(v1ind,1) = 1 + (1/exp(1))*(support1(supN1)/support1(v1ind)^2);
            end
        end
        if v1ind==1
            if support2(v2ind) >= support2(supN2)/exp(1) && v2ind<supN2
                conpdf(1,v2ind) = 1 + (1/exp(1))*(support2(supN2)/support2(v2ind)^2);
            end
        end
    end
end
conpdf = (conpdf./sum(sum(conpdf)))*(1 - 1/exp(1));
conpdf(supN1,1) = (1/2)*(1/exp(1));
conpdf(1,supN2) = (1/2)*(1/exp(1));

%% Nominal distribution
nomProb = zeros(supN1,supN2);
for v1ind = 1:supN1
    for v2ind = 1:supN2
        nomProb(v1ind,v2ind) = mvnpdf([support1(v1ind),support2(v2ind)],mu,cov);
    end
end
nomProb = nomProb./sum(sum(nomProb));

% Marginal Nominal Probability of v1
margNomProb1 = sum(nomProb,2);

%% Reverse nominal distribution
revNomProb = ones(supN1,supN2) - nomProb;
revNomProb = revNomProb./sum(sum(revNomProb));

%% Compute optimal (revenue maximizing) mechanism under nominal distribution
q1 = sdpvar(supN1, supN2,'full');
m1 = sdpvar(supN1, supN2,'full');
q2 = sdpvar(supN1, supN2,'full');
m2 = sdpvar(supN1, supN2,'full');
%Constraints in Vectorized Format
constraints = [];
%IR Constraints
sup1Rep = repmat(transpose(support1),1,supN2);
constraints = [constraints, q1(:,:).*sup1Rep(:,:) - m1(:,:) >= 0];
sup2Rep = repmat(support2,supN1,1);
constraints = [constraints, q2(:,:).*sup2Rep(:,:) - m2(:,:) >= 0];
%Allocation Feasibility
constraints = [constraints, q1(:) + q2(:) <= 1];
constraints = [constraints, q1(:) >= 0];
constraints = [constraints, q2(:) >= 0];
%IC Constraints
for w1Ind=1:supN1
    Qtemp = repmat(q1(w1Ind,:),supN1,1);
    Mtemp = repmat(m1(w1Ind,:),supN1,1);
    constraints = [constraints, q1(:,:).*sup1Rep(:,:) - m1(:,:) >= Qtemp(:,:).*sup1Rep(:,:) - Mtemp(:,:)];
end
for w2Ind=1:supN2
    Qtemp = repmat(q2(:,w2Ind),1,supN2);
    Mtemp = repmat(m2(:,w2Ind),1,supN2);
    constraints = [constraints, q2(:,:).*sup2Rep(:,:) - m2(:,:) >= Qtemp(:,:).*sup2Rep(:,:) - Mtemp(:,:)];
end
%Objective -- maximize expected revenue
objective = 0;
for v1ind=1:supN1
    for v2ind=1:supN2
        objective = objective + nomProb(v1ind,v2ind)*(m1(v1ind,v2ind) + m2(v1ind,v2ind));
    end
end
%Solve the optimal mechanism design problem
ops = sdpsettings('solver', 'gurobi', 'verbose', 1);
diagnosis = optimize(constraints, -objective, ops); %maximize
opt_objective = value(objective);
q1 = value(q1);
q2 = value(q2);
m1 = value(m1);
m2 = value(m2);

%% Proposed minimax regret mechanism, second-price auction without reserve price, and single-sample mechanism
m1prop = zeros(supN1, supN2);
m2prop = zeros(supN1, supN2);
m1SPA = zeros(supN1, supN2);
m2SPA = zeros(supN1, supN2);
m1SSM = zeros(supN1, supN2);
m2SSM = zeros(supN1, supN2);
for v1Ind=1:supN1
    for v2Ind=1:supN2
        if support1(v1Ind)>=support2(v2Ind)
            m1SPA(v1Ind,v2Ind) = support2(v2Ind);
            if support1(v1Ind)>=max(v1bar,v2bar)/exp(1)
                if support2(v2Ind)>=max(v1bar,v2bar)/exp(1)
                    m1prop(v1Ind,v2Ind) = support1(v1Ind) + support2(v2Ind)*log(support2(v2Ind)/max(v1bar,v2bar));
                else
                    m1prop(v1Ind,v2Ind) = support1(v1Ind) - max(v1bar,v2bar)/exp(1);
                end
            end
        end
        if support2(v2Ind)>support1(v1Ind) %tie-breaking included here -- in case of equalities bidder 1 always wins
            m2SPA(v1Ind,v2Ind) = support1(v1Ind);
            if support2(v2Ind)>=v2bar/exp(1)
                if support1(v1Ind)>=max(v1bar,v2bar)/exp(1)
                    m2prop(v1Ind,v2Ind) = support2(v2Ind) + support1(v1Ind)*log(support1(v1Ind)/max(v1bar,v2bar));
                else
                    m2prop(v1Ind,v2Ind) = support2(v2Ind) - max(v1bar,v2bar)/exp(1);
                end
            end
        end
        for art3Ind=1:supN1 %Creating an artifical bidder whose value will be reserve price %Only reasonable in case of symmetric bidders, ow, can fail fundamentally
            if support1(v1Ind)>=support2(v2Ind) && support1(v1Ind)>=support1(art3Ind)
                m1SSM(v1Ind,v2Ind) = m1SSM(v1Ind,v2Ind) + margNomProb1(art3Ind)*max(support2(v2Ind),support1(art3Ind));
            elseif support2(v2Ind)>support1(v1Ind) && support2(v2Ind)>=support1(art3Ind)
                m2SSM(v1Ind,v2Ind) = m2SSM(v1Ind,v2Ind) + margNomProb1(art3Ind)*max(support1(v1Ind),support1(art3Ind));
            end
        end
    end
end

%% Computing regret for each mechanism and scenario of values
regNom = zeros(supN1, supN2);
regProp = zeros(supN1, supN2);
regSPA = zeros(supN1, supN2);
regSSM = zeros(supN1, supN2);
for v1Ind=1:supN1
    for v2Ind=1:supN2
        regNom(v1Ind,v2Ind) = max(support1(v1Ind),support2(v2Ind)) - m1(v1Ind,v2Ind) - m2(v1Ind,v2Ind);
        regProp(v1Ind,v2Ind) = max(support1(v1Ind),support2(v2Ind)) - m1prop(v1Ind,v2Ind) - m2prop(v1Ind,v2Ind);
        regSPA(v1Ind,v2Ind) = max(support1(v1Ind),support2(v2Ind)) - m1SPA(v1Ind,v2Ind) - m2SPA(v1Ind,v2Ind);
        regSSM(v1Ind,v2Ind) = max(support1(v1Ind),support2(v2Ind)) - m1SSM(v1Ind,v2Ind) - m2SSM(v1Ind,v2Ind);
    end
end

for contInd=1:length(contLevels)
    contLevel = contLevels(contInd);

    %% True Distribution - Contamination with worst-case inspired distribution
    contpdf = (1-contLevel)*nomProb + contLevel*conpdf;

    %% SSM - with sample from true distribution
    % Marginal Nominal Probability of v1
    margTrueProb1 = sum(contpdf,2);
    m1SSMTrue = zeros(supN1, supN2);
    m2SSMTrue = zeros(supN1, supN2);
    for v1Ind=1:supN1
        for v2Ind=1:supN2
            for art3Ind=1:supN1 %Creating an artifical bidder whose value will be reserve price %Only reasonable in case of symmetric bidders, ow, can fail fundamentally
                if support1(v1Ind)>=support2(v2Ind) && support1(v1Ind)>=support1(art3Ind)
                    m1SSMTrue(v1Ind,v2Ind) = m1SSMTrue(v1Ind,v2Ind) + margTrueProb1(art3Ind)*max(support2(v2Ind),support1(art3Ind));
                elseif support2(v2Ind)>support1(v1Ind) && support2(v2Ind)>=support1(art3Ind)
                    m2SSMTrue(v1Ind,v2Ind) = m2SSMTrue(v1Ind,v2Ind) + margTrueProb1(art3Ind)*max(support1(v1Ind),support1(art3Ind));
                end
            end
        end
    end
    regSSMTrue = zeros(supN1, supN2);
    for v1Ind=1:supN1
        for v2Ind=1:supN2
            regSSMTrue(v1Ind,v2Ind) = max(support1(v1Ind),support2(v2Ind)) - m1SSMTrue(v1Ind,v2Ind) - m2SSMTrue(v1Ind,v2Ind);
        end
    end

    %% Compute the distributions of regret values for each mechanism
    regNom = round(regNom,3);
    regArrNom = regNom(:);
    regArrNom = unique(regArrNom);
    for k=1:length(regArrNom)
        [r c]=find(regNom==regArrNom(k));
        probRegNom(k) = 0;
        for l=1:length(r)
            probRegNom(k) = probRegNom(k) + contpdf(r(l),c(l));
        end
    end
    regProp = round(regProp,3);
    regArrProp = regProp(:);
    regArrProp = unique(regArrProp);
    for k=1:length(regArrProp)
        [r c]=find(regProp==regArrProp(k));
        probRegProp(k) = 0;
        for l=1:length(r)
            probRegProp(k) = probRegProp(k) + contpdf(r(l),c(l));
        end
    end
    regSPA = round(regSPA,3);
    regArrSPA = regSPA(:);
    regArrSPA = unique(regArrSPA);
    for k=1:length(regArrSPA)
        [r c]=find(regSPA==regArrSPA(k));
        probRegSPA(k) = 0;
        for l=1:length(r)
            probRegSPA(k) = probRegSPA(k) + contpdf(r(l),c(l));
        end
    end
    regSSM = round(regSSM,3);
    regArrSSM = regSSM(:);
    regArrSSM = unique(regArrSSM);
    for k=1:length(regArrSSM)
        [r c]=find(regSSM==regArrSSM(k));
        probRegSSM(k) = 0;
        for l=1:length(r)
            probRegSSM(k) = probRegSSM(k) + contpdf(r(l),c(l));
        end
    end
    regSSMTrue = round(regSSMTrue,3);
    regArrSSMTrue = regSSMTrue(:);
    regArrSSMTrue = unique(regArrSSMTrue);
    for k=1:length(regArrSSMTrue)
        [r c]=find(regSSMTrue==regArrSSMTrue(k));
        probRegSSMTrue(k) = 0;
        for l=1:length(r)
            probRegSSMTrue(k) = probRegSSMTrue(k) + contpdf(r(l),c(l));
        end
    end

    %% Standardize the regret distributions to pre-selected regret values
    regVec = [0:0.001:1];
    cdfRegNom = zeros(length(regVec),1);
    cdfRegProp = zeros(length(regVec),1);
    cdfRegSPA = zeros(length(regVec),1);
    cdfRegSSM = zeros(length(regVec),1);
    cdfRegSSMTrue = zeros(length(regVec),1);
    for l=1:length(regVec)
        for k=1:length(regArrNom)
            if regArrNom(k)<=regVec(l)
                cdfRegNom(l) = cdfRegNom(l)+probRegNom(k);
            end
        end
        for k=1:length(regArrProp)
            if regArrProp(k)<=regVec(l)
                cdfRegProp(l) = cdfRegProp(l)+probRegProp(k);
            end
        end
        for k=1:length(regArrSPA)
            if regArrSPA(k)<=regVec(l)
                cdfRegSPA(l) = cdfRegSPA(l)+probRegSPA(k);
            end
        end
        for k=1:length(regArrSSM)
            if regArrSSM(k)<=regVec(l)
                cdfRegSSM(l) = cdfRegSSM(l)+probRegSSM(k);
            end
        end
        for k=1:length(regArrSSMTrue)
            if regArrSSMTrue(k)<=regVec(l)
                cdfRegSSMTrue(l) = cdfRegSSMTrue(l)+probRegSSMTrue(k);
            end
        end
    end
    %Recompute pdf
    pdfRegNomStand(1) = cdfRegNom(1);
    pdfRegPropStand(1) = cdfRegProp(1);
    pdfRegSPAStand(1) = cdfRegSPA(1);
    pdfRegSSMStand(1) = cdfRegSSM(1);
    pdfRegSSMTrueStand(1) = cdfRegSSMTrue(1);
    for l=2:length(regVec)
        pdfRegNomStand(l) = cdfRegNom(l) - cdfRegNom(l-1);
        pdfRegPropStand(l) = cdfRegProp(l) - cdfRegProp(l-1);
        pdfRegSPAStand(l) = cdfRegSPA(l) - cdfRegSPA(l-1);
        pdfRegSSMStand(l) = cdfRegSSM(l) - cdfRegSSM(l-1);
        pdfRegSSMTrueStand(l) = cdfRegSSMTrue(l) - cdfRegSSMTrue(l-1);
    end

    %% Compute regret quantiles
    sumpdf = 0;
    count = 1;
    while sumpdf<percentile
        sumpdf = sumpdf + pdfRegNomStand(count);
        quantileInd = count;
        count = count + 1;
    end
    quantileNom(contInd) = regVec(quantileInd);

    sumpdf = 0;
    count = 1;
    while sumpdf<percentile
        sumpdf = sumpdf + pdfRegPropStand(count);
        quantileInd = count;
        count = count + 1;
    end
    quantileProp(contInd) = regVec(quantileInd);

    sumpdf = 0;
    count = 1;
    while sumpdf<percentile
        sumpdf = sumpdf + pdfRegSPAStand(count);
        quantileInd = count;
        count = count + 1;
    end
    quantileSPA(contInd) = regVec(quantileInd);

    sumpdf = 0;
    count = 1;
    while sumpdf<percentile
        sumpdf = sumpdf + pdfRegSSMStand(count);
        quantileInd = count;
        count = count + 1;
    end
    quantileSSM(contInd) = regVec(quantileInd);

    sumpdf = 0;
    count = 1;
    while sumpdf<percentile
        sumpdf = sumpdf + pdfRegSSMTrueStand(count);
        quantileInd = count;
        count = count + 1;
    end
    quantileSSMTrue(contInd) = regVec(quantileInd);
end

%% Plot quantile regrets
figure;
set(0,'DefaultTextFontname', 'CMU Serif')
set(0,'DefaultAxesFontName', 'CMU Serif')
plot(contLevels, quantileProp,'LineWidth',3)
hold on;
plot(contLevels, quantileNom,':','LineWidth',3)
hold on;
plot(contLevels, quantileSSMTrue,'-o','LineWidth',3)
hold on;
plot(contLevels, quantileSPA,'--','LineWidth',3)
set(findall(gcf,'-property','FontSize'),'FontSize',18);
legend('Minimax Regret Mechanism','Nominal Mechanism','Single Sample Mechanism','Second Price Auction','Location','southeast');
%title("$\overline{v} =[$" + v1bar + "," + v2bar + "$]$" + ", " + "$\mu=[$" + mean1 + "," + mean2 + "$]$" + ", " + "$\Sigma = [$" + var1 + "," + corr + ";" + corr + "," + var2 + "$]$" + ", ",'FontSize',18,'interpreter','latex')
xlabel('Contamination Level ($\epsilon$)','FontSize',24,'interpreter','latex');
ylabel(percentile*100 + "th Percentile Regret",'FontSize',24,'interpreter','latex');
ylim([0 1]);

%% Worst-case regrets
worstregNom = max(max(regNom));
worstregProp = max(max(regProp));
worstregSPA = max(max(regSPA));
worstregSSM = max(max(regSSM));




