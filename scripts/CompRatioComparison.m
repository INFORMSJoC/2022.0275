%% Regret Minimization and Separation - 1 Item 2 Bidders Case 
% We work with a believed-to-be-true distribution to compute an expected revenue maximizing mechanism
% Believed-to-be-true distribution is systematically contamined with a worst-case-inspired distribution
% Performance of the computed mechanism is compared with the performance of our proposed mechanism under the contamined distribution
% Comparison performed against second price auction (SPA) and single sample mechanism (SSM) as well

clear;

%% Parameters
%To generate the right-hand side figures in Figures 2,3 in main text and Figures 1,2,3 in online appendix, 
%the following parameters need to be adjusted to the values described in the paper
%Current values generate right-hand side figure of Figure 2 in main text
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
contLevel=[0:0.05:1]; %contamination level -- informally this is the difference between the training and test distributions

%% Worst-case inspired distribution to be used for contamination -- generate pdf directly
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

%% Evaluation on contamined distribution -- for different contamination levels
optRev = zeros(length(contLevel),1); %optimal revenues under full distributional knowledge
nomRev = zeros(length(contLevel),1); %revenues of nominal optimal mechanism
robRev = zeros(length(contLevel),1); %revenues of minimax regret mechanism
SPARev = zeros(length(contLevel),1); %revenues of SPA without reserve price
SSMRev = zeros(length(contLevel),1); %revenues of SSM
SSMTrueRev = zeros(length(contLevel),1); %revenues of SSM with sample from true distribution
for contInd=1:length(contLevel)
    %% True Distribution - Contamination with worst-case inspired distribution
    contpdf = (1-contLevel(contInd))*nomProb + contLevel(contInd)*conpdf;

    %% Optimal Mechanism with Full Distributional Knowledge
    q1Opt = sdpvar(supN1, supN2,'full');
    m1Opt = sdpvar(supN1, supN2,'full');
    q2Opt = sdpvar(supN1, supN2,'full');
    m2Opt = sdpvar(supN1, supN2,'full');
    %Constraints in Vectorized Format
    constraintsOpt = [];
    %IR Constraints
    sup1Rep = repmat(transpose(support1),1,supN2);
    constraintsOpt = [constraintsOpt, q1Opt(:,:).*sup1Rep(:,:) - m1Opt(:,:) >= 0];
    sup2Rep = repmat(support2,supN1,1);
    constraintsOpt = [constraintsOpt, q2Opt(:,:).*sup2Rep(:,:) - m2Opt(:,:) >= 0];
    %Allocation Feasibility
    constraintsOpt = [constraintsOpt, q1Opt(:) + q2Opt(:) <= 1];
    constraintsOpt = [constraintsOpt, q1Opt(:) >= 0];
    constraintsOpt = [constraintsOpt, q2Opt(:) >= 0];
    %IC Constraints
    for w1Ind=1:supN1
        Qtemp = repmat(q1Opt(w1Ind,:),supN1,1);
        Mtemp = repmat(m1Opt(w1Ind,:),supN1,1);
        constraintsOpt = [constraintsOpt, q1Opt(:,:).*sup1Rep(:,:) - m1Opt(:,:) >= Qtemp(:,:).*sup1Rep(:,:) - Mtemp(:,:)];
    end
    for w2Ind=1:supN2
        Qtemp = repmat(q2Opt(:,w2Ind),1,supN2);
        Mtemp = repmat(m2Opt(:,w2Ind),1,supN2);
        constraintsOpt = [constraintsOpt, q2Opt(:,:).*sup2Rep(:,:) - m2Opt(:,:) >= Qtemp(:,:).*sup2Rep(:,:) - Mtemp(:,:)];
    end
    %Objective -- maximize expected revenue
    objectiveOpt = 0;
    for v1ind=1:supN1
        for v2ind=1:supN2
            objectiveOpt = objectiveOpt + contpdf(v1ind,v2ind)*(m1Opt(v1ind,v2ind) + m2Opt(v1ind,v2ind));
        end
    end
    %Solve the optimal mechanism design problem
    ops = sdpsettings('solver', 'gurobi', 'verbose', 1); %Solver GUROBI - one can replace this another LP solver here that works with YALMIP
    diagnosis = optimize(constraintsOpt, -objectiveOpt, ops); %maximize
    opt_objectiveOpt = value(objectiveOpt);
    q1Opt = value(q1Opt);
    q2Opt = value(q2Opt);
    m1Opt = value(m1Opt);
    m2Opt = value(m2Opt);

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

    %% Performance Metrics
    optRev(contInd) = sum(sum(contpdf.*(m1Opt + m2Opt)));
    nomRev(contInd) = sum(sum(contpdf.*(m1 + m2)));
    robRev(contInd) = sum(sum(contpdf.*(m1prop + m2prop)));
    SPARev(contInd) = sum(sum(contpdf.*(m1SPA + m2SPA)));
    SSMRev(contInd) = sum(sum(contpdf.*(m1SSM + m2SSM)));
    SSMTrueRev(contInd) = sum(sum(contpdf.*(m1SSMTrue + m2SSMTrue)));
end

%% Plot results
figure;
set(0,'DefaultTextFontname', 'CMU Serif')
set(0,'DefaultAxesFontName', 'CMU Serif')
set(gcf,'color','w');
plot(contLevel, robRev./optRev,'LineWidth',3)
hold on;
plot(contLevel, nomRev./optRev,':','LineWidth',3)
hold on;
plot(contLevel, SSMTrueRev./optRev,'-*','LineWidth',3)
hold on;
plot(contLevel, SPARev./optRev,'--','LineWidth',3)
hold on;
set(findall(gcf,'-property','FontSize'),'FontSize',18);
legend('Minimax Regret Mechanism','Nominal Mechanism', 'Single Sample Mechanism','Second Price Auction','Location','southwest');
%title("$\overline{v} =[$" + v1bar + "," + v2bar + "$]$" + ", " + "$\mu=[$" + mean1 + "," + mean2 + "$]$" + ", " + "$\Sigma = [$" + var1 + "," + corr + ";" + corr + "," + var2 + "$]$" + ", ",'FontSize',18,'interpreter','latex')
xlabel('Contamination Level ($\epsilon$)','FontSize',24,'interpreter','latex');
ylabel('Revenues/Optimal Revenues','FontSize',24,'interpreter','latex');
ylim([0 1]);

