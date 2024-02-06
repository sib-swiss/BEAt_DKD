% ConsensusOPLS Script example
% 
% Please cite: 
% J. Boccard, D.N. Rutledge.
% A consensus OPLS-DA strategy for multiblock Omics data fusion
% Analytica Chimica Acta, 769, 30-39 (2013).
%
%
%% Data preparation

% Unit Variance scaling of each data matrix
collection(1)=matrix2saisir(zscore(RNA_mat));
collection(2)=matrix2saisir(zscore(Prot_mat));

BlockNames={'Genes', 'Proteins'};

collection(1).i=RNA_varNames;
collection(2).i=Prot_varNames;

collection(1).v=ObsNames;
collection(2).v=ObsNames;


%% Compute a consensusOPLS model with leave-one-out cross-validation

% Definition of parameters
LVsPred=1; % Number of predictive component(s)
LVsOrtho=1; % Maximum number of orthogonal components
CVfolds=14; % Number of cross-validation folds

% This is the main function to compute the consensusOPLS model
RVConsensusOPLSModelCV=RVConsensusOPLS(collection,Y,LVsPred,LVsOrtho,CVfolds,'nfold','mc','mc',0.75,'da',0);

disp('RVConsensusOPLS Results ');
disp(['R2: ' num2str(RVConsensusOPLSModelCV.koplsModel.R2Yhat(RVConsensusOPLSModelCV.cv.OrthoLVsOptimalNum+1))]);
disp(['Q2: ' num2str(RVConsensusOPLSModelCV.cv.Q2Yhat(RVConsensusOPLSModelCV.cv.OrthoLVsOptimalNum+1))]);

%% Plot the results

% Consensus Score plot
figure; plot(RVConsensusOPLSModelCV.koplsModel.T(:,1),RVConsensusOPLSModelCV.koplsModel.To(:,1),'k.');
title('ConsensusOPLS Score plot');
axis([-0.5 0.5 -1 1]);
xlabel('Predictive component')
ylabel('Orthogonal component')
text(RVConsensusOPLSModelCV.koplsModel.T(:,1), RVConsensusOPLSModelCV.koplsModel.To(:,1), ObsNames(:,1), 'HorizontalAlignment','left','VerticalAlignment','bottom','FontSize',8)

% Block contributions to the predictive component
figure; bar(RVConsensusOPLSModelCV.koplsModel.lambda(:,1));
set(gca,'xticklabel',BlockNames)
ylabel('Weight');
title('Block contributions to the predictive component');

% Block contributions to the first orthogonal component
figure; bar(RVConsensusOPLSModelCV.koplsModel.lambda(:,2));
set(gca,'xticklabel',BlockNames)
ylabel('Weight');
title('Block contributions to the first orthogonal component');

% Block contributions predictive vs. orthogonal
figure; scatter(RVConsensusOPLSModelCV.koplsModel.lambda(:,1),RVConsensusOPLSModelCV.koplsModel.lambda(:,2));
axis([0 0.8 0 0.8]);
xlabel('Predictive component')
ylabel('Orthogonal component')
text(RVConsensusOPLSModelCV.koplsModel.lambda(:,1), RVConsensusOPLSModelCV.koplsModel.lambda(:,2), BlockNames, 'HorizontalAlignment','left','VerticalAlignment','bottom','FontSize',8)
title('Block contributions');

% Loading plots (one for each table)
figure;
title('ConsensusOPLS Loading plot');
xlabel('Predictive component');
ylabel('Orthogonal component');
hold all;
for i=1:length(collection)
    scatter(RVConsensusOPLSModelCV.koplsModel.loadings{i,1},RVConsensusOPLSModelCV.koplsModel.loadings{i,2},'.','DisplayName',['Block ' num2str(i)]);
end
legend(gca,'show')

RVConsensusOPLSModelCV=RVConsensusOPLS(collection,Y,LVsPred,LVsOrtho,CVfolds,'nfold','mc','mc',0.75,'da',0);

nbruns=999;  % C'est le nombre de permutations
PermRes=RVConsensusOPLSPermReg(collection,Y,nbruns,1,3);

figure; hist(PermRes.Q2val,100);

% VIP calculation
VIP=MBVIP(collection,Y,RVConsensusOPLSModelCV)

%% save results

% Build results matrices
scores = table(ObsNames, RVConsensusOPLSModelCV.koplsModel.T, RVConsensusOPLSModelCV.koplsModel.To, 'VariableNames', {'ObsNames' 'Pred' 'Ortho'});
loadings_VIP_RNA = table(RNA_varNames, RVConsensusOPLSModelCV.koplsModel.loadings{1,1}, VIP{1,1}, 'VariableNames', {'RNA' 'Pred' 'VIPpred'});
loadings_VIP_Prot = table(Prot_varNames, RVConsensusOPLSModelCV.koplsModel.loadings{2,1}, VIP{1,2}, 'VariableNames', {'Prot' 'Pred' 'VIPpred'});
Desc = [transpose(RVConsensusOPLSModelCV.koplsModel.R2X) transpose(RVConsensusOPLSModelCV.koplsModel.R2Yhat) RVConsensusOPLSModelCV.cv.DQ2Yhat(1:2,1) RVConsensusOPLSModelCV.cv.Q2Yhat(1:2,1)];

% write csv files
writetable(scores, 'scores.txt');
writetable(loadings_VIP_RNA, 'loadings_RNA.txt');
writetable(loadings_VIP_Prot, 'loadings_Prot.txt');

csvwrite('Desc.csv', Desc);
csvwrite('BlocksContributions.csv', RVConsensusOPLSModelCV.koplsModel.lambda)
csvwrite('PermQ2.csv', transpose(PermRes.Q2val));


