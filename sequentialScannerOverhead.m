%% Overhead script for determining knockout and overexpression candidates by brute force sequential scanning.

% Initialise model and simulation environment
close all
clearvars
fprintf('Initiating environment...\n');
modelName='model100.mat';
load(modelName)
model=eval(extractBefore(modelName,'.mat'));
changeCobraSolver('ibm_cplex');

% %Setup model
% excs={'EX_ala_L(e)','EX_arg_L(e)','EX_asn_L(e)','EX_asp_L(e)','EX_cys_L(e)','EX_gln_L(e)','EX_glu_L(e)','EX_gly(e)','EX_his_L(e)','EX_ile_L(e)','EX_leu_L(e)','EX_lys_L(e)','EX_met_L(e)','EX_phe_L(e)','EX_pro_L(e)','EX_ser_L(e)','EX_thr_L(e)','EX_trp_L(e)','EX_tyr_L(e)','EX_val_L(e)'};
% model=changeRxnBounds(model,'EX_glc(e)',0,'l');
% constrs=[-4.47,-0.86,-26.81,-0.12,-0.13,0,-2.47,0,-0.11,-0.32,-0.32,-0.04,-0.07,-0.11,-0.26,-1.15,-0.60,0,-0.11,-0.31]';
% model=changeRxnBounds(model,excs,constrs,'l');
% model=changeRxnBounds(model,{'EX_o2(e)'},-100,'l');
% model=changeRxnBounds(model,{'EX_nh4(e)'},-5,'l');
% model=changeRxnBounds(model,{'EX_pi(e)'},-5,'l');
% model=changeRxnBounds(model,{'EX_so4(e)'},-5,'l');

% Reduce gene candidate set to the carbohydrate and amino acid metabolism
selSubSys={'Amino acids and related molecules','Carbohydrates and related molecules'};
sel=union(findRxnsFromSubSystem(model,selSubSys{1}),findRxnsFromSubSystem(model,selSubSys{2}));
genesel={};
for r=1:length(sel)
    grxn=findGenesFromRxns(model,sel(r));
    for g=1:length(grxn{1})
        genesel=[genesel;grxn{1}(g)];
    end
end
genesel=unique(genesel);

% The actual scanning loop
prodRxn='EX_nh4(e)';
MM=18;
optimVar='Flux';
nGenes=2;
upregFactor=1.5;
selectedPercentile=0.1;
muThreshold=0.75;
tic
res=sequentialScanner(model,genesel,nGenes,prodRxn,MM,upregFactor,selectedPercentile,muThreshold,optimVar);
disp(toc)