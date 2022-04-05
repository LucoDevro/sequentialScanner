%% Worker algorithm for the sequential scanner. Executes the actual scanning loop by looping over all gene candidates for a provided model.
%
% INPUTS:
%   model           COBRA model structure
%   genesel         cell array containing the candidate genes with names as in model.genes
%   n               number of affected genes
%   prodRxn         char containing the reaction to be optimised with name as in model.rxns
%   MM              molar mass of the product produced by prodRxn, used to calculate the product-biomass yield
%   upRegFactor     gene upregulation factor for use as in Wang et al. (Biochem Eng J doi: 10.1016/j.bej.2017.03.017)
%   solWT           wild-type solution structure
%   muThreshold     growth rate threshold
%   optimCol        column number within res to which the optimisation variable is assigned (4 for flux, 5 for yield)
%
% OUTPUTS:
%   res             cell array containing the phenotype ranking for this number of deletions, extended with the flux distribution,
%                   the growth rate, the product flux and the product yield resp. The 6th column is for ranking purposes.
%   model_c         cell array containing the model structure of each phenotype encountered in this scan
%
% Author: Lucas De Vrieze (14 Apr 2021)
%
function[res,model_c]=sequentialScannerWorker(model,genesel,prodRxn,MM,upRegFactor,solWT,muThreshold,optimCol)

    prodRxn=find(strcmp(model.rxns,prodRxn));
    muWT=solWT.f;
    
    % Initialising the loop
    res=cell(length(genesel)*2,6);
    jump=length(genesel);
    model_c=cell(length(genesel)*2,1);
    fprintf(append('\tCurrent gene:\t\t',repmat(' ',1,max(cellfun('length',genesel)))));
    name='';
    
    % The loop itself, scanning all candidate genes
    for g=1:length(genesel)
        nu=genesel{g};
        
        % String handling for verbose feedback
        fprintf(repmat('\b',1,length(name)+1));
        name=model.geneNames{contains(model.genes,nu)};
        fprintf(append(name,'-'));
        t=GetGeneNameFromSysGeneName(model,{nu});
        
        % Assess a deletion and collect results.
        res{g,1}=append(t{1},'-');
        model_del=deleteModelGenes(model,nu);
        model_c{g}=model_del;
        sol=MOMA(model,model_del,'max',false,'one',true,solWT);
        if isempty(sol.x) % An 'infeasible' solver feedback is penalised.
            res{g,6}=-Inf;
        else
            res{g,2}=sol.x; % Flux distribution
            res{g,3}=sol.f; % Growth rate
            res{g,4}=sol.x(prodRxn); % Product flux
            res{g,5}=sol.x(prodRxn)/sol.f/1000*MM;
            if sol.f<muWT*muThreshold % A growth rate below the threshold is penalised.
                res{g,6}=0;
            else
                res{g,6}=res{g,optimCol};
            end
        end
        
        % Assess an upregulation and collect results. Remove GECKO constraint if encountered.
        res{g+jump,1}=append(t{1},'+');
        fprintf('\b+');
        model_up=upregulateModelGenes(model,nu,upRegFactor,solWT);
        [~,geneInfo]=findRxnsFromGenes(model,t,true);
        formsplit=strsplit(geneInfo{2});
        GECKO=find(contains(formsplit,'_ab'));
        if ~isempty(GECKO)
            model_up=removeMetabolites(model_up,formsplit(GECKO),true);
        end
        model_c{g+jump}=model_up;
        sol_up=MOMA(model,model_up,'max',false,'one',true,solWT);
        if isempty(sol_up.x)
            res{g+jump,6}=-Inf;
        else
            res{g+jump,2}=sol_up.x;
            res{g+jump,3}=sol_up.f;
            res{g+jump,4}=sol_up.x(prodRxn);
            res{g+jump,5}=sol_up.x(prodRxn)/sol_up.f/1000*MM;
            if sol_up.f<muWT*muThreshold
                res{g+jump,6}=0;
            else
                res{g+jump,6}=res{g+jump,optimCol};
            end
        end
    end
    
    % Rank results
    [res,ind]=sortrows(res,6,'descend');
    model_c=model_c(ind);
    fprintf('\n')
end