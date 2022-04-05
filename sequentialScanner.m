%% Sequential scanner for both gene deletions and upregulations, adaptable for number of affected genes, loop selection percentile, upregulation factor & growth rate threshold
%
% INPUTS:
%   model           COBRA model structure
%   genesel         cell array containing the candidate genes with names as in model.genes
%   n               number of affected genes
%   prodRxn         char containing the reaction to be optimised with name as in model.rxns
%   MM              molar mass of the product produced by prodRxn, used to calculate the product-biomass yield
%   upRegFactor     gene upregulation factor for use as in Wang et al. (Biochem Eng J doi: 10.1016/j.bej.2017.03.017)
%   selPerc         selected percentile of the phenotype ranking when increasing the number of deletions
%   muThreshold     growth rate threshold
%   optimVar        optimisation variable, choice between product 'yield' and 'flux'
%
% OUTPUTS:
%   resPerLevel     cell array containing the phenotype ranking for each number of deletions extended with the flux distribution,
%                   the growth rate, the product flux and the product yield resp. The 6th column is for ranking purposes.
%
% Author: Lucas De Vrieze (14 Apr 2021)
%
function[resPerLevel]=sequentialScanner(model,genesel,n,prodRxn,MM,upRegFactor,selPerc,muThreshold,optimVar)
    % Initialisation
    j=1;
    resPerLevel=cell(n,1);
    solWT=optimizeCbModel(model,'max','one');
    if strcmpi(optimVar,'yield')
        optimCol=5;
    elseif strcmpi(optimVar,'flux')
        optimCol=4;
    else
        error('Optimisation variable not understood!')
    end
    
    while j<=n
        fprintf(append('Entering loop ',num2str(j),'...\n'));
        
        % Full first-level scan for the first gene
        if j==1
            [res,model_c]=sequentialScannerWorker(model,genesel,prodRxn,MM,upRegFactor,solWT,muThreshold,optimCol);
            resPerLevel{j}=res;
            cutoff=ceil(size(model_c,1)*selPerc);
            topModels=model_c(1:cutoff);
        
        % New scan for selected percentile of previous loop
        else
            resThisLevel=cell(0,5);
            model_c_ThisLevel=cell(0,1);
            for m=1:length(topModels)
                
                % Determine genes already affected in previous loops and remove them from the candidate list
                geneIndices=contains(topModels{m}.genes,'_')&~contains(topModels{m}.genes,'E_');
                genesUsed=model.genes(geneIndices);
                [~,genesAlreadyUsed,~]=intersect(genesel,genesUsed);
                geneselUnused=genesel;
                geneselUnused(genesAlreadyUsed)=[];
                
                % Some string handling for a proper gene representation in the results cell array
                geneIndicesDels=contains(topModels{m}.genes,'_d');
                geneNamesDels=append(model.geneNames(geneIndicesDels),'-');
                geneIndicesUps=xor(geneIndices,geneIndicesDels);
                geneNamesUps=append(model.geneNames(geneIndicesUps),'+');
                geneNames=[geneNamesDels;geneNamesUps];
                geneNames=strjoin(geneNames,' & ');
                
                % Calling the scanner worker and storing its results
                fprintf(append('Starting scan for gene ',num2str(m),':\t',geneNames,'...\n'));
                [res,model_c]=sequentialScannerWorker(topModels{m},geneselUnused,prodRxn,MM,upRegFactor,solWT,muThreshold,optimCol);
                res(:,1)=cellfun(@(x) strjoin({x,geneNames},' & '),res(:,1),'UniformOutput',false);
                
                % Store all results of this loop together
                resThisLevel=[resThisLevel;res];
                model_c_ThisLevel=[model_c_ThisLevel;model_c];
            end
            
            % Deduplicate phenotype combinations (e.g. 'geneA+ & geneB-' and 'geneB- & geneA+')
            allCombinations=cellfun(@(x) strjoin(sort(strsplit(x,' & ')),''),resThisLevel(:,1),'UniformOutput',false);
            [~,uniqInd,~]=unique(allCombinations);
            resThisLevel=resThisLevel(uniqInd,:);
            model_c_ThisLevel=model_c_ThisLevel(uniqInd);
            
            % Rank the results of this loop
            [resThisLevel,ind]=sortrows(resThisLevel,6,'descend');
            model_c_ThisLevel=model_c_ThisLevel(ind);
            resPerLevel{j}=resThisLevel;
            
            % Determine selected percentile for the next loop
            cutoff=ceil(size(model_c_ThisLevel,1)*selPerc);
            topModels=model_c_ThisLevel(1:cutoff);
        end
        j=j+1;
    end
end