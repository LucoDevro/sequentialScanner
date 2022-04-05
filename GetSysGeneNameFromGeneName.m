function[sysGeneNames]=GetSysGeneNameFromGeneName(model,GeneNames)
sysGeneNames=cell(length(GeneNames),1);
for geneInd=1:length(GeneNames)
    index=strcmpi(model.geneNames,GeneNames(geneInd));
    if ~any(index)
        error('Gene name not found!')
    else
        sysGeneNames{geneInd}=model.genes{index};
    end
end
end