function[geneNames]=GetGeneNameFromSysGeneName(model,sysGeneNames)
geneNames=cell(length(sysGeneNames),1);
for sysGeneInd=1:length(sysGeneNames)
    index=strcmpi(model.genes,sysGeneNames(sysGeneInd));
    if ~any(index)
        error('Systematic gene name not found!')
    else
        geneNames{sysGeneInd}=model.geneNames{index};
    end
end
end