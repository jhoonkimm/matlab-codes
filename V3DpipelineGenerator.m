function V3DpipelineGenerator(path,folder)
% V3DpipelineGenerator(path,folder)
%
% Generates pipelines for all trial files in a single folder.
%
% No flags
%% Pipeline


    trialCell = {'pws','fast','slow','split','post'};
    
    for m = 1:length(trialCell)
        V3DGeneratePipeline(path,folder,'trial',trialCell{m})
    end
