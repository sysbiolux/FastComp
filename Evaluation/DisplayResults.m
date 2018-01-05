function [FCResult, FCPureResult,MOResult,RandomResult] = DisplayResults(Percentages,SampleCount,FileName,Comps,model,excludeC,onlypred)

if nargin < 6
    excludeC = 0;
end
[imps,exps] = findImpsAndExpRxns(model);
exchangereactions = model.rxns(imps|exps);
for prec=1:numel(Percentages)
    FCResults = [];
    FCPureResults = [];
    MOResults = []; 
    RandomResults = [];
    samples = 0;
    FCTimes = [];
    MOTimes = [];    
    FCPureTimes = [];
    
    for r = 1:SampleCount
            try
            load([FileName num2str(Percentages(prec)) '-' num2str(r) '.mat'])
            
            %fprintf(['loaded ' FileName num2str(Percentages(prec)) '-' num2str(r) '.mat\n'])
            for i = 1:length(fieldnames(Predictions))
                
                nonexchange = ~ismember(Predictions.(['Replicate' num2str(i)]){5},exchangereactions);
                origdata = Predictions.(['Replicate' num2str(i)]){4}';
                origdata = origdata(nonexchange);
                preddata1 = Predictions.(['Replicate' num2str(i)]){1};
                preddata2 = Predictions.(['Replicate' num2str(i)]){2};
                preddata4 = Predictions.(['Replicate' num2str(i)]){3};
                [fcEMR,fcL_Score,fcR_Score] = Evaluate(origdata,preddata1,Comps,excludeC,0,onlypred);
                [fcfullEMR,fcfullL_Score,fcfullR_Score] = Evaluate(origdata,preddata2,Comps,excludeC,0,onlypred);
                [MOEMR,MOL_Score,MOR_Score] = Evaluate(origdata,preddata4,Comps,excludeC,0,onlypred);
                [randEMR,randL_Score,randR_Score] = Evaluate(origdata,'',Comps,0,1);    
                FCResults(end+1,:) = [fcEMR,fcL_Score,fcR_Score];
                FCPureResults(end+1,:) = [fcfullEMR,fcfullL_Score,fcfullR_Score];
                MOResults(end+1,:) = [MOEMR,MOL_Score,MOR_Score];
                RandomResults(end+1,:) = [randEMR,randL_Score,randR_Score];
                samples = samples + 1;                
            end
            %fprintf('Success\n')
            FCTimes = [FCTimes ; CompartResults(:,5)];
            FCPureTimes = [FCPureTimes ; ResultFC(:,5)];
            MOTimes = [MOTimes ; ResultMO(:,5)];
            catch MException
                %disp(MException);
                continue
            end       
    end    
    fprintf('For %i percent the following statistics can be computed for a total number of %i:\n',Percentages(prec),samples)
    fprintf('Algorithm\tEMR   \tL_Score\tR_Score       \tTime     \n with Cytosolic Reactions\n')
    fprintf('Pure FC  \t %f5.2 ± %f5.2 \t %f5.2 ± %f5.2\t %f5.2 ± %f5.2\t %f5.2 ± %f5.2\n',mean(FCResults(:,1)),std(FCResults(:,1)),mean(FCResults(:,2)),std(FCResults(:,2)),mean(FCResults(:,3)),std(FCResults(:,3)),mean(FCTimes),std(FCTimes))
    fprintf('FCComp  \t %f5.2 ± %f5.2 \t %f5.2 ± %f5.2\t %f5.2 ± %f5.2\t %f5.2 ± %f5.2\n',mean(FCPureResults(:,1)),std(FCPureResults(:,1)),mean(FCPureResults(:,2)),std(FCPureResults(:,2)),mean(FCPureResults(:,3)),std(FCPureResults(:,3)),mean(FCPureTimes),std(FCPureTimes))
    fprintf('MO  \t %f5.2 ± %f5.2 \t %f5.2 ± %f5.2\t %f5.2 ± %f5.2\t %f5.2 ± %f5.2\n',mean(MOResults(:,1)),std(MOResults(:,1)),mean(MOResults(:,2)),std(MOResults(:,2)),mean(MOResults(:,3)),std(MOResults(:,3)),mean(MOTimes),std(MOTimes))    
    fprintf('Random  \t %f5.2 ± %f5.2 \t %f5.2 ± %f5.2\t %f5.2 ± %f5.2\n',mean(RandomResults(:,1)),std(RandomResults(:,1)),mean(RandomResults(:,2)),std(RandomResults(:,2)),mean(RandomResults(:,3)),std(RandomResults(:,3)))
    FCResult = [mean(FCResults(:,1)),std(FCResults(:,1)),...
                mean(FCResults(:,2)),std(FCResults(:,2)),...
                mean(FCResults(:,3)),std(FCResults(:,3)),...
                mean(FCTimes),std(FCTimes)];
    FCPureResult = [mean(FCPureResults(:,1)),std(FCPureResults(:,1)),...
                    mean(FCPureResults(:,2)),std(FCPureResults(:,2)),...
                    mean(FCPureResults(:,3)),std(FCPureResults(:,3)),...
                    mean(FCPureTimes),std(FCPureTimes)];
    MOResult = [mean(MOResults(:,1)),std(MOResults(:,1)),...
                mean(MOResults(:,2)),std(MOResults(:,2)),...
                mean(MOResults(:,3)),std(MOResults(:,3)),...
                mean(MOTimes),std(MOTimes)];
    RandomResult = [mean(RandomResults(:,1)),std(RandomResults(:,1)),mean(RandomResults(:,2)),std(RandomResults(:,2)),mean(RandomResults(:,3)),std(RandomResults(:,3)),0,0];

end


end

function randomsample = createRandomSampleSet(origdata,samplecount)
randomsample = createRandomSample(origdata);

for i = 2:samplecount
    randomsample = [randomsample ; createRandomSample(origdata)];
end

end

function randdata = createRandomSample(origdata)

randdata = zeros(size(origdata));
for i=1:size(origdata,2)
    pos = numel(find(origdata(:,i)));
    sample = randsample(1:size(origdata,1),pos);
    randdata(sample,i) = 1;
end
end

