function loss = hammingLoss(origdata,preddata)

loss = sum(sum(origdata == preddata)) / (size(origdata,1) * size(origdata,2));