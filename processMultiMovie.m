outPath='E:\fly movies';
v=VideoReader(movID);
nFrames=v.NumberOfFrames;
batchSize=20000;
nCycles=ceil(nFrames/batchSize);

for i = 1:nCycles
    if batchSize*i < nFrames
        outputStruct = runAlignment(movID,outPath,batchSize*(i-1)+1,batchSize*i);
    else
        outputStruct = runAlignment(movID,outPath,batchSize*(i-1)+1,nFrames);
    end
    clearvars -except i batchSize movID outPath nCycles outputStruct nFrames
end