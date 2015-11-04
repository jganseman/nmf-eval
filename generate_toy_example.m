




function [adj,X] = generate_toy_example(vNum,recordNum,lenInterval)
    % this function is to generate data for the toy example
    
    % define random seed
    randn('state',0);
    rand('twister',0);
    
    % number of vertices in the graph
    %vNum = 5000;
    % number of record
    %recordNum = 400;
    
    % generate adjacenty matrix
    adj = zeros(vNum);
    for i=1:(vNum-1)
        adj(i,i+1) = 1;
        adj(i+1,i) = 1;
    end
    adj(vNum/2,vNum/2+1) = 0;
    
    % generate running record
    X = zeros(recordNum,vNum);
    for i=1:recordNum
        for k=1:1
            if i<=recordNum/2
                len      = round(randsample(lenInterval,1)/2);
                posStart = max(randsample(1:vNum/2,1)-len,1);
                %posStart = max(round(normrnd(vNum/recordNum*i,300,1))-len,1);
                posEnd   = min(posStart + len, vNum/2);
            else
                len      = round(randsample(lenInterval,1)/2);
                posStart = max(randsample(vNum/2+1:vNum,1)-len,vNum/2);
                %posStart = max(round(normrnd(vNum/recordNum*i,300,1))-len,vNum/2);
                posEnd   = min(posStart + len,vNum);
            end
            X(i,posStart:posEnd) = X(i,posStart:posEnd)+ ones(1,posEnd-posStart+1);
        end
    end
    
    X(X<=0)=0;
    
end