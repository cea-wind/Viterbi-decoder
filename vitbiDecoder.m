function [ decoded ] = vitbiDecoder( code,trel,tblen )
%The simplest implementation of viterbi decoder
%This function can be optimazed in two aspect.
%   Using C program
%   Adjust algorithm by the structure of trellis
%The BER is lower of this function than vitdec in most cases. The reasons may includes
%   We choose the maximum metric at last.
%   Traceback depth largest.
K = log2(trel.numInputSymbols);
N = log2(trel.numOutputSymbols);
codeReshape = reshape(code,N,[])';
savedStates= zeros(trel.numStates,size(codeReshape,1)+1);
decoded =  zeros(trel.numStates,size(codeReshape,1));
metric = zeros(1,trel.numStates);
% initial state is 0, some states can't reach if m < log2(trel.numStates)/K
for m = 1:min(log2(trel.numStates)/K,size(codeReshape,1))
    
    selected = reshape(trel.nextStates(1+savedStates(:,m),:)',1,[]);
    outputs =  reshape(trel.outputs(1+savedStates(:,m),:)',1,[]);
    inputs = repmat(0:2^K-1,1,trel.numStates)';
    outputs = dec2bin(outputs) - 48;
    metric_add = repmat(codeReshape(m,:),size(outputs,1),1).*(1-2*outputs);
    metric_add = sum(metric_add,2);
    selectedMetric = reshape(repmat(metric,2^K,1),1,[]) + metric_add';
    
    savedStates(:,1:m) = savedStates(floor((0:trel.numStates-1)/2^K)+1,1:m);
    decoded(:,1:m-1) = decoded(floor((0:trel.numStates-1)/2^K)+1,1:m-1);
    savedStates(:,m+1) = selected(1:trel.numStates);
    decoded(:,m) = inputs(1:trel.numStates);
    metric = selectedMetric(1:trel.numStates);
    
end

if(ceil(log2(trel.numStates)/K)==log2(trel.numStates)/K)
    [~,initSort] = sort(savedStates(:,m+1));
    savedStates = savedStates(initSort,:);
    decoded = decoded(initSort,:);
    metric = metric(initSort);
else
    savedStatesTemp = savedStates;
    metricFlag = zeros(size(metric));
    metricTemp = metricFlag;
    decodedTemp = decoded;
    for n = 1:K^(N*floor((log2(trel.numStates)/K)))
        nextState = trel.nextStates(1+savedStates(n,floor(log2(trel.numStates)/K)),:);
        outputs =  trel.outputs(1+savedStates(n,floor(log2(trel.numStates)/K)),:);
        outputs = dec2bin(outputs) - 48;
        metric_add = repmat(codeReshape(floor(log2(trel.numStates)/K)...
            ,:),size(outputs,1),1).*(1-2*outputs);
        for k = 1:length(nextState)
            if(metricFlag(nextState(k)+1)==0)
                savedStatesTemp(nextState(k)+1,:) = savedStates(n,:);
                metricTemp(nextState(k)+1) = metric(n) + metric_add(k);
                decodedTemp(nextState(k)+1,:) = decoded(n,:);
                decodedTemp(nextState(k)+1,floor(log2(trel.numStates)/K)+1) = k;
                metricFlag(nextState(k)+1) = 1;
            elseif(metric(n) + metric_add(k)> metricTemp(nextState(k)+1))
                savedStatesTemp(nextState(k)+1,:) = savedStates(n,:);
                metricTemp(nextState(k)+1) = metric(n) + metric_add(k);
                decodedTemp(nextState(k)+1,:) = decoded(n,:);
                decodedTemp(nextState(k)+1,floor(log2(trel.numStates)/K)+1) = k;
            end
        end
    end
    savedStates = savedStatesTemp;
    metric = metricTemp;
    decoded = decodedTemp;
    savedStates(:,floor(log2(trel.numStates)/K)+1) = 0:trel.numStates-1;
end

% next states , output and metric for every possible inputs
% the order of current states is 0,1,2,...
selected = reshape(trel.nextStates(1:trel.numStates,:)',1,[]);
currentStates = reshape(repmat(0:trel.numStates-1,2^K,1),1,[]);
[~,stateSort] = sort(selected);
currentStates = currentStates(stateSort);
outputs =  reshape(trel.outputs(1:trel.numStates,:)',1,[]);
inputs = repmat(0:2^K-1,1,trel.numStates)';
inputs = inputs(stateSort);
outputs = dec2bin(outputs) - 48;

% survival path, include every states.
for m = ceil(log2(trel.numStates)/K)+1: size(codeReshape,1)
    % In AWGN channel,LLR = C1(r*v) + C2.
    metric_add = repmat(codeReshape(m,:),size(outputs,1),1).*(1-2*outputs);
    metric_add = sum(metric_add,2);
    selectedMetric = reshape(repmat(metric,2^K,1),1,[]) + metric_add';
    % order selectedMetric by state, and choose the maximum metric
    selectedMetricsort = selectedMetric(stateSort);
    selectedMetricsort = reshape(selectedMetricsort,2^K,[]);
    [metric_update,survior] = max(selectedMetricsort);
    survior = survior + (0:trel.numStates-1)*2^K;
    
    savedStates(:,1:m) = savedStates(currentStates(survior)+1,1:m);
    decoded(:,1:m-1) = decoded(currentStates(survior)+1,1:m-1);
    savedStates(:,m+1) = 0:trel.numStates-1;
    decoded(:,m) = inputs(survior);
    metric = metric_update;
end
[~,maxMetricIndex] = max(metric);
decoded = [zeros(1,tblen) decoded(maxMetricIndex,1:end-tblen)];
end

