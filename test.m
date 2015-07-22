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

[~,initSort] = sort(savedStates(:,m+1));
savedStates = savedStates(initSort,:);
decoded = decoded(initSort,:);
metric = metric(initSort);

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
for m = log2(trel.numStates)/K+1: size(codeReshape,1)
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
decoded = decoded(1,:);