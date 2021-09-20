function stop = fminoutfun(x, optimValues, state)
	global iter
	global fmin_history
	stop = false;
	switch state
  		case 'init'
			fmin_history = struct();
		    % fmin_history.x = x;
		    % fmin_history.optimValues = optimValues;
		    % fmin_history.timerVal = tic;
			fmin_history.x = zeros(100,length(x));
			fmin_history.gradient = zeros(100,length(x));
			fmin_history.fval = zeros(100,1);
			fmin_history.timerVal = zeros(100,1);
			fmin_history.timerStart = tic;
			fprintf('iter value at init state = %i',iter)
		case 'iter'
			% ind = length(fmin_history)+1;
		    % fmin_history(iter).x = x;
		    % fmin_history(iter).optimValues = optimValues;
		    % fmin_history(iter).timerVal = toc(fmin_history(1).timerVal);
			fmin_history.x(iter,:) = x;
			fmin_history.fval(iter) = optimValues.fval;
			fmin_history.gradient(iter,:) = optimValues.gradient;
			fmin_history.timerVal(iter) = toc(fmin_history.timerStart);
			fprintf('iter value at iter state = %i',iter)
		case 'done'
			%
  	end
	%history.x = [history.x; x];
	%history.fval = [history.fval; optimValues.fval];
	%history.gradient = [history.gradient; optimValues.gradient'];
	%save(par.fhistory,'history') %try making history an output of fminunc after exitflag, then only need to save once
end
