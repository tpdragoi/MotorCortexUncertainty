function aligned = doAlign(data, trials, dtBins)
    [T,N,TT] = size(data);
    aligned  = nan(T,N,TT);
    for k = 1:numel(trials)
      tr  = trials(k);
      sh  = dtBins(k);
      X   = data(:,:,tr);      % T×N
      if sh >= 0
        aligned(1:T-sh,:,tr) = X(1+sh:end,:);
      else
        aligned(1-sh:T,:,tr) = X(1:end+sh,:);
      end
    end
end