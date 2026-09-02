function [first_20, last_20] = my_callback(f)
    % Declare persistent variables inside this function
    persistent iter_count xlog_first xlog_last

    if nargin == 0
        % Retrieval mode: return stored values
        first_20 = xlog_first;
        last_20 = xlog_last;
        return;
    end

    % Callback mode: store values
    stop = false; %#ok<NASGU> % IPOPT expects this, but we don't use it
    x_val = full(f.x);

    % Initialize on first call
    if isempty(iter_count)
        iter_count = 0;
        xlog_first = [];
        xlog_last = [];
    end

    iter_count = iter_count + 1;

    % Store first 20
    if iter_count <= 20
        xlog_first = [xlog_first; x_val'];
    end

    % Keep rolling log of last 20
    if size(xlog_last, 1) >= 20
        xlog_last(1,:) = []; % remove oldest
    end
    xlog_last = [xlog_last; x_val'];

    first_20 = [];
    last_20 = [];
end