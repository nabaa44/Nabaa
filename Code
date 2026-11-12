function Ainv = myinv(A)
    [n, m] = size(A);

    if n == m
        Ainv = gaussJordan(A);
    elseif n > m
        Ainv = gaussJordan(A' * A) * A';
    else
        Ainv = A' * gaussJordan(A * A');
    end
end

function B = gaussJordan(A)
    n = size(A, 1);
    M = [A, eye(n)];

    for i = 1:n
        [~, k] = max(abs(M(i:n, i)));
        k = k + i - 1;

        if abs(M(k, i)) < 1e-12
            error('Matrix is singular: no inverse exists.');
        end

        M([i k], :) = M([k i], :);

        M(i, :) = M(i, :) / M(i, i);

        for j = [1:i-1, i+1:n]
            M(j, :) = M(j, :) - M(j, i) * M(i, :);
        end
    end

    B = M(:, n+1:end);
end
