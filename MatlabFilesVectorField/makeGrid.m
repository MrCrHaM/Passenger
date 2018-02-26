function A = makeGrid(img)
    [row, col] = size(img);
    % Assume row based grid!
    A = zeros(row * col);
    for i = 0:1:row - 1
        for j = 0:1:col - 1
            if j ~= col - 1
                A(i * col + j + 1, i * col + (j + 1) + 1) = 1;
            end
            if j ~= 0
                A(i * col + j + 1, i * col + (j - 1) + 1) = 1;
            end
            if i ~= row - 1
                A((i + 1) * col + j + 1, i * col + j + 1) = 1;
            end
            if i ~= 0
                A((i - 1) * col + j + 1, i * col + j + 1) = 1;
            end
        end
    end
end