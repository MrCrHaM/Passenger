function graph = fully_connected_e_neighbour_graph(img, sigma, e)
    [rowN, colN] = size(img);
    graph = zeros(rowN * colN);
    % Manioulating directly as sparse matrix might be better
    for row = 1:1:rowN
        for col = 1:1:colN
            for rowDiff = -e:1:e
                for colDiff = -e:1:e
                    % Avoid self edges
                    if rowDiff == 0 && colDiff == 0
                        continue
                    end
                    if(row + rowDiff <= rowN && 0 < row + rowDiff && col + colDiff <= colN && 0 < col + colDiff)
                        if(rowDiff * rowDiff + colDiff * colDiff <= e * e)
                            row2 = row + rowDiff;
                            col2 = col + colDiff;
                            weight = exp(- (img(row, col) - img(row2, col2)) ^ 2 / (2 * sigma ^ 2));
                            graph((row - 1) * colN + (col - 1) + 1, (row2 - 1) * colN + (col2 - 1) + 1) = weight;
                        end
                    end
                end
            end
        end
    end
    graph = sparse(graph);
end