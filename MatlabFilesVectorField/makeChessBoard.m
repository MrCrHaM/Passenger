function img = makeChessBoard(WhiteFrameWidth, chessNodeWidth, chessElementNumberPerRow)
    size = WhiteFrameWidth * 2 + chessElementNumberPerRow * chessNodeWidth;
    img = zeros(size);
    for i = 0:1:size - 2 * WhiteFrameWidth - 1
        for j = 0:1:size - 2 * WhiteFrameWidth - 1
            if mod((fix(i / chessNodeWidth)) + (fix(j / chessNodeWidth)), 2) == 0
                img(WhiteFrameWidth + i + 1, WhiteFrameWidth + j + 1) = 200;
            end
        end
    end
end