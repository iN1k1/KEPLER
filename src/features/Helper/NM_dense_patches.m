function [ patches, centers ] = NM_dense_patches( image, cellWidth, cellHeight, stepRight, stepBottom )

[rows, cols, ~] = size(image);
centers = [];
numPatches = 1;
for r=1:stepBottom:rows
    toRow = (r+cellHeight)-1;
    if toRow>rows
        toRow = rows;
    end
    for c=1:stepRight:cols
        toCol = (c+cellWidth)-1;
        if toCol>cols
            toCol = cols;
        end
        patches{numPatches} = image(r:toRow, c:toCol, :);
        numPatches = numPatches + 1;
        centers = [centers; mean(r:toRow), mean(c:toCol)];
        if toCol == cols
            break;
        end
    end
    if toRow == rows
        break;
    end
    
end

end

