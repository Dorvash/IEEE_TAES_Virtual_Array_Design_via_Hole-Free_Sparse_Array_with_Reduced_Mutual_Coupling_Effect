function coeefMat = double_blocked(A,sz)

output_row_num = sz(1) + size(A,1) - 1;
output_col_num = sz(2) + size(A,2) - 1;

A_zero_padded = zeros(output_row_num, output_col_num);
start_row = output_row_num - size(A, 1) + 1;
start_col = 1;
A_zero_padded(start_row:output_row_num, start_col:size(A,2)) = A;

blk_idx = zeros(size(A_zero_padded,1));

T = cell(size(A_zero_padded,1),1);

for i = 1:size(A_zero_padded,1)

    r = [A_zero_padded(i,1) zeros(1,sz(2)-1)];

    T{i} = toeplitz(A_zero_padded(i,:),r);

    blk_idx(i:size(blk_idx,1), i) = size(A_zero_padded,1):-1:i;

end

blk_idx(:,sz(1)+1:end) = [];
[T_row_num, T_col_num] = size(T{1});

coeffMatCell = cell(size(blk_idx));

for i = 1:size(blk_idx, 1)
    for j = 1:size(blk_idx, 2)
        if blk_idx(i,j) == 0
            coeffMatCell{i,j} = zeros(T_row_num, T_col_num);
        else
            coeffMatCell{i,j} = T{blk_idx(i,j)};
        end
    end
end

coeefMat = cell2mat(coeffMatCell);

end