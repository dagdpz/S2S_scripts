function mat=cell2mat_expanded(cellarray,mode)
if nargin==1
   mode='beginning';
end
n_columns=max(cellfun(@(x) numel(x),cellarray));
n_rows=numel(cellarray);


if iscell(cellarray{1})
mat=NaN(n_rows,n_columns);
mat=num2cell(mat);
else
mat=NaN(n_rows,n_columns);
end

for row=1:n_rows
    n_cell_elements=numel(cellarray{row});
    if strcmp(mode,'end')
    mat(row,end-n_cell_elements+1:end)=cellarray{row}; 
    else
    mat(row,1:n_cell_elements)=cellarray{row}; 
    end
end

end