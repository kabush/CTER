function data = parse_confound(str_mat)

data = zeros(size(str_mat,1),1);
for i=1:size(str_mat,1)
    data(i,1) = str2double(str_mat(i,:));
end

%replace leading nans with 0
data(find(isnan(data)))=0;