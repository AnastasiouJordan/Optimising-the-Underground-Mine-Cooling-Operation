function S = V2S(V,fields)
n = length(fields);
for i = 1:n
    S.(fields{i}) = V(i,:);
end