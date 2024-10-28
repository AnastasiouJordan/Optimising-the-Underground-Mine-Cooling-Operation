function S = V2S(V,fields,S)
n = length(fields);
for i = 1:n
    S.(fields{i}) = V(i,:);
end