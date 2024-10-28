function V = S2V(S,fields)
n = length(fields);
for i = 1:n
     V(i,:) = S.(fields{i});
end