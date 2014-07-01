function result = testCall

parfor k=1:5

    for j =1:5
        a(j).xy = j;
    end
end


for j = 1:5
    result(j).xy = gplus(a(j).xy);
end

% result = a;