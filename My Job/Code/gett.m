count = 1;
tall = zeros(218,2);
for j = 4:12
    tall(count,1) = 2002;
    tall(count,2) = j;
    count = count+1;
end
for i = 2003:2019
    for j = 1:12
        tall(count,1) = i;
        tall(count,2) = j;
        count = count+1;
    end
end
 for j = 1:5
    tall(count,1) = 2020;
    tall(count,2) = j;
    count = count+1;
 end