clc;
den = 0.01;
n = 4000;
m = den * n ^ 2;
line = "%%MatrixMarket matrix coordinate real unsymmetric";
output = line + newline + num2str(n) + ' ' + num2str(n) + ' ' + num2str(m) + newline;

coor = zeros(m, 2);
val = zeros(m, 1);
U = randperm(n^2, m);

for k = 1:m
    i = floor(U(k) / n) + 1;
    j = mod(U(k), n) + 1;
    coor(k, 1) = i;
    coor(k, 2) = j;
    val(k, 1) = unifrnd(1, 100);
end

coor = sortrows(coor, 1);

for k = 1:m
    output = output + num2str(coor(k,1)) + ' ';
    output = output + num2str(coor(k,2)) + ' ';
    output = output + num2str(val(k,1)) + newline;
end


writelines(output, 'file.txt');




