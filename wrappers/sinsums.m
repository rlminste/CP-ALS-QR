function T = sinsums(d,n)
% create rank d ktensor representation of sin of sums tensor with d modes
% Inputs: 
%    - d: number of modes of the tensor
%    - n: dimension of tensor
% Outputs:
%    - T: ktensor of d modes, size n in each mode, rank d


lambda = ones(d,1);

x = linspace(0,2*pi,n)';
% choose alphas
a = linspace(.1,3.1,d);
diff = a(2)-a(1);

offs = cell(d,d);

for i = 1:d
    for j = 1:d
        if i == j
            offs{i,j} = sin(x);
        else
            offs{i,j} = sin(x+(j-i)*diff)/sin((j-i)*diff);
        end
    end
end

% form factor matrices
A = cell(1,d);
for i = 1:d
    y = offs{1,i};
    for j = 2:d
        y = [y,offs{j,i}];
    end
    A{i} = y;
end

T = ktensor(lambda,A);

end
