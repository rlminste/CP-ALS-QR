function err = normdiff(A,B)
% computes exact relative error between two ktensors A,B in a
% memory-efficient manner

d = length(B.U);
dims = size(A);
halves = floor(dims/2);

Asub = A;
Bsub = B;

sqerr = 0;

for i = 0:2^(d)-1
    for j = 1:d
        if bitget(i,j)
            Asub.U{j} = A.U{j}(1:halves(j),:);
            Bsub.U{j} = B.U{j}(1:halves(j),:);
        else
            Asub.U{j} = A.U{j}(halves(j)+1:end,:);
            Bsub.U{j} = B.U{j}(halves(j)+1:end,:);
        end
    end
    sqerr = sqerr + (norm(full(Asub)-full(Bsub)))^2;
end

err = sqrt(sqerr);