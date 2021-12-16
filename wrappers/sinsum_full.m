function T = sinsum_full(d,n)
% creates ktensor for representation of sin of sums using trig identities for 
% sum of angles
% Inputs:
%   - d: number of angles (modes)
%   - n: dimension of each angle (mode)
% Outputs:
%   - T: d-way ktensor with size n in each mode, and with rank 2^(d-1)



% discretize 
x = linspace(0,2*pi,n)';

% initialize factor matrices
A = cell(1,d);

% intialize terms for expanded sum
Ptotal = cell(1,2^(d-1));

count = 1;
for k = 0:floor((d-1)/2)
    S = nchoosek(1:d,2*k+1);
    for i = 1:size(S,1)
        if ismember(1,S(i,1))
            P = (-1)^k*sin(x);
        else
            P = (-1)^k*cos(x);
        end
        for j = 2:d
            if ismember(j,S(i,:))
                P = [P,sin(x)];
            else
                P = [P,cos(x)];
            end
        end
        Ptotal{count} = P;
        count = count+1;
    end
end
        
for i = 1:d
    for j = 1:2^(d-1)
        A{i}(:,j) = Ptotal{j}(:,i);
    end
end

T = ktensor(A);

end
