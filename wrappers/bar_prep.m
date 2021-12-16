function [sdata,ratios] = bar_prep(tals,tpinv,tqr,tsvd)
% preps iteration data to be plotted as a grouped/stacked bar chart
% inputs:
%   - tals, tpinv, tqr, tsvd: structs with 3,4,5-way timing results
% outputs:
%   - sdata: struct with padded data to be plotted as a stacked bar chart
%   - ratios: total slowdown ratios between CP-ALS-QR and CP-ALS for all
%             three dimensions


%% pad with zeros
tals.r3 = [tals.r3,zeros(4,2)];
tpinv.r3 = [tpinv.r3,zeros(4,2)];

tals.r4 = [tals.r4,zeros(4,2)];
tpinv.r4 = [tpinv.r4,zeros(4,2)];

tals.r5 = [tals.r5,zeros(3,2)];
tpinv.r5 = [tpinv.r5,zeros(3,2)];

%% 3way prep
nals3 = zeros(4,5);
npinv3 = zeros(4,5);
nqr3 = zeros(4,5);
nsvd3 = zeros(4,5);

% put zeros in place to correspond with CP-ALS-QR/SVD 
% sum last columns into other category
nals3(:,1:2) = tals.r3(:,1:2);
nals3(:,3) = tals.r3(:,6);
nals3(:,4) = tals.r3(:,7);
nals3(:,5) = sum(tals.r3(:,3:5),2);

npinv3(:,1:2) = tpinv.r3(:,1:2);
npinv3(:,3) = tpinv.r3(:,6);
npinv3(:,4) = tpinv.r3(:,7);
npinv3(:,5) = sum(tpinv.r3(:,3:5),2);

% sum last columns into other category
nqr3(:,1:4) = tqr.r3(:,1:4);
nqr3(:,5) = sum(tqr.r3(:,5:7),2);

nsvd3(:,1:4) = tsvd.r3(:,1:4);
nsvd3(:,5) = sum(tsvd.r3(:,5:7),2);

%% 4way prep
nals4 = zeros(4,5);
npinv4 = zeros(4,5);
nqr4 = zeros(4,5);
nsvd4 = zeros(4,5);

% put zeros in place to correspond with CP-ALS-QR/SVD 
% sum last columns into other category
nals4(:,1:2) = tals.r4(:,1:2);
nals4(:,3) = tals.r4(:,6);
nals4(:,4) = tals.r4(:,7);
nals4(:,5) = sum(tals.r4(:,3:5),2);

npinv4(:,1:2) = tpinv.r4(:,1:2);
npinv4(:,3) = tpinv.r4(:,6);
npinv4(:,4) = tpinv.r4(:,7);
npinv4(:,5) = sum(tpinv.r4(:,3:5),2);

% sum last columns into other category
nqr4(:,1:4) = tqr.r4(:,1:4);
nqr4(:,5) = sum(tqr.r4(:,5:7),2);

nsvd4(:,1:4) = tsvd.r4(:,1:4);
nsvd4(:,5) = sum(tsvd.r4(:,5:7),2);

%% 5way prep
nals5 = zeros(3,5);
npinv5 = zeros(3,5);
nqr5 = zeros(3,5);
nsvd5 = zeros(3,5);

% put zeros in place to correspond with CP-ALS-QR/SVD 
% sum last columns into other category
nals5(:,1:2) = tals.r5(:,1:2);
nals5(:,3) = tals.r5(:,6);
nals5(:,4) = tals.r5(:,7);
nals5(:,5) = sum(tals.r5(:,3:5),2);

npinv5(:,1:2) = tpinv.r5(:,1:2);
npinv5(:,3) = tpinv.r5(:,6);
npinv5(:,4) = tpinv.r5(:,7);
npinv5(:,5) = sum(tpinv.r5(:,3:5),2);

% sum last columns into other category
nqr5(:,1:4) = tqr.r5(:,1:4);
nqr5(:,5) = sum(tqr.r5(:,5:7),2);

nsvd5(:,1:4) = tsvd.r5(:,1:4);
nsvd5(:,5) = sum(tsvd.r5(:,5:7),2);


%% creating stacked data

sdata3 = [];
sdata4 = [];
sdata5 = [];
for i = 1:4
    sdata3 = [sdata3; nals3(i,:); npinv3(i,:); nqr3(i,:); nsvd3(i,:); zeros(1,5)];
    sdata4 = [sdata4; nals4(i,:); npinv4(i,:); nqr4(i,:); nsvd4(i,:); zeros(1,5)];
end

for i = 1:3
    sdata5 = [sdata5; nals5(i,:); npinv5(i,:); nqr5(i,:); nsvd5(i,:); zeros(1,5)];
end

sdata = struct;
sdata.N3 = sdata3;
sdata.N4 = sdata4;
sdata.N5 = sdata5;

%% computing totals for ratios
tot_als3 = sum(tals_r3,2);
tot_qr3 = sum(tqr_r3,2);

tot_als4 = sum(tals_r4,2);
tot_qr4 = sum(tqr_r4,2);

tot_als5 = sum(tals_r5,2);
tot_qr5 = sum(tqr_r5,2);

% ratios
ratios = struct;
ratios.N3 = tot_qr3./tot_als3;
ratios.N4 = tot_qr4./tot_als4;
ratios.N5 = tot_qr5./tot_als5;

end
