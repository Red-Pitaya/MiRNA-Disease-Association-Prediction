clear
clc

load .\data\diseasesim.txt
load .\data\mirsim.txt

A = load ('.\data\Known_miRNA-disease_association_number.txt');
%disp(diseasesim)

row=max(A(:,1));
col=max(A(:,2));
data=zeros(row,col);
for i=1:length(A(:,1))
    data(A(i,1),A(i,2))=1;
    data(A(i,2),A(i,1))=1;
end
%disp(data)
%save data.mat data

dataT = data';

%关联矩阵
M0 = zeros(878);
for i = 1:878
    if(i<=495)
       for j = 1:495
            M0(i,j) = mirsim(i,j);
       end
       for j = 496:878
            M0(i,j) = data(i,j-495);
       end
    else
       for j = 1:495
            M0(i,j) = dataT(i-495,j);
       end
       for j = 496:878
            M0(i,j) = diseasesim(i-495,j-495);
       end
    end
end
%disp(M1)

[~,n] = size(M0);

%KATZ算法
A = zeros(n);
for i =1:n
    A(:,i) = M0(:,i)/sum(M0(:,i));
end
beta = 0.1;
katz1 = beta*A;
katz2 = beta^2*A^2;
katz3 = beta^3*A^3;
katz4 = beta^4*A^4;
skatz = katz1+katz2+katz3+katz4;

r1 = 281;
%数据归一化
MaxValue=max(skatz(:,r1));
MinValue=min(skatz(:,r1));
vlabel = (skatz-MinValue)/(MaxValue-MinValue);

datalabel = vlabel(1:495,1);
for i=1:495
    if(datalabel(i,1)>0.2)
        datalabel(i,1)=1;
    else
        datalabel(i,1)=0;
    end
end
%disp(datalabel) %得到的标签
label = data(:,r1); %原始标签

%ROC
plotroc(label',datalabel')
%AUC
pos_num = sum(label == 1);
neg_num = sum(label == 0);
m = length(label);
[~, index] = sort(datalabel);
label = label(index);
PX = zeros(m+1,1);
PY = zeros(m+1,1);
Auc = 0;
PX(1) = 1; PY(1) = 1;
for i = 2:m
    TP = sum(label(i:m)==1);
    FP = sum(label(i:m)==0);
    PX(i) = FP/neg_num;
    PY(i) = TP/pos_num;
    Auc = Auc + (PY(i)+PY(i-1))*(PX(i-1)-PX(i))/2;     % 梯形面积：（上底+下底）*高/2
end
PX(m+1) = 0;
PY(m+1) = 0;
Auc = Auc + PY(m)*PX(m)/2


