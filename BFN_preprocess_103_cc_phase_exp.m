% List of known 103 cell-cycle genes along with their phase annotation
List_103_cc_phase = importdata('list_103_cell_cycle_genes_phase_annotation.xlsx');
List_103_cc_phase = List_103_cc_phase.Sheet1;
% S.cerevisiae whole genome expression data, obtained from 
% http://genome-www.stanford.edu/cellcycle/data/rawdata/combined.txt
Raw_data = importdata('combined.txt');
Data = Raw_data.data;
% Select yeast cultures synchronized with alpha factor arrest; that's 18
% columns from 7th to 24th
Data_alpha = Data(:,7:24);
m = size(Data_alpha,2);
Text = Raw_data.textdata;
%Remove header from table
Gene_names = Text(2:end,1);
%Remove empty rows
Gene_names(sum(isnan(Data_alpha),2)==m,:) = [];
Data_alpha(sum(isnan(Data_alpha),2)==m,:) = [];
%Select gene's expression profiles according to list of 103 genes 
[Cell_cycle_idx, Names_loc] = ismember(List_103_cc_phase(:,1), Gene_names);
X = Data_alpha(Names_loc,:);
%Input missing values
One_missed_val = find(sum(isnan(X),2)==1);
Two_missed_val = find(sum(isnan(X),2)==2);
Three_missed_val = find(sum(isnan(X),2)==3);
Four_missed_val_and_more = find(sum(isnan(X),2)>=4);
for i = 1:length(One_missed_val(:,1))
    xx = find(isnan(X(One_missed_val(i),:))==1);
    x = 1:m;
    x(xx) = [];
    y = X(One_missed_val(i),:);
    y(xx) = [];
    yy = interp1(x,y,xx,'spline');
    X(One_missed_val(i),xx) = interp1(x,y,xx,'spline');
end
for i = 1:length(Two_missed_val(:,1))
    xx = find(isnan(X(Two_missed_val(i),:))==1);
    x = 1:m;
    x(xx) = [];
    y = X(Two_missed_val(i),:);
    y(xx) = [];
    yy = interp1(x,y,xx,'spline');
    X(Two_missed_val(i),xx) = yy;
end
for i = 1:length(Three_missed_val(:,1))
    xx = find(isnan(X(Three_missed_val(i),:))==1);
    x = 1:m;
    x(xx) = [];
    y = X(Three_missed_val(i),:);
    y(xx) = [];
    yy = interp1(x,y,xx,'spline');
    X(Three_missed_val(i),xx) = yy;
end
X(Four_missed_val_and_more,:) = [];

% ecdf discretization to P matrix
n = size(X,1);
Ecdf = zeros(n,m);
[XS, IS] = sort(X,2);

for i = 1:length(XS(:,1))
        XSU = unique(XS(i,:));
        XSUC = hist(XS(i,:),XSU);
        k = XSUC(1);
        Ecdf(i,IS(i,1:k)) = 1/m;
        for j = 2:length(XSUC)
            if XSU(j) < median(XS(i,:))
                Ecdf(i,IS(i,k+1:XSUC(j)+k)) = (k+1)/m;
            else
                Ecdf(i,IS(i,k+1:XSUC(j)+k)) = (k+XSUC(j))/m;
            end
            k = XSUC(j)+k;
        end
end

PX = zeros(n,m,2);
for i = 1:n
    for j = 1:m
        PX(i,j,2) = Ecdf(i,j); %P(X=1)
        PX(i,j,1) = 1 - Ecdf(i,j);%P(X=0)
    end
end
dlmwrite('px0_103_cc_phase.csv', PX(:,:,1), 'precision', 16);
dlmwrite('px1_103_cc_phase.csv', PX(:,:,2), 'precision', 16);