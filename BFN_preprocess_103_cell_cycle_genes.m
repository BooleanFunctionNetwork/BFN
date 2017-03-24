% List of known 103 cell-cycle genes according to Spellman et al.
Cell_cycle_103_list = importdata('list_103_cell_cycle_genes.xlsx');
Cell_cycle_103 = Cell_cycle_103_list.Sheet1;
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
%Select gene names of 103 cell cycle genes and their corresponding expression
%data from whole genome expression data syncronized with alpha factor arrest
[C, Cell_cycle_idx, Names_idx] = intersect(Cell_cycle_103, Gene_names);
Data_alpha_cc = Data_alpha(Names_idx,:);
Gene_names_cc = Gene_names(Names_idx);
%The list of all S.cerevisiae genes known to be Transcription Factors
%according to SGD (http://www.yeastgenome.org/) and having targets
Tf_sgd = importdata('tf_sgd.xlsx');
Tf_sgd_with_targets = Tf_sgd.tf_with_targets;
% Select genes which are TFs accroding to SGD
[C, Gene_names_cc_idx, Tf_sgd_with_targets_idx] = intersect(Gene_names_cc, Tf_sgd_with_targets(:,1));
Names_regulators = Gene_names_cc(Gene_names_cc_idx);
Data_regulators = Data_alpha_cc(Gene_names_cc_idx,:);
Names_targets = Gene_names_cc(setdiff(1:size(Gene_names_cc,1),Gene_names_cc_idx));
Data_targets = Data_alpha_cc(setdiff(1:size(Data_alpha_cc,1),Gene_names_cc_idx),:);
% Order gene names in vector Names and their corresponding expression in matrix
% X so that 185 TFs are followed by 6016 potential regulatory targets
Names = [Names_regulators; Names_targets];
X = [Data_regulators; Data_targets];
% Input NA values
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
Names(Four_missed_val_and_more,:) = [];

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
dlmwrite('px0_103_cell_cycle.csv', PX(:,:,1), 'precision', 16);
dlmwrite('px1_103_cell_cycle.csv', PX(:,:,2), 'precision', 16);
T = cell2table(Names);
writetable(T,'gene_names_103_cell_cycle.csv');
csvwrite('number_of_regulators_103_cell_cycle.csv',length(C));