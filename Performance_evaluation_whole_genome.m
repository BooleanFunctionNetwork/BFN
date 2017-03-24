%Computes number of Positives, True Positives, False positives, Negatives, True
%negatives and False Negatives that BFN produces while using as reference
%golden-standard GRN#1, GRN#2, GRN#3 from Ma et al.(2014)
N = importdata('gene_names_whole_genome.csv');
Names = N(2:end);
Gene_no = size (Names,1);
Tf_no = importdata('number_of_regulators_whole_genome.csv');
Result = importdata('S1_Table.xlsx');
D = Result.data.Sheet1;
Result_no_GRN1 = D(:,1:2);
Result_no_GRN2 = D(:,1:2);
Result_no_GRN3 = D(:,1:2);
Gold_stand = importdata('goldstandard_grns.xlsx');
Gold_stand = Gold_stand(2:end,:);
GRN1 = Gold_stand(:,1:2);
GRN2 = Gold_stand(:,3:4);
GRN3 = Gold_stand(:,5:6);
[P, TP, FP, N, TN, FN] = deal(zeros(3,1));
%Limit GRN1 network in accordance with list of genes we operate with
[grn1_idx_source, loc_grn1_source] = ismember(GRN1(:,1),Names);
[grn1_idx_target, loc_grn1_target] = ismember(GRN1(:,2),Names);
GRN1_no = [loc_grn1_source loc_grn1_target];
GRN1_no_limited = GRN1_no(grn1_idx_source>0 & grn1_idx_target>0,:);
%Limit BFN result on whole genome in accordance with list of genes in GRN1
GRN1_unique_source = unique(GRN1_no_limited(:,1));
c1 = setdiff(1:Tf_no,GRN1_unique_source);
GRN1_unique_target = unique(GRN1_no_limited(:,2));
c2 = setdiff(1:Gene_no,GRN1_unique_target);
for i = 1:length(c1)
    q = find(Result_no_GRN1(:,1)==c1(i));
    Result_no_GRN1(q,:) = [];
end
for i=1:length(c2)
    q = find(Result_no_GRN1(:,2)==c2(i));
    Result_no_GRN1(q,:) = [];
end
GRN1_tp_list=intersect(GRN1_no_limited,Result_no_GRN1,'rows');
P(1) = size(Result_no_GRN1,1);
TP(1) = size(GRN1_tp_list,1);
FP(1) = P(1) - TP(1);
N(1) = size(GRN1_unique_source,1)*size(GRN1_unique_target,1)-P(1);
FN(1) = size(GRN1_no_limited,1)-TP(1);
TN(1) = N(1)-FN(1);

%Limit GRN2 network in accordance with list of genes we operate with
[grn2_idx_source, loc_grn2_source] = ismember(GRN2(:,1),Names);
[grn2_idx_target, loc_grn2_target] = ismember(GRN2(:,2),Names);
GRN2_no = [loc_grn2_source loc_grn2_target];
GRN2_no_limited = GRN2_no(grn2_idx_source>0 & grn2_idx_target>0,:);
%Limit BFN result on whole genome in accordance with list of genes in GRN2
GRN2_unique_source = unique(GRN2_no_limited(:,1));
c1 = setdiff(1:Tf_no,GRN2_unique_source);
GRN2_unique_target = unique(GRN2_no_limited(:,2));
c2 = setdiff(1:Gene_no,GRN2_unique_target);
for i = 1:length(c1)
    q = find(Result_no_GRN2(:,1)==c1(i));
    Result_no_GRN2(q,:) = [];
end
for i=1:length(c2)
    q = find(Result_no_GRN2(:,2)==c2(i));
    Result_no_GRN2(q,:) = [];
end
GRN2_tp_list=intersect(GRN2_no_limited,Result_no_GRN2,'rows');
P(2) = size(Result_no_GRN2,1);
TP(2) = size(GRN2_tp_list,1);
FP(2) = P(2) - TP(2);
N(2) = size(GRN2_unique_source,1)*size(GRN2_unique_target,1)-P(2);
FN(2) = size(GRN2_no_limited,1)-TP(2);
TN(2) = N(2)-FN(2);

%Limit GRN3 network in accordance with list of genes we operate with
[grn3_idx_source, loc_grn3_source] = ismember(GRN3(:,1),Names);
[grn3_idx_target, loc_grn3_target] = ismember(GRN3(:,2),Names);
GRN3_no = [loc_grn3_source loc_grn3_target];
GRN3_no_limited = GRN3_no(grn3_idx_source>0 & grn3_idx_target>0,:);
%Limit BFN result on whole genome in accordance with list of genes in GRN3
GRN3_unique_source = unique(GRN3_no_limited(:,1));
c1 = setdiff(1:Tf_no,GRN3_unique_source);
GRN3_unique_target = unique(GRN3_no_limited(:,2));
c2 = setdiff(1:Gene_no,GRN3_unique_target);
for i = 1:length(c1)
    q = find(Result_no_GRN3(:,1)==c1(i));
    Result_no_GRN3(q,:) = [];
end
for i=1:length(c2)
    q = find(Result_no_GRN3(:,2)==c2(i));
    
    Result_no_GRN3(q,:) = [];
end
GRN3_tp_list=intersect(GRN3_no_limited,Result_no_GRN3,'rows');
P(3) = size(Result_no_GRN3,1);
TP(3) = size(GRN3_tp_list,1);
FP(3) = P(3) - TP(3);
N(3) = size(GRN3_unique_source,1)*size(GRN3_unique_target,1)-P(3);
FN(3) = size(GRN3_no_limited,1)-TP(3);
TN(3) = N(3)-FN(3);
GRN_names = {'GRN1'; 'GRN2'; 'GRN3'};
Performance_characteristics = [P TP FP N TN FN];
Result_table = [GRN_names, num2cell(Performance_characteristics)];
T = cell2table(Result_table,'VariableNames',{'Reference','P','TP','FP','N','TN','FN'});
writetable(T,'performance_whole_genome.csv');