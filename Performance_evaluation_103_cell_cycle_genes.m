N = importdata('gene_names_103_cell_cycle.csv');
Names = N(2:end);
Gene_no = size (Names,1);
Tf_no = importdata('number_of_regulators_103_cell_cycle.csv');
Result = importdata('test1_result_103_cell_cycle.csv');
Result = Result.data;
p_val = [1 0.5 0.4 0.3 0.2 0.1 0.095 0.09 0.085 0.08 0.075 0.07 0.065 ...
    0.06 0.055 0.05 0.045 0.04 0.035 0.03 0.025 0.02 0.015 0.01 0.0075 0.005];

SGD = importdata('sgd_regulatory_relations.xlsx');
SGD = SGD.Sheet1;
[P, TP, FP, N, TN, FN] = deal(zeros(26,1));
%Select all possible regulatory relations from SGD based on list of 103
%cell cycle genes
[sgd_idx_source, loc_sgd_source] = ismember(SGD(:,1),Names);
[sgd_idx_target, loc_sgd_target] = ismember(SGD(:,2),Names);
SGD_no = [loc_sgd_source loc_sgd_target];
SGD_no_limited = SGD_no(sgd_idx_source>0 & sgd_idx_target>0,:);
SGD_limited_wo_self_links = SGD_no_limited((SGD_no_limited(:,1)==SGD_no_limited(:,2))==0,:);
for i=1:26
    Result_no = Result(Result(:,6) <= p_val(i), 1:2);
    SGD_tp_list=intersect(SGD_limited_wo_self_links,Result_no,'rows');
    P(i) = size(Result_no,1);
    TP(i) = size(SGD_tp_list,1);
    FP(i) = P(i) - TP(i);
    N(i) = size(Result,1)-P(i);
    FN(i) = size(SGD_limited_wo_self_links,1)-TP(i);
    TN(i) = N(i)-FN(i);
end
Performance_characteristics = [P TP FP N TN FN p_val'];
Result_table = num2cell(Performance_characteristics);
T = cell2table(Result_table,'VariableNames',{'P','TP','FP','N','TN','FN','p_val'});
writetable(T,'performance_103_cell_cycle.csv');