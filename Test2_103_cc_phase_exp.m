%Select links to enter Test2: if for link x1->x2 there exist links x1->x3
%and x3->x2, which time delays are compatible, i.e. td13 + td32 <= td13
%then such link enters Tes2, otherwise it's recorded as Direct
test1_result_table = importdata('test1_result_103_cc_phase.csv');
test1_result = test1_result_table.data;
n = size(test1_result,1);
Gene_names = importdata('list_103_cell_cycle_genes_phase_annotation.xlsx');
Gene_names = Gene_names.Sheet1;
PX0 = csvread('px0_103_cc_phase.csv');
PX1 = csvread('px1_103_cc_phase.csv');
PX(:,:,1) = PX0;
PX(:,:,2) = PX1;
Result_test2 = zeros(size(test1_result,1),8);
Direct = [];

for i=1:n
    %Check the if the link has potential intermediates
    x_intermed = (1:n);
    x_intermed(x_intermed==test1_result(i,1)) = [];
    x_intermed(x_intermed==test1_result(i,2)) = [];
    x1 = repmat(test1_result(i,1),length(x_intermed),1);
    x2 = repmat(test1_result(i,2),length(x_intermed),1);
    INDEX = cat(2,x1,x_intermed',x2);
    [idx1, test1_link1_loc] = ismember(INDEX(:,1:2), test1_result(:,1:2),'rows');
    [idx2, test1_link2_loc] = ismember(INDEX(:,2:3), test1_result(:,1:2),'rows');
    INDEX_sub_links_exist = INDEX((idx1>0 & idx2>0),:);
    test1_link1_loc = test1_link1_loc(idx1>0 & idx2>0);
    test1_link2_loc = test1_link2_loc(idx1>0 & idx2>0);
    INDEX_sub_links_exist_and_td_conform = INDEX_sub_links_exist((test1_result(test1_link1_loc,4)+test1_result(test1_link2_loc,4)) <= test1_result(i,4),:);
    TD = test1_result(i,4);
    
    if isempty(INDEX_sub_links_exist_and_td_conform)==1
        Direct = [Direct; test1_result(i,:)];
    else
        %Test2 procedure: if there are several possible intermediate genes,
        %choose the one which maximizes likelihood ratio between linked and
        %non-linked models, further choose time-delay which maximizes this
        %ratio
        INDEX=INDEX_sub_links_exist_and_td_conform;
        R=zeros(length(INDEX(:,1)),TD+1);
        I=zeros(length(INDEX(:,1)),TD+1);
        p=zeros(length(INDEX(:,1)),TD+1);
        [R(:,1), I(:,1), p(:,1)] = test2_fun(PX,INDEX,TD,TD);
        for j=1:TD
            [R(:,j+1), I(:,j+1), p(:,j+1)] = test2_fun(PX,INDEX,TD,TD-j);
        end
        if length(R(:,1))>1
            [R_gmax, Ind_gmax] = max(R,[],1);
            idx = sub2ind(size(R), Ind_gmax, 1:size(R,2));
            F_gmax = I(idx);
            p_gmax = p(idx);
            [R_gmax_tmax, Ind_gmax_tmax] = max(R_gmax);
            Index_gmax_tmax = INDEX(Ind_gmax(Ind_gmax_tmax),:);  
        else
            R_gmax = R;
            F_gmax = I;
            p_gmax = p;
            [R_gmax_tmax, Ind_gmax_tmax] = max(R_gmax);
            Index_gmax_tmax = INDEX;
        end
        Result_test2(i,:) = [Index_gmax_tmax,F_gmax(Ind_gmax_tmax),R_gmax_tmax,TD,Ind_gmax_tmax-1,1-p_gmax(Ind_gmax_tmax)];
    end
end    
Result_test2(all(Result_test2==0,2),:) = [];
Result_test2_p05 = Result_test2(Result_test2(:,8)<0.05,:);
Test2 = [Result_test2_p05(:,1) Result_test2_p05(:,3)];
Test1 = [test1_result(:,1) test1_result(:,2)];
[idxt2, loct1] = ismember(Test2, Test1, 'rows');
locidx = loct1(loct1>0);
Test2_result = test1_result(locidx,:);
Result_test1_test2_positive = [Direct(Direct(:,3)==1,:); Test2_result(Test2_result(:,3)==1,:)];
Regulator_name = Gene_names(Result_test1_test2_positive(:,1),:);
Target_name = Gene_names(Result_test1_test2_positive(:,2),:);
Result_table = [Regulator_name, Target_name, num2cell(Result_test1_test2_positive)];
T = cell2table(Result_table,'VariableNames',{'Source','Source_Phase','Target','Target_Phase','Source_no','Target_no','Function','Time_delay','Lik_Ratio','p_value'});
writetable(T,'test2_result_103_cc_phase.csv');