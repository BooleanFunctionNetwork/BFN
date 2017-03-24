%Set the range of time-delays from 0 to 5
td = 0:5;
PX0 = csvread('px0_whole_genome.csv');
PX1 = csvread('px1_whole_genome.csv');
PX(:,:,1)=PX0;
PX(:,:,2)=PX1;
Gene_names = importdata('gene_names_whole_genome.csv');
Gene_names = Gene_names(2:end);
Tf_no = csvread('number_of_regulators_whole_genome.csv');
n = size(PX,1);
m = size(PX,2);
%Form the table of all possible gene pairs, which will enter the test
INDEX = zeros(n*Tf_no,2);
for i = 1:Tf_no
INDEX((i-1)*n+1:i*n,1) = i;
end
ind = 1:n;
INDEX(:,2) = repmat(ind',Tf_no,1);
INDEX(INDEX(:,1)==INDEX(:,2),:)=[];
%For each time delay we identify loglikelihoods of model with links and
%model without link with corresponding best fitting boolean function for
%linked model
tic
[L_nolink, L_link, F_link] = deal(zeros(size(INDEX,1), length(td)));
for i=1:length(td)
[L_nolink(:,i) L_link(:,i) F_link(:,i)] = test1_fun(PX,INDEX,td(i));
end
L_link_minus_nolink = L_link-L_nolink;
toc
%Choose the maximum log likelihood and corresponding boolean function
%over all possible time-delays
[L_link_max Ind_link_max] = max(L_link_minus_nolink,[],2);
time_delay_link_max = Ind_link_max-1;
idx = sub2ind(size(F_link), (1:size(F_link,1))',Ind_link_max);
F_link_max = F_link(idx);
%Compute p-values for each gene pair
p_val = 1-chi2cdf(2*L_link_max,1);
Result_test1 = [INDEX, F_link_max, time_delay_link_max, L_link_max, p_val];
%Select links with p < 0.05
Result_test1_p05 = Result_test1(Result_test1(:,6) <0.05,:);
Regulator_name = Gene_names(Result_test1_p05(:,1));
Target_name = Gene_names(Result_test1_p05(:,2));
Result_table = [Regulator_name, Target_name, num2cell(Result_test1_p05)];
T = cell2table(Result_table,'VariableNames',{'Source','Target','Source_No','Target_No','Function','Time_delay','Lik_Ratio','p_value'});
writetable(T,'test1_result_whole_genome.csv');