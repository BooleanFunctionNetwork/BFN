function [L_nolink L_link I_link L_link_] = test1_fun(PX,INDEX,tau)
eps = 0.005;
n = length(PX(:,1,1));
m = length(PX(1,:,1));
q0 = zeros(n,1);
q1 = zeros(n,1);
for i = 1:n
    q0(i) = sum(PX(i,:,1))/m;
    q1(i) = sum(PX(i,:,2))/m;
end
[p00, p11, p01, p10] = deal(zeros(size(INDEX,1),m-tau));
[q02, q12] = deal(zeros(size(INDEX,1), 1));

for i=1:length(INDEX(:,1))
     for j=1:m-tau
            
            p00(i,j)=PX(INDEX(i,1),j,1)*PX(INDEX(i,2),j+tau,1);
            p11(i,j)=PX(INDEX(i,1),j,2)*PX(INDEX(i,2),j+tau,2);
            p01(i,j)=PX(INDEX(i,1),j,1)*PX(INDEX(i,2),j+tau,2);
            p10(i,j)=PX(INDEX(i,1),j,2)*PX(INDEX(i,2),j+tau,1);  
     end
    q02(i)=q0(INDEX(i,2));
    q12(i)=q1(INDEX(i,2));
end

L_nolink=log(prod((p00+p01+p10+p11),2));

[Lfun_1, Lfun_2, Lfun_3, Lfun_4, Lfun_5, Lfun_6] = deal(zeros(size(INDEX,1),m-tau));
   
for i=1:length(INDEX(:,1))
    Lfun_1(i,:)=p00(i,:)*(1-eps)/q02(i)+p11(i,:)*(1-eps)/q12(i)+p01(i,:)*eps/q12(i)+p10(i,:)*eps/q02(i);
    Lfun_2(i,:)=p01(i,:)*(1-eps)/q12(i)+p10(i,:)*(1-eps)/q02(i)+p00(i,:)*eps/q02(i)+p11(i,:)*eps/q12(i);
    Lfun_3(i,:)=p10(i,:)+p11(i,:)+p00(i,:)*(1-eps)/q02(i) + p01(i,:)*eps/q12(i);
    Lfun_4(i,:)=p10(i,:)+p11(i,:)+p01(i,:)*(1-eps)/q12(i) + p00(i,:)*eps/q02(i);
    Lfun_5(i,:)=p00(i,:)+p01(i,:)+p10(i,:)*(1-eps)/q02(i) + p11(i,:)*eps/q12(i);
    Lfun_6(i,:)=p00(i,:)+p01(i,:)+p11(i,:)*(1-eps)/q12(i) + p10(i,:)*eps/q02(i);
end
Lfun1=sum(log(Lfun_1),2);
Lfun2=sum(log(Lfun_2),2);
Lfun3=sum(log(Lfun_3),2);
Lfun4=sum(log(Lfun_4),2);
Lfun5=sum(log(Lfun_5),2);
Lfun6=sum(log(Lfun_6),2); 
        
L_link_=cat(2, Lfun1, Lfun2, Lfun3, Lfun4, Lfun5, Lfun6);
[L_link I_link]=max(L_link_,[],2);
