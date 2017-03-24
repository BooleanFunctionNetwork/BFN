function [R, I, p_chi2] = test2_fun(PX, INDEX, tau, tau_prime)

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

for i = 1:length(INDEX(:,1))
     for j = 1:m-tau
            p00(i,j) = PX(INDEX(i,2),j+tau_prime,1)*PX(INDEX(i,3),j+tau,1);
            p11(i,j) = PX(INDEX(i,2),j+tau_prime,2)*PX(INDEX(i,3),j+tau,2);
            p01(i,j) = PX(INDEX(i,2),j+tau_prime,1)*PX(INDEX(i,3),j+tau,2);
            p10(i,j) = PX(INDEX(i,2),j+tau_prime,2)*PX(INDEX(i,3),j+tau,1);    
     end
    q02(i) = q0(INDEX(i,3));
    q12(i) = q1(INDEX(i,3));
end
[Lfun1, Lfun2, Lfun3, Lfun4, Lfun5, Lfun6] = deal(zeros(size(INDEX,1),m-tau));
for i = 1:length(INDEX(:,1))
    Lfun1(i,:) = p00(i,:)*(1-eps)/q02(i)+p11(i,:)*(1-eps)/q12(i)+p01(i,:)*eps/q12(i)+p10(i,:)*eps/q02(i);
    Lfun2(i,:) = p01(i,:)*(1-eps)/q12(i)+p10(i,:)*(1-eps)/q02(i)+p00(i,:)*eps/q02(i)+p11(i,:)*eps/q12(i);
    Lfun3(i,:) = p10(i,:)+p11(i,:)+p00(i,:)*(1-eps)/q02(i) + p01(i,:)*eps/q12(i);
    Lfun4(i,:) = p10(i,:)+p11(i,:)+p01(i,:)*(1-eps)/q12(i) + p00(i,:)*eps/q02(i);
    Lfun5(i,:) = p00(i,:)+p01(i,:)+p10(i,:)*(1-eps)/q02(i) + p11(i,:)*eps/q12(i);
    Lfun6(i,:) = p00(i,:)+p01(i,:)+p11(i,:)*(1-eps)/q12(i) + p10(i,:)*eps/q02(i);
end
Lfun_1 = sum(log(Lfun1),2);
Lfun_2 = sum(log(Lfun2),2);
Lfun_3 = sum(log(Lfun3),2);
Lfun_4 = sum(log(Lfun4),2);
Lfun_5 = sum(log(Lfun5),2);
Lfun_6 = sum(log(Lfun6),2);
   
L0 = [Lfun_1, Lfun_2, Lfun_3, Lfun_4, Lfun_5, Lfun_6];
[L0_max, I0] = max(L0,[],2);

[p000, p001, p010, p011, p100, p101, p110, p111] = deal(zeros(size(INDEX,1),m-tau));   

for i = 1:length(INDEX(:,1))
    for j = 1:m-tau
    
        p000(i,j) = PX(INDEX(i,1),j,1)*PX(INDEX(i,2),j+tau_prime,1)*PX(INDEX(i,3),j+tau,1);
        p001(i,j) = PX(INDEX(i,1),j,1)*PX(INDEX(i,2),j+tau_prime,1)*PX(INDEX(i,3),j+tau,2);
        p010(i,j) = PX(INDEX(i,1),j,1)*PX(INDEX(i,2),j+tau_prime,2)*PX(INDEX(i,3),j+tau,1);
        p011(i,j) = PX(INDEX(i,1),j,1)*PX(INDEX(i,2),j+tau_prime,2)*PX(INDEX(i,3),j+tau,2);
        p100(i,j) = PX(INDEX(i,1),j,2)*PX(INDEX(i,2),j+tau_prime,1)*PX(INDEX(i,3),j+tau,1);
        p101(i,j) = PX(INDEX(i,1),j,2)*PX(INDEX(i,2),j+tau_prime,1)*PX(INDEX(i,3),j+tau,2);
        p110(i,j) = PX(INDEX(i,1),j,2)*PX(INDEX(i,2),j+tau_prime,2)*PX(INDEX(i,3),j+tau,1);
        p111(i,j) = PX(INDEX(i,1),j,2)*PX(INDEX(i,2),j+tau_prime,2)*PX(INDEX(i,3),j+tau,2);
       
        
    end
end
[Lf1, Lf2, Lf3, Lf4, Lf5, Lf6, Lf7, Lf8, Lf9, Lf10, Lf11, Lf12,...
    Lf13, Lf14, Lf15, Lf16, Lf17, Lf18, Lf19, Lf20, Lf21, Lf22,...
    Lf23, Lf24, Lf25, Lf26, Lf27, Lf28, Lf29, Lf30, Lf31, Lf32,...
    Lf33, Lf34, Lf35, Lf36, Lf37, Lf38, Lf39, Lf40, Lf41, Lf42]...
    = deal(zeros(size(INDEX,1),m-tau));
for i = 1:size(INDEX,1)
        Lf1(i,:) = (p000(i,:)/q02(i)+p010(i,:)/q02(i)+p100(i,:)/q02(i)+p111(i,:)/q12(i))*(1-eps)+(p001(i,:)/q12(i)+p011(i,:)/q12(i)+p101(i,:)/q12(i)+p110(i,:)/q02(i))*eps;
        Lf2(i,:) = (p000(i,:)/q02(i)+p010(i,:)/q02(i)+p101(i,:)/q12(i)+p110(i,:)/q02(i))*(1-eps)+(p001(i,:)/q12(i)+p011(i,:)/q12(i)+p100(i,:)/q02(i)+p111(i,:)/q12(i))*eps;
        Lf3(i,:) = (p000(i,:)/q02(i)+p011(i,:)/q12(i)+p100(i,:)/q02(i)+p110(i,:)/q02(i))*(1-eps)+(p001(i,:)/q12(i)+p010(i,:)/q02(i)+p101(i,:)/q12(i)+p111(i,:)/q12(i))*eps;
        Lf4(i,:) = (p000(i,:)/q02(i)+p011(i,:)/q12(i)+p101(i,:)/q12(i)+p110(i,:)/q02(i))*(1-eps)+(p001(i,:)/q12(i)+p010(i,:)/q02(i)+p100(i,:)/q02(i)+p111(i,:)/q12(i))*eps;
        Lf5(i,:) = (p000(i,:)/q02(i)+p011(i,:)/q12(i)+p101(i,:)/q12(i)+p111(i,:)/q12(i))*(1-eps)+(p001(i,:)/q12(i)+p010(i,:)/q02(i)+p100(i,:)/q02(i)+p110(i,:)/q02(i))*eps;
        Lf6(i,:) = (p001(i,:)/q12(i)+p010(i,:)/q02(i)+p100(i,:)/q02(i)+p110(i,:)/q02(i))*(1-eps)+(p000(i,:)/q02(i)+p011(i,:)/q12(i)+p101(i,:)/q12(i)+p111(i,:)/q12(i))*eps;
        Lf7(i,:) = (p001(i,:)/q12(i)+p010(i,:)/q02(i)+p100(i,:)/q02(i)+p111(i,:)/q12(i))*(1-eps)+(p000(i,:)/q02(i)+p011(i,:)/q12(i)+p101(i,:)/q12(i)+p110(i,:)/q02(i))*eps;
        Lf8(i,:) = (p001(i,:)/q12(i)+p010(i,:)/q02(i)+p101(i,:)/q12(i)+p111(i,:)/q12(i))*(1-eps)+(p000(i,:)/q02(i)+p011(i,:)/q12(i)+p100(i,:)/q02(i)+p110(i,:)/q02(i))*eps;
        Lf9(i,:) = (p001(i,:)/q12(i)+p011(i,:)/q12(i)+p100(i,:)/q02(i)+p111(i,:)/q12(i))*(1-eps)+(p000(i,:)/q02(i)+p010(i,:)/q02(i)+p101(i,:)/q12(i)+p110(i,:)/q02(i))*eps;
        Lf10(i,:)= (p001(i,:)/q12(i)+p011(i,:)/q12(i)+p101(i,:)/q12(i)+p110(i,:)/q02(i))*(1-eps)+(p000(i,:)/q02(i)+p010(i,:)/q02(i)+p100(i,:)/q02(i)+p111(i,:)/q12(i))*eps;
        Lf11(i,:)= p000(i,:)+p001(i,:)+(p010(i,:)/q02(i)+p100(i,:)/q02(i)+p110(i,:)/q02(i))*(1-eps)+(p011(i,:)/q12(i)+p101(i,:)/q12(i)+p111(i,:)/q12(i))*eps;
        Lf12(i,:)= p000(i,:)+p001(i,:)+(p010(i,:)/q02(i)+p100(i,:)/q02(i)+p111(i,:)/q12(i))*(1-eps)+(p011(i,:)/q12(i)+p101(i,:)/q12(i)+p110(i,:)/q02(i))*eps;
        Lf13(i,:)= p000(i,:)+p001(i,:)+(p010(i,:)/q02(i)+p101(i,:)/q12(i)+p110(i,:)/q02(i))*(1-eps)+(p011(i,:)/q12(i)+p100(i,:)/q02(i)+p111(i,:)/q12(i))*eps;
        Lf14(i,:)= p000(i,:)+p001(i,:)+(p010(i,:)/q02(i)+p101(i,:)/q12(i)+p111(i,:)/q12(i))*(1-eps)+(p011(i,:)/q12(i)+p100(i,:)/q02(i)+p110(i,:)/q02(i))*eps;
        Lf15(i,:)= p000(i,:)+p001(i,:)+(p011(i,:)/q12(i)+p100(i,:)/q02(i)+p110(i,:)/q02(i))*(1-eps)+(p010(i,:)/q02(i)+p101(i,:)/q12(i)+p111(i,:)/q12(i))*eps;
        Lf16(i,:)= p000(i,:)+p001(i,:)+(p011(i,:)/q12(i)+p100(i,:)/q02(i)+p111(i,:)/q12(i))*(1-eps)+(p010(i,:)/q02(i)+p101(i,:)/q12(i)+p110(i,:)/q02(i))*eps;
        Lf17(i,:)= p000(i,:)+p001(i,:)+(p011(i,:)/q12(i)+p101(i,:)/q12(i)+p110(i,:)/q02(i))*(1-eps)+(p010(i,:)/q02(i)+p100(i,:)/q02(i)+p111(i,:)/q12(i))*eps;
        Lf18(i,:)= p000(i,:)+p001(i,:)+(p011(i,:)/q12(i)+p101(i,:)/q12(i)+p111(i,:)/q12(i))*(1-eps)+(p010(i,:)/q02(i)+p100(i,:)/q02(i)+p110(i,:)/q02(i))*eps;
        Lf19(i,:)= p010(i,:)+p011(i,:)+(p000(i,:)/q02(i)+p100(i,:)/q02(i)+p110(i,:)/q02(i))*(1-eps)+(p001(i,:)/q12(i)+p101(i,:)/q12(i)+p111(i,:)/q12(i))*eps;
        Lf20(i,:)= p010(i,:)+p011(i,:)+(p000(i,:)/q02(i)+p100(i,:)/q02(i)+p111(i,:)/q12(i))*(1-eps)+(p001(i,:)/q12(i)+p101(i,:)/q12(i)+p110(i,:)/q02(i))*eps;
        Lf21(i,:)= p010(i,:)+p011(i,:)+(p000(i,:)/q02(i)+p101(i,:)/q12(i)+p110(i,:)/q02(i))*(1-eps)+(p001(i,:)/q12(i)+p100(i,:)/q02(i)+p111(i,:)/q12(i))*eps;
        Lf22(i,:)= p010(i,:)+p011(i,:)+(p000(i,:)/q02(i)+p101(i,:)/q12(i)+p111(i,:)/q12(i))*(1-eps)+(p001(i,:)/q12(i)+p100(i,:)/q02(i)+p110(i,:)/q02(i))*eps;
        Lf23(i,:)= p010(i,:)+p011(i,:)+(p001(i,:)/q12(i)+p100(i,:)/q02(i)+p110(i,:)/q02(i))*(1-eps)+(p000(i,:)/q02(i)+p101(i,:)/q12(i)+p111(i,:))*eps;
        Lf24(i,:)= p010(i,:)+p011(i,:)+(p001(i,:)/q12(i)+p100(i,:)/q02(i)+p111(i,:)/q12(i))*(1-eps)+(p000(i,:)/q02(i)+p101(i,:)/q12(i)+p110(i,:)/q02(i))*eps;
        Lf25(i,:)= p010(i,:)+p011(i,:)+(p001(i,:)/q12(i)+p101(i,:)/q12(i)+p110(i,:)/q02(i))*(1-eps)+(p000(i,:)/q02(i)+p100(i,:)/q02(i)+p111(i,:))*eps;
        Lf26(i,:)= p010(i,:)+p011(i,:)+(p001(i,:)/q12(i)+p101(i,:)/q12(i)+p111(i,:)/q12(i))*(1-eps)+(p000(i,:)/q02(i)+p100(i,:)/q02(i)+p110(i,:)/q02(i))*eps;
        Lf27(i,:)= p100(i,:)+p101(i,:)+(p000(i,:)/q02(i)+p010(i,:)/q02(i)+p110(i,:)/q02(i))*(1-eps)+(p001(i,:)/q12(i)+p011(i,:)/q12(i)+p111(i,:)/q12(i))*eps;
        Lf28(i,:)= p100(i,:)+p101(i,:)+(p000(i,:)/q02(i)+p010(i,:)/q02(i)+p111(i,:)/q12(i))*(1-eps)+(p001(i,:)/q12(i)+p011(i,:)/q12(i)+p110(i,:)/q02(i))*eps;
        Lf29(i,:)= p100(i,:)+p101(i,:)+(p000(i,:)/q02(i)+p011(i,:)/q12(i)+p110(i,:)/q02(i))*(1-eps)+(p001(i,:)/q12(i)+p010(i,:)/q02(i)+p111(i,:)/q12(i))*eps;
        Lf30(i,:)= p100(i,:)+p101(i,:)+(p000(i,:)/q02(i)+p011(i,:)/q12(i)+p111(i,:)/q12(i))*(1-eps)+(p001(i,:)/q12(i)+p010(i,:)/q02(i)+p110(i,:)/q02(i))*eps;
        Lf31(i,:)= p100(i,:)+p101(i,:)+(p001(i,:)/q12(i)+p010(i,:)/q02(i)+p110(i,:)/q02(i))*(1-eps)+(p000(i,:)/q02(i)+p011(i,:)/q12(i)+p111(i,:)/q12(i))*eps;
        Lf32(i,:)= p100(i,:)+p101(i,:)+(p001(i,:)/q12(i)+p010(i,:)/q02(i)+p111(i,:)/q12(i))*(1-eps)+(p000(i,:)/q02(i)+p011(i,:)/q12(i)+p110(i,:)/q02(i))*eps;
        Lf33(i,:)= p100(i,:)+p101(i,:)+(p001(i,:)/q12(i)+p011(i,:)/q12(i)+p110(i,:)/q02(i))*(1-eps)+(p000(i,:)/q02(i)+p010(i,:)/q02(i)+p111(i,:)/q12(i))*eps;
        Lf34(i,:)= p100(i,:)+p101(i,:)+(p001(i,:)/q12(i)+p011(i,:)/q12(i)+p111(i,:)/q12(i))*(1-eps)+(p000(i,:)/q02(i)+p010(i,:)/q02(i)+p110(i,:)/q02(i))*eps;
        Lf35(i,:)= p110(i,:)+p111(i,:)+(p000(i,:)/q02(i)+p010(i,:)/q02(i)+p100(i,:)/q02(i))*(1-eps)+(p001(i,:)/q12(i)+p011(i,:)/q12(i)+p101(i,:)/q12(i))*eps;
        Lf36(i,:)= p110(i,:)+p111(i,:)+(p000(i,:)/q02(i)+p010(i,:)/q02(i)+p101(i,:)/q12(i))*(1-eps)+(p001(i,:)/q12(i)+p011(i,:)/q12(i)+p100(i,:)/q02(i))*eps;
        Lf37(i,:)= p110(i,:)+p111(i,:)+(p000(i,:)/q02(i)+p011(i,:)/q12(i)+p100(i,:)/q02(i))*(1-eps)+(p001(i,:)/q12(i)+p010(i,:)/q02(i)+p101(i,:)/q12(i))*eps;
        Lf38(i,:)= p110(i,:)+p111(i,:)+(p000(i,:)/q02(i)+p011(i,:)/q12(i)+p101(i,:)/q12(i))*(1-eps)+(p001(i,:)/q12(i)+p010(i,:)/q02(i)+p100(i,:)/q02(i))*eps;
        Lf39(i,:)= p110(i,:)+p111(i,:)+(p001(i,:)/q12(i)+p010(i,:)/q02(i)+p100(i,:)/q02(i))*(1-eps)+(p000(i,:)/q02(i)+p011(i,:)/q12(i)+p101(i,:)/q12(i))*eps;
        Lf40(i,:)= p110(i,:)+p111(i,:)+(p001(i,:)/q12(i)+p010(i,:)/q02(i)+p101(i,:)/q12(i))*(1-eps)+(p000(i,:)/q02(i)+p011(i,:)/q12(i)+p100(i,:)/q02(i))*eps;
        Lf41(i,:)= p110(i,:)+p111(i,:)+(p001(i,:)/q12(i)+p011(i,:)/q12(i)+p100(i,:)/q02(i))*(1-eps)+(p000(i,:)/q02(i)+p010(i,:)/q02(i)+p101(i,:)/q12(i))*eps;
        Lf42(i,:)= p110(i,:)+p111(i,:)+(p001(i,:)/q12(i)+p011(i,:)/q12(i)+p101(i,:)/q12(i))*(1-eps)+(p000(i,:)/q02(i)+p010(i,:)/q02(i)+p100(i,:)/q02(i))*eps;
end
Lf_1 = sum(log(Lf1),2);
Lf_2 = sum(log(Lf2),2);
Lf_3 = sum(log(Lf3),2);
Lf_4 = sum(log(Lf4),2);
Lf_5 = sum(log(Lf5),2);
Lf_6 = sum(log(Lf6),2);
Lf_7 = sum(log(Lf7),2);
Lf_8 = sum(log(Lf8),2);
Lf_9 = sum(log(Lf9),2);
Lf_10 = sum(log(Lf10),2);
Lf_11 = sum(log(Lf11),2);
Lf_12 = sum(log(Lf12),2);
Lf_13 = sum(log(Lf13),2);
Lf_14 = sum(log(Lf14),2);
Lf_15 = sum(log(Lf15),2);
Lf_16 = sum(log(Lf16),2);
Lf_17 = sum(log(Lf17),2);
Lf_18 = sum(log(Lf18),2);
Lf_19 = sum(log(Lf19),2);
Lf_20 = sum(log(Lf20),2);
Lf_21 = sum(log(Lf21),2);
Lf_22 = sum(log(Lf22),2);
Lf_23 = sum(log(Lf23),2);
Lf_24 = sum(log(Lf24),2);
Lf_25 = sum(log(Lf25),2);
Lf_26 = sum(log(Lf26),2);
Lf_27 = sum(log(Lf27),2);
Lf_28 = sum(log(Lf28),2);
Lf_29 = sum(log(Lf29),2);
Lf_30 = sum(log(Lf30),2);
Lf_31 = sum(log(Lf31),2);
Lf_32 = sum(log(Lf32),2);
Lf_33 = sum(log(Lf33),2);
Lf_34 = sum(log(Lf34),2);
Lf_35 = sum(log(Lf35),2);
Lf_36 = sum(log(Lf36),2);
Lf_37 = sum(log(Lf37),2);
Lf_38 = sum(log(Lf38),2);
Lf_39 = sum(log(Lf39),2);
Lf_40 = sum(log(Lf40),2);
Lf_41 = sum(log(Lf41),2);
Lf_42 = sum(log(Lf42),2);

L1=cat(2,Lf_1,Lf_2,Lf_3,Lf_4,Lf_5,Lf_6,Lf_7,Lf_8,Lf_9,Lf_10,Lf_11,Lf_12,Lf_13,Lf_14,Lf_15,Lf_16,Lf_17,Lf_18,Lf_19,Lf_20,Lf_21,Lf_22,Lf_23,Lf_24,Lf_25,Lf_26,Lf_27,Lf_28,Lf_29,Lf_30,Lf_31,Lf_32,Lf_33,Lf_34,Lf_35,Lf_36,Lf_37,Lf_38,Lf_39,Lf_40,Lf_41,Lf_42);
[L1_max I1]=max(L1,[],2);
        
R = L1_max-L0_max;
p_chi2=chi2cdf(2*R,1);
I=zeros(length(R),1);

for i =1:length(R)
    if R(i)>=0
        I(i)=I1(i);
    else
        I(i)=0;
    end
end