close all

a=2.^linspace(log2(0.25),log2(8),61);
s=4;
s=linspace(0.1,1,61);
z=linspace(-2,2,1001);

sumd=[];
for si=1:length(s)
yd=exp(-(abs(z)/s(si)).^a(ai))/sqrt(s(si));
sumd(si)=sqrt(sum(yd.*z.^2)/sum(yd));
end

plot(s,sumd)


