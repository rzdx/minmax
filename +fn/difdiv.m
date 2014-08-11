function [ O ] = difdiv(A)
O=zeros(length(A)-1,1);
for i=1:length(A)-1
    O(i)=mean(A{i+1}./(A{i}+realmin));
end
end