clf; hold on;
hold on; plot(G(1,:),G(2,:),'.');

for i = 450:500
%hold on; plot(rslt(i).F(1,:),rslt(i).F(2,:),'*g');
%hold on; plot(rslt(i).F1(1,:),rslt(i).F1(2,:),'.');
  for j = 1:(size(G,2))
    H=circle(G(:,j),rslt(j).R1,30,'r');
  end;
end;

