clf; hold on;
hold on; plot(state(1).G(1,:),state(1).G(2,:),'.');
for i = 100:2000
  for j = 1:length(state(i).l_F)
    H=circle(state(i).F(:,j),state(i).r_F(j),30,'r');
  end;
end;

