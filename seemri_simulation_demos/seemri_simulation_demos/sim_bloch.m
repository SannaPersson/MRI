
dt = 50e-9;

B0 = 0.010;

Bv = [0;0;B0];

gamma = 2*pi*42.58e6;

tp = 10e-6;
B1 = pi/2/gamma/tp

B1env = @(t) B1.*(t >= 0).*(t <= tp);

Ms = zeros(3,11);
Ms(3,:) = 1;

w_rf = gamma*B0;

Gz = 1e-3;
Gs = Gz*((1:size(Ms,2))-(size(Ms,2)+1)/2);

bh = plot3([0 0 0 0 1 1 1 1]*(size(Ms,2)+1), ...
           [-1 -1 1 1 -1 -1 1 1], ...
           [-1 1 -1 1 -1 1 -1 1], '.');
hold on
mh = plot3([(1:size(Ms,2)); (1:size(Ms,2))+Ms(1,:)], ...
           [zeros(1, size(Ms,2)); Ms(2,:)], ...
           [zeros(1, size(Ms,2)); Ms(3,:)]);
t = 0;
th = title(num2str(t, 't = %g'));
axis equal
hold off
drawnow

for tind = 1:1000;  
  Beff = [ones(1,size(Ms,2))*B1env(t); ...
          zeros(1,size(Ms,2)); ...
          Gs+(B0 - w_rf/gamma)];
  dalpha = gamma*sqrt(sum(Beff.^2))*dt;
  absBeff = sqrt(sum(Beff.^2));
  absBeff(absBeff==0) = 1;
  Ms = rotvec(Ms, Beff./([1;1;1]*absBeff), -dalpha);
  
  arrayfun(@(h, Mx, My, Mz, x) set(h, 'XData', [x Mx+x], ...
                                      'YData', [0 My], ...
                                      'ZData', [0 Mz]), ...
           mh', Ms(1,:), Ms(2,:), Ms(3,:), 1:size(Ms,2));
  set(th, 'String', num2str(t, 't = %g'));
  drawnow
  
  t = t + dt;
end

