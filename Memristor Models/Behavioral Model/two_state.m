clear;
clc;

f = 10;
tstop = 10/f;
dt = tstop*1e-6;
time = 0:dt:tstop-dt;
steps = length(time);
Vamp = 2;

t0 = 3e5;
t1 = 3e5;
v0 = 0.05;
v1 = 0.05;

states = 2;
init_state = 2;

P = zeros(states,length(time));
the_state = zeros(1,length(time));
P(init_state,1) = 1;
the_state(1) = init_state;

Ron = 1e3;
Roff = 10e3;

I{1} = @(V) V/Ron;
I{2} = @(V) V/Roff;

r_set = @(V) (V>0).*exp(abs(V)/v0)/t0*dt;
r_reset = @(V) (V<0).*exp(abs(V)/v1)/t1*dt;

consv_rate = 1.5e6*dt;
W = @(V,state) [-(state==2)*consv_rate-r_reset(V) (state==1)*consv_rate+r_set(V);...
    (state==2)*consv_rate+r_reset(V) -(state==1)*consv_rate-r_set(V)];

the_bar = waitbar(0,'Please wait...');
for i = 2:length(time)
    V(i-1) = Vamp*sin(2*pi*f*time(i-1));
    yy1 = W(V(i-1),the_state(i-1))*P(:,i-1);
    yy2 = yy1 + dt/2*W(V(i-1),the_state(i-1))*yy1;
    yy3 = yy2 + dt/2*W(V(i-1),the_state(i-1))*yy2;
    yy4 = yy3 + dt*W(V(i-1),the_state(i-1))*yy3;
    P(:,i) = P(:,i-1) + dt/6*(yy1 + 2*yy2 + 2*yy3 + yy4);
    the_rand = rand;
    new_state = 1*(the_rand<P(1,i)) + 2*(the_rand>=P(1,i));
    if new_state~=the_state(i-1)
        the_state(i) = new_state;
        P(:,i) = zeros(states,1);
        P(new_state,i) = 1;
    else
        the_state(i) = new_state;
    end
    if mod(i,steps/10)==0
        waitbar(i/steps,the_bar,['Simulation Completed: ' num2str(i*100/steps) '%']);
    end
end
close(the_bar);
V(length(time)) = 0;

figure(1);
fig = gcf;
fig.Color = 'white';
fig.Position = [fig.Position(1) fig.Position(2) 516 400];

subplot(3,1,1);
plot(time,V,'--');
ax = gca;
ax.Children(1).LineWidth = 1.5;
ax.XTick = [];
ax.YTick = [-1.5 0 1.5];
ax.YLim = [-1.55 1.55];
ax.FontSize = 10;
ylabel('Applied Voltage (V)','FontSize',14);

subplot(3,1,2);
P(P<1e-12) = 1e-12;
plot(time,P','LineWidth',2);
ax = gca;
ax.YScale = 'log';
ax.YGrid = 'on';
ax.YMinorGrid = 'off';
ax.YTick = 10.^[-12 -9 -6 -3 0];
ax.XTick = [];
ax.YLim = [1e-12 2];
ax.FontSize = 10;
ylabel('Probability','FontSize',14);
legend({'P_{ON}','P_{OFF}'},'FontSize',14);

subplot(3,1,3);
plot(time,the_state,'o-');
ax = gca;
ax.FontSize = 10;
ax.YTick = 1:2;
ax.YTickLabel = {'P_{ON}','P_{OFF}'};
ax.Children(1).LineWidth = 1;
ylim([0.5 2.5]);
ylabel('State','FontSize',14);
xlabel('Time (s)','FontSize',14);

III_prop = zeros(1,length(V));
for iii = 1:length(V)
    III_prop(iii) = I{the_state(iii)}(V(iii));
end

mean_III_prop = mean(reshape(III_prop,round(1/dt/f),round(length(III_prop)/(1/dt/f))),2);

figure(2);
fig = gcf;
fig.Color = 'white';
fig.Position = [fig.Position(1) fig.Position(2) 516 fig.Position(4)];
plot(V,III_prop,'o-',V(1:round(1/dt/f)),mean_III_prop)
ax = gca;
ax.FontSize = 18;
ax.Children(1).LineWidth = 2;
xlabel('Applied Voltage (V)','FontSize',20);
ylabel('I_{MEM}','FontSize',24);
ylim([-3e-3 3e-3]);