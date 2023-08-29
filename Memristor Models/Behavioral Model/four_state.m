clear;
clc;

redondo = importdata('sin_10Hz.csv');
redondo_time = redondo.data(:,1);
redondo_state = round((redondo.data(:,4)*1e7)+1);
redondo_input = redondo.data(:,2);

NumPeriods = 100;
f = 10;
tstop = NumPeriods/f;
dt = tstop*1e-6;
time = 0:dt:tstop-dt;
steps = length(time);

states = 4;
init_state = 1;

P = zeros(states,length(time));
the_state = zeros(1,length(time));
P(init_state,1) = 1;
the_state(1) = init_state;

Kb = 1.38e-23;
q = 1.6e-19;
eps0 = 8.85e-12;
area = 6.4e-9;
d = 4e-8;
scl1 = 5;
scl2 = 7;
scl3 = 9;
AA1 = 1e6;
AA2 = 2e6;
AA3 = 5e7;
Ub = 0.9;
epsr = 8;
T0 = 300;
ktemp = 1;
epsi = epsr*eps0;

ron = 6.1e3;

p_th_1_2 = 0.9e-12;
p_th_2_3 = p_th_1_2 + 1e-11;
p_th_3_4 = p_th_2_3 + 5e-9;
p_th_4_3 = -1.5e-7;
p_th_3_2 = p_th_4_3 - 0.8e-11;
p_th_2_1 = p_th_3_2 - 0.5e-11;
consv_rate = 3e8*dt;

I1 = @(V) area*AA1*ktemp*(T0^2)*exp(-q*Ub/(Kb*T0))*exp(sqrt(abs(V)))*(scl1+q/(Kb*T0)*(sqrt(q/(d*4*pi*epsi))));
I2 = @(V) area*AA2*ktemp*(T0^2)*exp(-q*Ub/(Kb*T0))*exp(sqrt(abs(V)))*(scl2+q/(Kb*T0)*(sqrt(q/(d*4*pi*epsi))));
I3 = @(V) area*AA3*ktemp*(T0^2)*exp(-q*Ub/(Kb*T0))*exp(sqrt(abs(V)))*(scl3+q/(Kb*T0)*(sqrt(q/(d*4*pi*epsi))));
I4 = @(V) abs(V)/ron;

t1_p = @(V) p_th_1_2./(V.*I1(V));
t2_p = @(V) (p_th_2_3 - p_th_1_2)./(V.*I2(V));
t3_p = @(V) (p_th_3_4 - p_th_2_3)./(V.*I3(V));
t3_n = @(V) p_th_4_3./(V.*I4(V));
t2_n = @(V) (p_th_3_2 - p_th_4_3)./(V.*I3(V));
t1_n = @(V) (p_th_2_1 - p_th_3_2)./(V.*I2(V));

r1 = @(V) (V>0)./t1_p(V);
r2 = @(V) (V>0)./t2_p(V);
r3 = @(V) (V>0)./t3_p(V);
r4 = @(V) (V<0)./t3_n(V);
r5 = @(V) (V<0)./t2_n(V);
r6 = @(V) (V<0)./t1_n(V);

Wv = @(V) [-r1(V) r6(V) 0 0;...
    r1(V) -r2(V)-r6(V) r5(V) 0;...
    0 r2(V) -r3(V)-r5(V) r4(V);...
    0 0 r3(V) -r4(V)];

Ws = @(state) [-((state==2)+(state==3)+(state==4))*consv_rate (state==1)*consv_rate (state==1)*consv_rate (state==1)*consv_rate;...
    (state==2)*consv_rate -((state==1)+(state==3)+(state==4))*consv_rate (state==2)*consv_rate (state==2)*consv_rate;...
    (state==3)*consv_rate (state==3)*consv_rate -((state==1)+(state==2)+(state==4))*consv_rate (state==3)*consv_rate;...
    (state==4)*consv_rate (state==4)*consv_rate (state==4)*consv_rate -((state==1)+(state==2)+(state==3))*consv_rate];

W = @(V,state) Wv(V) + Ws(state);

the_bar = waitbar(0,'Please wait...');
for i = 2:length(time)
    V(i-1) = 3*sin(2*pi*f*time(i-1));
    yy1 = W(V(i-1),the_state(i-1))*P(:,i-1);
    yy2 = yy1 + dt/2*W(V(i-1),the_state(i-1))*yy1;
    yy3 = yy2 + dt/2*W(V(i-1),the_state(i-1))*yy2;
    yy4 = yy3 + dt*W(V(i-1),the_state(i-1))*yy3;
    P(:,i) = P(:,i-1) + dt/6*(yy1 + 2*yy2 + 2*yy3 + yy4);
    the_rand = rand;
    new_state = 1*(the_rand<P(1,i)) + ...
        2*(the_rand>=P(1,i) && the_rand<(P(1,i)+P(2,i)))...
        + 3*(the_rand>(P(1,i)+P(2,i)) && the_rand<(P(1,i)+P(2,i)+P(3,i)))...
        + 4*(the_rand>(P(1,i)+P(2,i)+P(3,i)));
    if new_state~=the_state(i-1)
        the_state(i) = new_state;
        P(:,i) = zeros(states,1);
        P(new_state,i) = 1;
    else
        the_state(i) = new_state;
    end
    if mod(i,steps/10)==0
        waitbar(i/steps,the_bar,...
            ['Simulation Completed: ' num2str(i*100/steps) '%']);
    end
end
close(the_bar);
V(length(time)) = 0;

I{1} = @(V) area*AA1*ktemp*(T0^2)*exp(-q*Ub/(Kb*T0))*exp(sqrt(abs(V)))*(scl1+q/(Kb*T0)*(sqrt(q/(d*4*pi*epsi))));
I{2} = @(V) area*AA2*ktemp*(T0^2)*exp(-q*Ub/(Kb*T0))*exp(sqrt(abs(V)))*(scl2+q/(Kb*T0)*(sqrt(q/(d*4*pi*epsi))));
I{3} = @(V) area*AA3*ktemp*(T0^2)*exp(-q*Ub/(Kb*T0))*exp(sqrt(abs(V)))*(scl3+q/(Kb*T0)*(sqrt(q/(d*4*pi*epsi))));
I{4} = @(V) abs(V)/ron;

III_prop = zeros(1,length(V));
for iii = 1:length(V)
    III_prop(iii) = I{the_state(iii)}(V(iii));
end

mean_III_prop = mean(reshape(III_prop,round(1/dt/f),round(length(III_prop)/(1/dt/f))),2);
III_redondo = zeros(1,length(redondo_time));
for iii = 1:length(redondo_time)
    III_redondo(iii) = I{redondo_state(iii)}(redondo_input(iii));
end

figure(2);
fig = gcf;
fig.Color = 'white';
fig.Position = [fig.Position(1) fig.Position(2) 516 fig.Position(4)];
plot(V,abs(III_prop),'o-',redondo_input,abs(III_redondo))
ax = gca;
ax.FontSize = 18;
ax.YScale = 'log';
xlabel('Applied Voltage (V)','FontSize',20);
ylabel('I_{MEM}','FontSize',24);
ylim([1e-12 3e-3]);

switching = [0 (the_state(2:end)-1).^2-(the_state(1:end-1)-1).^2];
the_time = 0:dt:tstop/NumPeriods-dt;
the_switching = reshape(switching,round(length(time)/NumPeriods),NumPeriods);
[m_t_1_2,~] = find(the_switching==1);
the_mean_time_1_2 = mean(the_time(m_t_1_2));
[m_t_2_3,~] = find(the_switching==3);
the_mean_time_2_3 = mean(the_time(m_t_2_3));
[m_t_3_4,~] = find(the_switching==5);
the_mean_time_3_4 = mean(the_time(m_t_3_4));
[m_t_4_3,~] = find(the_switching==-5);
the_mean_time_4_3 = mean(the_time(m_t_4_3));
[m_t_3_2,~] = find(the_switching==-3);
the_mean_time_3_2 = mean(the_time(m_t_3_2));
[m_t_2_1,~] = find(the_switching==-1);
the_mean_time_2_1 = mean(the_time(m_t_2_1));
mean_sw = ones(1,round(length(time)/NumPeriods));
mean_sw(the_time>mean(the_time(m_t_1_2))) = 2;
mean_sw(the_time>mean(the_time(m_t_2_3))) = 3;
mean_sw(the_time>mean(the_time(m_t_3_4))) = 4;
mean_sw(the_time>mean(the_time(m_t_4_3))) = 3;
mean_sw(the_time>mean(the_time(m_t_3_2))) = 2;
mean_sw(the_time>mean(the_time(m_t_2_1))) = 1;


III_mean_time = zeros(1,length(the_time));
for iii = 1:length(the_time)
    III_mean_time(iii) = I{mean_sw(iii)}(V(iii));
end

hold on;
plot(V(1:length(the_time)),III_mean_time,'LineWidth',2);
hold off;
ax.Children(1).LineWidth = 2;
ax.Children(2).LineWidth = 2;