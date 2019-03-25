global f_var;

y0 = initialValues();
options = odeset('RelTol',1e-9,'AbsTol',1e-9);

t = 0:1:1800;

condition = 3;

sim_RP        = zeros(length(t),condition);
sim_ShP       = zeros(length(t),condition);
sim_PI3K_act  = zeros(length(t),condition);
sim_Raf_act   = zeros(length(t),condition);
sim_MEKPP     = zeros(length(t),condition);
sim_ERKPP     = zeros(length(t),condition);
sim_Akt_PI_PP = zeros(length(t),condition);

for i = 1:condition
    if i == 1
        y0(f_var.HRG) = 330.0; % 10nM
    elseif i == 2
        y0(f_var.HRG) = 33.0;  % 1nM
    elseif i == 3
        y0(f_var.HRG) = 3.3;   % 0.1nM
    end

    [T,Y] = ode15s(@diffeq,[0 1800],y0,options);

    sim_RP(:,i) = interp1(T,2*(Y(:,f_var.RP)+Y(:,f_var.R_PI3K)+Y(:,f_var.R_PI3K_act)+...
                  Y(:,f_var.R_ShGS)+Y(:,f_var.R_ShP)+Y(:,f_var.R_Shc))./y0(f_var.R).*100.0,t);
    sim_ShP(:,i) = interp1(T,Y(:,f_var.ShP)./y0(f_var.Shc).*100.0,t);
    sim_PI3K_act(:,i) = interp1(T,Y(:,f_var.PI3K_act)./y0(f_var.PI3K).*100.0,t);
    sim_Raf_act(:,i) = interp1(T,Y(:,f_var.Raf_act)./y0(f_var.Raf).*100.0,t);
    sim_MEKPP(:,i) = interp1(T,Y(:,f_var.MEKPP)./y0(f_var.MEK).*100.0,t);
    sim_ERKPP(:,i) = interp1(T,Y(:,f_var.ERKPP)./y0(f_var.ERK).*100.0,t);
    sim_Akt_PI_PP(:,i) = interp1(T,Y(:,f_var.Akt_PI_PP)./y0(f_var.Akt).*100.0,t);

end