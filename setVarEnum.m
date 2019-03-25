function setVarEnum()
  global f_var;

  f_var.Akt = 0;
  f_var.Akt_PIP3 = 0;
  f_var.Akt_PI_P = 0;
  f_var.Akt_PI_PP = 0;
  f_var.ERK = 0;
  f_var.ERKP = 0;
  f_var.ERKPP = 0;
  f_var.GS = 0;
  f_var.HRG = 0;
  f_var.internalization = 0;
  f_var.MEK = 0;
  f_var.MEKP = 0;
  f_var.MEKPP = 0;
  f_var.PI = 0;
  f_var.PI3K = 0;
  f_var.PI3K_act = 0;
  f_var.PIP3 = 0;
  f_var.R = 0;
  f_var.RP = 0;
  f_var.R_HRG = 0;
  f_var.R_HRG2 = 0;
  f_var.R_PI3K = 0;
  f_var.R_PI3K_act = 0;
  f_var.R_ShGS = 0;
  f_var.R_ShP = 0;
  f_var.R_Shc = 0;
  f_var.Raf = 0;
  f_var.Raf_act = 0;
  f_var.RasGDP = 0;
  f_var.RasGTP = 0;
  f_var.ShGS = 0;
  f_var.ShP = 0;
  f_var.Shc = 0;
  f_var.E = 0;
  f_var.MKP3 = 0;
  f_var.PP2A = 0;

  variable = fieldnames(f_var);
  lenVars = length(variable);
  for i=1:lenVars
    f_varName = char(variable(i));
    f_var.(f_varName) = i;
  end

  f_var.size = lenVars;

end