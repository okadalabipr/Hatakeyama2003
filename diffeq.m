function dydt = diffeq(t,y)
  global f_var;

  dydt = zeros(f_var.size,1);

  v = reaction(y);

  dydt(f_var.Akt) = -v(29);
  dydt(f_var.Akt_PIP3) = v(29) -v(30) + v(31);
  dydt(f_var.Akt_PI_P) = v(30) - v(31) - v(32) + v(33);
  dydt(f_var.Akt_PI_PP) = v(32) - v(33);
  dydt(f_var.ERK) = -v(19) + v(20);
  dydt(f_var.ERKP) = v(19) - v(20) - v(21) + v(22);
  dydt(f_var.ERKPP) = v(21) - v(22);
  dydt(f_var.GS) = -v(7) + v(9);
  dydt(f_var.HRG) = -v(1);
  dydt(f_var.internalization) = v(34);
  dydt(f_var.MEK) = -v(15) + v(16);
  dydt(f_var.MEKP) = v(15) - v(16) - v(17) + v(18);
  dydt(f_var.MEKPP) = v(17) - v(18);
  dydt(f_var.PI) = -v(27) + v(28);
  dydt(f_var.PI3K) = -v(23) + v(26);
  dydt(f_var.PI3K_act) = v(25) - v(26);
  dydt(f_var.PIP3) = v(27) - v(28) - v(29);
  dydt(f_var.R) = -v(1);
  dydt(f_var.RP) = v(3) - v(4) - v(5) + v(8) - v(23) + v(25) - v(34);
  dydt(f_var.R_HRG) = v(1) - 2*v(2);
  dydt(f_var.R_HRG2) = v(2) - v(3) + v(4);
  dydt(f_var.R_PI3K) = v(23) - v(24);
  dydt(f_var.R_PI3K_act) = v(24) - v(25);
  dydt(f_var.R_ShGS) = v(7) - v(8);
  dydt(f_var.R_ShP) = v(6) - v(7);
  dydt(f_var.R_Shc) = v(5) - v(6);
  dydt(f_var.Raf) = -v(13) + v(14);
  dydt(f_var.Raf_act) = v(13) - v(14);
  dydt(f_var.RasGDP) = -v(11) + v(12);
  dydt(f_var.RasGTP) = v(11) - v(12);
  dydt(f_var.ShGS) = v(8) - v(9);
  dydt(f_var.ShP) = v(9) - v(10);
  dydt(f_var.Shc) = -v(5) + v(10);

end
