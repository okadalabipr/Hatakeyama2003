function y0 = initialValues()
    global f_var;

    y0 = zeros(f_var.size,1);

    y0(f_var.R) = 80.0;
    y0(f_var.Shc) = 1000.0;
    y0(f_var.GS) = 10.0;
    y0(f_var.RasGDP) = 120.0;
    y0(f_var.Raf) = 100.0;
    y0(f_var.MEK) = 120.0;
    y0(f_var.ERK) = 1000.0;
    y0(f_var.PI3K) = 10.0;
    y0(f_var.PI) = 800.0;
    y0(f_var.Akt) = 10.0;
    y0(f_var.E) = 7.0;
    y0(f_var.MKP3) = 2.4;
    y0(f_var.PP2A) = 11.4;

end