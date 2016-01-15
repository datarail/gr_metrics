function t_out = evaluate_GRvalue(t_in)
t_out = t_in;
t_out.GRvalue = 2.^( log2(t_in.cell_count./t_in.cell_count__time0) ./ ...
    log2(t_in.cell_count__ctrl./t_in.cell_count__time0) ) -1;