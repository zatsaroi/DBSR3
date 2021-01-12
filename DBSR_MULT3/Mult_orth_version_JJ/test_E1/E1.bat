cp ../test_mult_2conf_jj.exe .
test_mult_2conf_jj  4p5_5s.c  4p5_5p.c  E1  tab=coef_E1
  
rm mult_bnk*
dbsr_mult3 4p5_5s.c 4p5_5p.c E1  
mult_jj c1=4p5_5s.c c2=4p5_5p.c  tab=coef_E1_mult   jort=-1

