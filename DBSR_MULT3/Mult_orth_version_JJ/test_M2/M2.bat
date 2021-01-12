cp ../test_mult_2conf_jj.exe .
test_mult_2conf_jj  4p5_5s.c  4p5_5p.c  M2  tab=coef_M2
  
rm mult_bnk*
dbsr_mult3a 4p5_5s.c 4p5_5p.c M2  
mult_jj c1=4p5_5s.c c2=4p5_5p.c  tab=coef_M2_mult   jort=-1   bnk=mult_bnk_M2

