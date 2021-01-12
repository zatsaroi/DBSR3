cp ../test_mult_2conf_jj.exe .
test_mult_2conf_jj  1.c  2.c  E2  tab=coef_E2
  
rm mult_bnk*
dbsr_mult3 1.c 2.c E2  
mult_jj c1=1.c c2=2.c  tab=coef_E2_mult   jort=-1   bnk=mult_bnk_E2

