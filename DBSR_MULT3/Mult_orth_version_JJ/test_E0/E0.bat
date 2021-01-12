cp ../test_mult_2conf_jj.exe .
test_mult_2conf_jj  1.c  2.c  E0  tab=coef_E0
  
rm mult_bnk*
dbsr_mult3a 1.c 2.c E0  
mult_jj c1=1.c c2=2.c  tab=coef_E0_mult   jort=-1   bnk=mult_bnk_E0



