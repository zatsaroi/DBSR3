cp ../test_mult_2conf_jj.exe .
test_mult_2conf_jj  1.c  2.c  M1  tab=coef_M1
  
rm mult_bnk*
dbsr_mult3a 1.c 2.c M1  
mult_jj c1=1.c c2=2.c  tab=coef_M1_mult   jort=-1   bnk=mult_bnk_M1

