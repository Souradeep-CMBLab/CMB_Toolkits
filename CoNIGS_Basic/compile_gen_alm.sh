#!bin/bash
g++ -c gen_alm_img_L_1_final_new_3.c 
g++ -o gen_alm_img_L_1_final_new_3.exe gen_alm_img_L_1_final_new_3.o
g++ -c gen_alm_real_L_1_final_new_3.c                                         
g++ -o gen_alm_real_L_1_final_new_3.exe gen_alm_real_L_1_final_new_3.o
rm *.o
