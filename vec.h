#ifndef __it_vec_h_nto
#define __it_vec_h_nto

void vec_rand_bin( vec v, int n );

/*
void porteuses_2_init (int sub_dim, vec v1, vec v2);
void porteuses_3_init (int sub_dim, vec v1, vec v2, vec v3);*/

void vec_fill_ones(vec v, int dim_SP);
void vec_full_ones(vec v);
void vec_show(vec v);
void vec_norm_patern(vec p);
void split_HFBF(vec all, vec BF, vec HF);

#endif
