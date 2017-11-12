#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <string.h>
#include <time.h>

typedef struct sort_param {
    int n;
    int m;
    int* array;
    int* tmp;
} sort_param;

typedef struct merge {
    int l_size;
    int r_size;
    int* left;
    int* right;
    int* target;
    int* tmp;
    int m;
} merge;

int cmp (const void * a, const void * b) {
    return ( *(int*)a - *(int*)b );
}

void merger(merge* chunk) {
    int l = 0;
    int r = 0;
    int index = 0;

    while (l < chunk->l_size && r < chunk->r_size) {
        if (chunk->left[l] < chunk->right[r]) {
            chunk->tmp[index] = chunk->left[l];
            l++;
        } else {
            chunk->tmp[index] = chunk->right[r];
            r++;
        }
        index++;
    }

    while (l < chunk->l_size) {
        chunk->tmp[index] = chunk->left[l];
        l++;
        index++;
    }
    while(r < chunk->r_size) {
        chunk->tmp[index] = chunk->right[r];
        r++;
        index++;
    }

    if (chunk->tmp != chunk->target) {
        memcpy(chunk->target, chunk->tmp, (chunk->l_size + chunk->r_size) * sizeof(int));
    }
}

void merge_sort(sort_param* par) {
    if (par->n <= par->m) {
        qsort(par->array, par->n, sizeof(int), cmp);
    } else {

        int mid = par->n / 2;
        
        sort_param left_par;
        left_par.n = mid;
        left_par.m = par->m;
        left_par.array = par->array;
        left_par.tmp = par->tmp;

	    sort_param right_par;
        right_par.n = par->n - mid;
        right_par.m = par->m;
        right_par.array = par->array + mid;
        right_par.tmp = par->tmp + mid;

        merge_sort(&left_par);
        merge_sort(&right_par);

        merge curr_merge; 
        curr_merge.l_size = mid;
        curr_merge.r_size = par->n - mid;
        curr_merge.left = par->array;
        curr_merge.right = par->array + mid;
        curr_merge.target = par->array;
        curr_merge.tmp = par->tmp;
        curr_merge.m =par->m;

        merger(&curr_merge);
    }
}

int main(int argc, char** argv) {
    if (argc != 4) {
        printf("Wrong number of arguments");
        return 0;
    }

    sort_param par;
    par.n = atoi(argv[1]);
    par.m = atoi(argv[2]);
    int P = atoi(argv[3]);

    int* ar = (int*) malloc(par.n * sizeof(int));
    int* tmp = (int*) malloc(par.n * sizeof(int));

    FILE* data = fopen("data.txt", "w");
    srand(time(NULL));
    for (int i = 0; i < par.n; i++) {
        ar[i] = rand();
        fprintf(data, "%d ", ar[i]);
    }
    fprintf(data, "\n");

    par.array = ar;
    par.tmp = tmp;

    free(ar);
    free(tmp);

    float start = omp_get_wtime();
    merge_sort(&par);
    float end = omp_get_wtime();
    float work_time = end - start;

    FILE* stats = fopen("stats.txt", "w");
    fprintf(stats, "%f %d %d %d\n", work_time, par.n, par.m, P);
    fclose(stats);

    for (int i = 0; i < par.n; i++) {
        fprintf(data, "%d ", par.array[i]);
    }
    fprintf(data, "\n");
    fclose(data);

    return 0;
}
