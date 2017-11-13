#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <string.h>
#include <pthread.h>

typedef struct sort_param {
    int n;
    int m;
    int* array;
    int* tmp;
    int num_thr;
    pthread_t* threads;
} sort_param;

typedef struct merge {
    int l_size;
    int r_size;
    int* left;
    int* right;
    int* dest;
    int* tmp;
    int m;
    int num_thr;
    pthread_t* threads;
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

    if (chunk->tmp != chunk->dest) {
        memcpy(chunk->dest, chunk->tmp, (chunk->l_size + chunk->r_size) * sizeof(int));
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
        curr_merge.dest = par->array;
        curr_merge.tmp = par->tmp;
        curr_merge.m =par->m;

        merger(&curr_merge);
    }
}

int bin_search(int* ar, int l, int r, int val) {
    int mid = l + (r - l) / 2;
    if (r - l == 0) {
        return r;
    } else
    if (ar[mid] == val) {
        return mid;
    } else
    if (ar[mid] > val) {
        return bin_search(ar, l, mid, val);
    } else {
        return bin_search(ar, mid + 1, r, val);
    }
}

void* parallel_merge(void* ptr_chunk) {
    merge* chunk = (merge*) ptr_chunk;
    if (chunk->num_thr <= 1 || chunk->l_size < chunk->m || chunk->r_size < chunk->m) {
        merger(chunk);
    } else {
        int left_value = chunk->left[chunk->l_size / 2];
        int right_pos = bin_search(chunk->right, 0, chunk->r_size, left_value);

        merge left_chunk;
        left_chunk.l_size = chunk->l_size / 2;
        left_chunk.r_size = right_pos;
        left_chunk.left = chunk->left;
        left_chunk.right = chunk->right;
        left_chunk.dest = chunk->dest;
        left_chunk.tmp = chunk->tmp;
        left_chunk.m = chunk->m;
        left_chunk.num_thr = (chunk->num_thr - 2) / 2;
        left_chunk.threads = chunk->threads + 2;

        merge right_chunk;
        right_chunk.l_size = chunk->l_size - chunk->l_size / 2;
        right_chunk.r_size = chunk->r_size - right_pos;
        right_chunk.left = chunk->left + chunk->l_size / 2;
        right_chunk.right = chunk->right + right_pos;
        right_chunk.dest = chunk->dest + chunk->l_size / 2 + right_pos;
        right_chunk.tmp = chunk->tmp + chunk->l_size / 2 + right_pos;
        right_chunk.m = chunk->m;
        right_chunk.num_thr = chunk->num_thr - 2 -left_chunk.num_thr;
        right_chunk.threads = chunk->threads + 2 + left_chunk.num_thr;

        pthread_t left_thread = *(chunk->threads);
        pthread_t right_thread = *(chunk->threads + 1);
        pthread_create(&left_thread, NULL, parallel_merge, &left_chunk);
        pthread_create(&right_thread, NULL, parallel_merge, &right_chunk);

        pthread_join(left_thread, NULL);
        pthread_join(right_thread, NULL);
    }
    return NULL;
}

void* parallel_merge_sort(void * ptr_par) {
    sort_param* par = (sort_param*)ptr_par;
    if (par->num_thr <= 1) {
        merge_sort(par);
    } else if (par->n <= par->m) {
        qsort(par->array, par->n, sizeof(int), cmp);
    } else {
        int mid = par->n / 2;

        sort_param left_par;
        left_par.n = mid;
        left_par.m = par->m;
        left_par.array = par->array;
        left_par.tmp = par->tmp;
        left_par.threads = par->threads + 2;
        left_par.num_thr = (par->num_thr - 2) / 2;

        sort_param right_par;
        right_par.n = par->n - mid;
        right_par.m = par->m;
        right_par.array = par->array + mid;
        right_par.tmp = par->tmp + mid;
        right_par.threads = par->threads + 2 + left_par.num_thr;
        right_par.num_thr = par->num_thr - 2 - left_par.num_thr;

        pthread_t left_thread = *(par->threads);
        pthread_t right_thread = *(par->threads + 1);
        pthread_create(&left_thread, NULL, parallel_merge_sort, &left_par);
        pthread_create(&right_thread, NULL, parallel_merge_sort, &right_par);

        pthread_join(left_thread, NULL);
        pthread_join(right_thread, NULL);

        merge curr_merge;
        curr_merge.l_size = mid;
        curr_merge.r_size = par->n - mid;
        curr_merge.left = par->array;
        curr_merge.right = par->array + mid;
        curr_merge.dest = par->array;
        curr_merge.tmp = par->tmp;
        curr_merge.m =par->m;

        parallel_merge(&curr_merge);
        memcpy(par->array, curr_merge.dest, sizeof(int) * par->n);
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
    par.num_thr = atoi(argv[3]);

    pthread_t threads[par.num_thr];
    par.threads = threads;

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
    parallel_merge_sort(&par);
    float end = omp_get_wtime();
    float work_time = end - start;

    FILE* stats = fopen("stats.txt", "a");
    fseek(stats, 0, SEEK_END);
    fprintf(stats, "%f %d %d %d\n", work_time, par.n, par.m, par.num_thr);
    fclose(stats);

    for (int i = 0; i < par.n; i++) {
        fprintf(data, "%d ", par.array[i]);
    }
    fprintf(data, "\n");
    fclose(data);


    return 0;
}