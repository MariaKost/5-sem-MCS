#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <string.h>
#include <time.h>

typedef struct sort_par {
    int n;
    int m;
    int* array;
    int* tmp;
} sort_par;

typedef struct merge {
    int l_size;
    int r_size;
    int* left;
    int* right;
    int* dest;
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
    if (chunk->tmp != chunk->dest) {
        memcpy(chunk->dest, chunk->tmp, (chunk->l_size + chunk->r_size) * sizeof(int));
    }
}

void merge_sort(sort_par* par) {
    if (par->n <= par->m) {
        qsort(par->array, par->n, sizeof(int), cmp);
    } else {
        int mid = par->n / 2;

       sort_par left_par;
        left_par.n = mid;
        left_par.m = par->m;
        left_par.array = par->array;
        left_par.tmp = par->tmp;

        sort_par right_par;
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

int binary_search(int* array, int l, int r, int val) {
    int mid = l + (r - l) / 2;
    if (r - l == 0)
        return r;
    else
        if (array[mid] == val)
            return mid;
        else
            if (array[mid] > val)
                return binary_search(array, l, mid, val);
            else
                return binary_search(array, mid + 1, r, val);

}


void mergerr(int* array, int* second_array, int left1, int right1, int left2, int right2, int index) {
    int l = left1;
    int r = left2;
    while (l < right1 && r < right2) {
        if (array[l] < array[r]) {
            second_array[index] = array[l];
            l++;
        } else {
            second_array[index] = array[r];
            r++;
        }
        index++;
    }
    while (l < right1) {
        second_array[index] = array[l];
        l++;
        index++;
    }
    while (r < right2) {
        second_array[index] = array[r];
        r++;
        index++;
    }
}

void paralell_merge_sort(int* array, int* second_array, int left, int right, int m) {
    if (right - left <= m) {
        qsort(&array[left], right - left, sizeof(int), cmp);
    } else {

        int mid = (right + left) / 2;

#pragma omp parallel
        {
#pragma omp single
            {
#pragma omp task
                {
                    paralell_merge_sort(array, second_array, left, mid, m);
                }
#pragma omp task
                {
                    paralell_merge_sort(array, second_array, mid, right, m);
                }
#pragma omp taskwait
            }
        }

        int left_mid = (left + mid) / 2;
        int right_mid = binary_search(array, mid, right, array[left_mid]);

        second_array[left_mid + right_mid - mid] = array[left_mid];

#pragma omp parallel
        {
#pragma omp single
            {
#pragma omp task
                {
                    mergerr(array, second_array, left, left_mid, mid, right_mid, left);
                }
#pragma omp task
                {
                    mergerr(array, second_array, left_mid + 1, mid, right_mid, right, left_mid + right_mid - mid + 1);
                }
#pragma omp taskwait
            }
        }
        memcpy(&array[left], &second_array[left], sizeof(int) * (right - left));
    }
}

int main(int argc, char** argv) {
    if (argc != 4) {
        printf("Wrong number of arguments");
        return 0;
    }

    int n = atoi(argv[1]);
    int m = atoi(argv[2]);
    int P = atoi(argv[3]);

    srand(time(NULL));
    omp_set_num_threads(P);

    int* ar = (int*) malloc(n * sizeof(int));
    int* tmp = (int*) malloc(n * sizeof(int));

    FILE* data = fopen("data.txt", "w");
    fprintf(data, "Not sorted: \n");

    for (int i = 0; i < n; i++) {
        ar[i] = rand();
        fprintf(data, "%d ", ar[i]);
    }
    fprintf(data, "\n");

    memcpy(tmp, ar, n * sizeof(int));

    sort_par par = {
            .n = n,
            .m = m,
            .array = ar,
            .tmp = tmp,
    };

    double start = omp_get_wtime();
    if (P == 1) {
        merge_sort(&par);
    } else {
        paralell_merge_sort(ar, tmp, 0, n, m);
    }
    double end = omp_get_wtime();
    double work_time = end - start;
    FILE* stats = fopen("stats.txt", "a");
    fprintf(stats, "%.5fs %d %d %d\n", work_time, n, m, P);
    fclose(stats);

    fprintf(data, "\nSorted: \n");
    for (int i = 0; i < n; i++) {
        fprintf(data, "%d ", par.array[i]);
    }
    fprintf(data, "\n");
    fclose(data);

    free(ar);
    free(tmp);
    return 0;
}
