#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <omp.h>

typedef struct one_particle {
    int finish_in_b;
    int life_time;
} one_particle;

typedef struct result {
    float p;
    float average_life_time;
    float work_time;
} result;

typedef struct rand_walk_param {
    int a;
    int b;
    int x;
    int N;
    int P;
    float p;
} rand_walk_param;

unsigned int seed_for_threads[128];

int new_step(rand_walk_param* par) {
    int thread = omp_get_thread_num();
    float probability = (float) rand_r(&seed_for_threads[thread]) / RAND_MAX;
    if (probability > par->p)
        return -1;
    else
        return 1;
}

one_particle one_par_motion(rand_walk_param* par) {
    int position = par->x;
    int finish_in_b = 0;
    int life_time = 0;

    while (true) {
        if (position <= par->a) {
            finish_in_b = 0;
            break;
        } else
            if (position >= par->b) {
                finish_in_b = 1;
                break;
            }
        position += new_step(par);
        life_time++;
    }

    one_particle res;
    res.finish_in_b = finish_in_b;
    res.life_time = life_time;

    return res;
}

void random_walk_paralell(rand_walk_param* par, result* res) {

    int num_finished_in_b = 0;
    int sum_life_time = 0;

    float t_start = omp_get_wtime();
#pragma omp parallel for schedule(dynamic) reduction(+: num_finished_in_b, sum_life_time)
    for (int i = 0; i < par->N; ++i) {
        one_particle particle_result = one_par_motion(par);
        num_finished_in_b += particle_result.finish_in_b;
        sum_life_time += particle_result.life_time;
    }
    float t_finish = omp_get_wtime();

    res->average_life_time = (float) sum_life_time / par->N;
    res->p = (float) num_finished_in_b / par->N;
    res->work_time = t_finish - t_start;
}

void random_walk(rand_walk_param* par, result* res) {

    int num_finished_in_b = 0;
    int sum_life_time = 0;
    one_particle particle_result;

    float t_start = omp_get_wtime();
    for (int i = 0; i < par->N; ++i) {
        particle_result = one_par_motion(par);
        num_finished_in_b += particle_result.finish_in_b;
        sum_life_time += particle_result.life_time;
    }
    float t_finish = omp_get_wtime();

    res->average_life_time = (float) sum_life_time / par->N;
    res->p = (float) num_finished_in_b / par->N;
    res->work_time = t_finish - t_start;
}


int main(int argc, char** argv) {
    if (argc != 7) {
        printf("Wrong number of arguments\n");
        return 0;
    }

    rand_walk_param par;
    par.a = atoi(argv[1]);
    par.b = atoi(argv[2]);
    par.x = atoi(argv[3]);
    par.N = atoi(argv[4]);
    par.p = atof(argv[5]);
    par.P = atoi(argv[6]);

    omp_set_num_threads(par.P);
    srand(time(NULL));

    for (int i = 0; i < par.P; i++) {
        seed_for_threads[i] = rand();
    }
    result res;

    if (par.P < 2) {
        random_walk(&par, &res);
    } else {
        random_walk_paralell(&par, &res);
    }

    FILE* file = fopen("stats.txt", "a");
    if (file == NULL) {
        printf("Error opening file!\n");
        exit(1);
    }

    fseek(file, 0, SEEK_END);
    fprintf(file, "%.2f %.1f %.5fs %d %d %d %d %.2f %d\n", res.p, res.average_life_time, res.work_time, par.a, par.b, par.x, par.N, par.p, par.P);
    fclose(file);

    return 0;
}
