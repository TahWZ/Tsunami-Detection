#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>
#include <pthread.h>
#include <mpi.h>
#include <signal.h>
#include <time.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <arpa/inet.h>
#include <unistd.h>
//##########(Time related)##########
#define TIME_INTERVAL 2 //The time intervals
//##########(Event related)##########
#define ALTI_TOLERANCE 300 //The tolerance range between node and altimeter
#define ALTI_RANGE 1200 //The value range the altimeter produces (threshold, threshold + ALTI_RANGE)
//##########( FIFO related )##########
#define QUEUE_SIZE 8 //Size of the FIFO queue data strcuture
//##########( Thread related )##########
#define BS_NUM_THREAD 2 //Number of threads used by the base station
#define SENTINAL -1 //The sentinal value for base station shut down (Only for manual termination)
//##########( Local imports )##########
#include "shared.c" //Contains shared structures and functions
#include "nodes.c"  //Contains functions used to simulate the sensor nodes

/* Structures */
struct alti_thread_args
/* Arguments for altimeter thread */
{
    bool* flag;
    int *intervals;
    int threshold;
    int rows;
    int cols;
};

struct time_thread_args
/* Arguments for base station shutdown thread */
{
    bool* flag;
    int time;
};

/* Shared global arrays */
int coordsR[QUEUE_SIZE], coordsC[QUEUE_SIZE];
double seaHeight[QUEUE_SIZE];
char alti_report_time[QUEUE_SIZE][50];

/* Deadlock intializer */
pthread_mutex_t g_Mutex = PTHREAD_MUTEX_INITIALIZER;

/* Function prototypes */
void *timeEnd(void *arg);
void *manualEnd(void *arg);
void *startAltimeter(void *arg);

int main(int argc, char *argv[])
{
    void *res;
    int size, my_rank, provided;
    MPI_Comm commNodes;
    MPI_Datatype report_cont;
    int nrows, ncols, node_count;
    pthread_t hThread[BS_NUM_THREAD];
    int threshold;
    int i,k;
    FILE *fp;

    double start, end, time_taken;
    struct report_contents report_full;

    /* Start up initial MPI environment */
    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    /* Process command line arguments*/
    if (argc < 4){
        if (my_rank == 0)
            printf("Insufficient parameters, at least 3 required:\n@Parameters (m, n, threshold, duration(Optional))\n");
        MPI_Finalize();
        return 0;
    } else {
        nrows = atoi(argv[1]);
        ncols = atoi(argv[2]);
        threshold = atoi(argv[3]);
        node_count = nrows * ncols;
        /* Check if dimension is valid */
        if ((nrows * ncols) != size - 1)
        {
            if (my_rank == 0)
                printf("ERROR: nrows*ncols)=%d * %d = %d != %d\n", nrows, ncols, nrows * ncols, size);
            MPI_Finalize();
            return 0;
        }
    }
    /* Section responsible for resetting the report file and opening it */
    fp = fopen("report.txt", "w"); //Clear the file
    fclose(fp);
    MPI_File report_file;
    int access_mode = MPI_MODE_CREATE | MPI_MODE_WRONLY;
    char file_buf[256];
    MPI_File_open(MPI_COMM_WORLD, "report.txt", access_mode, MPI_INFO_NULL, &report_file);

    /* Creation of MPI Datatype which contains information for the report */
    int blocklengths[9] = {5, 2, 1, 1, 50, 50, 5, 1, 1};
    MPI_Datatype types[9] = {MPI_DOUBLE, MPI_INT, MPI_INT, MPI_INT, MPI_CHAR, MPI_CHAR, MPI_INT, MPI_DOUBLE, MPI_INT};
    MPI_Aint offsets[9];
    MPI_Aint base_address;

    MPI_Get_address(&report_full, &base_address);
    MPI_Get_address(&report_full.node_report, &offsets[0]);
    MPI_Get_address(&report_full.coords, &offsets[1]);
    MPI_Get_address(&report_full.iteration, &offsets[2]);
    MPI_Get_address(&report_full.tolerance, &offsets[3]);
    MPI_Get_address(&report_full.curr_datetime[0], &offsets[4]);
    MPI_Get_address(&report_full.ip_address[0], &offsets[5]);
    MPI_Get_address(&report_full.node_ranks, &offsets[6]);
    MPI_Get_address(&report_full.comm_time, &offsets[7]);
    MPI_Get_address(&report_full.match_count, &offsets[8]);

    for (i = 0; i< 9; i++){
        offsets[i] = MPI_Aint_diff(offsets[i], base_address);
    }

    MPI_Type_create_struct(9, blocklengths, offsets, types, &report_cont);
    MPI_Type_commit(&report_cont);
    /* End of creation */

    /* Create communicators containing processes simulating the sensor nodes */
    MPI_Comm_split(MPI_COMM_WORLD, my_rank == 0, 0, &commNodes);
    MPI_Barrier(MPI_COMM_WORLD);
    if (my_rank == 0)
    {
        int intervals = 0;
        bool alti_term_flag = true; //Termination flag for altimeter
        bool shutdown = true; //Shutdown flag for base station

        /* Arguments for altimeter */
        struct alti_thread_args arguments;
        arguments.rows = nrows;
        arguments.cols = ncols;
        arguments.threshold = threshold;
        arguments.flag = &alti_term_flag;
        arguments.intervals = &intervals;

        /* Initialization of POSIX threads simulating the altimeter */
        pthread_create(&hThread[0], NULL, startAltimeter, (void *)&arguments);
        if (argc == 5 && atoi(argv[4])>0){
            /* If time argument passed, base station ends after the stated time has passed */
            struct time_thread_args thread_time;
            thread_time.flag = &shutdown;
            thread_time.time = atoi(argv[4]);
            printf("Duration parameter provided, base station will run for %d seconds.\n", atoi(argv[4]));
            pthread_create(&hThread[1], NULL, timeEnd, (void *)&thread_time);
        }else{
            /* If time argument not passed, base station required manual input for termination */
            printf("Duration parameter not given (Manual termination).\n");
            pthread_create(&hThread[1], NULL, manualEnd, (void *)&shutdown);
        }

        int recv_flag; //Flag which checks if any communication requests exist for the base station
        while (true)
        {
            /* Check for any communication requests */
            MPI_Iprobe(MPI_ANY_SOURCE, MSG_REPORT, MPI_COMM_WORLD, &recv_flag, MPI_STATUS_IGNORE);
            if (recv_flag)
            {
                MPI_Recv(&report_full, 1, report_cont, MPI_ANY_SOURCE, MSG_REPORT, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                /* Mutex lock to ensure that the shared global array is not being accessed whilst the code checks for matches between the altimeter and the nodes */
                report_full.comm_time = MPI_Wtime() - report_full.comm_time;
                pthread_mutex_lock(&g_Mutex);
                /* Loop to compare values in the shared global array with report */
                for (k = 0; k < QUEUE_SIZE; k++)
                {
                    /* If global array not filled entirely */
                    if (k > intervals){
                        break;
                    } 
                    /* Section responsible for logging the reports whenever a match is found between the altimeter and a node */
                    else if (coordsR[k]==report_full.coords[0] && coordsC[k]==report_full.coords[1]  &&  abs(seaHeight[k] - report_full.node_report[0]) <= ALTI_TOLERANCE) {
                        /* Log information */
                        sprintf(file_buf,"Iteration: %d\nAlert time: %sAlert type: Match\n\n", report_full.iteration, report_full.curr_datetime);
                        MPI_File_write(report_file, file_buf, strlen(file_buf), MPI_CHAR,MPI_STATUS_IGNORE);
                        /* Report information (Reporting node) */
                        sprintf(file_buf,"Reporting Node     Coord     Height (m)     IPv4\n");
                        MPI_File_write(report_file, file_buf, strlen(file_buf), MPI_CHAR,MPI_STATUS_IGNORE);
                        sprintf(file_buf,"%-19d(%d,%d)     %-15.3f%s\n", report_full.node_ranks[0], report_full.coords[0], report_full.coords[1],report_full.node_report[0],report_full.ip_address);
                        MPI_File_write(report_file, file_buf, strlen(file_buf), MPI_CHAR,MPI_STATUS_IGNORE);
                        /* Report information (Adjacent nodes) */
                        sprintf(file_buf,"\nAdjacent Nodes     Coord     Height (m)\n");
                        MPI_File_write(report_file, file_buf, strlen(file_buf), MPI_CHAR,MPI_STATUS_IGNORE);
                        if (report_full.node_ranks[1] != -2){
                            sprintf(file_buf,"%-19d(%d,%d)     %-15.3f\n", report_full.node_ranks[1],report_full.coords[0]-1, report_full.coords[1],report_full.node_report[1]);
                            MPI_File_write(report_file, file_buf, strlen(file_buf), MPI_CHAR,MPI_STATUS_IGNORE);
                        }
                        if (report_full.node_ranks[2] != -2){
                            sprintf(file_buf,"%-19d(%d,%d)     %-15.3f\n", report_full.node_ranks[2],report_full.coords[0]+1, report_full.coords[1],report_full.node_report[2]);
                            MPI_File_write(report_file, file_buf, strlen(file_buf), MPI_CHAR,MPI_STATUS_IGNORE);
                        }
                        if (report_full.node_ranks[3] != -2){
                            sprintf(file_buf,"%-19d(%d,%d)     %-15.3f\n", report_full.node_ranks[3],report_full.coords[0], report_full.coords[1]-1,report_full.node_report[3]);
                            MPI_File_write(report_file, file_buf, strlen(file_buf), MPI_CHAR,MPI_STATUS_IGNORE);
                        }
                        if (report_full.node_ranks[4] != -2){
                            sprintf(file_buf,"%-19d(%d,%d)     %-15.3f\n", report_full.node_ranks[4],report_full.coords[0], report_full.coords[1]+1,report_full.node_report[4]);
                            MPI_File_write(report_file, file_buf, strlen(file_buf), MPI_CHAR,MPI_STATUS_IGNORE);
                        }
                        /* Satelite alimeter details */
                        sprintf(file_buf,"\nSatelite altimeter reporting time: %s", alti_report_time[k]);
                        MPI_File_write(report_file, file_buf, strlen(file_buf), MPI_CHAR,MPI_STATUS_IGNORE);
                        sprintf(file_buf,"Satelite altimeter reporting height (m): %.3f\n", seaHeight[k]);
                        MPI_File_write(report_file, file_buf, strlen(file_buf), MPI_CHAR,MPI_STATUS_IGNORE);
                        sprintf(file_buf,"Satelite altimeter reporting Coords: (%d,%d)\n", coordsR[k],coordsC[k]);
                        MPI_File_write(report_file, file_buf, strlen(file_buf), MPI_CHAR,MPI_STATUS_IGNORE);

                        /* Additional details */
                        sprintf(file_buf,"\nCommunication Time (seconds): %.8f\n", report_full.comm_time);
                        MPI_File_write(report_file, file_buf, strlen(file_buf), MPI_CHAR,MPI_STATUS_IGNORE);
                        sprintf(file_buf,"Total Messages send between reporting node and base station: 1\n");
                        MPI_File_write(report_file, file_buf, strlen(file_buf), MPI_CHAR,MPI_STATUS_IGNORE);
                        sprintf(file_buf,"Number of adjacent matches to reporing node: %d\n", report_full.match_count);
                        MPI_File_write(report_file, file_buf, strlen(file_buf), MPI_CHAR,MPI_STATUS_IGNORE);
                        sprintf(file_buf,"Max. tolerance range between nodes readings (m): %d\n", NODE_TOLERANCE);
                        MPI_File_write(report_file, file_buf, strlen(file_buf), MPI_CHAR,MPI_STATUS_IGNORE);
                        sprintf(file_buf,"Max. tolerance range between satelite altimeter and reporting nodes readings (m): %d\n", ALTI_TOLERANCE);
                        MPI_File_write(report_file, file_buf, strlen(file_buf), MPI_CHAR,MPI_STATUS_IGNORE);
                        /* Line skip */
                        sprintf(file_buf,"-------------------------------------------------------------------------------------\n");
                        MPI_File_write(report_file, file_buf, strlen(file_buf), MPI_CHAR,MPI_STATUS_IGNORE);
                        break;
                    }
                }
                /* Mutex unlocks when the checks are finished */
                pthread_mutex_unlock(&g_Mutex);
            } else if (!shutdown) {
                break;
            }
        }
        MPI_Request term[node_count];
        /* Termination of the nodes */
        bool term_msg = false;
        for (i = 1; i <= node_count; i++)
        {   
            MPI_Isend(&term_msg, 1, MPI_C_BOOL, i, MSG_TERM_NODES, MPI_COMM_WORLD, &term[i-1]);
        }
        /* Termination of the altimeter */
        alti_term_flag = false;
        for (i = 0; i < BS_NUM_THREAD; i++)
        {
            pthread_join(hThread[i], &res); //Checks if threads are all terminated
        }
        printf("Altimeter terminated successfully\n");
        /* Waits for sensor nodes to terminate */
        MPI_Waitall(node_count, term, MPI_STATUSES_IGNORE);
        printf("Termination message successfully sent to sensor nodes\n");
    }
    else
    {
        /* Function which simulates the sensor node (In nodes.c) */
        start_nodes(report_cont,MPI_COMM_WORLD, commNodes, nrows, ncols, threshold, TIME_INTERVAL);
    }
    /* File is closed after specified runtime is finished */
    MPI_File_close(&report_file);
    MPI_Barrier(MPI_COMM_WORLD);
    if (my_rank==0){
        printf("Sensor nodes terminated successfully\n");
        printf("Base station will now shut down\n");
    }
    MPI_Barrier(MPI_COMM_WORLD);
    pthread_mutex_destroy(&g_Mutex);
    MPI_Finalize();
    return 0;
}

/* Function definitions */
void *manualEnd(void *arg)
{
    /* Thread function used to allows users to manually terminate the base station */
    int *end = (int *)arg;
    char check[5];
    printf("Enter -1 to shut down the base station: \n");
    while (1){
        fgets(check, 5, stdin);
        /* Check if user input is the sentinal value */
        if (atoi(check) == SENTINAL)
            break;
    }
    printf("Sentinal value entered, base station will now begin shut down sequence\n");
    *end = false;
}

void *timeEnd(void *arg)
{
    /* Thread function used to terminate the base station after a set time */
    struct time_thread_args *args = (struct time_thread_args *)arg;
    delay(args->time);
    printf("Time limit reached, base station will now begin shut down sequence\n");
    *args->flag = false;
}

void *startAltimeter(void *arg)
{
    /* Thread function which simulates the altimeter */
    int i;
    time_t clock_node;
    struct alti_thread_args *args = (struct alti_thread_args *)arg;
    double start, end, time_taken;
    double height; 
    int coordr, coordc; //Represents the random height value generated at a certain coordinate
    int seed = time(0) + pthread_self();
    while (*args->flag)
    {
        start = MPI_Wtime();
        /* Mutex lock to ensure that the shared global array is not being accessed whilst it is being updated */
        pthread_mutex_lock(&g_Mutex);
        height = args -> threshold + (rand_r(&seed) % ((int)(ALTI_RANGE*pow(10,PRECISION))))/pow(10, PRECISION);
        coordr = rand_r(&seed) % args->rows;
        coordc = rand_r(&seed) % args->cols;
        /* First In First Out (FIFO) implementation of generating coordinates */
        if (*args->intervals >= QUEUE_SIZE)
        {
            for (i = 0; i < QUEUE_SIZE-1; i++)
            {
                /* Values are pushed to the left of the array */
                coordsR[i] = coordsR[i + 1];    
                coordsC[i] = coordsC[i + 1];
                seaHeight[i] = seaHeight[i + 1];
                sprintf(alti_report_time[i], "%s", alti_report_time[i+1]);
            }
            /* New value is then inserted at the last index */
            coordsR[QUEUE_SIZE-1] = coordr;
            coordsC[QUEUE_SIZE-1] = coordc;
            seaHeight[QUEUE_SIZE-1] = height;
            clock_node = time(NULL);
            sprintf(alti_report_time[QUEUE_SIZE-1], "%s", asctime(localtime(&clock_node)));
            *args->intervals += 1;
        }
        else
        {
            /* Normal insertion into the array (for when ) */
            coordsR[*args->intervals] = coordr;
            coordsC[*args->intervals] = coordc;
            seaHeight[*args->intervals] = height;
            clock_node = time(NULL);
            sprintf(alti_report_time[*args->intervals], "%s", asctime(localtime(&clock_node)));
            *args->intervals += 1;
        }
        /* Mutex is unlocked after update */
        pthread_mutex_unlock(&g_Mutex);
        end = MPI_Wtime();
        time_taken = (end - start);
        delay(TIME_INTERVAL - time_taken);
    }
}