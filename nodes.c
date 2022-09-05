#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>
#include <mpi.h>
#include <time.h>
#include <netdb.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <arpa/inet.h>
//##########( Local imports )##########
//#include "shared.c" //The shared functions (Added to base_station.c)
//##########( Topology related )##########
#define SHIFT_ROW 0 //Row (MPI_Cart_shift)
#define SHIFT_COL 1 //Column (MPI_Cart_shift)
#define DISP 1 //Displacement
//##########( Random related )##########
#define NODE_RAND_MAX 7000 //The upper bound for the node's random value
#define NODE_RAND_MIN 5500 //The lower bound for the node's random value
//##########( Event related )##########
//#define NODE_TOLERANCE 300 //Tolerance range between nodes (Added to shared.c)
//##########( SMA related )##########
#define K 5 //The k value

/* Structures */
struct n_thread_args { 
    double * value;
    bool * term_flag;
    MPI_Comm comm;
};

/* Function prototypes */
void *sendFunc(void *arg);
double SMA(int k, double SMA_val, double new_val, double* SMA_array,int SMA_count);

int start_nodes(MPI_Datatype report_type,MPI_Comm commMain, MPI_Comm commNodes, int nrows, int ncols, int threshold ,int t_interval)
{
    /* Function which simulates the sensor nodes */
    int i;
    time_t clock_node = time(NULL);
    int ndims = 2, size, my_rank, reorder, my_cart_rank, ierr;
    MPI_Comm comm2D;
    int nbr_ranks[4]; //Top, Bottom, Left, Right
    int dims[ndims], coord[ndims];
    int wrap_around[ndims];

    double start, end, time_taken; //, comm_time;

    MPI_Comm_size(commNodes, &size);
    MPI_Comm_rank(commNodes, &my_rank);
    /* set dimension size */
    dims[0] = nrows;
    dims[1] = ncols;
    /*************************************************************/
    /* create cartesian topology for processes */
    /*************************************************************/
    MPI_Dims_create(size, ndims, dims);
    /* create cartesian mapping */
    wrap_around[0] = wrap_around[1] = 0; /* periodic shift is false. */
    reorder = 1;
    ierr = 0;
    ierr = MPI_Cart_create(commNodes, ndims, dims, wrap_around, reorder, &comm2D);
    if (ierr != 0)
        printf("ERROR[%d] creating CART\n", ierr);
    fflush(stdout);
    MPI_Cart_coords(comm2D, my_rank, ndims, coord);
    MPI_Cart_rank(comm2D, coord, &my_cart_rank);
    /* find adjacent processes in the cartesian communicator group */
    MPI_Cart_shift(comm2D, SHIFT_ROW, DISP, &nbr_ranks[0], &nbr_ranks[1]);
    MPI_Cart_shift(comm2D, SHIFT_COL, DISP, &nbr_ranks[2], &nbr_ranks[3]);

    /*************************************************************/
    /* communication between nodes */
    /*************************************************************/
    int adj_count = 4; // The number of adjacent nodes (4 for 2D)
    pthread_t hThread;
    MPI_Request send_request[adj_count];
    MPI_Request receive_request[adj_count];
    MPI_Request term;
    MPI_Status dummy_status;

    unsigned int seed;

    report_format full_report;
    /* For height calculation */
    double recvVal[4] = {0, 0, 0, 0}; 
    double SMA_array[K];
    double SMA_val, randomVal;
    int match_count, SMA_count = 0;
    bool term_check = true; //Condition to check if sensor node should terminate
    bool thread_term_flag = true; //For terminating the thread
    void *res;

    struct report_contents report_full;

    /* Thread arguments (for handling communication requests) */
    struct n_thread_args thread_args;
    thread_args.value = &SMA_val; //Reference to the SMA value
    thread_args.comm = comm2D; //The MPI communicator
    thread_args.term_flag = &thread_term_flag; //The termination flag
    pthread_create(&hThread, NULL, sendFunc, (void *)&thread_args);

    seed = time(NULL) + my_rank * my_cart_rank;
    /* Termination message receiver */
    MPI_Irecv(&term_check, 1, MPI_C_BOOL, 0, MSG_TERM_NODES, commMain, &term);
    /* While the termination message is not sent by base station */
    while (term_check)
    {
        match_count = 0;
        start = MPI_Wtime();
        /* Calculates the height */
        randomVal = (rand_r(&seed) % (NODE_RAND_MAX - NODE_RAND_MIN)) + NODE_RAND_MIN;
        SMA_val = SMA(K, SMA_val, randomVal, SMA_array, SMA_count);
        SMA_val = roundf(SMA_val * pow(10, PRECISION)) / pow(10, PRECISION);

        MPI_Barrier(comm2D); //To avoid race condtion
        /* Generation of contents for the report whenever the node height value exceeds the threshold */
        if (SMA_val > threshold)
        {
            /* Send communication requests to receive adjacent nodes height values */
            for (i = 0; i < 4; i++)
            {
                MPI_Isend(&my_cart_rank, 1, MPI_INT, nbr_ranks[i], MSG_NRANK, comm2D, &send_request[i]);
                MPI_Irecv(&recvVal[i], 1, MPI_DOUBLE, nbr_ranks[i], MSG_NVAL, comm2D, &receive_request[i]);
            }
            MPI_Waitall(adj_count, receive_request, MPI_STATUSES_IGNORE);
            /* Check if the value difference between nodes are within the tolerance range */
            for (i = 0; i < adj_count; i++)
            {
                if (abs(recvVal[i] - SMA_val) < NODE_TOLERANCE)
                    match_count += 1;
            }
            /* If more than two adjacent nodes match */
            if (match_count >= 2)
            {
                full_report.node_report[0] = SMA_val;
                full_report.match_count = match_count;
                /* Information about self */
                full_report.coords[0] = coord[0];
                full_report.coords[1] = coord[1];
                full_report.node_ranks[0] = my_cart_rank;
                /* Information about neighbour nodes */
                for (i = 0; i < adj_count; i++){
                    full_report.node_ranks[i + 1] = nbr_ranks[i];
                    full_report.node_report[i + 1] = recvVal[i];
                }
                full_report.iteration = SMA_count + 1;
                full_report.tolerance = NODE_TOLERANCE;
                clock_node = time(NULL);
                sprintf(full_report.curr_datetime, "%s", asctime(localtime(&clock_node)));

                char host[256];
                char *IP;
                struct hostent *host_entry;
                int hostname;
                hostname = gethostname(host, sizeof(host)); //find the host name
                host_entry = gethostbyname(host); //find host information
                IP = inet_ntoa(*((struct in_addr *)host_entry->h_addr_list[0]));
                sprintf(full_report.ip_address, "%s", IP);
                full_report.comm_time = MPI_Wtime();
                MPI_Send(&full_report, 1, report_type, 0, MSG_REPORT, MPI_COMM_WORLD);
            }
        }

        MPI_Barrier(comm2D); //To avoid race condition
        SMA_count += 1;
        end = MPI_Wtime();
        time_taken = (end - start);
        delay(t_interval - time_taken);
    }
    /* Terminates the thread */
    thread_term_flag = false; 
    pthread_join(hThread, &res);
    MPI_Barrier(comm2D);
    MPI_Comm_free(&comm2D);
    return 0;
}

/* Function Definitions */
void *sendFunc(void *arg)
{
    /* Thread function which handles the communication between nodes */
    struct n_thread_args *args = (struct n_thread_args *)arg;
    int n_rank;
    int recv_flag;
    MPI_Status receive_status;
    /* While termination flag not set by main thread */
    while (*args->term_flag)
    {
        /* Probes to check for any communication requests */
        MPI_Iprobe(MPI_ANY_SOURCE, 1, args->comm, &recv_flag, MPI_STATUS_IGNORE);
        if (recv_flag)
        {
            /* Receive rank of requesting node and sends them the height value of this sensor node */
            MPI_Recv(&n_rank, 1, MPI_INT, MPI_ANY_SOURCE, MSG_NRANK, args->comm, &receive_status);
            MPI_Send(args->value, 1, MPI_DOUBLE, n_rank, MSG_NVAL, args->comm);
        }
    }
}

double SMA(int k, double SMA_val, double new_val, double* SMA_array,int SMA_count)
{
    /* Function used to calculate the SMA value for the height */
    int i;
    double res;
    if (SMA_count == 0){ 
        /* Base case (for first entry) */
        res = new_val;
        SMA_array[0] = res;
    } else if (k > SMA_count){ 
        /* When current entry is less than k (Cumulative moving average) */
        res = SMA_val+(new_val-SMA_val)/(SMA_count+1);
        SMA_array[SMA_count] = new_val;
    } else {
        /* Simple moving average */
        res = SMA_val+(new_val-SMA_array[0])/k;
        for(i=0;i<k-1;i++){
            SMA_array[i] = SMA_array[i+1];
        }
        SMA_array[k-1] = new_val;
    }
    return res;
}