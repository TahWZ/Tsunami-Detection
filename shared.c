#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>

//##########( Shared structures )##########

typedef struct report_contents
/* Report structure */
{
    double node_report[5]; //Contains (Current node's value, match count, adjacent nodes values)
    int coords[2];
    int iteration; 
    int tolerance;
    char curr_datetime[50];
    char ip_address[50];
    int node_ranks[5]; //Contains (Current node's rank, adjacent nodes ranks)
    double comm_time;
    int match_count;
} report_format;

//##########( Shared functions )##########

int delay(double seconds)
/* Function used to make the application wait (Used for setting time intervals)*/
{
    struct timespec remaining;
    struct timespec request= { (
        int)(seconds),
        (((int)seconds) % 1) * 1000000000
    };

    return nanosleep(&request , &remaining);
}

//###########( Shared constants )#############
#define NODE_TOLERANCE 300 //Tolerance range between nodes
#define PRECISION 3 //The decimal places for results produced

//###########( MPI tags )#############
#define MSG_TERM_NODES 0 //Used when sending termination message
#define MSG_NRANK 1 //Used for sending node rank
#define MSG_NVAL 2 //Used for sending node's height value
#define MSG_REPORT 3 //Used for sending report