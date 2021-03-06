/// \file
/// Performance timer functions.
///
/// Use the timer functionality to collect timing and number of calls
/// information for chosen computations (such as force calls) and
/// communication (such as sends, receives, reductions).  Timing results
/// are reported at the end of the run showing overall timings and
/// statistics of timings across ranks.
///
/// A new timer can be added as follows:
/// -# add new handle to the TimerHandle in performanceTimers.h
/// -# provide a corresponding name in timerName
///
/// Note that the order of the handles and names must be the
/// same. This order also determines the order in which the timers are
/// printed. Names can contain leading spaces to show a hierarchical
/// ordering.  Timers with zero calls are omitted from the report.
///
/// Raw timer data is obtained from the getTime() and getTick()
/// functions.  The supplied portable versions of these functions can be
/// replaced with platform specific versions for improved accuracy or
/// lower latency.
/// \see TimerHandle
/// \see getTime
/// \see getTick
///


#include "performanceTimers.h"

#include <stdio.h>
#include <sys/time.h>
#include <string.h>
#include <stdint.h>
#include <inttypes.h>
#include <math.h>

#include "performanceTimers.h"
#include "mytype.h"
#include "parallel.h"
#include "yamlOutput.h"
#include "read_netstats.h"

#ifdef RUN_PAPI
#include "papi.h"
#endif

static uint64_t getTime(void);
static double getTick(void);
static void timerStats(void);
// aggregate netstats
static void aggregateNetStats(void);

static uint64_t total_sent_bytes;
static uint64_t total_recv_bytes;

/// You must add timer name in same order as enum in .h file.
/// Leading spaces can be specified to show a hierarchy of timers.
char* timerName[numberOfTimers] = {
   "total",
   "loop",
   "timestep",
   "  position",
   "  velocity",
   "  redistribute",
   "    atomHalo",
   "  force",
   "   force 1",
   "   force 2",
   "   force 3",
   "    eamHalo",
   "commHalo",
   "commReduce"
};

/// Timer data collected.  Also facilitates computing averages and
/// statistics.
typedef struct TimersSt
{
   uint64_t start;     //!< call start time
   uint64_t total;     //!< current total time
   uint64_t count;     //!< current call count
   uint64_t elapsed;   //!< lap time
   struct net_event **device_list; //!< net_device list
   uint64_t net_tx_total;     //!< total tx'd
   uint64_t net_rx_total;     //!< total rx'd
   uint64_t net_total;        //!< total bytes
   uint64_t net_total_ranks;  //!< total across ranks
   uint64_t net_tx_bytes_start[10]; //!< network bytes tx'd
   uint64_t net_rx_bytes_start[10]; //!< network bytes rx'd
   uint64_t net_tx_bytes[10]; //!< network bytes tx'd
   uint64_t net_rx_bytes[10]; //!< network bytes rx'd

   int minRank;        //!< rank with min value
   int maxRank;        //!< rank with max value

   double minValue;    //!< min over ranks
   double maxValue;    //!< max over ranks
   double average;     //!< average over ranks
   double stdev;       //!< stdev across ranks

} Timers;

/// Global timing data collected.
typedef struct TimerGlobalSt
{
   double atomRate;       //!< average time (us) per atom per rank
   double atomAllRate;    //!< average time (us) per atom
   double atomsPerUSec;   //!< average atoms per time (us)
} TimerGlobal;

static Timers perfTimer[numberOfTimers];
static TimerGlobal perfGlobal;
static int papi_error = 0;

#ifdef RUN_PAPI
#define NUM_EVENTS 3
static int Events[NUM_EVENTS] = {PAPI_TOT_INS, PAPI_L3_TCM, PAPI_DP_OPS};
long_long values[NUM_EVENTS];
long_long all_procs[NUM_EVENTS];
#endif

static struct net_event **device_list = NULL;
extern int num_devices;

void profileStart(const enum TimerHandle handle)
{
   int ii, retval;
   perfTimer[handle].start = getTime();

#ifdef RUN_PAPI
   if (handle == totalTimer) {
     int retval = PAPI_start_counters(Events, NUM_EVENTS);
     if (retval != PAPI_OK) {
       printf("Error starting counters %s\n",PAPI_strerror(retval));
       printf("\n");
       papi_error = 1;
     }
   }
#endif

#ifdef PRINT_NET_STATS
   if (!netstats_init_complete || !device_list) {
      init(&device_list);
   }
   retval = read_stats(&device_list);
   if (!device_list || retval) {
      fprintf(screenOut, "Unable to init network counters\n");
   } else {
      for (ii = 0; ii < num_devices; ++ii) {
         perfTimer[handle].net_tx_bytes_start[ii] = device_list[ii]->tx_bytes;
         perfTimer[handle].net_rx_bytes_start[ii] = device_list[ii]->rx_bytes;
      }
   }
#endif
}

void profileStop(const enum TimerHandle handle)
{
   int ii, retval;
   perfTimer[handle].count += 1;
   uint64_t delta = getTime() - perfTimer[handle].start;
   perfTimer[handle].total += delta;
   perfTimer[handle].elapsed += delta;

#ifdef RUN_PAPI
   if (handle == totalTimer) {
     if (PAPI_accum_counters(values, NUM_EVENTS) != PAPI_OK)
       papi_error = 1;
   }
#endif

#ifdef PRINT_NET_STATS
   retval = read_stats(&device_list);
   if (retval) {
     fprintf(screenOut, "Unable to read netstats\n");
   } else {
     for (ii = 0; ii < num_devices; ++ii) {
        perfTimer[handle].net_tx_bytes[ii] = device_list[ii]->tx_bytes - perfTimer[handle].net_tx_bytes_start[ii];
        perfTimer[handle].net_tx_bytes[ii] = device_list[ii]->tx_bytes - perfTimer[handle].net_tx_bytes_start[ii];
     }
   }
#endif
}

/// \details
/// Return elapsed time (in seconds) since last call with this handle
/// and clear for next lap.
double getElapsedTime(const enum TimerHandle handle)
{
   double etime = getTick() * (double)perfTimer[handle].elapsed;
   perfTimer[handle].elapsed = 0;

   return etime;
}

/// \details
/// The report contains two blocks.  The upper block is performance
/// information for the printRank.  The lower block is statistical
/// information over all ranks.
void printPerformanceResults(int nGlobalAtoms, int printRate)
{
   // Aggregate net stats
   aggregateNetStats();
   // Collect timer statistics overall and across ranks
   timerStats();

   if (!printRank())
      return;

   uint64_t time = getTime();
   uint64_t time_new = getTime();

   while (time == time_new) {time_new = getTime();}


   // only print timers with non-zero values.
   double tick = getTick();
   double loopTime = perfTimer[loopTimer].total*tick;

   fprintf(screenOut, "\n\nTimings for Rank %d\n", getMyRank());
   fprintf(screenOut, "Timer resolution %"PRIu64" ms \n", (time_new-time));
   fprintf(screenOut, "        Timer        # Calls    Avg/Call (s)   Total (s)    %% Loop\n");
   fprintf(screenOut, "___________________________________________________________________\n");
   for (int ii=0; ii<numberOfTimers; ++ii)
   {
      double totalTime = perfTimer[ii].total*tick;
      if (perfTimer[ii].count > 0)
         fprintf(screenOut, "%-16s%12"PRIu64"     %12.4g      %8.4f    %8.2f\n",
                 timerName[ii],
                 perfTimer[ii].count,
                 totalTime/(double)perfTimer[ii].count,
                 totalTime,
                 totalTime/loopTime*100.0);
   }

   fprintf(screenOut, "\nTiming Statistics Across %d Ranks:\n", getNRanks());
   fprintf(screenOut, "        Timer        Rank: Min(s)       Rank: Max(s)      Avg(s)    Stdev(s)\n");
   fprintf(screenOut, "_____________________________________________________________________________\n");

   for (int ii = 0; ii < numberOfTimers; ++ii)
   {
      if (perfTimer[ii].count > 0)
         fprintf(screenOut, "%-16s%6d:%10.4f  %6d:%10.4f  %10.4f  %10.4f\n",
            timerName[ii],
            perfTimer[ii].minRank, perfTimer[ii].minValue*tick,
            perfTimer[ii].maxRank, perfTimer[ii].maxValue*tick,
            perfTimer[ii].average*tick, perfTimer[ii].stdev*tick);
   }

   fprintf(screenOut, "\nNetwork Statistics Across 1 Rank\n");
   fprintf(screenOut, "%-24s%10d\n","Bytes Sent 1 Rank", bytes_sent);
   fprintf(screenOut, "%-24s%10d\n","Bytes Received 1 Rank", bytes_recvd);

   fprintf(screenOut, "\nNetwork Statistics Across%d Ranks:\n", getNRanks());
   fprintf(screenOut, "%-24s%10"PRIu64"\n","Bytes Sent", total_sent_bytes);
   fprintf(screenOut, "%-24s%10"PRIu64"\n","Bytes Received", total_recv_bytes);

#ifdef RUN_PAPI
   fprintf(screenOut, "\nPapi statistics\n");
   fprintf(screenOut, "Papi total instruction count: %llu\n", values[0]);
   fprintf(screenOut, "Papi total l3 misses: %llu\n", values[1]);
   fprintf(screenOut, "Papi total dp operations: %llu\n", values[2]);
#endif

#ifdef PRINT_NET_STATS
   fprintf(screenOut, "\nNetwork Results:\n");
   fprintf(screenOut, "Performance Results For Rank %d:\n", getMyRank());
   for (int ii = 0; ii < numberOfTimers; ii++) {
      if (perfTimer[ii].count > 0)
      {
         double totalTime = perfTimer[ii].total*tick;
         fprintf(screenOut, "  Timer: %s\n", timerName[ii]);
/*         for (int jj = 0; jj < num_devices; ++jj) {*/
/*           fprintf(screenOut, "    NetDevice: %d\n", jj);//perfTimer[ii].device_list[jj]->device_name);*/
/*           fprintf(screenOut, "      TotalTime: %f\n", totalTime);*/
/*           fprintf(screenOut, "      RXCount: %llu bytes\n", perfTimer[ii].net_rx_bytes[jj]);*/
/*           fprintf(screenOut, "      TXCount: %llu bytes\n", perfTimer[ii].net_tx_bytes[jj]);*/
/*           double rate = ((double)perfTimer[ii].net_rx_bytes[jj]) / totalTime;*/
/*           fprintf(screenOut, "      AverageRXRate: %f bytes/second\n", rate);*/
/*           rate = ((double)perfTimer[ii].net_tx_bytes[jj]) / totalTime;*/
/*           fprintf(screenOut, "      AverageTXRate: %f bytes/second\n", rate);*/
/*         }*/
        fprintf(screenOut, "    NetDevice All\n");
        fprintf(screenOut, "      RXCount: %llu\n", perfTimer[ii].net_rx_total);
        fprintf(screenOut, "      TXCount: %llu\n", perfTimer[ii].net_tx_total);
        fprintf(screenOut, "      TransferCount: %llu\n", perfTimer[ii].net_total);
        fprintf(screenOut, "      AverageRate: %f\n", ((double)perfTimer[ii].net_total) / totalTime);
      }
      fprintf(screenOut, "    NetDevice All Ranks\n");
      fprintf(screenOut, "      TransferCount: %llu\n", perfTimer[ii].net_total_ranks);
      fprintf(screenOut, "      AverageRate: %f\n",
            ((double) perfTimer[ii].net_total_ranks) / (perfTimer[ii].average * (double)getNRanks()));
   }
#endif

   double atomsPerTask = nGlobalAtoms/(real_t)getNRanks();
   perfGlobal.atomRate = perfTimer[timestepTimer].average * tick * 1e6 /
      (atomsPerTask * perfTimer[timestepTimer].count * printRate);
   perfGlobal.atomAllRate = perfTimer[timestepTimer].average * tick * 1e6 /
      (nGlobalAtoms * perfTimer[timestepTimer].count * printRate);
   perfGlobal.atomsPerUSec = 1.0 / perfGlobal.atomAllRate;

   fprintf(screenOut, "\n---------------------------------------------------\n");
   fprintf(screenOut, " Average atom update rate:     %6.2f us/atom/task\n", perfGlobal.atomRate);
   fprintf(screenOut, "---------------------------------------------------\n\n");

   fprintf(screenOut, "\n---------------------------------------------------\n");
   fprintf(screenOut, " Average all atom update rate: %6.2f us/atom\n", perfGlobal.atomAllRate);
   fprintf(screenOut, "---------------------------------------------------\n\n");

   fprintf(screenOut, "\n---------------------------------------------------\n");
   fprintf(screenOut, " Average atom rate:            %6.2f atoms/us\n", perfGlobal.atomsPerUSec);
   fprintf(screenOut, "---------------------------------------------------\n\n");
}

void printPerformanceResultsYaml(FILE* file)
{
   // Aggregate net stats
   aggregateNetStats();

   if (! printRank())
      return;

   double tick = getTick();
   double loopTime = perfTimer[loopTimer].total*tick;

   fprintf(file,"\nPerformance Results:\n");
   fprintf(file, "  TotalRanks: %d\n", getNRanks());
   fprintf(file, "  ReportingTimeUnits: seconds\n");
   fprintf(file, "Performance Results For Rank %d:\n", getMyRank());
   for (int ii = 0; ii < numberOfTimers; ii++)
   {
      if (perfTimer[ii].count > 0)
      {
         double totalTime = perfTimer[ii].total*tick;
         fprintf(file, "  Timer: %s\n", timerName[ii]);
         fprintf(file, "    CallCount:  %"PRIu64"\n", perfTimer[ii].count);
         fprintf(file, "    AvgPerCall: %8.4f\n", totalTime/(double)perfTimer[ii].count);
         fprintf(file, "    Total:      %8.4f\n", totalTime);
         fprintf(file, "    PercentLoop: %8.2f\n", totalTime/loopTime*100);
      }
   }

   fprintf(file, "Performance Results Across Ranks:\n");
   for (int ii = 0; ii < numberOfTimers; ii++)
   {
      if (perfTimer[ii].count > 0)
      {
         fprintf(file, "  Timer: %s\n", timerName[ii]);
         fprintf(file, "    MinRank: %d\n", perfTimer[ii].minRank);
         fprintf(file, "    MinTime: %8.4f\n", perfTimer[ii].minValue*tick);
         fprintf(file, "    MaxRank: %d\n", perfTimer[ii].maxRank);
         fprintf(file, "    MaxTime: %8.4f\n", perfTimer[ii].maxValue*tick);
         fprintf(file, "    AvgTime: %8.4f\n", perfTimer[ii].average*tick);
         fprintf(file, "    StdevTime: %8.4f\n", perfTimer[ii].stdev*tick);
      }
   }

   fprintf(file,"Performance Global Update Rates:\n");
   fprintf(file, "  AtomUpdateRate:\n");
   fprintf(file, "    AverageRate: %6.2f\n", perfGlobal.atomRate);
   fprintf(file, "    Units: us/atom/task\n");
   fprintf(file, "  AllAtomUpdateRate:\n");
   fprintf(file, "    AverageRate: %6.2f\n", perfGlobal.atomAllRate);
   fprintf(file, "    Units: us/atom\n");
   fprintf(file, "  AtomRate:\n");
   fprintf(file, "    AverageRate: %6.2f\n", perfGlobal.atomsPerUSec);
   fprintf(file, "    Units: atoms/us\n");

   fprintf(file, "Network Results:\n");
   fprintf(file, "  TotalRanks: %d\n", getNRanks());
   fprintf(file, "  ReportingTransferUnits: bytes\n");
   fprintf(file, "Performance Results For Rank %d:\n", getMyRank());

 #ifdef RUN_PAPI
    fprintf(screenOut, "\nPapi statistics\n");
    fprintf(screenOut, "Papi total instruction count: %llu\n", values[0]);
    fprintf(screenOut, "Papi total l3 misses: %llu\n", values[1]);
    fprintf(screenOut, "Papi total dp operations: %llu\n", values[2]);
 #endif

#ifdef PRINT_NET_STATS
   for (int ii = 0; ii < numberOfTimers; ii++) {
      if (perfTimer[ii].count > 0)
      {
         double totalTime = perfTimer[ii].total*tick;
         fprintf(file, "  Timer: %s\n", timerName[ii]);
/*         for (int jj = 0; jj < num_devices; ++jj) {*/
/*           fprintf(file, "    NetDevice: %d\n", jj);//perfTimer[ii].device_list[jj]->device_name);*/
/*           fprintf(file, "      RXCount: %llu\n", perfTimer[ii].net_rx_bytes[jj]);*/
/*           fprintf(file, "      TXCount: %llu\n", perfTimer[ii].net_tx_bytes[jj]);*/
/*           double rate = ((double)perfTimer[ii].net_rx_bytes[jj]) / totalTime;*/
/*           fprintf(file, "      AverageRXRate: %f\n", rate);*/
/*           rate = ((double)perfTimer[ii].net_tx_bytes[jj]) / totalTime;*/
/*           fprintf(file, "      AverageTXRate: %f\n", rate);*/
/*        }*/
        fprintf(file, "    NetDevice All\n");
        fprintf(file, "      RXCount: %"PRIu64"\n", perfTimer[ii].net_rx_total);
        fprintf(file, "      TXCount: %"PRIu64"\n", perfTimer[ii].net_tx_total);
        fprintf(file, "      TransferCount: %"PRIu64"\n", perfTimer[ii].net_total);
        fprintf(file, "      AverageRate: %f\n", (perfTimer[ii].net_total) / totalTime);
      }
      fprintf(file, "    NetDevice All Ranks\n");
      fprintf(file, "      TransferCount: %"PRIu64"\n", perfTimer[ii].net_total_ranks);
      fprintf(file, "      AverageRate: %f\n",
            ((double) perfTimer[ii].net_total_ranks) / (perfTimer[ii].average * (double)getNRanks()));
   }
#endif
   fprintf(file, "\n");
}

/// Returns current time as a 64-bit integer.  This portable version
/// returns the number of microseconds since mindight, Jamuary 1, 1970.
/// Hence, timing data will have a resolution of 1 microsecond.
/// Platforms with access to calls with lower latency or higher
/// resolution (such as a cycle counter) may wish to replace this
/// implementation and change the conversion factor in getTick as
/// appropriate.
/// \see getTick for the conversion factor between the integer time
/// units of this function and seconds.
static uint64_t getTime(void)
{
   struct timeval ptime;
   uint64_t t = 0;
   gettimeofday(&ptime, (struct timezone *)NULL);
   t = ((uint64_t)1000000)*(uint64_t)ptime.tv_sec + (uint64_t)ptime.tv_usec;

   return t;
}

/// Returns the factor for converting the integer time reported by
/// getTime into seconds.  The portable getTime returns values in units
/// of microseconds so the conversion is simply 1e-6.
/// \see getTime
static double getTick(void)
{
   double seconds_per_cycle = 1.0e-6;
   return seconds_per_cycle;
}

/// Collect timer statistics across ranks.
void timerStats(void)
{
   double sendBuf[numberOfTimers], recvBuf[numberOfTimers];

   // Determine average of each timer across ranks
   for (int ii = 0; ii < numberOfTimers; ii++)
      sendBuf[ii] = (double)perfTimer[ii].total;
   addDoubleParallel(sendBuf, recvBuf, numberOfTimers);

   for (int ii = 0; ii < numberOfTimers; ii++)
      perfTimer[ii].average = recvBuf[ii] / (double)getNRanks();

   // Determine total transfer across all ranks
   for (int ii = 0; ii < numberOfTimers; ii++)
      sendBuf[ii] = (double)perfTimer[ii].net_total;
   addDoubleParallel(sendBuf, recvBuf, numberOfTimers);

   for (int ii = 0; ii < numberOfTimers; ii++)
      perfTimer[ii].net_total_ranks = recvBuf[ii];

   // totalbytes sent/recvd
   int tmp_recv;
   addIntParallel(&bytes_sent, &tmp_recv, 1);
   total_sent_bytes = tmp_recv;
   addIntParallel(&bytes_recvd, &tmp_recv, 1);
   total_recv_bytes = tmp_recv;

   // Determine min and max across ranks and which rank
   RankReduceData reduceSendBuf[numberOfTimers], reduceRecvBuf[numberOfTimers];
   for (int ii = 0; ii < numberOfTimers; ii++)
   {
      reduceSendBuf[ii].val = (double)perfTimer[ii].total;
      reduceSendBuf[ii].rank = getMyRank();
   }
   minRankDoubleParallel(reduceSendBuf, reduceRecvBuf, numberOfTimers);
   for (int ii = 0; ii < numberOfTimers; ii++)
   {
      perfTimer[ii].minValue = reduceRecvBuf[ii].val;
      perfTimer[ii].minRank = reduceRecvBuf[ii].rank;
   }
   maxRankDoubleParallel(reduceSendBuf, reduceRecvBuf, numberOfTimers);
   for (int ii = 0; ii < numberOfTimers; ii++)
   {
      perfTimer[ii].maxValue = reduceRecvBuf[ii].val;
      perfTimer[ii].maxRank = reduceRecvBuf[ii].rank;
   }

   // Determine standard deviation
   for (int ii = 0; ii < numberOfTimers; ii++)
   {
      double temp = (double)perfTimer[ii].total - perfTimer[ii].average;
      sendBuf[ii] = temp * temp;
   }
   addDoubleParallel(sendBuf, recvBuf, numberOfTimers);
   for (int ii = 0; ii < numberOfTimers; ii++)
   {
      perfTimer[ii].stdev = sqrt(recvBuf[ii] / (double) getNRanks());
   }
}

// aggregate netstats
static void aggregateNetStats(void)
{
  for (int ii = 0; ii < numberOfTimers; ii++) {
    if (perfTimer[ii].count > 0)
    {
       perfTimer[ii].net_rx_total = 0;
       perfTimer[ii].net_tx_total = 0;
       perfTimer[ii].net_total = 0;
       for (int jj = 0; jj < num_devices; ++jj) {
           perfTimer[ii].net_rx_total += perfTimer[ii].net_rx_bytes[jj];
           perfTimer[ii].net_tx_total += perfTimer[ii].net_tx_bytes[jj];
           perfTimer[ii].net_total = perfTimer[ii].net_total + perfTimer[ii].net_rx_bytes[jj] + perfTimer[ii].net_tx_bytes[jj];
       }
    }
  }
}
